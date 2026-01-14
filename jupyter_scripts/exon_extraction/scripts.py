#!/usr/bin/env python3
import os
import re
import sys
from Bio import SeqIO

def parse_gtf(gtf_file):
    """
    解析 GTF 文件中所有 exon 记录，返回字典：
      {transcript_id: [(exon_number, start, end, strand), ...]}
    如果某个 exon 没有 exon_number，则 exon_number 设为 None。
    """
    transcripts = {}
    with open(gtf_file) as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature != "exon":
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            strand = fields[6]
            attr_field = fields[8]
            m = re.search(r'transcript_id "([^"]+)"', attr_field)
            if not m:
                continue
            transcript_id = m.group(1)
            m2 = re.search(r'exon_number "([^"]+)"', attr_field)
            exon_number = int(m2.group(1)) if m2 else None

            exon_info = (exon_number, start, end, strand)
            transcripts.setdefault(transcript_id, []).append(exon_info)
    # 对每个转录本中的 exon 进行排序（优先依据 exon_number，否则依据起始位置）
    for tid in transcripts:
        if all(x[0] is not None for x in transcripts[tid]):
            transcripts[tid].sort(key=lambda x: x[0])
        else:
            transcripts[tid].sort(key=lambda x: x[1])
    return transcripts

def partition_aligned_sequence_with_coords(aligned_seq, exon_lengths):
    """
    根据原始 exon 长度（exon_lengths 列表），遍历包含 gap 的比对序列，
    按累计非 gap 字符数将比对序列分段，返回两个列表：
      - segments：每个 exon 对应的比对序列片段（包含 gap）
      - seg_coords：每个片段在对齐序列中的起始和结束坐标（1-based），即 (start, end)
    """
    segments = []
    seg_coords = []
    # 预先计算每个 exon 的累计目标值，例如 [L1, L1+L2, L1+L2+L3, ...]
    cumulative_targets = []
    cum = 0
    for length in exon_lengths:
        cum += length
        cumulative_targets.append(cum)
    raw_count = 0      # 累计非 gap 字符数（跨所有 exon）
    current_segment = []
    seg_start = None   # 当前片段在对齐序列中的起始位置（1-based）
    exon_index = 0     # 当前正在处理第几个 exon（从 0 开始）
    for i, char in enumerate(aligned_seq, start=1):
        if seg_start is None:
            seg_start = i
        current_segment.append(char)
        if char != '-':
            raw_count += 1
        # 判断是否达到当前 exon 的累计目标
        if exon_index < len(cumulative_targets) and raw_count == cumulative_targets[exon_index]:
            seg_end = i
            segments.append("".join(current_segment))
            seg_coords.append((seg_start, seg_end))
            # 准备下一个 exon
            current_segment = []
            seg_start = None
            exon_index += 1
    return segments, seg_coords

def main(mafft_file):
    # 读取 MAFFT 对齐结果（FASTA格式）
    records = list(SeqIO.parse(mafft_file, "fasta"))
    if not records:
        sys.stderr.write("未能读取到任何序列！\n")
        sys.exit(1)
    # 假定所有序列经过 MAFFT 对齐后长度相同
    alignment_length = len(records[0].seq)
    
    # 用于缓存各个物种对应的 GTF 解析结果： {species: {transcript_id: exon_info_list}}
    species_to_transcripts = {}
    # 收集每条记录的处理结果，便于最后制作对齐坐标汇总表
    summary_data = []

    for record in records:
        header = record.id  # 例如 "Bicyclus_anynana_Bany005764.1"
        # 按 "_" 分割，假定 transcriptID 为最后一部分，其余部分构成 species 名称
        parts = header.split("_")
        if len(parts) < 2:
            sys.stderr.write(f"错误的 header 格式: {header}\n")
            continue
        species = "_".join(parts[:2])
        transcript_id = "_".join(parts[2:])

        # 加载对应物种的 GTF 文件（只加载一次）
        if species not in species_to_transcripts:
            gtf_file = f"{species}.exon.gtf"
            if not os.path.exists(gtf_file):
                sys.stderr.write(f"未找到 GTF 文件：{gtf_file}\n")
                species_to_transcripts[species] = None
            else:
                transcripts = parse_gtf(gtf_file)
                species_to_transcripts[species] = transcripts

        if species_to_transcripts[species] is None:
            sys.stderr.write(f"跳过 {header}，因为找不到对应的 GTF 文件。\n")
            continue

        transcripts = species_to_transcripts[species]
        if transcript_id not in transcripts:
            sys.stderr.write(f"在 {species}.exon.gtf 中未找到转录本 {transcript_id} 的 exon 信息\n")
            continue

        exons = transcripts[transcript_id]
        # 计算每个 exon 的原始长度（注意 GTF 坐标为 1-based）
        exon_lengths = [exon[2] - exon[1] + 1 for exon in exons]
        aligned_seq = str(record.seq)
        segments, coords = partition_aligned_sequence_with_coords(aligned_seq, exon_lengths)

        # 输出当前记录的详细信息（展示每个 exon 在对齐序列中的坐标）
        print(f">{header}")
        print(f"对齐序列长度： {alignment_length}")
        for i, (start, end) in enumerate(coords, start=1):
            print(f"exon{i}: 对齐坐标 {start} - {end}，片段： {segments[i-1]}")
        print()

        summary_data.append({
            "header": header,
            "species": species,
            "transcript": transcript_id,
            "coords": coords,
        })

    # --- 输出在统一对齐坐标系下不同物种（转录本）间外显子相对位置的汇总表 ---
    print("不同物种/转录本在统一对齐坐标系下的 exon 坐标：")
    # 确定所有记录中 exon 数目的最大值（通常同源转录本的 exon 数目应相同）
    max_exons = max(len(item["coords"]) for item in summary_data)
    # 表头
    header_line = f"{'Species_Transcript':30}"
    for i in range(1, max_exons+1):
        header_line += f"Exon{i:>15}"
    print(header_line)
    # 表体
    for item in summary_data:
        row = f"{item['header']:30}"
        for (start, end) in item["coords"]:
            row += f"{(str(start)+'-'+str(end)):>15}"
        print(row)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("用法: python script.py <mafft_alignment.fasta>\n")
        sys.exit(1)
    main(sys.argv[1])

