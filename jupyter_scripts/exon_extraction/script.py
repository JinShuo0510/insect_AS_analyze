#!/usr/bin/env python3
import os
import re
import sys
from Bio import SeqIO

def parse_gtf(gtf_file):
    """
    解析 GTF 文件中所有的 exon 记录，
    返回一个字典：{transcript_id: [(exon_number, start, end, strand), ...]}
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
            # 提取 transcript_id ，形如 transcript_id "XXX";
            m = re.search(r'transcript_id "([^"]+)"', attr_field)
            if not m:
                continue
            transcript_id = m.group(1)
            # 尝试提取 exon_number
            m2 = re.search(r'exon_number "([^"]+)"', attr_field)
            exon_number = int(m2.group(1)) if m2 else None

            exon_info = (exon_number, start, end, strand)
            transcripts.setdefault(transcript_id, []).append(exon_info)
    # 对每个 transcript 的 exon 按照 exon_number 排序（若 exon_number 可用，否则按 start 升序）
    for tid in transcripts:
        if all(x[0] is not None for x in transcripts[tid]):
            transcripts[tid].sort(key=lambda x: x[0])
        else:
            transcripts[tid].sort(key=lambda x: x[1])
    return transcripts

def partition_aligned_sequence(aligned_seq, exon_lengths):
    """
    给定比对序列 aligned_seq（包含 gap '-'）以及各个 exon 的原始长度 exon_lengths（列表），
    根据累计非 gap 字符数将 aligned_seq 分段，返回分段后的列表，每段对应一个 exon。
    """
    segments = []
    current_segment = []
    raw_count = 0  # 计数非 gap 字符个数
    exon_index = 0
    # target 为当前 exon 截取时非 gap字符的累计目标值
    target = exon_lengths[exon_index] if exon_lengths else 0

    for char in aligned_seq:
        # 如果已经分完所有 exon，则忽略后面可能存在的多余 gap（例如末端 gap）
        if exon_index >= len(exon_lengths):
            break
        current_segment.append(char)
        if char != '-':
            raw_count += 1
            if raw_count == target:
                # 到达当前 exon 的边界，保存当前 segment
                segments.append("".join(current_segment))
                current_segment = []
                exon_index += 1
                if exon_index < len(exon_lengths):
                    target += exon_lengths[exon_index]
    return segments

def main(mafft_file):
    # 读取 mafft fasta 文件
    records = list(SeqIO.parse(mafft_file, "fasta"))
    # 用于缓存各个物种对应的 GTF 解析结果，格式：{species: {transcript_id: exon_info_list}}
    species_to_transcripts = {}

    for record in records:
        header = record.id  # 如 "Bicyclus_anynana_Bany005764.1"
        # 解析 species 和 transcript_id：
        # 假设 transcript_id 为最后一个 "_" 后面的部分，species 为其余部分
        parts = header.split("_")
        if len(parts) < 2:
            sys.stderr.write(f"错误的 header 格式: {header}\n")
            continue
        print(parts[1:])
        species = "_".join(parts[:2])
        transcript_id = "_".join(parts[2:])

        # 加载该物种对应的 GTF 文件（只加载一次）
        if species not in species_to_transcripts:
            gtf_file = f"{species}.exon.gtf"
            if not os.path.exists(gtf_file):
                sys.stderr.write(f"未找到 GTF 文件：{gtf_file}\n")
                species_to_transcripts[species] = None
            else:
                transcripts = parse_gtf(gtf_file)
                species_to_transcripts[species] = transcripts

        # 如果未能加载 GTF，则跳过
        if species_to_transcripts[species] is None:
            sys.stderr.write(f"跳过 {header}，因为找不到对应的 GTF 文件。\n")
            continue

        transcripts = species_to_transcripts[species]
        if transcript_id not in transcripts:
            sys.stderr.write(f"在 {species}.exon.gtf 中未找到转录本 {transcript_id} 的 exon 信息\n")
            continue

        exons = transcripts[transcript_id]
        # 计算各 exon 的原始长度（注意：GTF 中坐标为 1-based，长度计算公式为 end - start + 1）
        exon_lengths = [exon[2] - exon[1] + 1 for exon in exons]
        aligned_seq = str(record.seq)
        segments = partition_aligned_sequence(aligned_seq, exon_lengths)

        # 输出结果：先输出序列名，再逐条输出 exon 片段
        print(f">{header}")
        for i, seg in enumerate(segments, 1):
            print(f"exon{i}: {seg}")
        print()  # 空行分隔不同记录

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("用法: python script.py <mafft_alignment.fasta>\n")
        sys.exit(1)
    main(sys.argv[1])

