#!/usr/bin/env python3
import os
import re
import sys
import itertools
from Bio import SeqIO

def parse_gtf(gtf_file):
    """
    解析 GTF 文件中所有 exon 记录，返回字典：
      {transcript_id: [(exon_number, start, end, strand), ...]}
    若某个 exon 没有 exon_number，则 exon_number 设为 None。
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
    # 对每个转录本的 exon 按 exon_number（若可用，否则按 start）排序
    for tid in transcripts:
        if all(x[0] is not None for x in transcripts[tid]):
            transcripts[tid].sort(key=lambda x: x[0])
        else:
            transcripts[tid].sort(key=lambda x: x[1])
    return transcripts

def partition_aligned_sequence_with_coords(aligned_seq, exon_lengths):
    """
    根据原始 exon 长度（exon_lengths 列表），遍历带 gap 的比对序列，
    按累计非 gap 字符数将比对序列分段，返回两个列表：
      - segments：每个 exon 对应的原始比对片段（含 gap，可供参考）
      - seg_coords：每个 exon 经过修剪后在对齐序列中的坐标 (start, end)，
                    修剪时将片段开头和结尾处的 gap 忽略掉（即取第一个、最后一个非 gap 字符的位置）。
    """
    segments = []
    seg_coords = []
    # 计算累计目标值，如 [L1, L1+L2, L1+L2+L3, ...]
    cumulative_targets = []
    cum = 0
    for length in exon_lengths:
        cum += length
        cumulative_targets.append(cum)
    raw_count = 0      # 跨所有 exon 的累计非 gap 字符数
    current_segment_chars = []
    seg_start = None   # 当前片段在对齐序列中的起始位置（1-based）
    exon_index = 0     # 当前处理的 exon 编号（从 0 开始）
    for i, char in enumerate(aligned_seq, start=1):
        if seg_start is None:
            seg_start = i
        current_segment_chars.append(char)
        if char != '-':
            raw_count += 1
        if exon_index < len(cumulative_targets) and raw_count == cumulative_targets[exon_index]:
            seg_end = i
            segment_str = "".join(current_segment_chars)
            # 修剪：去除片段首尾的 gap
            left_offset = 0
            while left_offset < len(segment_str) and segment_str[left_offset] == '-':
                left_offset += 1
            right_offset = len(segment_str) - 1
            while right_offset >= 0 and segment_str[right_offset] == '-':
                right_offset -= 1
            adjusted_start = seg_start + left_offset if left_offset < len(segment_str) else seg_start
            adjusted_end   = seg_start + right_offset if right_offset >= 0 else seg_end
            segments.append(segment_str)
            seg_coords.append((adjusted_start, adjusted_end))
            # 重置，准备处理下一个 exon
            current_segment_chars = []
            seg_start = None
            exon_index += 1
    return segments, seg_coords

def compute_intersection_union(intervals):
    """
    给定多个区间（列表，每个区间为 (start, end)），计算交集与并集。
      - 交集：start = max(starts), end = min(ends)，若无重叠则长度为 0；
      - 并集：start = min(starts), end = max(ends)。
    返回：(inter_start, inter_end, union_start, union_end, overlap_ratio)
    """
    starts = [s for s, _ in intervals]
    ends   = [e for _, e in intervals]
    union_start = min(starts)
    union_end   = max(ends)
    inter_start = max(starts)
    inter_end   = min(ends)
    inter_len = max(0, inter_end - inter_start + 1)
    union_len = union_end - union_start + 1
    ratio = inter_len / union_len if union_len > 0 else 0
    return inter_start, inter_end, union_start, union_end, ratio

def main(mafft_file):
    # 读取 MAFFT 对齐结果（FASTA格式）
    records = list(SeqIO.parse(mafft_file, "fasta"))
    if not records:
        sys.stderr.write("未能读取到任何序列！\n")
        sys.exit(1)
    alignment_length = len(records[0].seq)
    
    # 用于缓存各物种的 GTF解析结果： {species: {transcript_id: exon_info_list}}
    species_to_transcripts = {}
    # 汇总每条记录的处理结果，后续用于统一对齐坐标和同源外显子组分析
    summary_data = []  # 每个元素为字典，包含 "header", "species", "transcript", "coords"
    
    for record in records:
        header = record.id  # 例如 "Bicyclus_anynana_Bany005764.1"
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
        # 计算每个 exon 的原始长度（GTF 坐标为 1-based，长度 = end - start + 1）
        exon_lengths = [exon[2] - exon[1] + 1 for exon in exons]
        aligned_seq = str(record.seq)
        segments, coords = partition_aligned_sequence_with_coords(aligned_seq, exon_lengths)
        # 输出当前记录的详细信息
        print(f">{header}")
        print(f"对齐序列长度： {alignment_length}")
        for i, (start, end) in enumerate(coords, start=1):
            print(f"exon{i}: 对齐坐标 {start} - {end}，原始片段： {segments[i-1]}")
        print()
        summary_data.append({
            "header": header,
            "species": species,
            "transcript": transcript_id,
            "coords": coords,
        })
    
    # --- 输出统一对齐坐标下各物种/转录本的 exon 坐标汇总表 ---
    print("各物种/转录本在统一对齐坐标系下的 exon 坐标：")
    max_exons = max(len(item["coords"]) for item in summary_data)
    header_line = f"{'Species_Transcript':30}"
    for i in range(1, max_exons+1):
        header_line += f"Exon{i:>15}"
    print(header_line)
    for item in summary_data:
        row = f"{item['header']:30}"
        for (start, end) in item["coords"]:
            row += f"{(str(start)+'-'+str(end)):>15}"
        print(row)
    
    # --- 同源外显子组筛选（基于参考物种的坐标匹配） ---
    #
    # 思路：选取一个参考物种，使用该物种中代表性转录本的各 exon 作为参考，
    #       对于其他物种，从其所有转录本中选取中点与参考 exon 最接近的外显子，
    #       然后计算这些候选外显子的交集与并集比值，若 > 0.6，则认为是一组直系同源外显子组。
    #
    # 先按物种分组
    species_records = {}
    for item in summary_data:
        species_records.setdefault(item["species"], []).append(item)
    species_list = sorted(species_records.keys())
    if not species_list:
        sys.stderr.write("未找到任何物种记录！\n")
        sys.exit(1)
    
    # 选择参考物种：这里取排序后第一个
    ref_species = species_list[0]
    # 从该物种中选择一个代表性转录本（例如取外显子数最多的）
    ref_transcript = max(species_records[ref_species], key=lambda rec: len(rec["coords"]))
    
    print("\n基于参考物种【{}】的同源外显子组筛选：".format(ref_species))
    homologous_groups = []  # 保存满足条件的组
    # 对参考转录本中每个外显子进行处理
    for idx, ref_interval in enumerate(ref_transcript["coords"]):
        ref_start, ref_end = ref_interval
        ref_mid = (ref_start + ref_end) / 2
        group_candidates = { ref_species: (ref_interval, ref_transcript["header"]) }
        # 对于其他物种，选择中点与参考中点最接近的外显子（遍历所有转录本中所有外显子）
        valid_group = True
        for sp in species_list:
            if sp == ref_species:
                continue
            best_diff = None
            best_candidate = None  # (interval, transcript_header)
            for rec in species_records[sp]:
                for interval in rec["coords"]:
                    mid = (interval[0] + interval[1]) / 2
                    diff = abs(mid - ref_mid)
                    if best_diff is None or diff < best_diff:
                        best_diff = diff
                        best_candidate = (interval, rec["header"])
            if best_candidate is None:
                sys.stderr.write(f"物种 {sp} 无外显子候选，跳过参考转录本的 exon{idx+1}\n")
                valid_group = False
                break
            group_candidates[sp] = best_candidate
        if not valid_group:
            continue
        # 计算这一组候选外显子的交集与并集
        intervals = [ group_candidates[sp][0] for sp in group_candidates ]
        inter_start, inter_end, union_start, union_end, ratio = compute_intersection_union(intervals)
        if ratio > 0.6:
            homologous_groups.append({
                "ref_exon_index": idx+1,
                "intersection": (inter_start, inter_end),
                "union": (union_start, union_end),
                "ratio": ratio,
                "group_candidates": group_candidates,
            })
            print("参考转录本 exon{}：".format(idx+1))
            print("  交集：{}-{}，并集：{}-{}，重叠比：{:.3f}".format(inter_start, inter_end, union_start, union_end, ratio))
            for sp in species_list:
                interval, hdr = group_candidates[sp]
                print("  {} : {}  (外显子坐标: {}-{})".format(sp, hdr, interval[0], interval[1]))
            print()
        else:
            print("参考转录本 exon{}：重叠比 {:.3f} 未达到阈值".format(idx+1, ratio))
    
    if not homologous_groups:
        print("未找到满足条件的同源外显子组。")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("用法: python script.py <mafft_alignment.fasta>\n")
        sys.exit(1)
    main(sys.argv[1])

