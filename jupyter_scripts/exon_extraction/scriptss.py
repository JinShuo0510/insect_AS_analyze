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
    # 汇总每条记录的处理结果，后续用于同源外显子组分析
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
        # 计算每个 exon 的原始长度（GTF 坐标为 1-based）
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
    
    # 输出统一对齐坐标下各物种/转录本的 exon 坐标汇总表（每条记录一行）
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
    
    # --- 同源外显子组筛选：放宽要求，同一物种中只需要有一个转录本满足即可 ---
    # 首先按物种分组（同一物种可能有多个转录本）
    species_records = {}
    for item in summary_data:
        species_records.setdefault(item["species"], []).append(item)
    species_list = sorted(species_records.keys())
    # 对于每个物种，计算该物种所能提供的 exon 组数（取其所有转录本中外显子数的最大值）
    species_exon_availability = {}
    for sp in species_list:
        exon_counts = [len(item["coords"]) for item in species_records[sp]]
        species_exon_availability[sp] = max(exon_counts)
    # 只有当所有物种都提供了该 exon 组时，才能考虑
    min_exons = min(species_exon_availability.values())
    
    print("\n同源外显子组筛选（放宽要求：同一物种下只需有一个转录本满足）：")
    print("{:>10} {:>20} {:>20} {:>15}".format("ExonGroup", "Intersection", "Union", "Overlap Ratio"))
    
    homologous_groups = []  # 用于保存满足条件的 exon 组信息
    # 针对 exon 组索引 0 到 min_exons-1 进行分析
    for i in range(min_exons):
        # 对于每个物种，收集该物种中所有转录本关于 exon i 的候选：
        # 每个候选记录为 (coordinate, transcript_header)
        candidates_per_species = {}
        valid = True
        for sp in species_list:
            cand = []
            for rec in species_records[sp]:
                if len(rec["coords"]) > i:
                    cand.append( (rec["coords"][i], rec["header"]) )
            if not cand:
                valid = False
                sys.stderr.write(f"物种 {sp} 无 exon{i+1} 信息，跳过 exon组{i+1}\n")
                break
            candidates_per_species[sp] = cand
        if not valid:
            continue
        
        # 遍历所有物种候选的笛卡尔积，寻找最佳组合（即使整体交集/并集比值最大）
        best_ratio = 0
        best_combo = None  # 记录最佳组合，格式为 { species: (coordinate, header) }
        # 为了保证顺序，按 species_list 的顺序构造候选列表
        candidate_lists = [ candidates_per_species[sp] for sp in species_list ]
        for combo in itertools.product(*candidate_lists):
            # combo 的顺序与 species_list对应，元素为 (coordinate, header)
            intervals = [ cand[0] for cand in combo ]  # 仅取坐标 (start, end)
            _, _, _, _, ratio = compute_intersection_union(intervals)
            if ratio > best_ratio:
                best_ratio = ratio
                best_combo = { sp: combo[j] for j, sp in enumerate(species_list) }
        # 若最佳组合的比值超过阈值，则认为该 exon 组为同源外显子组
        if best_ratio > 0.6:
            # 重新计算交集与并集信息
            best_intervals = [ best_combo[sp][0] for sp in species_list ]
            inter_start, inter_end, union_start, union_end, _ = compute_intersection_union(best_intervals)
            homologous_groups.append({
                "group": i+1,
                "intersection": (inter_start, inter_end),
                "union": (union_start, union_end),
                "ratio": best_ratio,
                "combo": best_combo
            })
            print("{:>10} {:>20} {:>20} {:>15.3f}".format("Exon"+str(i+1),
                                                          f"{inter_start}-{inter_end}",
                                                          f"{union_start}-{union_end}",
                                                          best_ratio))
            for sp in species_list:
                coord, hdr = best_combo[sp]
                print("{:>30} : {}  (Exon{}坐标: {}-{})".format(sp, hdr, i+1, coord[0], coord[1]))
            print()
        else:
            print(f"Exon组{i+1} 未达到阈值（最佳比值：{best_ratio:.3f}）")
    
    if not homologous_groups:
        print("未找到满足条件的同源外显子组。")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("用法: python script.py <mafft_alignment.fasta>\n")
        sys.exit(1)
    main(sys.argv[1])

