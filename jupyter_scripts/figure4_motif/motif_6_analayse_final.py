import itertools
from collections import defaultdict, Counter
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import pandas as pd
import os
import glob  # 用于查找文件
import time  # 用于计时
import re    # 用于从文件名提取组织和方向

SIGNIFICANCE_THRESHOLD = 0.05  # 常用的 FDR 阈值
BASES = 'ACGT'  # 定义有效碱基

def parse_fasta(filename):
    """解析 FASTA 文件，按 'upstream' 和 'downstream' 分组  
       header 格式示例: {GeneID}_upstream::chr:start-end 或 {GeneID}_downstream::chr:start-end"""
    sequences = defaultdict(list)
    current_label = None
    current_seq = []

    if not os.path.exists(filename):
        print(f"警告：找不到 FASTA 文件 '{filename}'，跳过。")
        return None  # 返回 None 表示文件有问题

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # 保存上一条序列 (如果存在)
                    if current_label and current_seq:
                        seq_str = "".join(current_seq).upper()
                        if any(base in BASES for base in seq_str):
                            sequences[current_label].append(seq_str)
                    # 解析新的 header  
                    # 例如：Gland_upstream::NC_085107.1:4047384-4047584
                    header = line[1:].strip()
                    gene_info = header.split("::")[0]
                    if gene_info.endswith('_upstream'):
                        current_label = 'upstream'
                    elif gene_info.endswith('_downstream'):
                        current_label = 'downstream'
                    else:
                        current_label = None  # 无法识别则忽略
                    current_seq = []  # 重置序列
                elif current_label:  # 仅在有效标签下添加序列
                    current_seq.append(line)
            # 添加最后一条序列
            if current_label and current_seq:
                seq_str = "".join(current_seq).upper()
                if any(base in BASES for base in seq_str):
                    sequences[current_label].append(seq_str)
    except Exception as e:
        print(f"错误：解析文件 {filename} 时出错: {e}")
        return None

    if not sequences.get('upstream') or not sequences.get('downstream'):
         print(f"警告：文件 {filename} 中上游或下游序列组为空，可能无法进行比较。")
    return dict(sequences)

def generate_hexamers(bases=BASES):
    """生成所有可能的六核苷酸序列"""
    return [''.join(p) for p in itertools.product(bases, repeat=6)]

def count_hexamers_in_list(sequences, all_hexamers_set):
    """统计列表中所有序列的 hexamer 计数和总有效位置"""
    hexamer_counts = Counter()
    total_valid_positions = 0
    valid_bases = set(BASES)
    skipped_hexamers = 0

    start_time = time.time()
    for seq in sequences:
        seq_len = len(seq)
        if seq_len < 6:
            continue  # 序列太短无法形成 hexamer
        for i in range(seq_len - 5):
            hexamer = seq[i:i+6]
            if all(base in valid_bases for base in hexamer):
                if hexamer in all_hexamers_set:
                    hexamer_counts[hexamer] += 1
                total_valid_positions += 1
            else:
                skipped_hexamers += 1

    elapsed = time.time() - start_time
    return hexamer_counts, total_valid_positions

def run_enrichment_analysis(target_label, background_label, target_counts, target_total_pos,
                            background_counts, background_total_pos, all_hexamers):
    """执行 Fisher 精确检验和多重检验校正"""
    results = []
    raw_p_values = []

    if target_total_pos == 0 or background_total_pos == 0:
        print(f"错误：{target_label} 或 {background_label} 组的总有效位置数为 0，无法进行分析。")
        return pd.DataFrame()

    for hexamer in all_hexamers:
        count_target = target_counts.get(hexamer, 0)
        count_background = background_counts.get(hexamer, 0)
        # 使用各自的总位置数
        table = [
            [count_target, count_background],
            [target_total_pos - count_target, background_total_pos - count_background]
        ]
        try:
            _, p_value = stats.fisher_exact(table, alternative='greater')
            raw_p_values.append(p_value)
            results.append({
                'hexamer': hexamer,
                f'count_{target_label}': count_target,
                f'count_{background_label}': count_background,
                'raw_p_value': p_value
            })
        except ValueError as e:
            print(f"警告：计算 hexamer '{hexamer}' 时出错: {e}。列联表: {table}。跳过。")
            raw_p_values.append(1.0)
            results.append({
                'hexamer': hexamer,
                f'count_{target_label}': count_target,
                f'count_{background_label}': count_background,
                'raw_p_value': 1.0
            })

    if not raw_p_values:
         print("警告：没有有效的 p 值进行多重检验校正。")
         return pd.DataFrame(results)

    try:
        reject, corrected_p_values, _, _ = multipletests(raw_p_values, method='fdr_bh', alpha=SIGNIFICANCE_THRESHOLD)
        for i, res in enumerate(results):
            res['corrected_p_value (q-value)'] = corrected_p_values[i]
            res['significant (FDR < 0.05)'] = reject[i]
    except Exception as e:
        print(f"错误：多重检验校正时出错: {e}")
        for res in results:
             res['corrected_p_value (q-value)'] = float('nan')
             res['significant (FDR < 0.05)'] = False

    return pd.DataFrame(results)

def extract_tissue_info(filename):
    """
    从文件名提取组织和方向信息  
    文件名格式：{tissue}_[up/down]_filtered_exons.fa  
    例如：Head_down_filtered_exons.fa  → tissue: Head, direction: down
    """
    basename = os.path.basename(filename)
    match = re.match(r'^(.*?)_((?:up)|(?:down))_filtered_exons\.fa$', basename)
    if match:
         tissue = match.group(1)
         file_direction = match.group(2)
         return tissue, file_direction
    else:
         print(f"警告：无法从文件名 '{filename}' 提取组织和方向信息。将使用完整文件名作为标识。")
         return basename, "NA"

# --- 主程序 ---
if __name__ == "__main__":
    start_overall_time = time.time()

    # 查找所有匹配的前景文件
    fasta_files = glob.glob('*_filtered_exons.fa')
    if not fasta_files:
        print("错误：在当前目录下未找到匹配 '*_filtered_exons.fa' 模式的文件。")
        exit()

    print(f"找到 {len(fasta_files)} 个匹配的前景 FASTA 文件。开始处理...")

    # 只生成一次 hexamers 列表
    all_hexamers_list = generate_hexamers()
    all_hexamers_set = set(all_hexamers_list)
    print(f"已生成 {len(all_hexamers_list)} 种 Hexamers。")

    all_enriched_list = []

    # 遍历每个前景文件进行分析
    for i, foreground_file in enumerate(fasta_files):
        print(f"\n--- [{i+1}/{len(fasta_files)}] 开始处理前景文件: {os.path.basename(foreground_file)} ---")
        file_start_time = time.time()
        tissue, file_direction = extract_tissue_info(foreground_file)

        # 构造对应的背景文件名，例如：extracted_new_intron_output_without_up_and_down_Gland_regions.fa
        background_file = f"extracted_new_intron_output_without_up_and_down_{tissue}_regions.fa"
        if not os.path.exists(background_file):
            print(f"警告：背景文件 {background_file} 不存在，跳过 {os.path.basename(foreground_file)} 的分析。")
            continue

        # 解析前景和背景 FASTA 文件
        foreground_seqs = parse_fasta(foreground_file)
        background_seqs = parse_fasta(background_file)

        if foreground_seqs is None or background_seqs is None:
            print(f"警告：解析 {os.path.basename(foreground_file)} 或背景文件失败，跳过。")
            continue

        # 对上游和下游区域分别进行分析
        for region in ['upstream', 'downstream']:
            fg_region_seqs = foreground_seqs.get(region, [])
            bg_region_seqs = background_seqs.get(region, [])

            if not fg_region_seqs:
                print(f"警告：前景文件 {os.path.basename(foreground_file)} 中 {region} 序列为空，跳过该区域分析。")
                continue
            if not bg_region_seqs:
                print(f"警告：背景文件 {background_file} 中 {region} 序列为空，跳过该区域分析。")
                continue

            # 计数 hexamer
            fg_counts, fg_total = count_hexamers_in_list(fg_region_seqs, all_hexamers_set)
            bg_counts, bg_total = count_hexamers_in_list(bg_region_seqs, all_hexamers_set)

            if fg_total == 0 or bg_total == 0:
                print(f"警告：文件 {os.path.basename(foreground_file)} 或背景文件 {background_file} 中 {region} 的有效 hexamer 位置数为零，跳过该区域分析。")
                continue

            # 执行富集分析：前景 vs 背景
            results = run_enrichment_analysis(
                target_label='foreground',
                background_label='background',
                target_counts=fg_counts,
                target_total_pos=fg_total,
                background_counts=bg_counts,
                background_total_pos=bg_total,
                all_hexamers=all_hexamers_list
            )

            if results.empty:
                print(f"警告：{region} 区域没有生成有效结果。")
                continue

            significant = results[results['significant (FDR < 0.05)']].copy()
            if not significant.empty:
                # 标记额外信息：组织、调控方向（up/down，用来表示上调或下调）以及基因区域（upstream/downstream）
                significant.loc[:, 'Tissue'] = tissue
                significant.loc[:, 'regulation'] = file_direction  # 此处 up/down 代表基因调控方向
                significant.loc[:, 'region'] = region       # 此处 upstream/downstream 代表基因组区域
                all_enriched_list.append(significant)
                print(f"  在 {region} 区域发现 {len(significant)} 个显著富集的 Hexamers。")
            else:
                print(f"  在 {region} 区域未发现显著富集的 Hexamers。")

        file_end_time = time.time()
        print(f"--- 处理前景文件 {os.path.basename(foreground_file)} 完成，用时 {file_end_time - file_start_time:.2f} 秒 ---")

    # --- 汇总并保存结果 ---
    print("\n--- 所有文件处理完毕，开始汇总结果 ---")
    output_file = 'all_tissues_enriched_in_foreground_vs_background.csv'

    if all_enriched_list:
        final_df = pd.concat(all_enriched_list, ignore_index=True)
        # 根据组织、调控方向和区域排序
        final_df = final_df.sort_values(by=['Tissue', 'regulation', 'region', 'corrected_p_value (q-value)'])
        cols = ['Tissue', 'regulation', 'region', 'hexamer'] + [col for col in final_df.columns if col not in ['Tissue', 'regulation', 'region', 'hexamer']]
        final_df = final_df[cols]
        try:
            final_df.to_csv(output_file, index=False)
            print(f"富集分析结果已保存到: {output_file}")
        except Exception as e:
            print(f"错误：保存富集结果文件时出错: {e}")
    else:
        print("在所有文件中均未发现显著富集的 Hexamers，不生成结果文件。")

    end_overall_time = time.time()
    print(f"\n--- 整体分析完成，总用时 {end_overall_time - start_overall_time:.2f} 秒 ---")

