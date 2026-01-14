import os # 用于文件系统操作，如列出目录内容
import re # 用于正则表达式匹配
import pandas as pd # 用于数据处理和CSV输出
from collections import defaultdict # 虽然最终没直接用，但有时处理这类问题有用

# --- 可配置参数 ---
target_motif = 'TCTCTC' # 你要搜索的目标基序
output_csv_file = 'aggregated_motif_TCTCTC_positions.csv' # 最终输出的汇总文件名
fasta_directory = '.' # 当前目录。如果文件在别处，修改此路径
# ---

def parse_fasta(filename):
    """
    解析 FASTA 文件，产生标题和序列。
    处理多行序列。
    """
    header = None # 当前序列的标题
    sequence_parts = [] # 存储当前序列的各部分
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip() # 去除行首尾空白
                if not line: # 跳过空行
                    continue
                if line.startswith('>'): # 如果是标题行
                    # 如果存在上一个序列，先产生它
                    if header:
                        yield header, "".join(sequence_parts)
                    # 开始处理新序列
                    header = line
                    sequence_parts = []
                elif header: # 只有在找到标题后才添加序列行
                    sequence_parts.append(line.upper()) # 将序列存储为大写

            # 产生文件中的最后一个序列
            if header:
                yield header, "".join(sequence_parts)
    except FileNotFoundError:
        print(f"错误：文件未在 {filename} 找到")
        return # 文件未找到则停止执行
    except Exception as e:
        print(f"读取文件 {filename} 时发生错误: {e}")
        return

def process_single_fasta(fasta_filename, motif, tissue, direction, results_list):
    """
    处理单个FASTA文件，查找基序并将结果添加到全局列表。
    """
    motif_upper = motif.upper()
    motif_len = len(motif_upper)
    local_motif_found_count = 0

    # 用于解析标题的正则表达式: >{Gene_ID}_[region]::{chr}:{start}-{end}
    header_pattern = re.compile(r"^>(.+)_(upstream|downstream|exon)::([^:]+):(\d+)-(\d+)$")

    print(f"  处理文件: {os.path.basename(fasta_filename)} (Tissue: {tissue}, Direction: {direction})")

    for header, sequence in parse_fasta(fasta_filename):
        match = header_pattern.match(header)
        if not match:
            print(f"    警告：无法解析标题: {header} (在文件 {os.path.basename(fasta_filename)} 中)")
            continue

        gene_id, region, chrom, start_str, end_str = match.groups()
        try:
            start = int(start_str)
            end = int(end_str)
        except ValueError:
            print(f"    警告：无法解析坐标 '{start_str}-{end_str}' 为整数 (标题: {header})")
            continue

        seq_len = len(sequence)

        start_index = 0
        while True:
            pos = sequence.find(motif_upper, start_index)
            if pos == -1:
                break

            local_motif_found_count += 1
            relative_pos = None

            if region == "upstream":
                relative_pos = pos - seq_len
            elif region == "downstream":
                relative_pos = pos + 1
            elif region == "exon":
                relative_pos = pos + 1

            results_list.append({
                "Tissue": tissue, # 新增列
                "Direction": direction, # 新增列
                "Gene_ID": gene_id,
                "Chr": chrom,
                "Start (Exon)": start,
                "End (Exon)": end,
                "Region": region,
                "Motif_Start_Pos_In_Region": pos,
                "RelativePositionSpliceSite": relative_pos,
                "Count": 1
            })
            start_index = pos + 1

    print(f"    文件 {os.path.basename(fasta_filename)} 处理完成，找到 {local_motif_found_count} 个基序。")

# --- 主程序执行 ---

all_results = [] # 用于存储所有文件结果的列表

# 定义文件名匹配模式
# 模式1: {Tissue}_[up|down]_filtered_exons.fa
pattern1 = re.compile(r"^(\w+?)_(up|down)_filtered_exons\.fa$")
# 模式2: {Tissue}_sequences_not_in_csv.fa
pattern2 = re.compile(r"^(\w+?)_sequences_not_in_csv\.fa$")

print(f"开始在目录 '{fasta_directory}' 中搜索 FASTA 文件并分析基序 '{target_motif}'...")

# 遍历指定目录下的所有文件
processed_files_count = 0
try:
    for filename in os.listdir(fasta_directory):
        # 构建完整的文件路径
        filepath = os.path.join(fasta_directory, filename)

        # 确保是文件而不是目录
        if not os.path.isfile(filepath):
            continue

        tissue = None
        direction = None
        matched = False

        # 尝试匹配模式1
        match1 = pattern1.match(filename)
        if match1:
            tissue, direction = match1.groups()
            matched = True

        # 如果模式1不匹配，尝试匹配模式2
        if not matched:
            match2 = pattern2.match(filename)
            if match2:
                tissue = match2.group(1)
                direction = "NA" # 对于这个模式，方向是NA
                matched = True

        # 如果文件名匹配了任一模式
        if matched:
            process_single_fasta(filepath, target_motif, tissue, direction, all_results)
            processed_files_count += 1

except FileNotFoundError:
    print(f"错误：指定的目录 '{fasta_directory}' 不存在。")
    exit()
except Exception as e:
    print(f"遍历文件或处理文件时发生错误: {e}")
    exit()

print(f"\n总共处理了 {processed_files_count} 个匹配的 FASTA 文件。")

# 将所有结果转换为DataFrame
if all_results:
    final_df = pd.DataFrame(all_results)

    # 定义期望的列顺序
    output_columns = [
        "Tissue", "Direction", "Gene_ID", "Chr", "Start (Exon)", "End (Exon)",
        "Region", "RelativePositionSpliceSite", "Motif_Start_Pos_In_Region", "Count"
    ]
    # 确保所有列都存在，如果某列在某些文件中没有结果，则填充NaN，这里因为我们每次都添加完整字典，应该不会缺列
    final_df = final_df[output_columns]

    # 保存汇总的CSV文件
    try:
        final_df.to_csv(output_csv_file, index=False)
        print(f"\n所有结果已汇总并保存到: {output_csv_file}")
        print(f"总共找到 {len(final_df)} 个基序实例。")
    except Exception as e:
        print(f"保存汇总 CSV 文件时出错: {e}")
else:
    print("\n在所有处理的文件中未找到目标基序或没有处理任何文件。")
