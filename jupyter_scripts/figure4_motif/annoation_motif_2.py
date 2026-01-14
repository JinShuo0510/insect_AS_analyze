import pandas as pd
import math
import re
from collections import defaultdict

def parse_pwm_file(pwm_filepath):
    """
    解析包含多个 PWM 条目的文件。
    每个条目以 'RBP\t' 开头，包含元数据和频率矩阵。
    """
    pwms = []
    current_pwm = None
    metadata = {}
    matrix_lines = []

    with open(pwm_filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: # 空行通常分隔条目，处理上一个条目
                if current_pwm and matrix_lines:
                    try:
                        # 将矩阵行解析为浮点数列表
                        matrix = [[float(v) for v in row.split()]
                                  for row in matrix_lines]
                        # 验证矩阵格式 (4列)
                        if not all(len(row) == 4 for row in matrix):
                             print(f"Warning: Skipping motif {metadata.get('Motif', 'Unknown')} due to incorrect matrix column count.")
                        else:
                            current_pwm['matrix'] = matrix
                            current_pwm['length'] = len(matrix)
                            pwms.append(current_pwm)
                    except ValueError:
                         print(f"Warning: Skipping motif {metadata.get('Motif', 'Unknown')} due to non-numeric matrix values.")

                current_pwm = None
                metadata = {}
                matrix_lines = []
                continue

            if line.startswith('RBP\t'):
                 # 如果遇到新的RBP行，并且上一个条目有数据，先处理上一个条目
                if current_pwm and matrix_lines:
                    try:
                        matrix = [[float(v) for v in row.split()] for row in matrix_lines]
                        if not all(len(row) == 4 for row in matrix):
                            print(f"Warning: Skipping motif {metadata.get('Motif', 'Unknown')} due to incorrect matrix column count.")
                        else:
                            current_pwm['matrix'] = matrix
                            current_pwm['length'] = len(matrix)
                            pwms.append(current_pwm)
                    except ValueError:
                        print(f"Warning: Skipping motif {metadata.get('Motif', 'Unknown')} due to non-numeric matrix values.")

                # 开始新的条目
                current_pwm = {}
                metadata = {}
                matrix_lines = []
                key_value = line.split('\t', 1)
                if len(key_value) == 2:
                     metadata['RBP'] = key_value[1]
                current_pwm = {'metadata': metadata}

            elif current_pwm and line.startswith(('RBP Name\t', 'Gene\t', 'Motif\t', 'Family\t', 'Species\t')):
                key_value = line.split('\t', 1)
                if len(key_value) == 2:
                    metadata[key_value[0].split(' ')[0]] = key_value[1] # 使用空格前的部分作为key

            elif current_pwm and line.startswith('Pos\t'):
                # 矩阵头行，忽略
                pass
            elif current_pwm and re.match(r'^\d+\s+', line):
                 # 矩阵数据行，提取A C G U值部分
                 parts = line.split(maxsplit=1) # 按第一个空白分割
                 if len(parts) == 2:
                     matrix_lines.append(parts[1])

        # 处理文件末尾最后一个条目
        if current_pwm and matrix_lines:
            try:
                matrix = [[float(v) for v in row.split()] for row in matrix_lines]
                if not all(len(row) == 4 for row in matrix):
                    print(f"Warning: Skipping motif {metadata.get('Motif', 'Unknown')} due to incorrect matrix column count.")
                else:
                    current_pwm['matrix'] = matrix
                    current_pwm['length'] = len(matrix)
                    pwms.append(current_pwm)
            except ValueError:
                 print(f"Warning: Skipping motif {metadata.get('Motif', 'Unknown')} due to non-numeric matrix values.")


    print(f"Parsed {len(pwms)} PWMs successfully.")
    return pwms

def calculate_max_pwm_prob(pwm_matrix):
    """计算给定PWM可能产生的最大概率得分"""
    max_prob = 1.0
    for position in pwm_matrix:
        # position 是 [pA, pC, pG, pU]
        if not position: # 跳过空行（理论上不应出现）
            continue
        # 找到当前位置的最大概率
        max_pos_prob = max(position) if position else 0
        if max_pos_prob == 0: # 如果某位置最大概率为0，则整个序列概率为0
             # 实践中可能需要处理伪计数，但这里严格按定义
            return 0.0
        max_prob *= max_pos_prob
    return max_prob

def calculate_hexamer_prob_from_pwm(hexamer, pwm, base_map={'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}):
    """
    计算PWM生成给定hexamer的最大概率（考虑所有可能的6nt对齐窗口）。
    """
    hexamer = hexamer.upper().replace('T', 'U') # 标准化为大写U
    if len(hexamer) != 6:
        return 0.0

    pwm_matrix = pwm['matrix']
    pwm_len = pwm['length']

    if pwm_len < 6:
        return 0.0 # PWM太短，无法生成hexamer

    max_prob_for_hexamer = 0.0

    # 遍历PWM中所有可能的6nt窗口起始位置
    for start_pos in range(pwm_len - 5):
        current_prob = 1.0
        possible = True
        for i in range(6):
            base = hexamer[i]
            base_index = base_map.get(base)
            if base_index is None: # 非法碱基
                current_prob = 0.0
                possible = False
                break

            pwm_row = pwm_matrix[start_pos + i]
            if not pwm_row or len(pwm_row) != 4: # 检查行有效性
                 current_prob = 0.0
                 possible = False
                 break

            prob_base = pwm_row[base_index]
            if prob_base == 0: # 如果某碱基概率为0，则此窗口概率为0
                current_prob = 0.0
                possible = False
                break
            current_prob *= prob_base

        if possible:
            max_prob_for_hexamer = max(max_prob_for_hexamer, current_prob)

    return max_prob_for_hexamer


# --- 主程序 ---
pwm_filepath = 'PWM.txt' # <-- 修改为你的PWM文件名
results_filepath = 'all_tissues_enriched_in_foreground_vs_background.csv' # <-- 修改为你的结果文件名
output_filepath = 'all_tissues_enriched_in_foreground_vs_background.csv_annotated_results.csv' # <-- 输出文件名

# 1. 解析PWM文件
print(f"Parsing PWM file: {pwm_filepath}...")
pwms = parse_pwm_file(pwm_filepath)
if not pwms:
    print("Error: No PWMs were parsed. Please check the PWM file format.")
    exit()

# 2. 预计算每个PWM的最大概率和阈值
print("Pre-calculating maximum probabilities and thresholds for PWMs...")
pwm_thresholds = {}
for pwm in pwms:
    motif_id = pwm['metadata'].get('Motif', 'UnknownMotif')
    if motif_id == 'UnknownMotif' or 'matrix' not in pwm:
        print(f"Warning: Skipping PWM with missing Motif ID or matrix.")
        continue
    max_p = calculate_max_pwm_prob(pwm['matrix'])
    if max_p > 0: # 只有当最大概率大于0时，阈值才有意义
        threshold = 0.5 * max_p
        pwm_thresholds[motif_id] = threshold
        # 存储 RBP Name 以便后续查找
        pwm['threshold'] = threshold
        pwm['max_p'] = max_p
    else:
         pwm['threshold'] = 0 # 如果最大概率为0，则没有序列能匹配
         pwm['max_p'] = 0

print(f"Calculated thresholds for {len(pwm_thresholds)} PWMs.")

# 3. 读取结果文件
print(f"Reading results file: {results_filepath}...")
try:
    results_df = pd.read_csv(results_filepath)
except FileNotFoundError:
    print(f"Error: Results file not found at {results_filepath}")
    exit()

# 4. 准备注释列
matched_rbp_names = []
matched_motif_ids = []

# 5. 遍历结果文件中的每个hexamer并进行注释
print("Annotating hexamers...")
total_hexamers = len(results_df)
for index, row in results_df.iterrows():
    hexamer = row['hexamer']
    current_matches_rbp = []
    current_matches_motif = []

    if pd.isna(hexamer): # 跳过空的hexamer
        matched_rbp_names.append("")
        matched_motif_ids.append("")
        continue

    for pwm in pwms:
         # 确保PWM有阈值和矩阵
         if 'threshold' not in pwm or 'matrix' not in pwm:
             continue

         threshold = pwm['threshold']
         if threshold == 0: # 如果阈值为0，直接跳过
            continue

         prob_hex = calculate_hexamer_prob_from_pwm(str(hexamer), pwm) # 确保hexamer是字符串

         if prob_hex >= threshold:
             rbp_name = pwm['metadata'].get('RBP', 'UnknownRBP') # 使用 'RBP' 作为键
             motif_id = pwm['metadata'].get('Motif', 'UnknownMotif')
             if rbp_name not in current_matches_rbp: # 避免重复添加
                 current_matches_rbp.append(rbp_name)
             if motif_id not in current_matches_motif:
                 current_matches_motif.append(motif_id)

    # 将匹配结果添加到列表
    matched_rbp_names.append(";".join(current_matches_rbp) if current_matches_rbp else "No_Match")
    matched_motif_ids.append(";".join(current_matches_motif) if current_matches_motif else "No_Match")

    # 打印进度
    if (index + 1) % 100 == 0 or (index + 1) == total_hexamers:
        print(f"Processed {index + 1}/{total_hexamers} hexamers...")

# 6. 将注释结果添加到DataFrame
results_df['Matched_RBP_Name'] = matched_rbp_names
results_df['Matched_Motif_ID'] = matched_motif_ids

# 7. 保存带有注释的结果文件
print(f"Saving annotated results to: {output_filepath}...")
results_df.to_csv(output_filepath, index=False)

print("Annotation complete.")
