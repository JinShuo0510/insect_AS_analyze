import pandas as pd
import numpy as np

def process_cell(cell):
    """
    对单个单元格处理：
      - 如果 cell 中有逗号，则分割成多个 log2(RPKM) 值，
      - 转换为 float 后用 2**(value) 得到表达水平，
      - 对多个重复值取平均；
      - 如果 cell 中只有一个数字，同样转换。
    返回该组织的表达水平（线性尺度）。
    """
    if pd.isna(cell):
        return 0.0
    # 将 cell 转为字符串，按逗号拆分
    parts = str(cell).split(',')
    try:
        # 将每个 log2 值转换为浮点数，并转为线性表达量 2^(log2(RPKM))
#        linear_vals = [2**float(x) for x in parts]
        linear_vals = [float(x) for x in parts]
    except Exception as e:
        print(f"转换出错，cell = {cell}")
        linear_vals = [0.0]
    return np.mean(linear_vals)

def compute_tau_for_row(row):
    """
    对一行（一个基因）的所有组织计算 Tau 值。
    row 为一个 Series，其索引为各组织的名称（第一列之外）。
    """
    # 对每个组织，计算线性表达水平
    expr_values = row.apply(process_cell).values  # 得到各组织表达值的数组
    n = len(expr_values)
    if n == 0:
        return np.nan
    max_expr = np.max(expr_values)
    # 若所有组织均为 0 则定义 tau = 0（或也可设为 NA）
    if max_expr == 0:
        return 0.0
    # 计算每个组织的归一化表达水平
    norm_expr = expr_values / max_expr
    # 当只有一个组织时避免除以 0
    if n == 1:
        return 0.0
    tau = np.sum(1 - norm_expr) / (n - 1)
    return tau

# ---------------------------
# 主程序部分

# 1. 读入数据文件（假设文件名为 "input.txt"，按制表符分隔）
df = pd.read_csv("GeneExpression_GroupedData_origin.tsv", sep="\t")

# 2. 假定第一列为 GeneID，后续列为不同组织
gene_ids = df.iloc[:, 0]
data = df.iloc[:, 1:]

# 3. 对每个基因（每行）计算 tau
tau_values = data.apply(compute_tau_for_row, axis=1)

# 4. 生成输出 DataFrame（包含基因ID和 tau 值）
df_tau = pd.DataFrame({
    "GeneID": gene_ids,
    "Tau": tau_values
})

# 5. 输出结果到文件（例如 "tau_output.txt"，制表符分隔）
df_tau.to_csv("tau_output.txt", sep="\t", index=False)

print("Tau 值计算完毕，结果保存在 tau_output.txt 文件中。")

