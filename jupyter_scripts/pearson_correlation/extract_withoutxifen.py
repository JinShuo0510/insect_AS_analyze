import csv
from collections import defaultdict

def normalize_field(text):
    """字段规范化函数：去前后空格、替换内部空格为下划线"""
    return text.strip().replace(' ', '_')

def process_data():
    run_group_mapping = {}
    with open('../2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 规范化每个字段
            age = normalize_field(row['Age'])
            tissue = normalize_field(row['Tissue.1'])
            
            # 构建分组名称
            group = f"{age}__{tissue}"
            run_group_mapping[row['Run']] = group

    with open('Bombyx_mori_cpm_gene_mat.tsv') as fin, \
         open('GeneExpression_GroupedData_gene.tsv', 'w') as fout:
        
        sample_ids = fin.readline().strip().split()
        
        group_indices = defaultdict(list)
        for idx, sample in enumerate(sample_ids):
            if group := run_group_mapping.get(sample):
                group_indices[group].append(idx)
        
        # 写入规范化的表头
        headers = ["GeneID"] + list(group_indices.keys())
        fout.write("\t".join(headers) + "\n")
        
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 1:
                continue
            
            gene_id = normalize_field(parts[0])  # 基因ID也做规范化
            try:
                values = [float(x) for x in parts[1:]]
            except ValueError as e:
                print(f"数值转换错误（行：{gene_id}）：{str(e)}")
                continue
            
            output = [gene_id]
            for group, indices in group_indices.items():
                try:
                    group_values = [f"{values[i]:.4f}" for i in indices]
                except IndexError:
                    group_values = ["NA"]
                output.append(",".join(group_values))
            
            fout.write("\t".join(output) + "\n")

if __name__ == "__main__":
    process_data()
    print("处理完成，结果已保存到 GeneExpression_GroupedData_withoutxifen.tsv")

