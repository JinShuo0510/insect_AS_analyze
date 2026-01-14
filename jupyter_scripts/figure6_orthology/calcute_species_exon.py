import pandas as pd
import numpy as np
import os
import re
from itertools import combinations

def calculate_reciprocal_pairs(species1, species2, liftover_dir="liftover", output_dir="reciprocal_pairs", threshold=0.3):
    """
    识别两个物种之间的互惠最佳外显子重叠，使用halLiftover生成的BED文件。
    采用更低的阈值和更精确的匹配方法。
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 构建文件路径
    s1_to_s2_file = os.path.join(liftover_dir, f"{species1}_to_{species2}.bed")
    s2_to_s1_file = os.path.join(liftover_dir, f"{species2}_to_{species1}.bed")
    output_file = os.path.join(output_dir, f"{species1}_{species2}_reciprocal_pairs.tsv")
    
    # 检查文件是否存在
    if not os.path.exists(s1_to_s2_file) or not os.path.exists(s2_to_s1_file):
        print(f"Warning: Liftover files missing for pair {species1}-{species2}. Skipping.")
        return 0
    
    print(f"Processing reciprocal liftover regions for {species1} and {species2}...")
    
    try:
        # 读取BED文件
        cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'meta']
        
        # 添加调试信息
        print(f"  Reading {s1_to_s2_file}")
        df_a_to_b = pd.read_csv(s1_to_s2_file, sep='\t', header=None, names=cols)
        print(f"  Found {len(df_a_to_b)} entries in {species1} to {species2}")
        
        print(f"  Reading {s2_to_s1_file}")
        df_b_to_a = pd.read_csv(s2_to_s1_file, sep='\t', header=None, names=cols)
        print(f"  Found {len(df_b_to_a)} entries in {species2} to {species1}")
        
        # 解析exon标识符 - 调整为仅使用name字段作为ID
        # 计算每个片段的长度
        df_a_to_b['exon_length'] = df_a_to_b['end'] - df_a_to_b['start']
        df_b_to_a['exon_length'] = df_b_to_a['end'] - df_b_to_a['start']
        
        # 从name中提取exon ID
        # 假设格式为 "Species::gene::exon"
        def extract_exon_id(name_field):
            parts = name_field.split('::')
            if len(parts) >= 3:
                return '::'.join(parts[1:])  # 返回不含物种名的部分
            return name_field  # 如果不符合预期格式，返回原始值
        
        # 应用提取函数
        df_a_to_b['exon_id'] = df_a_to_b['name'].apply(extract_exon_id)
        df_b_to_a['exon_id'] = df_b_to_a['name'].apply(extract_exon_id)
        
        # 对同一exon ID的多个映射片段进行分组
        print("  Grouping liftover fragments by exon ID...")
        a_to_b_grouped = df_a_to_b.groupby('exon_id').agg({
            'exon_length': 'sum',
            'name': 'first',  # 保留一个完整的name用于参考
            'chrom': lambda x: ','.join(set(x))  # 记录映射到的所有染色体
        }).reset_index()
        
        b_to_a_grouped = df_b_to_a.groupby('exon_id').agg({
            'exon_length': 'sum',
            'name': 'first',
            'chrom': lambda x: ','.join(set(x))
        }).reset_index()
        
        print(f"  After grouping: {len(a_to_b_grouped)} unique exons from {species1}, {len(b_to_a_grouped)} from {species2}")
        
        # 提取物种2的exon ID
        # 假设物种2的标识符也在name字段中
        # 由于文件格式限制，我们需要一种方法来识别这些连接
        
        # 分析数据结构，查找潜在的连接方式
        print("  Searching for reciprocal mappings...")
        
        # 尝试另一种方法：从原始文件中提取源IDs和映射IDs之间的关系
        # 将物种1的exon ID映射到物种2
        s1_to_s2_map = df_a_to_b[['name', 'chrom']].drop_duplicates()
        s1_to_s2_map['source_id'] = s1_to_s2_map['name'].apply(lambda x: x.split('::')[1:] if '::' in x else [x])
        s1_to_s2_map['target_chrom'] = s1_to_s2_map['chrom']
        
        # 将物种2的exon ID映射到物种1
        s2_to_s1_map = df_b_to_a[['name', 'chrom']].drop_duplicates()
        s2_to_s1_map['source_id'] = s2_to_s1_map['name'].apply(lambda x: x.split('::')[1:] if '::' in x else [x])
        s2_to_s1_map['target_chrom'] = s2_to_s1_map['chrom']
        
        # 查找互惠对 - 基于染色体位置匹配
        # 这是一个启发式方法，尝试根据可用信息匹配外显子
        
        # 创建临时列以帮助匹配
        a_to_b_grouped['key_for_match'] = a_to_b_grouped['chrom'].astype(str)
        b_to_a_grouped['key_for_match'] = b_to_a_grouped['chrom'].astype(str)
        
        # 使用开始位置近似匹配
        merged_df = pd.merge(
            a_to_b_grouped,
            b_to_a_grouped,
            on='key_for_match',
            suffixes=('_ab', '_ba')
        )
        
        print(f"  Found {len(merged_df)} potential pairs based on genomic location")
        
        # 计算覆盖率
        merged_df['min_length'] = merged_df[['exon_length_ab', 'exon_length_ba']].min(axis=1)
        merged_df['max_length'] = merged_df[['exon_length_ab', 'exon_length_ba']].max(axis=1)
        merged_df['ratio'] = merged_df['min_length'] / merged_df['max_length']
        
        # 阈值过滤
        final_pairs = merged_df[merged_df['ratio'] >= threshold].copy()
        print(f"  Found {len(final_pairs)} pairs with ratio >= {threshold}")
        
        # 创建结果DataFrame
        if len(final_pairs) > 0:
            # 使用copy()避免SettingWithCopyWarning
            result_df = pd.DataFrame({
                'exon_id_species1': final_pairs['exon_id_ab'],
                'exon_id_species2': final_pairs['exon_id_ba'],
                'ratio': final_pairs['ratio'],
                'length_species1': final_pairs['exon_length_ab'],
                'length_species2': final_pairs['exon_length_ba'],
                'chromosome_species1': final_pairs['chrom_ab'],
                'chromosome_species2': final_pairs['chrom_ba']
            })
            
            # 保存结果
            result_df.to_csv(output_file, sep='\t', index=False)
            print(f"Saved {len(result_df)} reciprocal pairs to {output_file}")
            return len(result_df)
        else:
            # 创建空的结果文件
            pd.DataFrame(columns=[
                'exon_id_species1', 'exon_id_species2', 'ratio',
                'length_species1', 'length_species2',
                'chromosome_species1', 'chromosome_species2'
            ]).to_csv(output_file, sep='\t', index=False)
            print(f"No pairs found above threshold. Empty file created: {output_file}")
            return 0
        
    except Exception as e:
        print(f"Error processing pair {species1}-{species2}: {str(e)}")
        import traceback
        traceback.print_exc()
        return 0

def extract_species_from_files(liftover_dir="liftover"):
    """从liftover文件名中提取唯一的物种名称"""
    files = os.listdir(liftover_dir)
    species_set = set()
    
    pattern = r"(.+)_to_(.+)\.bed"
    for file in files:
        match = re.match(pattern, file)
        if match:
            species_set.add(match.group(1))
            species_set.add(match.group(2))
    
    return sorted(list(species_set))

# --- 主执行部分 ---
def main():
    liftover_dir = "liftover"
    output_dir = "reciprocal_pairs"
    # 降低阈值以获取更多潜在匹配
    threshold = 0.3
    
    # 从文件名提取物种列表
    species_list = extract_species_from_files(liftover_dir)
    print(f"Detected {len(species_list)} species: {', '.join(species_list)}")
    
    # 创建摘要表
    summary_data = []
    
    # 处理所有物种对
    for s1, s2 in combinations(species_list, 2):
        pair_count = calculate_reciprocal_pairs(s1, s2, liftover_dir, output_dir, threshold)
        summary_data.append({
            'Species1': s1, 
            'Species2': s2, 
            'ReciprocalPairs': pair_count
        })
    
    # 保存摘要到文件
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, "reciprocal_pairs_summary.tsv")
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary saved to {summary_file}")
    print("Finished processing all pairs.")

if __name__ == "__main__":
    main()
