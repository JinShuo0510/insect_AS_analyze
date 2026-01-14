import csv
import re
from Bio import SeqIO

# 读取CSV文件
exon_data = []
with open('Bombyx_mori_merge_exon_filter.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        exon_data.append(row)

# 组织不同组织的数据
tissues = set(row['tissue'] for row in exon_data)

# 为每个组织处理fasta文件
for tissue in tissues:
    # 读取对应组织的fasta文件
    fasta_file = f'extracted_new_intron_{tissue}_regions.fa'
    
    # 创建输出文件
    output_file = f'{tissue}_filtered_exons.fa'
    
    # 筛选该组织的外显子数据
    tissue_exons = [row for row in exon_data if row['tissue'] == tissue]
    
    # 从fasta文件中获取序列
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    
    # 过滤出匹配的序列
    filtered_sequences = []
    
    for exon in tissue_exons:
        chr_id = exon['chr']
        start = int(exon['start'])
        end = int(exon['end'])
        gene = exon['gene']
        direction = exon['direction']
        
        # 构建可能的序列标识符模式
        pattern = f"{gene}_(upstream|downstream)::{chr_id}:"
        
        for seq in sequences:
            if re.search(pattern, seq.id):
                # 提取序列标识符中的区域信息
                region_match = re.search(r'::(.+):(\d+)-(\d+)', seq.id)
                if region_match:
                    seq_chr = region_match.group(1)
                    seq_start = int(region_match.group(2))
                    seq_end = int(region_match.group(3))
                    
                    # 检查是否与外显子区域重叠
                    if seq_chr == chr_id and (
                        (seq_start <= start <= seq_end) or 
                        (seq_start <= end <= seq_end) or
                        (start <= seq_start and end >= seq_end)
                    ):
                        filtered_sequences.append(seq)
    
    # 写入过滤后的序列到输出文件
    if filtered_sequences:
        with open(output_file, 'w') as output:
            SeqIO.write(filtered_sequences, output, 'fasta')
        print(f"已为{tissue}组织创建过滤后的外显子文件: {output_file}")
    else:
        print(f"在{tissue}组织中未找到匹配的外显子序列")
