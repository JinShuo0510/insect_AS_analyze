import csv
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
    
    # 存储不匹配的序列
    unmatched_sequences = []
    
    # 创建一个集合来存储CSV文件中的区域信息，便于快速查找
    csv_regions = set()
    for row in exon_data:
        if row['tissue'] == tissue:
            chr_id = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            csv_regions.add((chr_id, start, end))
    
    try:
        # 从fasta文件中获取序列
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))
        
        # 检查每个序列是否在CSV中
        for seq in sequences:
            # 提取区域信息
            region_match = re.search(r'::(.+):(\d+)-(\d+)', seq.id)
            if region_match:
                seq_chr = region_match.group(1)
                seq_start = int(region_match.group(2))
                seq_end = int(region_match.group(3))
                
                # 检查此序列是否在CSV中
                found_in_csv = False
                for csv_chr, csv_start, csv_end in csv_regions:
                    # 检查区域是否重叠
                    if (seq_chr == csv_chr and 
                        max(seq_start, csv_start) <= min(seq_end, csv_end)):
                        found_in_csv = True
                        break
                
                # 如果不在CSV中，添加到未匹配序列列表
                if not found_in_csv:
                    unmatched_sequences.append(seq)
        
        # 写入不在CSV中的序列到新文件
        if unmatched_sequences:
            output_file = f'{tissue}_sequences_not_in_csv.fa'
            SeqIO.write(unmatched_sequences, output_file, 'fasta')
            print(f"已为{tissue}组织创建不在CSV中的序列文件: {output_file}, 共{len(unmatched_sequences)}条序列")
        else:
            print(f"{tissue}组织中所有序列都在CSV文件中，没有生成新文件")
            
    except FileNotFoundError:
        print(f"警告：找不到文件 {fasta_file}，跳过{tissue}组织的处理")
