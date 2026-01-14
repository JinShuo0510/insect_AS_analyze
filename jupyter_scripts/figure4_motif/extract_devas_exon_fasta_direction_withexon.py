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
    fasta_file = f'extracted_new_intron_contain_{tissue}_regions.fa'
    
    # 根据direction分类存储序列
    up_sequences = []
    down_sequences = []
    
    # 筛选该组织的外显子数据
    tissue_exons = [row for row in exon_data if row['tissue'] == tissue]
    
    try:
        # 从fasta文件中获取序列
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))
        
        # 创建一个字典以存储exon序列及其关联的upstream/downstream序列
        exon_dict = {}  # 键为gene_chr_start_end，值为包含exon序列信息的字典
        all_sequences = {}  # 存储所有序列，键为序列ID
        
        # 先将所有序列放入字典便于查找
        for seq in sequences:
            all_sequences[seq.id] = seq
            
            # 检查是否为exon序列
            if "exon" in seq.id:
                region_match = re.search(r'::(.+):(\d+)-(\d+)', seq.id)
                if region_match:
                    seq_chr = region_match.group(1)
                    seq_start = int(region_match.group(2))
                    seq_end = int(region_match.group(3))
                    
                    # 从ID中提取基因名
                    gene_match = re.search(r'^(.+?)_exon', seq.id)
                    if gene_match:
                        gene = gene_match.group(1)
                        
                        # 创建唯一键
                        key = f"{gene}_{seq_chr}_{seq_start}_{seq_end}"
                        
                        exon_dict[key] = {
                            'exon_id': seq.id,
                            'exon_seq': seq,
                            'gene': gene,
                            'chr': seq_chr,
                            'start': seq_start,
                            'end': seq_end,
                            'upstream_id': None,
                            'downstream_id': None
                        }
        
        # 匹配upstream和downstream序列
        for seq_id in all_sequences:
            if "upstream" in seq_id:
                # 尝试找出关联的exon
                gene_match = re.search(r'^(.+?)_upstream', seq_id)
                region_match = re.search(r'::(.+):(\d+)-(\d+)', seq_id)
                
                if gene_match and region_match:
                    gene = gene_match.group(1)
                    seq_chr = region_match.group(1)
                    # 对于upstream和downstream，我们不需要使用它们的坐标进行匹配
                    # 只需通过基因名和染色体匹配相应的exon
                    
                    # 寻找匹配的exon
                    for key, data in exon_dict.items():
                        if data['gene'] == gene and data['chr'] == seq_chr:
                            exon_dict[key]['upstream_id'] = seq_id
                            break
                            
            elif "downstream" in seq_id:
                gene_match = re.search(r'^(.+?)_downstream', seq_id)
                region_match = re.search(r'::(.+):(\d+)-(\d+)', seq_id)
                
                if gene_match and region_match:
                    gene = gene_match.group(1)
                    seq_chr = region_match.group(1)
                    
                    # 寻找匹配的exon
                    for key, data in exon_dict.items():
                        if data['gene'] == gene and data['chr'] == seq_chr:
                            exon_dict[key]['downstream_id'] = seq_id
                            break
        
        # 与CSV中的数据匹配并生成输出序列
        for exon in tissue_exons:
            csv_gene = exon['gene']
            csv_chr = exon['chr']
            csv_start = int(exon['start'])
            csv_end = int(exon['end'])
            direction = exon['direction']
            
            # 寻找匹配的exon记录
            for key, data in exon_dict.items():
                if (data['gene'] == csv_gene and 
                    data['chr'] == csv_chr and 
                    ((data['start'] <= csv_start <= data['end']) or 
                     (data['start'] <= csv_end <= data['end']) or
                     (csv_start <= data['start'] and csv_end >= data['end']))):
                    
                    # 创建一个新的序列集合
                    matched_sequences = []
                    
                    # 添加exon序列
                    exon_seq = data['exon_seq']
                    new_exon_id = f"{csv_gene}_exon::{csv_chr}:{csv_start}-{csv_end}"
                    matched_sequences.append(SeqRecord(
                        exon_seq.seq,
                        id=new_exon_id,
                        description=""
                    ))
                    
                    # 添加upstream序列（如果存在）
                    if data['upstream_id']:
                        upstream_seq = all_sequences[data['upstream_id']]
                        new_upstream_id = f"{csv_gene}_upstream::{csv_chr}:{csv_start}-{csv_end}"
                        matched_sequences.append(SeqRecord(
                            upstream_seq.seq,
                            id=new_upstream_id,
                            description=""
                        ))
                    
                    # 添加downstream序列（如果存在）
                    if data['downstream_id']:
                        downstream_seq = all_sequences[data['downstream_id']]
                        new_downstream_id = f"{csv_gene}_downstream::{csv_chr}:{csv_start}-{csv_end}"
                        matched_sequences.append(SeqRecord(
                            downstream_seq.seq,
                            id=new_downstream_id,
                            description=""
                        ))
                    
                    # 根据direction添加到相应的列表
                    if direction == "up":
                        up_sequences.extend(matched_sequences)
                    elif direction == "down":
                        down_sequences.extend(matched_sequences)
                    
                    break  # 找到匹配后退出内循环
        
        # 写入过滤后的序列到输出文件
        if up_sequences:
            up_output_file = f'{tissue}_up_filtered_exons.fa'
            with open(up_output_file, 'w') as output:
                SeqIO.write(up_sequences, output, 'fasta')
            print(f"已为{tissue}组织创建上调(up)外显子文件: {up_output_file}")
        else:
            print(f"在{tissue}组织中未找到上调(up)外显子序列")
            
        if down_sequences:
            down_output_file = f'{tissue}_down_filtered_exons.fa'
            with open(down_output_file, 'w') as output:
                SeqIO.write(down_sequences, output, 'fasta')
            print(f"已为{tissue}组织创建下调(down)外显子文件: {down_output_file}")
        else:
            print(f"在{tissue}组织中未找到下调(down)外显子序列")
            
    except FileNotFoundError:
        print(f"警告：找不到文件 {fasta_file}，跳过{tissue}组织的处理")

