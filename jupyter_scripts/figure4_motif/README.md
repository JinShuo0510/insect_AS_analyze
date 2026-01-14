# Figure4 Motif序列分析

## 功能概述
分析可变剪接外显子上下游200bp区域的motif富集,识别剪接调控元件。

## 主要脚本

### 序列提取
1. **extract_devas_exon_fasta.py**
   - 提取DevAS外显子序列

2. **extract_devas_exon_fasta_direction_withexon.py**
   - 提取带方向的外显子及上下游序列
   - 推荐使用此脚本

3. **extract_devas_exon_fasta_direction_withexon_unfilter.py**
   - 未过滤版本(包含所有外显子)

### Motif分析
4. **motif_6_analayse_final.py**
   - **主要分析脚本**
   - 执行motif富集分析
   - FDR校正
   - 识别显著富集的motif

### 注释和汇总
5. **annoation_motif.py** / **annoation_motif_2.py**
   - 对识别的motif进行功能注释

6. **aggregate_motifs.py**
   - 聚合多个组织/时期的motif结果

## 分析流程
```bash
# 1. 提取序列
python extract_devas_exon_fasta_direction_withexon.py

# 2. Motif分析
python motif_6_analayse_final.py

# 3. 注释
python annoation_motif_2.py

# 4. 汇总
python aggregate_motifs.py
```

## 关键参数
- 上下游区域: 200bp
- FDR阈值: 0.05
- 最小motif长度: 6bp

## 输出结果
- 显著富集的motif列表
- Motif位置信息
- 富集分数和p值
