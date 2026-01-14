# DevAS主分析脚本集合

## 功能概述
识别和分析发育过程中的可变剪接(Developmental Alternative Splicing, DevAS)事件。

## 核心脚本分类

### 1. DevAS识别 (5个)
- `find_devas_all.R`: 识别所有DevAS事件(主脚本)
- `find_devas_all_full.R`: 完整版DevAS识别
- `find_devas.R`: 基础版本
- `find_fold_change_morethan_05.R`: 筛选fold change > 0.5的事件

### 2. 数据整合 (4个)
- `4_all_species_combined_tables.R`: 合并所有物种数据
- `species_combined_tables.R`: 物种数据合并
- `tran_long_format.R`: 转换为长数据格式
- `tran_orthogourp_transcript_to_geneid.R`: 转录本ID转基因ID

### 3. 数据提取 (6个)
- `extract_species_devas.R`: 提取物种特异DevAS
- `extract_species_common.R`: 提取物种共有AS事件
- `extract_species_common_xifen.R`: 细分版本
- `extract_devAS_gene.R`: 提取DevAS相关基因
- `extract_devas_exon.R`: 提取DevAS外显子
- `as_gene_extract.R`: AS基因提取

### 4. 过滤筛选 (4个)
- `filter_devAS_merge.R`: 合并并过滤DevAS
- `filter_as_events_at_least_2_age.R`: 筛选≥2个发育期的AS
- `filter_species.R`: 物种过滤
- `deletena.R`: 删除NA值

### 5. 统计分析 (4个)
- `MDS_analyse.R`: 多维尺度分析
- `pearson_psi_every.R`: 每个样本PSI相关性
- `pearson_psi.R`: PSI相关性分析
- `count.R` / `count_morethan2.R`: 统计计数

### 6. 可视化 (6个)
- `zhuxingtu.R`: 柱形图
- `plot_quxiantu.R`: 曲线图
- `sandiantu_3N.R`: 3N散点图
- `sandiantu_IDR.R`: IDR散点图
- `sandiantu_phastcons.R`: PhastCons保守性散点图
- `draw_plot-1.R` / `draw_plot-2.R`: 通用绘图

### 7. 其他工具 (7个)
- `analyse_test.R`: 测试分析
- `delete_results_NULL.R`: 删除空结果

## 典型分析流程

```bash
# 步骤1: 识别DevAS事件
Rscript find_devas_all.R

# 步骤2: 过滤和筛选
Rscript filter_devAS_merge.R
Rscript filter_as_events_at_least_2_age.R

# 步骤3: 多物种整合
Rscript 4_all_species_combined_tables.R

# 步骤4: 统计分析
Rscript MDS_analyse.R
Rscript pearson_psi_every.R

# 步骤5: 可视化
Rscript zhuxingtu.R
Rscript plot_quxiantu.R
```

## 关键概念

### DevAS定义
在发育过程中PSI值发生显著变化的可变剪接事件。

### 判定标准
- PSI变化阈值: ≥0.2
- 至少在2个发育期表达
- 通过统计检验(FDR < 0.05)

## 输入数据
- PSI矩阵文件
- 样本元数据
- 同源组信息(可选)

## 输出结果
- DevAS事件列表
- PSI变化矩阵
- 统计检验结果
- 可视化图表
