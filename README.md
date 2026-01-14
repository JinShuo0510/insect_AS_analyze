
# 昆虫可变剪接分析脚本集合

## 项目概述

本项目整理了昆虫可变剪接(Alternative Splicing)分析的所有脚本,包括差异表达分析、DevAS识别、序列分析和数据可视化。

**脚本总数**: 179个  

---

## Docker支持

本项目的DEG/DET分析流程提供Docker镜像,开箱即用:

**Docker Hub**: `kim0418/nondeg-det-pipeline:latest`

### 快速使用
```bash
# 拉取镜像
docker pull kim0418/nondeg-det-pipeline:latest

# 运行DEG/DET分析
docker run --rm -v /your/data:/data \
    kim0418/nondeg-det-pipeline:latest \
    Rscript /pipeline/analyse_DET_fifo_high_threshold_addpair_groupmean_cpm.R \
    -s Bombyx_mori --gtf-dir /data/Genome \
    --gene-dir /data/counts/gene --tx-dir /data/counts/tx \
    -m /data/samples.csv -o /data/results
```

详细使用文档请参考: [nonDEG_DET_pipeline](https://github.com/JinShuo0510/nonDEG_DET_pipeline)

---

## 目录结构

```
insect_as_analyze/
├── jupyter_scripts/              # Python/Jupyter分析脚本
│   ├── DEG_DET_analysis/        # 差异表达基因/转录本分析 (16个)
│   ├── exon_extraction/         # 外显子提取和序列分析 (10个)
│   ├── figure3_analysis/        # Figure3相关分析 (6个)
│   ├── figure4_motif/           # Motif序列分析 (7个)
│   ├── figure6_orthology/       # 同源性分析 (9个)
│   ├── figure7_matching/        # 匹配分析 (7个)
│   └── pearson_correlation/     # Pearson相关性分析 (11个)
│
└── rstudio_scripts/              # R绘图脚本
    ├── DevAS_main/              # DevAS主要分析 (36个)
    ├── figure_plots/            # 各图表绘制
    │   ├── figure2/             (3个)
    │   ├── figure3/             (8个)
    │   ├── figure4/             (18个)
    │   ├── figure5/             (11个)
    │   ├── figure6/             (19个)
    │   ├── figure7/             (1个)
    │   └── figure8/             (4个)
    ├── nonDEG_DET_analysis/     # 非差异基因DET分析 (6个)
    └── orthogroup_analysis/     # 同源组分析 (7个)
```

---

## 主要功能

### 1. 差异表达分析
- **位置**: `jupyter_scripts/DEG_DET_analysis/`
- **功能**: 使用edgeR进行DEG和DET识别
- **主要脚本**: `analyse_DET_fifo_high_threshold_addpair_groupmean_cpm.R`

### 2. DevAS识别
- **位置**: `rstudio_scripts/DevAS_main/`
- **功能**: 识别发育过程中的可变剪接事件
- **主要脚本**: `find_devas_all.R`

### 3. Motif分析
- **位置**: `jupyter_scripts/figure4_motif/`
- **功能**: 外显子上下游序列motif富集分析
- **主要脚本**: `motif_6_analayse_final.py`

### 4. 同源性分析
- **位置**: `jupyter_scripts/figure6_orthology/`
- **功能**: 跨物种外显子同源性分析
- **主要脚本**: `analyze_reciprocal_overlaps_youhua.py`

### 5. 数据可视化
- **位置**: `rstudio_scripts/figure_plots/`
- **功能**: 论文各图表绘制
- **包含**: Figure 2-8的绘图脚本

---

## 快速开始

### 场景1: DEG/DET分析
```bash
cd jupyter_scripts/DEG_DET_analysis/
Rscript analyse_DET_fifo_high_threshold_addpair_groupmean_cpm.R
```

### 场景2: DevAS识别
```bash
cd rstudio_scripts/DevAS_main/
Rscript find_devas_all.R
```

### 场景3: 绘图
```bash
cd rstudio_scripts/figure_plots/figure3/
Rscript plot_zhuxingtu.R
```

---

## 关键参数

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| CPM阈值 | ≥1 | 表达量过滤 |
| PSI变化 | ≥0.2 | AS显著性 |
| FDR | <0.05 | 统计检验 |

---

## 环境要求

### R环境
```R
# 基础包
install.packages(c("dplyr", "tidyr", "ggplot2"))

# Bioconductor
BiocManager::install(c("edgeR", "limma"))
```

### Python环境
```bash
pip install pandas numpy scipy statsmodels
```

---

## 输入数据

1. **FeatureCounts结果**: 基因和转录本level的counts文件
2. **样本元数据**: 包含Age、Tissue、分组信息的CSV文件
3. **GTF注释文件**: 基因组注释
4. **PSI矩阵**: 可变剪接PSI值(可选)

---

## 输出结果

- DEG/DET列表和统计表
- DevAS事件表
- 各类图表(PDF/PNG)
- 统计汇总文件

---

## 注意事项

1. 确保数据路径正确
2. 检查样本元数据格式
3. 大数据集建议使用并行处理
4. 查看各子目录README了解详细用法

---


**版本**: v1.0  
**日期**: 2026-01-14
