# DEG/DET分析脚本

## 功能概述
使用edgeR进行差异表达基因(DEG)和差异表达转录本(DET)的鉴定分析。

## 主要脚本

### 核心分析脚本
1. **analyse_DET_fifo_high_threshold_addpair_groupmean_cpm.R**
   - 高阈值DET分析主脚本
   - 使用CPM标准化
   - 并行处理提高效率

2. **edgeR_DTU_devDynamic_v2.6.R**
   - edgeR差异转录本使用(DTU)分析最新版本
   - 支持发育动态分析

### 辅助分析
- `isoform_var.R`: 单物种转录本变异分析
- `isoform_var_all.R`: 多物种转录本变异分析
- `count_det.py`: 统计Age维度的DET数量
- `count_det_tissue.py`: 统计Tissue维度的DET数量

## 使用方法
```bash
# 运行主分析
Rscript analyse_DET_fifo_high_threshold_addpair_groupmean_cpm.R

# 统计结果
python count_det.py
python count_det_tissue.py
```

## 输出结果
- DEG列表: 差异表达基因
- DET列表: 差异表达转录本  
- DTx_nonDEG: 非差异基因的差异转录本
- CPM矩阵: 标准化表达矩阵
