#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 并行 + 向量化 DTU 流程  - 高阈值版（含 Age、Tissue、Pair 三类结果）
# - 输出 DTx_nonDEG_orthogroups_matched_dt_updated.tsv
# - multi_match_count & total_match_count 已计算
# Date: 2025-05-13
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(data.table)
})
options(stringsAsFactors = FALSE)

# 1. 定义输入目录
age_dir    <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh_addpair/Age/Bombyx_mori"
tissue_dir <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh_addpair/Tissue/Bombyx_mori"
pair_dir   <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh_addpair/pair/Bombyx_mori"

# 2. 列出所有 *_DTx_nonDEG.tsv 文件
age_files    <- list.files(age_dir,    pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)
tissue_files <- list.files(tissue_dir, pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)
pair_files   <- list.files(pair_dir,   pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)

# 3. 定义函数：读取单个文件并格式化
read_and_format <- function(files, type_label) {
  dt_list <- lapply(files, function(fp) {
    dt <- fread(fp)
    # 解析文件名获取比较组，去掉后缀
    cmp   <- sub("_DTx_nonDEG\\.tsv$", "", basename(fp))
    parts <- tstrsplit(cmp, "_vs_")
    dt[, `:=`(
      Compare_type  = type_label,
      Compare1      = parts[[1]],
      Compare2      = parts[[2]],
      Transcript_id = tx_id,
      Gene_id       = gene_id
    )]
    # 保留所需列
    dt[, .(Compare_type, Compare1, Compare2,
           Transcript_id, Gene_id,
           logFC, AveExpr, FDR)]
  })
  rbindlist(dt_list)
}

# 4. 读取并合并 Age、Tissue、Pair 数据
dt_age    <- read_and_format(age_files,    "Age")
dt_tissue <- read_and_format(tissue_files, "Tissue")
dt_pair   <- read_and_format(pair_files,   "Pair")
combined_dt <- rbindlist(list(dt_age, dt_tissue, dt_pair))

# 5. 读取 Lepidoptera_Orthogroups.tsv 并计算匹配数
orth_dt <- fread("Lepidoptera_Orthogroups.tsv", check.names = FALSE)
species_cols <- setdiff(names(orth_dt), c("Orthogroup", "Bombyx_mori.protein"))

orth_dt[, `:=`(
  one2one = Reduce(`+`, lapply(.SD, function(col) lengths(strsplit(col, ",\\s*")) == 1)),
  multi_match_count = Reduce(`+`, lapply(.SD, function(col) lengths(strsplit(col, ",\\s*")) > 1))
), .SDcols = species_cols]
orth_dt[, total_match_count := one2one + multi_match_count]

# 6. 拆分 Bombyx_mori.protein，生成映射表
bombyx_map <- orth_dt[
  , .(
      Transcript_id      = unlist(strsplit(`Bombyx_mori.protein`, ",\\s*")),
      one2one, multi_match_count, total_match_count
    ),
  by = Orthogroup
]

# 7. 与 combined_dt 左连接，将 NA 填充为 0
setkey(bombyx_map, Transcript_id)
setkey(combined_dt, Transcript_id)
result_dt <- bombyx_map[combined_dt]
result_dt[is.na(one2one),           one2one           := 0]
result_dt[is.na(multi_match_count), multi_match_count := 0]
result_dt[is.na(total_match_count), total_match_count := 0]

# 8. 重命名列并输出最终结果
setnames(result_dt,
         c("one2one", "multi_match_count"),
         c("1:1_match_count", "multi_match_count"))

fwrite(result_dt,
       file = "DTx_nonDEG_orthogroups_matched_dt_updated_addpair.tsv",
       sep  = "\t")

