# 完整 R 脚本：data.table 加速遍历、合并并计算 multi_match_count 和 total_match_count

# 1. 加载包与全局选项
library(data.table)
options(stringsAsFactors = FALSE)

# 2. 定义输入目录
age_dir    <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/Age/Bombyx_mori"
tissue_dir <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/Tissue/Bombyx_mori"

# 3. 列出所有 *_DTx_nonDEG.tsv 文件
age_files    <- list.files(age_dir,    pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)
tissue_files <- list.files(tissue_dir, pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)

# 4. 定义函数：读取单个文件并格式化
read_and_format <- function(files, type_label) {
  dt_list <- lapply(files, function(fp) {
    dt <- fread(fp)
    # 解析文件名获取比较组
    cmp   <- sub("_DTx_nonDEG\\.tsv$", "", basename(fp))
    parts <- tstrsplit(cmp, "_vs_")
    dt[, `:=`(
      Compare_type     = type_label,
      Compare1         = parts[[1]],
      Compare2         = parts[[2]],
      Transcript_id    = tx_id,
      Gene_id          = gene_id
    )]
    # 只保留关键信息列
    dt[, .(Compare_type, Compare1, Compare2,
           Transcript_id, Gene_id,
           logFC, AveExpr, FDR)]
  })
  rbindlist(dt_list)
}

# 5. 读取并合并所有 Age 与 Tissue 数据
dt_age    <- read_and_format(age_files,    "Age")
dt_tissue <- read_and_format(tissue_files, "Tissue")
combined_dt <- rbindlist(list(dt_age, dt_tissue))

# 6. 读取 Lepidoptera_Orthogroups.tsv
orth_dt <- fread("Lepidoptera_Orthogroups.tsv", check.names = FALSE)
species_cols <- setdiff(names(orth_dt), c("Orthogroup", "Bombyx_mori.protein"))

# 7. 计算 one2one, multi_match_count, total_match_count
orth_dt[, `:=`(
  one2one = Reduce(
    `+`,
    lapply(.SD, function(col)
      lengths(strsplit(col, ",\\s*")) == 1)
  ),
  multi_match_count = Reduce(
    `+`,
    lapply(.SD, function(col)
      lengths(strsplit(col, ",\\s*")) > 1)
  )
), .SDcols = species_cols]

# 新增 total_match_count = one2one + multi_match_count
orth_dt[, total_match_count := one2one + multi_match_count]

# 8. 拆分 Bombyx_mori.protein，生成 Transcript→(one2one, multi, total) 映射表
bombyx_map <- orth_dt[
  , .(
      Transcript_id      = unlist(strsplit(`Bombyx_mori.protein`, ",\\s*")),
      one2one, multi_match_count, total_match_count
    ),
  by = Orthogroup
]

# 9. 与 combined_dt 一次性左连接，并将 NA 填充为 0
setkey(bombyx_map, Transcript_id)
setkey(combined_dt, Transcript_id)
result_dt <- bombyx_map[combined_dt]

result_dt[is.na(one2one),            one2one            := 0]
result_dt[is.na(multi_match_count),  multi_match_count  := 0]
result_dt[is.na(total_match_count),  total_match_count  := 0]

# 10. 重命名列，输出最终结果
setnames(result_dt,
         c("one2one", "multi_match_count"),
         c("1:1_match_count", "multi_match_count"))

fwrite(result_dt,
       file = "DTx_nonDEG_orthogroups_matched_dt_updated.tsv",
       sep = "\t")

# 脚本运行完毕

