# 加载必要包
library(dplyr)
library(purrr)
library(stringr)

# 1. 定义目录
age_dir    <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/Age/Bombyx_mori"
tissue_dir <- "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/Tissue/Bombyx_mori"

# 2. 获取所有 *_DTx_nonDEG.tsv 文件路径
age_files    <- list.files(age_dir,    pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)
tissue_files <- list.files(tissue_dir, pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)

# 3. 定义单文件处理函数
process_file <- function(file_path, type_label) {
  df <- read.delim(file_path, stringsAsFactors = FALSE)
  fname <- basename(file_path)
  comp  <- sub("_DTx_nonDEG\\.tsv$", "", fname)
  parts <- str_split(comp, "_vs_")[[1]]
  df %>%
    mutate(
      Compare_type = type_label,
      Compare1     = parts[1],
      Compare2     = parts[2],
      Transcript_id = tx_id,
      Gene_id       = gene_id
    ) %>%
    select(Compare_type, Compare1, Compare2, Transcript_id, Gene_id,
           logFC, AveExpr, FDR)
}

# 4. 逐目录读取并合并
df_age    <- map_df(age_files,    process_file, type_label = "Age")
df_tissue <- map_df(tissue_files, process_file, type_label = "Tissue")
combined_df <- bind_rows(df_age, df_tissue)

# 5. 读取orthogroups表格（保留原始列名）
orth <- read.delim("Lepidoptera_Orthogroups.tsv",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

# 6. 统计各物种列（除Bombyx_mori）中1:1与多对多匹配数
species_cols <- setdiff(names(orth), c("Orthogroup", "Bombyx_mori.protein"))

count_matches <- function(tx) {
  row <- orth[grep(paste0("\\b", tx, "\\b"), orth$`Bombyx_mori.protein`), ]
  if (nrow(row)==0) return(c("1:1_match_count"=0, "match_count"=0))
  row <- row[1, ]
  one2one <- 0; many <- 0
  for (col in species_cols) {
    val <- row[[col]]
    if (is.na(val) || val=="") next
    n <- length(str_split(val, ",\\s*")[[1]])
    if (n==1) one2one <- one2one + 1 else many <- many + 1
  }
  c("1:1_match_count"=one2one, "match_count"=many)
}

# 7. 应用匹配统计并拼回主表
match_df <- map_dfr(combined_df$Transcript_id, count_matches)
result_df <- bind_cols(combined_df, match_df)

# 8. 写出结果
write.table(result_df,
            file = "DTx_nonDEG_orthogroups_matched.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

