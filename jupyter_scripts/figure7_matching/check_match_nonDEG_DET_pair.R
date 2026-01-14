#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 批量提取 Le_bombyx_09 转录本在多个比较结果中的信息，包括 Age、Tissue 与 Pair 三类比较
# Date: 2025-05-13
# ─────────────────────────────────────────────────────────────────────────────

# 加载必要包（假如需要并行或其他包可在此添加）
suppressPackageStartupMessages({
  # library(data.table)
  # library(dplyr)
})

# 1. 读取转录本列表（去掉第一行物种名）
tran_ids <- readLines("Le_bombyx_09_tran_id.txt")[-1]

# 2. 定义三个目录：Age、Tissue、Pair
dirs <- list(
  Age    = "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/Age/Bombyx_mori",
  Tissue = "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/Tissue/Bombyx_mori",
  Pair   = "/data/home/jinshuo/jupyter-notebook/DEG_and_DET/edgeR_DTU_voom_highThresh/pair/Bombyx_mori"
)

# 3. 用于存放各比较结果的列表
results <- list()

# 4. 遍历各类型目录并筛选
for (type in names(dirs)) {
  files <- list.files(dirs[[type]], pattern = "_DTx_nonDEG\\.tsv$", full.names = TRUE)
  
  for (f in files) {
    fname <- basename(f)
    comp  <- sub("_DTx_nonDEG\\.tsv$", "", fname)
    parts <- strsplit(comp, "_vs_")[[1]]
    comp1 <- parts[1]
    comp2 <- parts[2]
    
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    # 只保留关注的转录本
    df_mat <- df[df$tx_id %in% tran_ids, ]
    
    if (nrow(df_mat) > 0) {
      df_out <- data.frame(
        `Compare type` = type,
        Compare1       = comp1,
        Compare2       = comp2,
        Transcript_id  = df_mat$tx_id,
        Gene_id        = df_mat$gene_id,
        logFC          = df_mat$logFC,
        AveExpr        = df_mat$AveExpr,
        FDR            = df_mat$FDR,
        stringsAsFactors = FALSE,
        check.names      = FALSE
      )
      results[[length(results) + 1]] <- df_out
    }
  }
}

# 5. 合并所有结果并写出到文件
final_df <- do.call(rbind, results)
write.table(
  final_df,
  file      = "matched_Le_bombyx_results.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

