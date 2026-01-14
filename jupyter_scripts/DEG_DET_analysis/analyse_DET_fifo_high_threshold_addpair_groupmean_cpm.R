#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 并行 + 向量化 DTU 流程 - 高阈值版（完整脚本）
# • 统一严格预过滤：基因/转录本 CPM > 1 且 ≥3 个样本表达
# • 比较：Age、Tissue、Age–Tissue (pair) 全部两两组合
# • half‑rule：比较两组内各 ≥⌈n/2⌉ 样本 CPM>1
# • 额外导出：
#     (1) 基因 CPM 矩阵（按 Age / Tissue / Pair）
#     (2) 每个比较的 **完整** 统计表：*_allGeneStats.tsv 与 *_allTxStats.tsv
#     (3) 常规子集：DEG / DETx / DTx_nonDEG
# Date: 2025‑06‑11
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(BiocParallel)
})

# ─── 参数 ────────────────────────────────────────────────────────────────────
species_list      <- c(
    "Bombyx_mori"
)
gtf_dir           <- "/data/home/jinshuo/transcript/Genome"
counts_gene_dir   <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
counts_tx_dir     <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file         <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root       <- "edgeR_DTU_voom_highThresh_groupmean_bo"
libsize_thresh    <- 1e6
workers           <- 4
cpm_thresh        <- 1
min_samples_expr  <- 3

# ─── 元数据 ──────────────────────────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run", "Age", "Tissue.1"), c("Sample", "Age", "Tissue"))
meta[, `:=`(Age    = gsub("\\s+", "_", trimws(Age)),
            Tissue = gsub("\\s+", "_", trimws(Tissue)))]

# ─── tx2gene 构建 ────────────────────────────────────────────────────────────
build_tx2gene <- function(sp) {
  f <- sprintf("%s_tx2gene.tsv", sp)
  if (file.exists(f))
    return(fread(f, col.names = c("tx_id", "gene_id"), data.table = FALSE))
  cmd <- sprintf("grep -v '^#' %s", shQuote(file.path(gtf_dir, sprintf("%s.gtf", sp))))
  gtf <- fread(cmd = cmd, sep = "\t", header = FALSE, select = c(3, 9),
               col.names = c("feature", "attr"))
  gtf <- gtf[feature == "transcript"]
  ext <- function(x, k) sub(sprintf(".*%s \"([^\"]+)\".*", k), "\\1", x)
  tx2gene <- unique(data.frame(tx_id  = ext(gtf$attr, "transcript_id"),
                               gene_id = ext(gtf$attr, "gene_id"),
                               stringsAsFactors = FALSE))
  fwrite(tx2gene, f, sep = "\t", quote = FALSE, col.names = FALSE)
  tx2gene
}

# ─── voom 辅助（修正：应用 contrasts.fit） ─────────────────────────────────
run_voom <- function(y, meta_df, fac, tx2gene = NULL) {
  lv   <- levels(factor(meta_df[[fac]]))
  dsgn <- model.matrix(~0 + factor(meta_df[[fac]], levels = lv))
  colnames(dsgn) <- make.names(lv)
  v    <- voom(y, dsgn, plot = FALSE)
  fit  <- lmFit(v, dsgn)
  # 构建两两对比
  cmb  <- combn(seq_along(lv), 2, simplify = FALSE)
  cnms <- sapply(cmb, function(i) paste(lv[i], collapse = "_vs_"))
  ctrs <- sapply(cmb, function(i) paste0(colnames(dsgn)[i[1]], "-", colnames(dsgn)[i[2]]))
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = ctrs, levels = dsgn)))
  # 提取结果
  res <- lapply(seq_along(cnms), function(i) {
    dt <- as.data.table(topTable(fit2, coef = i, number = Inf,
                                 sort.by = "none", adjust.method = "BH"),
                         keep.rownames = "id")
    if (is.null(tx2gene)) {
      setnames(dt, "id", "gene_id")
    } else {
      setnames(dt, "id", "tx_id")
      dt[, gene_id := tx2gene$gene_id[match(tx_id, tx2gene$tx_id)]]
    }
    dt[, FDR := adj.P.Val]
    dt
  })
  names(res) <- cnms
  res
}

# ─── half-rule 判断 ─────────────────────────────────────────────────────────
half_ok <- function(mat, idx) rowSums(mat[, idx, drop = FALSE] > 1) >= ceiling(length(idx)/2)

# ─── 主函数 ─────────────────────────────────────────────────────────────────
process_species <- function(sp) {
  message("==> ", sp)

  # 读取计数矩阵
  g_cnt <- fread(file.path(counts_gene_dir, sprintf("%s_combined_counts.txt", sp)), data.table = FALSE)
  t_cnt <- fread(file.path(counts_tx_dir,   sprintf("%s_combined_counts.txt", sp)), data.table = FALSE)
  rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[, -1]
  rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[, -1]

  # 样本过滤
  sam <- intersect(colnames(g_cnt), meta$Sample)
  g_cnt <- g_cnt[, sam]; t_cnt <- t_cnt[, sam]
  sam <- setdiff(sam, names(which(colSums(g_cnt) < libsize_thresh)))
  g_cnt <- g_cnt[, sam]; t_cnt <- t_cnt[, sam]
  meta_sp <- meta[Sample %in% sam]
  meta_sp[, pair := paste(Age, Tissue, sep = "_")]

  # 基因 CPM 计算
  yg <- DGEList(counts = g_cnt)
  yg <- yg[rowSums(cpm(yg) > cpm_thresh) >= min_samples_expr, , keep.lib.sizes = FALSE]
  yg <- yg[filterByExpr(yg, group = meta_sp$Age), , keep.lib.sizes = FALSE]
  yg <- calcNormFactors(yg)
  cpm_gene <- cpm(yg)

  # 转录本 CPM 计算
  tx2gene <- build_tx2gene(sp)
  yt <- DGEList(counts = t_cnt,
                 genes = data.frame(gene_id = tx2gene$gene_id[match(rownames(t_cnt), tx2gene$tx_id)]))
  yt <- yt[rowSums(cpm(yt) > cpm_thresh) >= min_samples_expr, , keep.lib.sizes = FALSE]
  yt <- yt[filterByExpr(yt, group = meta_sp$Age), , keep.lib.sizes = FALSE]
  yt <- calcNormFactors(yt)
  cpm_tx <- cpm(yt)

  # 输出基因 CPM 矩阵
  out_combo <- file.path(output_root, "combined_gene_CPM", sp)
  dir.create(out_combo, recursive = TRUE, showWarnings = FALSE)
  save_matrix <- function(level_vec, tag) {
    dt <- data.table(gene_id = rownames(cpm_gene))
    for (lv in unique(level_vec)) {
      dt[, (lv) := rowMeans(cpm_gene[, level_vec == lv, drop = FALSE])]
    }
    fwrite(dt, file.path(out_combo, sprintf("gene_CPM_by_%s.tsv", tag)), sep = "\t", quote = FALSE)
  }
  save_matrix(meta_sp$Age,    "Age")
  save_matrix(meta_sp$Tissue, "Tissue")
  save_matrix(meta_sp$pair,   "Pair")

  # 差异分析
  gene_age <- run_voom(yg,    meta_sp, "Age")
  gene_tis <- run_voom(yg,    meta_sp, "Tissue")
  gene_pr  <- run_voom(yg,    meta_sp, "pair")
  tx_age   <- run_voom(yt,    meta_sp, "Age",    tx2gene)
  tx_tis   <- run_voom(yt,    meta_sp, "Tissue", tx2gene)
  tx_pr    <- run_voom(yt,    meta_sp, "pair",   tx2gene)

  facets <- list(
    Age    = list(gene_age,   tx_age),
    Tissue = list(gene_tis,   tx_tis),
    pair   = list(gene_pr,    tx_pr)
  )

  for (fac in names(facets)) {
    g_list <- facets[[fac]][[1]]
    t_list <- facets[[fac]][[2]]
    for (con in names(t_list)) {
      comps <- strsplit(con, "_vs_")[[1]]
      idx1  <- which(meta_sp[[fac]] == comps[1])
      idx2  <- which(meta_sp[[fac]] == comps[2])

      # Gene 统计
      dt_g <- g_list[[con]]
      ok_g <- half_ok(cpm_gene, idx1) & half_ok(cpm_gene, idx2)
      dt_g <- dt_g[ok_g]
      dt_g[, MeanCPM_A := rowMeans(cpm_gene[gene_id, idx1, drop = FALSE])]
      dt_g[, MeanCPM_B := rowMeans(cpm_gene[gene_id, idx2, drop = FALSE])]

      # Tx 统计
      dt_tx <- t_list[[con]]
      ok_t  <- half_ok(cpm_tx, idx1) & half_ok(cpm_tx, idx2)
      dt_tx <- dt_tx[ok_t]
      dt_tx[, MeanCPM_A := rowMeans(cpm_tx[tx_id, idx1, drop = FALSE])]
      dt_tx[, MeanCPM_B := rowMeans(cpm_tx[tx_id, idx2, drop = FALSE])]

      # 完整统计输出
      outdir <- file.path(output_root, fac, sp)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      fwrite(dt_g,  file.path(outdir, sprintf("%s_allGeneStats.tsv", con)), sep = "\t", quote = FALSE)
      fwrite(dt_tx, file.path(outdir, sprintf("%s_allTxStats.tsv",  con)), sep = "\t", quote = FALSE)

      # 子集输出
      deg        <- dt_g[FDR < 0.05 & abs(logFC) >= 1]
      detx       <- dt_tx[FDR < 0.05 & abs(logFC) >= 1]
      dtx_nondeg <- detx[gene_id %in% dt_g[FDR >= 0.05, gene_id]]
      fwrite(deg,        file.path(outdir, sprintf("%s_DEG.tsv",        con)), sep = "\t", quote = FALSE)
      fwrite(detx,       file.path(outdir, sprintf("%s_DETx.tsv",       con)), sep = "\t", quote = FALSE)
      fwrite(dtx_nondeg, file.path(outdir, sprintf("%s_DTx_nonDEG.tsv", con)), sep = "\t", quote = FALSE)

      message(sprintf("    %s %s | DEG:%d | DETx:%d | DTx_nonDEG:%d",
                      fac, con, nrow(deg), nrow(detx), nrow(dtx_nondeg)))
    }
  }
}

# 并行执行
bplapply(species_list, process_species,
         BPPARAM = MulticoreParam(workers = workers))
message("★ 全部物种处理完成 - 高阈值版（含完整统计输出）")

