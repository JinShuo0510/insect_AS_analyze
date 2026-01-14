#!/usr/bin/env Rscript
# ───────────────────────────────────────────────────────────────
# devDTx_nonDEG (≥2 Age levels) — v2.4
# 将 Pupa 和 Egg 归并为 Whole_body，并包含详细日志和样本库大小过滤
# Date: 2025-05-06
# ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(BiocParallel)
})

## ─── 参数设置 ───────────────────────────────────────────────────────────
species_list <- c(
  "Acyrthosiphon_pisum","Aedes_aegypti","Apis_mellifera",
  "Blattella_germanica","Bombyx_mori","Drosophila_mojavensis",
  "Gryllus_bimaculatus","Helicoverpa_armigera","Tribolium_castaneum"
)
gtf_dir        <- "/data/home/jinshuo/transcript/Genome"
gene_cnt_dir   <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
tx_cnt_dir     <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file      <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root    <- "edgeR_DTU_devDynamic"

workers         <- 4
fdr_cut         <- 0.05
logfc_cut       <- 1
min_samples_age <- 2
min_age_levels  <- 2
libsize_thresh  <- 1e6

## ─── 读取并清洗样本注释 ───────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, `:=`(
  Age    = gsub("\\s+","_", trimws(Age)),
  Tissue = gsub("\\s+","_", trimws(Tissue))
)]
# 将 Pupa 和 Egg 统一视为 Whole_body 组织
meta[Tissue %in% c("Pupa","Egg"), Tissue := "Whole_body"]

## ─── tx2gene 缓存环境与构建函数 ────────────────────────────────────────
tx2gene_cache <- new.env(parent = emptyenv())
build_tx2gene <- function(sp){
  if (exists(sp, envir = tx2gene_cache, inherits = FALSE)) {
    return(get(sp, envir = tx2gene_cache))
  }
  gtf_path <- file.path(gtf_dir, paste0(sp, ".gtf"))
  gtf <- fread(cmd = paste("grep -v '^#'", shQuote(gtf_path)),
               sep = "\t", header = FALSE, select = c(3L,9L),
               col.names = c("feature","attr"))
  gtf <- gtf[feature == "transcript"]
  extract <- function(x, key) sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  dt <- unique(data.table(
    tx_id   = extract(gtf$attr, "transcript_id"),
    gene_id = extract(gtf$attr, "gene_id")
  ))
  assign(sp, dt, envir = tx2gene_cache)
  return(dt)
}

## ─── 动态检验函数 (voom + limma) ────────────────────────────────────────
analyse_dynamic <- function(y_dge, meta_sub, n_level){
  Age <- factor(meta_sub$Age, levels = unique(meta_sub$Age))
  design_matrix <- model.matrix(~0 + Age)
  colnames(design_matrix) <- make.names(levels(Age))

  v   <- voom(y_dge, design_matrix, plot = FALSE)
  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)

  if (n_level == 2) {
    res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    dt  <- as.data.table(res, keep.rownames = "id")[, .(id, logFC, FDR = adj.P.Val)]
  } else {
    res <- topTable(fit, number = Inf, sort.by = "none")
    dt  <- as.data.table(res, keep.rownames = "id")[, .(id, F_stat = F, FDR = adj.P.Val)]
  }
  return(dt)
}

## ─── 物种主流程 ─────────────────────────────────────────────────────────
process_species <- function(sp){
  tryCatch({
    message(">>> Processing species: ", sp)
    g_cnt <- fread(file.path(gene_cnt_dir, paste0(sp, "_combined_counts.txt")), data.table = FALSE)
    rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[,-1]
    t_cnt <- fread(file.path(tx_cnt_dir, paste0(sp, "_combined_counts.txt")), data.table = FALSE)
    rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[,-1]
    tx2gene <- build_tx2gene(sp)

    for (tis in unique(meta$Tissue)) {
      message("  - Tissue: ", tis)
      idx0 <- meta[Tissue == tis & Sample %in% colnames(g_cnt)]$Sample
      message("    Original samples: ", length(idx0))
      if (length(idx0) < min_samples_age) { message("    Skipped: < min_samples_age"); next }

      lib_sizes <- colSums(g_cnt[, idx0, drop = FALSE])
      keep_samp <- names(lib_sizes)[is.finite(lib_sizes) & lib_sizes >= libsize_thresh]
      message("    After libsize filter: ", length(keep_samp))
      if (length(keep_samp) < min_samples_age) { message("    Skipped: < min_samples_age after libsize"); next }

      g_mat <- g_cnt[, keep_samp, drop = FALSE]
      t_mat <- t_cnt[, keep_samp, drop = FALSE]
      meta_t <- meta[Sample %in% keep_samp]

      n_level <- length(unique(meta_t$Age))
      message("    Age levels: ", n_level, " (", paste(unique(meta_t$Age), collapse=","), ")")
      if (n_level < min_age_levels) { message("    Skipped: < min_age_levels"); next }

      y_g <- DGEList(counts = g_mat)
      keep_g <- filterByExpr(y_g, group = meta_t$Age)
      message("    Genes before/after filter: ", nrow(y_g), "/", sum(keep_g))
      y_g <- y_g[keep_g, , keep.lib.sizes = FALSE]; y_g <- calcNormFactors(y_g)
      g_res <- analyse_dynamic(y_g, meta_t, n_level)

      y_tx <- DGEList(
        counts = t_mat,
        genes  = data.frame(gene_id = tx2gene$gene_id[match(rownames(t_mat), tx2gene$tx_id)])
      )
      keep_tx <- filterByExpr(y_tx, group = meta_t$Age)
      message("    Transcripts before/after filter: ", nrow(y_tx), "/", sum(keep_tx))
      y_tx <- y_tx[keep_tx, , keep.lib.sizes = FALSE]; y_tx <- calcNormFactors(y_tx)
      tx_res <- analyse_dynamic(y_tx, meta_t, n_level)
      tx_res[, gene_id := y_tx$genes$gene_id[match(id, rownames(y_tx))]]

      if (n_level == 2) {
        DEG <- g_res[FDR < fdr_cut & abs(logFC) >= logfc_cut, id]
        DTx <- tx_res[FDR < fdr_cut & abs(logFC) >= logfc_cut]
      } else {
        DEG <- g_res[FDR < fdr_cut, id]
        DTx <- tx_res[FDR < fdr_cut]
      }
      message("    DEG/DTx counts: ", length(DEG), "/", nrow(DTx))
      devDTx_nonDEG <- DTx[!gene_id %in% DEG]
      message("    devDTx_nonDEG: ", nrow(devDTx_nonDEG))

      outdir <- file.path(output_root, sp, tis)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      fwrite(devDTx_nonDEG,
             file = file.path(outdir, "devDTx_nonDEG.tsv"),
             sep = "\t", quote = FALSE)
    }
  }, error = function(e){ message("!!! Error in species ", sp, ": ", e$message) })
}

## ─── 并行执行 ─────────────────────────────────────────────────────────
bplapply(species_list, process_species,
         BPPARAM = MulticoreParam(workers = workers))
message("★ devDTx_nonDEG v2.4 completed — Pupa/Egg merged to Whole_body")

