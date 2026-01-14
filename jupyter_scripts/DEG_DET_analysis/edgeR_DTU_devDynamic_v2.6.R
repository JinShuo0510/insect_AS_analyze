#!/usr/bin/env Rscript
# ───────────────────────────────────────────────────────────────
# edgeR_DTU_devDynamic_v2.6.R
# - 动态剪接 (devDTx_nonDEG) 分析，支持：
#   • Egg/Pupa -> Whole_body
#   • 四重过滤检查 + 日志
#   • 输出 devDTx_nonDEG 和 剩余 DTx_nonDEG，含外显子坐标
#   • FDR 排序，组织到文件
# Date: 2025-05-06
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(BiocParallel)
})

species_list <- c(
  "Acyrthosiphon_pisum","Aedes_aegypti","Apis_mellifera",
  "Blattella_germanica","Bombyx_mori","Drosophila_mojavensis",
  "Gryllus_bimaculatus","Helicoverpa_armigera","Tribolium_castaneum"
)

gtf_dir      <- "/data/home/jinshuo/transcript/Genome"
gene_cnt_dir <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
tx_cnt_dir   <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file    <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root  <- "edgeR_DTU_devDynamic"

workers        <- 4
fdr_cut        <- 0.05
logfc_cut      <- 1
min_samples_age<- 2
min_age_levels <- 2
libsize_thresh <- 1e6

# ───────────────────────────────────────────────────────────────
# 读取元数据，标准化组织名
# ───────────────────────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, `:=`(
  Age    = gsub("\\s+","_", trimws(Age)),
  Tissue = gsub("\\s+","_", trimws(Tissue))
)]
meta[Tissue %in% c("Pupa","Egg"), Tissue := "Whole_body"]

# ───────────────────────────────────────────────────────────────
# tx2gene / tx_exons 提取工具
# ───────────────────────────────────────────────────────────────
tx2gene_cache <- new.env(parent = emptyenv())
build_tx2gene <- function(sp){
  if (exists(sp, envir = tx2gene_cache)) return(get(sp, envir = tx2gene_cache))
  fn <- file.path(gtf_dir, paste0(sp, ".gtf"))
  gtf <- fread(cmd = paste("grep -v '^#'", shQuote(fn)), sep = "\t", header = FALSE,
               select = c(3L,9L), col.names = c("feature","attr"))
  gtf <- gtf[feature == "transcript"]
  extract <- function(x, k) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", x)
  dt <- unique(data.table(tx_id = extract(gtf$attr, "transcript_id"),
                          gene_id = extract(gtf$attr, "gene_id")))
  assign(sp, dt, envir = tx2gene_cache)
  dt
}

build_tx_exons <- function(sp){
  fn <- file.path(gtf_dir, paste0(sp, ".gtf"))
  gtf <- fread(cmd = paste("grep -v '^#'", shQuote(fn)), sep = "\t", header = FALSE,
               col.names = c("chr","src","feature","start","end","score","strand","phase","attr"))
  gtf <- gtf[feature == "exon"]
  extract <- function(x, k) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", x)
  data.table(
    tx_id = extract(gtf$attr, "transcript_id"),
    exon_chr = gtf$chr,
    exon_start = gtf$start,
    exon_end = gtf$end
  )
}

# ───────────────────────────────────────────────────────────────
# voom+limma 动态检验
# ───────────────────────────────────────────────────────────────
analyse_dynamic <- function(y_dge, meta_sub, n_level){
  Age <- factor(meta_sub$Age, levels = unique(meta_sub$Age))
  design_matrix <- model.matrix(~0 + Age)
  colnames(design_matrix) <- make.names(levels(Age))
  v <- voom(y_dge, design_matrix, plot = FALSE)
  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)
  if (n_level == 2) {
    res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    as.data.table(res, keep.rownames = "id")[, .(id, logFC, FDR = adj.P.Val)]
  } else {
    res <- topTable(fit, number = Inf, sort.by = "none")
    as.data.table(res, keep.rownames = "id")[, .(id, F_stat = F, FDR = adj.P.Val)]
  }
}

# ───────────────────────────────────────────────────────────────
# 每个物种主流程
# ───────────────────────────────────────────────────────────────
process_species <- function(sp){
  tryCatch({
    message(">>> ", sp)
    g_cnt <- fread(file.path(gene_cnt_dir, paste0(sp, "_combined_counts.txt")), data.table = FALSE)
    rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[,-1]
    t_cnt <- fread(file.path(tx_cnt_dir, paste0(sp, "_combined_counts.txt")), data.table = FALSE)
    rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[,-1]
    tx2gene <- build_tx2gene(sp)
    coords <- build_tx_exons(sp)

    for (tis in unique(meta$Tissue)) {
      meta_sub <- meta[Tissue == tis & Sample %in% colnames(g_cnt)]
      if (nrow(meta_sub) < min_samples_age) next
      libs <- colSums(g_cnt[, meta_sub$Sample, drop = FALSE])
      keep_samples <- names(libs)[libs >= libsize_thresh & is.finite(libs)]
      meta_sub <- meta_sub[Sample %in% keep_samples]
      if (nrow(meta_sub) < min_samples_age) next
      n_level <- length(unique(meta_sub$Age))
      if (n_level < min_age_levels) next

      g_mat <- g_cnt[, meta_sub$Sample, drop = FALSE]
      t_mat <- t_cnt[, meta_sub$Sample, drop = FALSE]

      y_g <- DGEList(counts = g_mat)
      y_g <- y_g[filterByExpr(y_g, group = meta_sub$Age), , keep.lib.sizes = FALSE]
      y_g <- calcNormFactors(y_g)
      g_res <- analyse_dynamic(y_g, meta_sub, n_level)

      y_tx <- DGEList(counts = t_mat,
                      genes = data.frame(gene_id = tx2gene$gene_id[match(rownames(t_mat), tx2gene$tx_id)]))
      y_tx <- y_tx[filterByExpr(y_tx, group = meta_sub$Age), , keep.lib.sizes = FALSE]
      y_tx <- calcNormFactors(y_tx)
      tx_res <- analyse_dynamic(y_tx, meta_sub, n_level)
      tx_res[, gene_id := y_tx$genes$gene_id[match(id, rownames(y_tx))]]

      stable_genes <- g_res[FDR >= fdr_cut, id]
      tx_nondeg_all <- tx_res[gene_id %in% stable_genes]

      if (n_level == 2) {
        sig    <- tx_nondeg_all[FDR < fdr_cut & abs(logFC) >= logfc_cut]
        nonsig <- tx_nondeg_all[!(FDR < fdr_cut & abs(logFC) >= logfc_cut)]
      } else {
        sig    <- tx_nondeg_all[FDR < fdr_cut]
        nonsig <- tx_nondeg_all[FDR >= fdr_cut]
      }

      dev_out <- merge(sig, coords, by.x = "id", by.y = "tx_id", allow.cartesian = TRUE)
      setorder(dev_out, FDR)

      other_out <- merge(nonsig, coords, by.x = "id", by.y = "tx_id", allow.cartesian = TRUE)
      setorder(other_out, FDR)

      outdir <- file.path(output_root, sp, tis)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

      fwrite(dev_out,
             file = file.path(outdir, "devDTx_nonDEG_with_exons.tsv"),
             sep = "\t", quote = FALSE, row.names = FALSE)
      fwrite(other_out,
             file = file.path(outdir, "DTx_nonDEG_with_exons.tsv"),
             sep = "\t", quote = FALSE, row.names = FALSE)

      message(sprintf("    tissue %-12s | devDTx_nonDEG: %4d | other DTx_nonDEG: %4d",
                      tis, nrow(dev_out), nrow(other_out)))
    }
  }, error = function(e){
    message("!!! Error in ", sp, ": ", e$message)
  })
}

# ───────────────────────────────────────────────────────────────
# 并行执行
# ───────────────────────────────────────────────────────────────
bplapply(species_list, process_species,
         BPPARAM = MulticoreParam(workers = workers))

message("★ devDTx_nonDEG v2.6 分析完成")

