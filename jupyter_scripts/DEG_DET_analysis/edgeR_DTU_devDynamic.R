#!/usr/bin/env Rscript
# ───────────────────────────────────────────────────────────────
# 计算“发育动态”的 DTx+nonDEG（devDTx_nonDEG）
# 逻辑：每个物种 × 每个组织内，Age 为唯一分组因子
# Date: 2025‑05‑06
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(BiocParallel)
})

## ─── 可调参数 ──────────────────────────────────────────────────
species_list    <- c("Acyrthosiphon_pisum","Aedes_aegypti","Apis_mellifera",
                     "Blattella_germanica","Bombyx_mori","Drosophila_mojavensis",
                     "Gryllus_bimaculatus","Helicoverpa_armigera","Tribolium_castaneum")
gtf_dir         <- "/data/home/jinshuo/transcript/Genome"
gene_cnt_dir    <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
tx_cnt_dir      <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file       <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root     <- "edgeR_DTU_devDynamic"
workers         <- 4
fdr_cut         <- 0.05      # 与前脚本保持一致
logfc_cut       <- 1
min_samples_age <- 2         # Age 内至少 2 个样本才分析
min_age_levels  <- 3         # 至少 3 个时期才判“动态”

## ─── 样本信息 ────────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, `:=`(Age = gsub("\\s+","_", trimws(Age)),
            Tissue = gsub("\\s+","_", trimws(Tissue)))]

## ─── 工具函数 ────────────────────────────────────────────────
build_tx2gene <- local({
  cache <- list()
  function(sp){
    if(!is.null(cache[[sp]])) return(cache[[sp]])
    fn <- file.path(gtf_dir, paste0(sp,".gtf"))
    gtf <- fread(cmd = paste("grep -v '^#'", shQuote(fn)),
                 sep="\t", header=FALSE, select=c(3L,9L),
                 col.names=c("feature","attr"))
    gtf <- gtf[feature=="transcript"]
    extr <- function(x,k) sub(paste0('.*',k,' "([^"]+)".*'),"\\1",x)
    dt <- unique(data.table(tx_id = extr(gtf$attr,"transcript_id"),
                            gene_id= extr(gtf$attr,"gene_id")))
    cache[[sp]] <<- dt
    dt
  }
})

analyse_age <- function(cnt_mat, meta_sub, fac_levels, design_prefix){
  # fac_levels 已保证有序
  grp      <- factor(meta_sub$Age, levels = fac_levels)
  design   <- model.matrix(~0+grp)
  colnames(design) <- make.names(levels(grp))
  v        <- voom(cnt_mat, design, plot = FALSE)
  fit      <- lmFit(v, design)

  ## --- 全局 F 检验：是否随 Age 有显著差异 ---
  fitF     <- eBayes(fit)
  ttF      <- topTableF(fitF, number=Inf, sort.by="none")
  F.tbl    <- data.table(id = rownames(ttF),
                         F_stat = ttF$F,
                         P.Value = ttF$P.Value,
                         FDR = p.adjust(ttF$P.Value, method="BH"))
  F.tbl
}

## ─── 物种 × 组织 主流程 ─────────────────────────────────────
process_species <- function(sp){

  message(">>> ", sp)
  g_cnt <- fread(file.path(gene_cnt_dir, paste0(sp,"_combined_counts.txt")), data.table=FALSE)
  rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[,-1]
  t_cnt <- fread(file.path(tx_cnt_dir, paste0(sp,"_combined_counts.txt")), data.table=FALSE)
  rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[,-1]
  tx2gene <- build_tx2gene(sp)

  ## 每个组织单独分析
  tissues_sp <- intersect(meta$Tissue, colnames(g_cnt))
  for(tis in unique(meta$Tissue)){
    idx_t <- meta[Tissue==tis & Sample %in% colnames(g_cnt)]$Sample
    if(length(idx_t) < min_samples_age) next        # 样本太少

    meta_t <- meta[Sample %in% idx_t]
    fac_lv <- unique(meta_t$Age)
    if(length(fac_lv) < min_age_levels) next        # Age 水平不足

    g_mat  <- g_cnt[, idx_t, drop=FALSE]
    tx_mat <- t_cnt[, idx_t, drop=FALSE]

    ## --- gene ---
    y_g <- DGEList(counts=g_mat)
    keep_g <- filterByExpr(y_g, group=meta_t$Age)
    y_g <- y_g[keep_g,,keep.lib.sizes=FALSE];  y_g <- calcNormFactors(y_g)
    g_res <- analyse_age(y_g, meta_t, fac_lv, "gene")

    ## --- tx ---
    y_tx <- DGEList(counts=tx_mat,
                    genes=data.frame(gene_id = tx2gene$gene_id[match(rownames(tx_mat), tx2gene$tx_id)]))
    keep_tx <- filterByExpr(y_tx, group=meta_t$Age)
    y_tx <- y_tx[keep_tx,,keep.lib.sizes=FALSE];  y_tx <- calcNormFactors(y_tx)
    tx_res <- analyse_age(y_tx, meta_t, fac_lv, "tx")
    tx_res[, gene_id := y_tx$genes$gene_id[match(id, rownames(y_tx))]]

    ## --- 分类 ---
    DEG_age  <- g_res[FDR < fdr_cut]
    DTx_age  <- tx_res[FDR < fdr_cut]
    devDTx_nonDEG <- DTx_age[ !gene_id %in% DEG_age$id &
                               abs(F_stat) >= logfc_cut ]  # logFC proxy：F_stat 大体相关

    ## --- 输出 ---
    outdir <- file.path(output_root, sp, tis)
    dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
    fwrite(devDTx_nonDEG,
           file=file.path(outdir, "devDTx_nonDEG.tsv"),
           sep="\t", quote=FALSE, row.names=FALSE)

    ## 可选：汇总代码
    summary_dt <- devDTx_nonDEG[, .(
        tissues = paste(unique(tis), collapse=","),
        max_F   = max(F_stat)
      ), by=.(gene_id, id)]
    fwrite(summary_dt,
           file=file.path(outdir, "devDTx_summary.tsv"),
           sep="\t", quote=FALSE, row.names=FALSE)

    message("    tissue ", tis,
            " | devDTx_nonDEG: ", nrow(devDTx_nonDEG))
  }
}

## ─── 并行执行 ────────────────────────────────────────────────
bplapply(species_list, process_species,
         BPPARAM = MulticoreParam(workers = workers))

message("★ devDTx_nonDEG 分析完成")

