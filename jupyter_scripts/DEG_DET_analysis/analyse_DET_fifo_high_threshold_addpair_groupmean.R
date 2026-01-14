#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 并行 + 向量化 DTU 流程  ‑ 高阈值版
# • 基因与转录本层面统一严格预过滤：CPM > 1 且 ≥3 个样本表达
# • Age、Tissue 以及 Age–Tissue 组合 (pair) 全部两两比较
# • 组内过滤采用“半数规则”：在每个比较组内 ≥⌈n/2⌉ 个样本 CPM > 1 即保留
# • 输出文件包含组内线性 CPM 平均值 MeanCPM_A / MeanCPM_B
# Date: 2025‑05‑12
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(limma)
    library(BiocParallel)
})

# ─── 参数 ────────────────────────────────────────────────────────────────────

#species_list      <- c(
#    "Bombyx_mori"
#)
species_list <- readLines("species_list.txt")
gtf_dir           <- "/data/home/jinshuo/transcript/Genome"
counts_gene_dir   <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
counts_tx_dir     <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file         <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root       <- "edgeR_DTU_voom_highThresh_groupmean"
libsize_thresh    <- 1e6
workers           <- 4

# 统一严格过滤阈值
cpm_thresh        <- 1     # CPM > 1
min_samples_expr  <- 3     # ≥3 个样本表达

# ─── 样本注释 ────────────────────────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, Age    := gsub("\\s+","_", trimws(Age))]
meta[, Tissue := gsub("\\s+","_", trimws(Tissue))]

# ─── tx2gene 构建 ────────────────────────────────────────────────────────────
build_tx2gene <- function(sp) {
    f <- paste0(sp, "_tx2gene.tsv")
    if (file.exists(f)) return(fread(f, col.names=c("tx_id","gene_id"), data.table=FALSE))

    cmd <- paste("grep -v '^#'", shQuote(file.path(gtf_dir, paste0(sp, ".gtf"))))
    gtf <- fread(cmd=cmd, sep="\t", header=FALSE, select=c(3,9),
                 col.names=c("feature","attr"))
    gtf <- gtf[feature=="transcript"]
    ext <- function(x,k) sub(paste0('.*',k,' "([^"]+)".*'),'\\1',x)
    tx2gene <- unique(data.frame(
        tx_id   = ext(gtf$attr,"transcript_id"),
        gene_id = ext(gtf$attr,"gene_id"), stringsAsFactors=FALSE))
    fwrite(tx2gene, f, sep="\t", quote=FALSE, col.names=FALSE)
    tx2gene
}

# ─── voom 差异分析 ───────────────────────────────────────────────────────────
analyse_voom <- function(y_dge, meta_sp, fac, tx2gene=NULL) {
    lv   <- levels(factor(meta_sp[[fac]]))
    dsgn <- model.matrix(~0 + factor(meta_sp[[fac]], levels=lv))
    colnames(dsgn) <- make.names(lv)

    v   <- voom(y_dge, dsgn, plot=FALSE)
    fit <- lmFit(v, dsgn)

    cmb  <- combn(seq_along(lv), 2, simplify=FALSE)
    cnms <- sapply(cmb, \(i) paste(lv[i], collapse="_vs_"))
    ctrs <- sapply(cmb, \(i) paste0(colnames(dsgn)[i[1]],"-",colnames(dsgn)[i[2]]))
    fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts=ctrs, levels=dsgn)))

    res <- lapply(seq_along(cnms), function(i){
        dt <- as.data.table(topTable(fit2, coef=i, number=Inf, sort.by="none",
                                     adjust.method="BH"), keep.rownames="id")
        if (is.null(tx2gene)) {
            setnames(dt,"id","gene_id")
        } else {
            setnames(dt,"id","tx_id")
            dt[, gene_id := tx2gene$gene_id[match(tx_id, tx2gene$tx_id)]]
        }
        dt[, FDR := adj.P.Val]
        dt
    })
    names(res) <- cnms
    res
}

# ─── 主流程 ────────────────────────────────────────────────────────────────
process_species <- function(sp){
    message(">>> 处理物种：", sp)

    # 读入 count
    g_cnt <- fread(file.path(counts_gene_dir, paste0(sp,"_combined_counts.txt")), data.table=FALSE)
    t_cnt <- fread(file.path(counts_tx_dir,   paste0(sp,"_combined_counts.txt")), data.table=FALSE)
    rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[,-1]
    rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[,-1]

    # 样本过滤
    sam0 <- intersect(colnames(g_cnt), meta$Sample)
    g_cnt <- g_cnt[, sam0, drop=FALSE]; t_cnt <- t_cnt[, sam0, drop=FALSE]
    keep  <- setdiff(sam0, names(which(colSums(g_cnt) < libsize_thresh)))
    g_cnt <- g_cnt[, keep, drop=FALSE]; t_cnt <- t_cnt[, keep, drop=FALSE]
    meta_sp <- meta[Sample %in% keep]
    meta_sp[, pair := paste(Age, Tissue, sep="_")]

    # ── 基因层面预过滤 ──
    y_gene <- DGEList(counts=g_cnt)
    keep_g <- rowSums(cpm(y_gene) > cpm_thresh) >= min_samples_expr
    y_gene <- y_gene[keep_g,,keep.lib.sizes=FALSE]
    y_gene <- y_gene[filterByExpr(y_gene, group=meta_sp$Age),,keep.lib.sizes=FALSE]
    y_gene <- calcNormFactors(y_gene)
    cpm_gene_mat <- cpm(y_gene)

    # ── 转录本层面预过滤 ──
    tx2gene <- build_tx2gene(sp)
    y_tx <- DGEList(counts=t_cnt,
                    genes=data.frame(gene_id=tx2gene$gene_id[match(rownames(t_cnt),tx2gene$tx_id)]))
    keep_tx <- rowSums(cpm(y_tx) > cpm_thresh) >= min_samples_expr
    y_tx <- y_tx[keep_tx,,keep.lib.sizes=FALSE]
    y_tx <- y_tx[filterByExpr(y_tx, group=meta_sp$Age),,keep.lib.sizes=FALSE]
    y_tx <- calcNormFactors(y_tx)
    cpm_tx_mat <- cpm(y_tx)

    # 差异分析
    gene_age    <- analyse_voom(y_gene, meta_sp,"Age")
    tx_age      <- analyse_voom(y_tx,   meta_sp,"Age",    tx2gene)
    gene_tissue <- analyse_voom(y_gene, meta_sp,"Tissue")
    tx_tissue   <- analyse_voom(y_tx,   meta_sp,"Tissue", tx2gene)
    gene_pair   <- analyse_voom(y_gene, meta_sp,"pair")
    tx_pair     <- analyse_voom(y_tx,   meta_sp,"pair",   tx2gene)

    # half‑rule helper
    half_ok <- function(mat, idx) rowSums(mat[,idx,drop=FALSE] > 1) >= ceiling(length(idx)/2)

    # 输出
    for (fac in c("Age","Tissue","pair")){
        g_res <- get(paste0("gene_",tolower(fac)))
        t_res <- get(paste0("tx_",tolower(fac)))
        for (con in names(t_res)){
            grps <- strsplit(con,"_vs_")[[1]]
            idx1 <- which(meta_sp[[fac]]==grps[1])
            idx2 <- which(meta_sp[[fac]]==grps[2])

            # 基因 half‑rule + 组均 CPM
            dt_g <- g_res[[con]]
            ok_g <- half_ok(cpm_gene_mat,idx1) & half_ok(cpm_gene_mat,idx2)
            dt_g <- dt_g[ok_g]
            dt_g[, MeanCPM_A := rowMeans(cpm_gene_mat[gene_id, idx1, drop=FALSE])]
            dt_g[, MeanCPM_B := rowMeans(cpm_gene_mat[gene_id, idx2, drop=FALSE])]

            # 转录本 half‑rule + 组均 CPM
            dt_tx <- t_res[[con]]
            ok_t <- half_ok(cpm_tx_mat,idx1) & half_ok(cpm_tx_mat,idx2)
            dt_tx <- dt_tx[ok_t]
            dt_tx[, MeanCPM_A := rowMeans(cpm_tx_mat[tx_id, idx1, drop=FALSE])]
            dt_tx[, MeanCPM_B := rowMeans(cpm_tx_mat[tx_id, idx2, drop=FALSE])]

            # DEG / DETx / DTx_nonDEG
            deg   <- dt_g[FDR<0.05 & abs(logFC)>=1]
            detx  <- dt_tx[FDR<0.05 & abs(logFC)>=1]
            dtx_nondeg <- detx[gene_id %in% dt_g[FDR>=0.05,gene_id]]

            outdir <- file.path(output_root, fac, sp)
            dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
            fwrite(deg,         file.path(outdir, paste0(con,"_DEG.tsv")),        sep="\t",quote=FALSE)
            fwrite(detx,        file.path(outdir, paste0(con,"_DETx.tsv")),       sep="\t",quote=FALSE)
            fwrite(dtx_nondeg,  file.path(outdir, paste0(con,"_DTx_nonDEG.tsv")), sep="\t",quote=FALSE)

            message(sprintf("    %s %s | DEG:%d | DETx:%d | DTx_nonDEG:%d",
                            fac, con, nrow(deg), nrow(detx), nrow(dtx_nondeg)))
        }
    }
}

# 并行执行
bplapply(species_list, process_species,
         BPPARAM=MulticoreParam(workers=workers))

message("★ 全部物种处理完成 - 高阈值版")

