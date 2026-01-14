#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 并行 + 向量化 DTU 流程  - 高阈值版（共享单一规模因子 + 变量持久化）
# • 基因与转录本共用同一套 TMM norm.factors，确保 CPM 一致
# • 样本注释 (meta) 行顺序显式对齐 count 矩阵列顺序
# • Age、Tissue 以及 Age–Tissue 组合 (pair) 全部两两比较
# • 组内过滤采用“半数规则”：在每个比较组内 ≥⌈n/2⌉ 个样本 CPM > 1 即保留
# • 输出文件包含组内线性 CPM 平均值 MeanCPM_A / MeanCPM_B
# • 保存关键 R 对象至工作目录，便于日后追溯
# Date: 2025-05-16
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(limma)
    library(BiocParallel)
})

# ─── 参数 ────────────────────────────────────────────────────────────────────
#species_list      <- c("Bombyx_mori")
species_list <- readLines("species_list.txt")
gtf_dir           <- "/data/home/jinshuo/transcript/Genome"
counts_gene_dir   <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
counts_tx_dir     <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file         <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root       <- "edgeR_DTU_voom_highThresh_groupmean"
libsize_thresh    <- 1e4
workers           <- 4

cpm_thresh       <- 1     # CPM > 1
min_samples_expr <- 1     # ≥3 个样本表达

# ─── 样本注释 ────────────────────────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, Age    := gsub("\\s+","_", trimws(Age))]
meta[, Tissue := gsub("\\s+","_", trimws(Tissue))]
setkey(meta, Sample)

# ─── tx2gene 构建 ────────────────────────────────────────────────────────────
build_tx2gene <- function(sp) {
    f <- paste0(sp, "_tx2gene.tsv")
    if (file.exists(f)) return(fread(f, col.names=c("tx_id","gene_id"), data.table=FALSE))
    cmd <- paste("grep -v '^#'", shQuote(file.path(gtf_dir, paste0(sp, ".gtf"))))
    gtf <- fread(cmd=cmd, sep="\t", header=FALSE, select=c(3,9),
                 col.names=c("feature","attr"))
    gtf <- gtf[feature=="transcript"]
    ext <- function(x,k) sub(paste0(".*",k," \"([^\"]+)\".*"),'\\1',x)
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
    cnms <- sapply(cmb, function(i) paste(lv[i], collapse="_vs_"))
    ctrs <- sapply(cmb, function(i) paste0(colnames(dsgn)[i[1]],"-",colnames(dsgn)[i[2]]))
    fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts=ctrs, levels=dsgn)))
    res <- lapply(seq_along(cnms), function(i){
        dt <- as.data.table(topTable(fit2, coef=i, number=Inf, sort.by="none", adjust.method="BH"), keep.rownames="id")
        if (is.null(tx2gene)) setnames(dt, "id", "gene_id") else {
            setnames(dt, "id", "tx_id")
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
    g_cnt <- fread(file.path(counts_gene_dir, paste0(sp, "_combined_counts.txt")), data.table=FALSE)
    t_cnt <- fread(file.path(counts_tx_dir,   paste0(sp, "_combined_counts.txt")), data.table=FALSE)
    rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[, -1]
    rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[, -1]

    # 样本过滤：基于库大小过滤不合格样本
    sam0 <- intersect(colnames(g_cnt), meta$Sample)
    g_cnt <- g_cnt[, sam0, drop=FALSE]
    t_cnt <- t_cnt[, sam0, drop=FALSE]
    
    # 计算库大小（基因和转录本）
    lib_g <- colSums(g_cnt)
    lib_t <- colSums(t_cnt)
    
    # 标记库大小为0和低于阈值的样本
    bad_g0    <- names(lib_g)[lib_g == 0]
    bad_t0    <- names(lib_t)[lib_t == 0]
    bad_low_g <- names(lib_g)[lib_g < libsize_thresh]
    bad_low_t <- names(lib_t)[lib_t < libsize_thresh]
    bad_samples <- unique(c(bad_g0, bad_t0, bad_low_g, bad_low_t))
    if (length(bad_samples) > 0) {
        message(sprintf("移除不合格样本（库大小0或小于%g）：%s", libsize_thresh, paste(bad_samples, collapse=", ")))
    }
    # 保留高质量样本
    keep <- setdiff(sam0, bad_samples)
    g_cnt <- g_cnt[, keep, drop=FALSE]
    t_cnt <- t_cnt[, keep, drop=FALSE]

    # 对齐 meta 排序
    meta_sp <- meta[colnames(g_cnt)]
    meta_sp[, pair := paste(Age, Tissue, sep="_")]

    # 统一 TMM 规模因子
    temp_dge <- calcNormFactors(DGEList(counts = t_cnt))
    nf <- temp_dge$samples$norm.factors
    y_gene <- DGEList(counts = g_cnt); y_gene$samples$norm.factors <- nf
    tx2gene <- build_tx2gene(sp)
    y_tx   <- DGEList(counts = t_cnt, genes = data.frame(gene_id = tx2gene$gene_id[match(rownames(t_cnt), tx2gene$tx_id)]))
    y_tx$samples$norm.factors <- nf

    # 预过滤
    keep_g <- rowSums(cpm(y_gene) > cpm_thresh) >= min_samples_expr
    y_gene <- y_gene[keep_g,, keep.lib.sizes=FALSE]
    y_gene <- y_gene[filterByExpr(y_gene, group=meta_sp$Age),, keep.lib.sizes=FALSE]
    cpm_gene_mat <- cpm(y_gene)

    keep_tx <- rowSums(cpm(y_tx) > cpm_thresh) >= min_samples_expr
    y_tx <- y_tx[keep_tx,, keep.lib.sizes=FALSE]
    y_tx <- y_tx[filterByExpr(y_tx, group=meta_sp$Age),, keep.lib.sizes=FALSE]
    cpm_tx_mat <- cpm(y_tx)

    # 差异分析
    gene_age    <- analyse_voom(y_gene, meta_sp, "Age")
    tx_age      <- analyse_voom(y_tx,   meta_sp, "Age",    tx2gene)
    gene_tissue <- analyse_voom(y_gene, meta_sp, "Tissue")
    tx_tissue   <- analyse_voom(y_tx,   meta_sp, "Tissue", tx2gene)
    gene_pair   <- analyse_voom(y_gene, meta_sp, "pair")
    tx_pair     <- analyse_voom(y_tx,   meta_sp, "pair",   tx2gene)

#    half_ok <- function(mat, idx) rowSums(mat[, idx, drop=FALSE] > 1) >= ceiling(length(idx)/2)
    half_ok <- function(mat, idx) rowSums(mat[, idx, drop=FALSE] > 1) >= 1
    for (fac in c("Age","Tissue","pair")){
        g_res <- get(paste0("gene_", tolower(fac)))
        t_res <- get(paste0("tx_",   tolower(fac)))
        for (con in names(t_res)){
            grps <- strsplit(con, "_vs_")[[1]]
            idx1 <- which(meta_sp[[fac]]==grps[1]); idx2 <- which(meta_sp[[fac]]==grps[2])
            dt_g <- g_res[[con]]
            ok_g <- half_ok(cpm_gene_mat, idx1) & half_ok(cpm_gene_mat, idx2)
            dt_g <- dt_g[ok_g]
            dt_g[, MeanCPM_A := rowMeans(cpm_gene_mat[match(gene_id, rownames(cpm_gene_mat)), idx1, drop=FALSE])]
            dt_g[, MeanCPM_B := rowMeans(cpm_gene_mat[match(gene_id, rownames(cpm_gene_mat)), idx2, drop=FALSE])]
            dt_tx <- t_res[[con]]
            ok_t <- half_ok(cpm_tx_mat, idx1) & half_ok(cpm_tx_mat, idx2)
            dt_tx <- dt_tx[ok_t]
            dt_tx[, MeanCPM_A := rowMeans(cpm_tx_mat[match(tx_id, rownames(cpm_tx_mat)), idx1, drop=FALSE])]
            dt_tx[, MeanCPM_B := rowMeans(cpm_tx_mat[match(tx_id, rownames(cpm_tx_mat)), idx2, drop=FALSE])]
            gene_tx_count <- table(tx2gene$gene_id)
            dt_tx <- dt_tx[gene_tx_count[gene_id] > 1]
            deg        <- dt_g[FDR<0.05 & abs(logFC)>=1]
            detx       <- dt_tx[FDR<0.05 & abs(logFC)>=1]
            dtx_nondeg <- detx[gene_id %in% dt_g[FDR>=0.05, gene_id]]

            outdir <- file.path(output_root, fac, sp)
            dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
            fwrite(deg,        file.path(outdir, paste0(con,"_DEG.tsv")),        sep="\t", quote=FALSE)
            fwrite(detx,       file.path(outdir, paste0(con,"_DETx.tsv")),       sep="\t", quote=FALSE)
            fwrite(dtx_nondeg, file.path(outdir, paste0(con,"_DTx_nonDEG.tsv")), sep="\t", quote=FALSE)
            message(sprintf("    %s %s | DEG:%d | DETx:%d | DTx_nonDEG:%d",
                            fac, con, nrow(deg), nrow(detx), nrow(dtx_nondeg)))
        }
    }

    # 保存关键对象以便溯源
    save(
        y_gene, y_tx,
        cpm_gene_mat, cpm_tx_mat,
        gene_age, tx_age,
        gene_tissue, tx_tissue,
        gene_pair, tx_pair,
        file = paste0(sp, "_debug_vars.RData")
    )
}

# 并行执行
bplapply(species_list, process_species,
         BPPARAM=MulticoreParam(workers=workers))

message("★ 全部物种处理完成 - 高阈值版 (共享 norm.factors + 持久化)")

