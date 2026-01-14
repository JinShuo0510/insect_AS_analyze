#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Gene-wise isoform diversity across tissues (CPM > 1 version)
# • 对每个 Tissue 统计 gene 的表达转录本数 IsoNum
# • 计算 IsoNum 在各组织间的方差 IsoVar
# • 输出长表 (gene × tissue) + 宽表 (gene × IsoVar)
# • 生成方差分布直方图与密度曲线
# Date: 2025-05-21
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(ggplot2)
    library(BiocParallel)
})

# ─── 参数 ────────────────────────────────────────────────────────────────────
species_list      <- c("Bombyx_mori")
counts_tx_dir     <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
gtf_dir           <- "/data/home/jinshuo/transcript/Genome"
meta_file         <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
out_root          <- "IsoformDiversity_CPM1"
libsize_thresh    <- 1e4
cpm_thresh        <- 1        # CPM > 1
workers           <- 4

# ─── 公共函数 ────────────────────────────────────────────────────────────────
build_tx2gene <- function(sp){
    f <- file.path(gtf_dir, sprintf("%s_tx2gene.tsv", sp))
    if (file.exists(f)) return(fread(f, col.names = c("tx_id","gene_id")))
    cmd  <- sprintf("grep -v '^#' %s", shQuote(file.path(gtf_dir, paste0(sp, ".gtf"))))
    gtf  <- fread(cmd = cmd, sep = "\t", select = c(3,9), col.names = c("feature","attr"))
    gtf  <- gtf[feature == "transcript"]
    ext  <- function(x,k) sub(paste0(".*",k," \"([^\"]+)\".*"), "\\1", x)
    tbl  <- unique(data.table(tx_id = ext(gtf$attr,"transcript_id"),
                              gene_id = ext(gtf$attr,"gene_id")))
    fwrite(tbl, f, sep = "\t", quote = FALSE, col.names = FALSE)
    tbl
}

get_libsize_bad <- function(mat, thresh){
    libs <- colSums(mat)
    bad0 <- names(libs)[libs == 0]
    badL <- names(libs)[libs < thresh]
    unique(c(bad0, badL))
}

# ─── 主流程 ────────────────────────────────────────────────────────────────
run_species <- function(sp){
    message(">>> 物种：", sp)
    dir.create(file.path(out_root, sp), recursive = TRUE, showWarnings = FALSE)

    # 载入 meta
    meta <- fread(meta_file)
    setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
    meta[, Age    := gsub("\\s+","_", trimws(Age))]
    meta[, Tissue := gsub("\\s+","_", trimws(Tissue))]
    setkey(meta, Sample)

    # 载入转录本计数
    cnt <- fread(file.path(counts_tx_dir, paste0(sp,"_combined_counts.txt")), data.table = FALSE)
    rownames(cnt) <- cnt[[1]]; cnt <- cnt[,-1]

    # 样本初筛
    sam0 <- intersect(colnames(cnt), meta$Sample)
    cnt  <- cnt[, sam0, drop = FALSE]

    # 去除库大小异常
    bad  <- get_libsize_bad(cnt, libsize_thresh)
    if(length(bad)){
        message("    剔除库大小不合格样本：", paste(bad, collapse = ", "))
        cnt <- cnt[, setdiff(colnames(cnt), bad), drop = FALSE]
    }
    meta <- meta[colnames(cnt)]
    stopifnot(identical(colnames(cnt), meta$Sample))  # 顺序锁定

    # TMM 与 CPM 计算
    dge <- calcNormFactors(DGEList(counts = cnt))
    cpm_mat <- cpm(dge)

    # 构建 tx2gene
    tx2gene <- build_tx2gene(sp)
    gene_of <- tx2gene$gene_id[match(rownames(cpm_mat), tx2gene$tx_id)]

    # —— Step 1. 逐组织统计 IsoNum ——
    tissues <- unique(meta$Tissue)
    iso_cnt_list <- lapply(tissues, function(tis){
        idx <- which(meta$Tissue == tis)
        present <- rowSums(cpm_mat[, idx, drop = FALSE] > cpm_thresh) > 0  # 至少一个样本表达
        data.table(gene_id = gene_of[present])[ , .(IsoNum = .N), by = gene_id][ , Tissue := tis]
    })
    iso_cnt_long <- rbindlist(iso_cnt_list, use.names = TRUE, fill = TRUE)
    fwrite(iso_cnt_long,
           file.path(out_root, sp, "gene_tissue_isoform_count.tsv"),
           sep = "\t", quote = FALSE)

    # —— Step 2. 计算组织间方差 ——
    iso_cnt_wide <- dcast(iso_cnt_long, gene_id ~ Tissue, value.var = "IsoNum", fill = 0)
    tissue_cols  <- setdiff(colnames(iso_cnt_wide), "gene_id")
    iso_cnt_wide[, IsoVar := apply(.SD, 1, var), .SDcols = tissue_cols]
    fwrite(iso_cnt_wide[, .(gene_id, IsoVar)],
           file.path(out_root, sp, "gene_isoform_variance.tsv"),
           sep = "\t", quote = FALSE)

    # —— Step 3. 方差分布可视化 ——
    p <- ggplot(iso_cnt_wide, aes(IsoVar)) +
        geom_histogram(bins = 50, alpha = 0.7) +
        geom_density(linewidth = 0.8) +
        labs(x = "Variance of expressed isoform number across tissues",
             y = "Gene count",
             title = sprintf("%s — Isoform diversity (CPM > 1)", sp)) +
        theme_bw()
    ggsave(file.path(out_root, sp, "IsoVar_distribution.png"),
           p, width = 6.5, height = 4.2)
    message("    完成：", sp,
            " | Genes: ", nrow(iso_cnt_wide),
            " | 图像已保存")
}

# 并行执行
bplapply(species_list, run_species,
         BPPARAM = MulticoreParam(workers = workers))

message("★ Isoform 多样性统计完成")

