#!/usr/bin/env Rscript
# ----------------------------------------
# 基因非显著 / 转录本显著 检索脚本（含样本库大小过滤）
# Date: 2025-05-06
# ----------------------------------------

suppressPackageStartupMessages({
    library(edgeR)      # ≥4.x
    library(data.table) # 高效 I/O
})

# ----------------------------------------
# 0. 工具函数
# ----------------------------------------

## 0.1 从 GTF 生成 tx2gene（纯 data.table 解析）
build_tx2gene <- function(sp,
                          gtf_dir = "/data/home/jinshuo/transcript/Genome",
                          out_dir = ".") {
    gtf_path     <- file.path(gtf_dir, paste0(sp, ".gtf"))
    tx2gene_file <- file.path(out_dir, paste0(sp, "_tx2gene.tsv"))
    if (file.exists(tx2gene_file)) {
        return(fread(tx2gene_file,
                     col.names = c("tx_id","gene_id"),
                     data.table = FALSE))
    }
    message("  • 构建 tx2gene：", sp)
    cmd <- paste("grep -v '^#'", shQuote(gtf_path))
    gtf <- fread(cmd = cmd,
                 sep = "\t", header = FALSE,
                 select = c(3L,9L),
                 col.names = c("feature","attr"))
    gtf <- gtf[grepl("transcript_id", attr, fixed=TRUE)]
    extract_field <- function(x, key) sub(paste0('.*',key,' "([^"]+)".*'), "\\1", x)
    tx2gene <- unique(data.frame(
        tx_id   = extract_field(gtf$attr, "transcript_id"),
        gene_id = extract_field(gtf$attr, "gene_id"),
        stringsAsFactors = FALSE
    ))
    fwrite(tx2gene,
           file      = tx2gene_file,
           sep       = "\t",
           quote     = FALSE,
           col.names = FALSE)
    tx2gene
}

## 0.2 生成所有两两对比名
pair_grid <- function(v) {
    cmb <- t(combn(unique(v), 2))
    apply(cmb, 1, function(x) paste(x, collapse="_vs_"))
}

## 0.3 差异显著过滤：FDR < 0.05 & |logFC| ≥ 1
sig_filter <- function(tab, fdr=0.05, lfc=1) {
    tab[tab$FDR < fdr & abs(tab$logFC) >= lfc, ]
}

# ----------------------------------------
# 1. 读取样本注释
# ----------------------------------------
meta <- fread("2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv")
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, Age    := gsub("\\s+","_", trimws(Age))]
meta[, Tissue := gsub("\\s+","_", trimws(Tissue))]

# ----------------------------------------
# 2. 主函数：对每个物种跑 Age 和 Tissue
# ----------------------------------------
run_species <- function(sp,
                        counts_gene_dir = "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene",
                        counts_tx_dir   = "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx",
                        output_root     = "edgeR_DTU",
                        libsize_thresh  = 1e6) {

    message(">>> 分析物种：", sp)

    # 2.1 读取 counts
    gene_file <- file.path(counts_gene_dir, paste0(sp,"_combined_counts.txt"))
    tx_file   <- file.path(counts_tx_dir,   paste0(sp,"_combined_counts.txt"))
    g_cnt     <- fread(gene_file, data.table=FALSE)
    rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[,-1]
    tx_cnt    <- fread(tx_file, data.table=FALSE)
    rownames(tx_cnt) <- tx_cnt[[1]]; tx_cnt <- tx_cnt[,-1]

    # 2.2 初始样本同步
    sam0      <- intersect(colnames(g_cnt), meta$Sample)
    g_cnt     <- g_cnt[, sam0, drop=FALSE]
    tx_cnt    <- tx_cnt[, sam0, drop=FALSE]

    # 2.3 库大小过滤（基于 gene-level counts）
    lib_sizes <- colSums(g_cnt)
    # 去除库大小为 0 的样本
    zero_lib <- names(lib_sizes)[lib_sizes == 0]
    if (length(zero_lib) > 0) {
        message("  • 移除零库样本：", paste(zero_lib, collapse=", "))
        keep <- setdiff(colnames(g_cnt), zero_lib)
        g_cnt  <- g_cnt[, keep, drop=FALSE]
        tx_cnt <- tx_cnt[, keep, drop=FALSE]
    }
    # 去除低库样本
    low_lib <- names(lib_sizes)[lib_sizes < libsize_thresh]
    if (length(low_lib) > 0) {
        message("  • 移除低库样本 (<", libsize_thresh,")：", paste(low_lib, collapse=", "))
        keep <- setdiff(colnames(g_cnt), low_lib)
        g_cnt  <- g_cnt[, keep, drop=FALSE]
        tx_cnt <- tx_cnt[, keep, drop=FALSE]
    }

    # 2.4 更新 meta_sp
    meta_sp <- meta[Sample %in% colnames(g_cnt)]

    # 2.5 构建 DGEList：gene-level
    y_gene <- DGEList(counts = g_cnt)
    keep_g <- filterByExpr(y_gene, group = meta_sp$Age)
    y_gene <- y_gene[keep_g, , keep.lib.sizes=FALSE]
    y_gene <- calcNormFactors(y_gene)

    # 2.6 构建 DGEList：tx-level
    tx2gene  <- build_tx2gene(sp)
    gene_map <- tx2gene$gene_id[match(rownames(tx_cnt), tx2gene$tx_id)]
    y_tx     <- DGEList(counts = tx_cnt, genes = data.frame(gene_id=gene_map))
    keep_tx  <- filterByExpr(y_tx, group = meta_sp$Age)
    y_tx     <- y_tx[keep_tx, , keep.lib.sizes=FALSE]
    y_tx     <- calcNormFactors(y_tx)

    # 内部函数：针对因子跑对比
    analyse_factor <- function(fac = c("Age","Tissue")) {
        fac      <- match.arg(fac)
        orig_lv  <- levels(factor(meta_sp[[fac]]))
        syn_lv   <- make.names(orig_lv)
        grp      <- factor(meta_sp[[fac]], levels = orig_lv)
        design   <- model.matrix(~0 + grp)
        colnames(design) <- syn_lv

        y_g  <- estimateDisp(y_gene, design); fit_g <- glmQLFit(y_g, design)
        y_t  <- estimateDisp(y_tx, design);   fit_t <- glmQLFit(y_t, design)

        for (p in pair_grid(orig_lv)) {
            parts <- strsplit(p, "_vs_")[[1]]
            s1 <- syn_lv[orig_lv == parts[1]]
            s2 <- syn_lv[orig_lv == parts[2]]
            contr <- makeContrasts(contrasts=paste0(s1,"-",s2), levels=design)

            # 基因
            res_g <- glmQLFTest(fit_g, contrast=contr)
            tab_g <- topTags(res_g, n=Inf)$table

            # 转录本
            res_t <- glmQLFTest(fit_t, contrast=contr)
            tab_t <- topTags(res_t, n=Inf)$table
            tab_t$gene_id <- y_t$genes$gene_id[match(rownames(tab_t), rownames(y_t))]

            # 筛：转录本显著 & 基因非显著
            sig_t  <- sig_filter(tab_t)
            nonDEG <- rownames(tab_g)[tab_g$FDR >= 0.05]
            final  <- sig_t[sig_t$gene_id %in% nonDEG, ]

            outdir <- file.path(output_root, fac, sp)
            dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
            fwrite(final,
                   file = file.path(outdir, paste0(p,"_DTx_nonDEG.tsv")),
                   sep = "\t", quote = FALSE, row.names = TRUE)

            message("    ", fac, " ", p, " => ", nrow(final), " 转录本")
        }
    }

    analyse_factor("Age")
    analyse_factor("Tissue")
}


# ----------------------------------------
# 3. 批量执行
# ----------------------------------------
#species_list <- c("Acyrthosiphon_pisum","Aedes_aegypti","Apis_mellifera",
#                  "Blattella_germanica","Bombyx_mori","Drosophila_mojavensis",
#                  "Gryllus_bimaculatus","Helicoverpa_armigera","Tribolium_castaneum")

species_list <- c("Apis_mellifera")
invisible(lapply(species_list, run_species))

message("★ 所有物种分析完成")

