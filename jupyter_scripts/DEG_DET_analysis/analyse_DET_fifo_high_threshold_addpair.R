#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# 并行 + 向量化 DTU 流程  - 高阈值版
# - 输出 DEG、DETx、DTx_nonDEG 三类结果
# - 额外使用更严格的转录本过滤：CPM > 1 且至少在 3 个样本中表达
# - 新增 Age–Tissue 组合（pair）比对
# Date: 2025-05-12
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(limma)
    library(BiocParallel)
})

# ─── 参数 ────────────────────────────────────────────────────────────────────
species_list      <- c(
    "Acyrthosiphon_pisum","Aedes_aegypti","Apis_mellifera",
    "Blattella_germanica","Bombyx_mori","Drosophila_mojavensis",
    "Gryllus_bimaculatus","Helicoverpa_armigera","Tribolium_castaneum"
)
gtf_dir           <- "/data/home/jinshuo/transcript/Genome"
counts_gene_dir   <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/gene"
counts_tx_dir     <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
meta_file         <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
output_root       <- "edgeR_DTU_voom_highThresh_addpair"
libsize_thresh    <- 1e6
workers           <- 4

# 转录本过滤阈值
cpm_thresh        <- 1          # CPM > 1
min_samples_expr  <- 3          # 至少 3 个样本

# ─── 样本注释 ────────────────────────────────────────────────────────────────
meta <- fread(meta_file)
setnames(meta, c("Run","Age","Tissue.1"), c("Sample","Age","Tissue"))
meta[, Age    := gsub("\\s+","_", trimws(Age))]
meta[, Tissue := gsub("\\s+","_", trimws(Tissue))]

# ─── tx2gene 生成 ────────────────────────────────────────────────────────────
build_tx2gene <- function(sp) {
    gtf_path     <- file.path(gtf_dir, paste0(sp, ".gtf"))
    tx2gene_file <- paste0(sp, "_tx2gene.tsv")
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
    gtf <- gtf[grepl("transcript_id", attr, fixed = TRUE)]
    extract_field <- function(x, key)
        sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
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

# ─── voom 分析 ───────────────────────────────────────────────────────────────
analyse_voom <- function(y_dge, meta_sp, fac, tx2gene = NULL) {
    orig_lv <- levels(factor(meta_sp[[fac]]))
    syn_lv  <- make.names(orig_lv)
    grp     <- factor(meta_sp[[fac]], levels = orig_lv)
    design  <- model.matrix(~0 + grp)
    colnames(design) <- syn_lv

    v   <- voom(y_dge, design, plot = FALSE)
    fit <- lmFit(v, design)

    combos <- combn(seq_along(syn_lv), 2, simplify = FALSE)
    contrast_names <- sapply(combos, function(idx)
        paste(orig_lv[idx], collapse = "_vs_"))
    contrast_expr  <- sapply(combos, function(idx)
        paste0(syn_lv[idx[1]], "-", syn_lv[idx[2]]))
    contMat <- makeContrasts(contrasts = contrast_expr, levels = design)

    fit2 <- eBayes(contrasts.fit(fit, contMat))

    res <- lapply(seq_along(contrast_names), function(i) {
        tt <- topTable(fit2, coef = i, number = Inf,
                       sort.by = "none", adjust.method = "BH")
        dt <- as.data.table(tt, keep.rownames = "id")
        if (is.null(tx2gene)) {
            setnames(dt, "id", "gene_id")
        } else {
            setnames(dt, "id", "tx_id")
            dt[, gene_id := tx2gene$gene_id[match(tx_id, tx2gene$tx_id)]]
        }
        dt[, FDR := adj.P.Val]
        dt
    })
    names(res) <- contrast_names
    res
}

# ─── 物种主流程 ──────────────────────────────────────────────────────────────
process_species <- function(sp) {
    message(">>> 处理物种：", sp)

    gfile <- file.path(counts_gene_dir, paste0(sp, "_combined_counts.txt"))
    tfile <- file.path(counts_tx_dir,   paste0(sp, "_combined_counts.txt"))

    g_cnt <- fread(gfile, data.table = FALSE); rownames(g_cnt) <- g_cnt[[1]]; g_cnt <- g_cnt[,-1]
    t_cnt <- fread(tfile, data.table = FALSE); rownames(t_cnt) <- t_cnt[[1]]; t_cnt <- t_cnt[,-1]

    sam0  <- intersect(colnames(g_cnt), meta$Sample)
    g_cnt <- g_cnt[, sam0, drop = FALSE]
    t_cnt <- t_cnt[, sam0, drop = FALSE]

    lib_sizes <- colSums(g_cnt)
    keep <- setdiff(sam0,
                    c(names(lib_sizes)[lib_sizes == 0],
                      names(lib_sizes)[lib_sizes < libsize_thresh]))
    g_cnt   <- g_cnt[, keep, drop = FALSE]
    t_cnt   <- t_cnt[, keep, drop = FALSE]
    meta_sp <- meta[Sample %in% keep]

    # ── 新增：生成 Age–Tissue 组合因子 ────────────────────────────────────────
    meta_sp[, pair := paste(Age, Tissue, sep = "_")]

    # ── 基因层面 ────────────────────────────────────────────────────────────
    y_gene <- DGEList(counts = g_cnt)
    y_gene <- y_gene[filterByExpr(y_gene, group = meta_sp$Age),
                     , keep.lib.sizes = FALSE]
    y_gene <- calcNormFactors(y_gene)

    # ── 转录本层面（更严格过滤）──────────────────────────────────────────────
    tx2gene <- build_tx2gene(sp)
    y_tx <- DGEList(counts = t_cnt,
                    genes  = data.frame(gene_id = tx2gene$gene_id[match(rownames(t_cnt),
                                                                         tx2gene$tx_id)]))
    cpm_tx <- cpm(y_tx)
    keep_tx <- rowSums(cpm_tx > cpm_thresh) >= min_samples_expr
    y_tx <- y_tx[keep_tx, , keep.lib.sizes = FALSE]
    y_tx <- y_tx[filterByExpr(y_tx, group = meta_sp$Age),
                 , keep.lib.sizes = FALSE]
    y_tx <- calcNormFactors(y_tx)

    # ── 差异分析 ────────────────────────────────────────────────────────────
    gene_age    <- analyse_voom(y_gene, meta_sp, "Age")
    tx_age      <- analyse_voom(y_tx,   meta_sp, "Age",    tx2gene)
    gene_tissue <- analyse_voom(y_gene, meta_sp, "Tissue")
    tx_tissue   <- analyse_voom(y_tx,   meta_sp, "Tissue", tx2gene)
    # ── 新增：pair 组合的差异分析 ───────────────────────────────────────────
    gene_pair   <- analyse_voom(y_gene, meta_sp, "pair")
    tx_pair     <- analyse_voom(y_tx,   meta_sp, "pair",   tx2gene)

    # ── 输出所有三种因子的结果 ─────────────────────────────────────────────
    for (fac in c("Age", "Tissue", "pair")) {
        g_res <- get(paste0("gene_", tolower(fac)))
        t_res <- get(paste0("tx_",   tolower(fac)))

        for (con in names(t_res)) {
            dt_g  <- g_res[[con]]
            dt_tx <- t_res[[con]]

            deg   <- dt_g[FDR < 0.05 & abs(logFC) >= 1]
            detx  <- dt_tx[FDR < 0.05 & abs(logFC) >= 1]

            nonDEG     <- dt_g[FDR >= 0.05, gene_id]
            dtx_nondeg <- detx[gene_id %in% nonDEG]

            outdir <- file.path(output_root, fac, sp)
            dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

            fwrite(deg,
                   file = file.path(outdir, paste0(con, "_DEG.tsv")),
                   sep  = "\t", quote = FALSE, row.names = FALSE)

            fwrite(detx,
                   file = file.path(outdir, paste0(con, "_DETx.tsv")),
                   sep  = "\t", quote = FALSE, row.names = FALSE)

            fwrite(dtx_nondeg,
                   file = file.path(outdir, paste0(con, "_DTx_nonDEG.tsv")),
                   sep  = "\t", quote = FALSE, row.names = FALSE)

            message(sprintf("    %s %s | DEG: %d | DETx: %d | DTx_nonDEG: %d",
                            fac, con, nrow(deg), nrow(detx), nrow(dtx_nondeg)))
        }
    }
}

# ─── 并行执行 ────────────────────────────────────────────────────────────────
bplapply(species_list, process_species,
         BPPARAM = MulticoreParam(workers = workers))

message("★ 全部物种处理完成 - 高阈值版")

