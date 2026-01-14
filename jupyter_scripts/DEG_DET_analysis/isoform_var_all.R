#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Isoform diversity across all species (meta-driven) – 带日志输出修正版
# Date: 2025-05-21
# ─────────────────────────────────────────────────────────────────────────────
# 依赖 ------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(data.table)
    library(edgeR)
    library(ggplot2)
    library(BiocParallel)
})

## ─── 路径参数 ────────────────────────────────────────────────────────────────
counts_tx_dir <- "/data/home/jinshuo/transcript/Featurecounts_combined_counts/tx"
gtf_dir       <- "/data/home/jinshuo/transcript/Genome"
meta_file     <- "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv"
out_root      <- "IsoformDiversity_CPM1_allSpecies"

## ─── 全局阈值 ────────────────────────────────────────────────────────────────
libsize_thresh <- 1e4
cpm_thresh     <- 1
workers        <- 6
eps            <- 1e-6

dir.create(out_root, showWarnings = FALSE)

## ─── 公共函数 ────────────────────────────────────────────────────────────────
# 生成或读取 tx2gene 对照表 ----------------------------------------------------
build_tx2gene <- function(sp){
    f_out <- file.path(gtf_dir, sprintf("%s_tx2gene.tsv", sp))
    if (file.exists(f_out))
        return(fread(f_out, col.names = c("tx_id","gene_id")))

    gtf_f <- file.path(gtf_dir, sprintf("%s.gtf", sp))
    if (!file.exists(gtf_f)){
        message(sprintf("[WARN] 缺失 GTF: %s", gtf_f))
        return(NULL)
    }
    cmd  <- sprintf("grep -v '^#' %s", shQuote(gtf_f))
    gtf  <- fread(cmd = cmd, sep = "\t", select = c(3,9),
                  col.names = c("feature","attr"))
    ext  <- function(x,k) sub(paste0(".*",k," \"([^\"]+)\".*"), "\\1", x)
    tbl  <- unique(data.table(
        tx_id   = ext(gtf$attr,"transcript_id"),
        gene_id = ext(gtf$attr,"gene_id")
    ))
    fwrite(tbl, f_out, sep = "\t", quote = FALSE, col.names = FALSE)
    tbl
}

# 统计零或低文库样本 ------------------------------------------------------------
count_bad <- function(mat, thresh){
    lib <- colSums(mat)
    unique(c(names(lib)[lib == 0], names(lib)[lib < thresh]))
}

# 核心计算：单物种 isoform 方差 ---------------------------------------------------
isoform_stats <- function(sp, meta){
    start_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] ▶ 开始处理物种: %s", start_ts, sp))

    res <- tryCatch({
        cnt_f <- file.path(counts_tx_dir, sprintf("%s_combined_counts.txt", sp))
        if (!file.exists(cnt_f)){
            message(sprintf("[WARN] 缺失 counts 文件：%s (跳过)", cnt_f))
            return(NULL)
        }
        cnt <- fread(cnt_f, data.table = FALSE)
        rownames(cnt) <- cnt[[1]]; cnt <- cnt[,-1]

        ## 样本过滤 ------------------------------------------------------------
        sam0 <- intersect(colnames(cnt), meta$Sample)
        cnt  <- cnt[, sam0, drop = FALSE]
        bad  <- count_bad(cnt, libsize_thresh)
        if (length(bad)) cnt <- cnt[, setdiff(colnames(cnt), bad), drop = FALSE]
        if (ncol(cnt) == 0){
            message(sprintf("[WARN] %s: 无可用样本，跳过", sp))
            return(NULL)
        }

        ## 精准对齐 meta ------------------------------------------------------
        meta_sp <- meta[Sample %in% colnames(cnt)]
        meta_sp <- meta_sp[match(colnames(cnt), meta_sp$Sample)]
        stopifnot(identical(meta_sp$Sample, colnames(cnt)))

        ## TMM + CPM -----------------------------------------------------------
        dge     <- calcNormFactors(DGEList(counts = cnt))
        cpm_mat <- cpm(dge)

        ## tx → gene ----------------------------------------------------------
        tx2gene <- build_tx2gene(sp)
        if (is.null(tx2gene)) return(NULL)
        gene_of <- tx2gene$gene_id[match(rownames(cpm_mat), tx2gene$tx_id)]

        ## 按组织统计 IsoNum ---------------------------------------------------
        tissues <- unique(meta_sp$Tissue)
        iso_cnt <- rbindlist(lapply(tissues, function(tis){
            idx <- which(meta_sp$Tissue == tis)
            expressed <- rowSums(cpm_mat[, idx, drop = FALSE] > cpm_thresh) > 0
            data.table(gene_id = gene_of[expressed])[
                , .(IsoNum = .N), by = gene_id
            ][, Tissue := tis]
        }))

        ## 计算方差并标注物种 ---------------------------------------------------
        wide <- dcast(iso_cnt, gene_id ~ Tissue, value.var = "IsoNum", fill = 0)
        wide[, IsoVar := apply(.SD, 1, var), .SDcols = setdiff(names(wide), "gene_id")]
        data.table(gene_id = wide$gene_id,
                   IsoVar  = wide$IsoVar,
                   Species = sp)
    }, error = function(e){
        message(sprintf("[ERROR] %s 处理失败: %s", sp, conditionMessage(e)))
        NULL
    })

    end_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    gene_n <- ifelse(is.null(res), 0, nrow(res))
    message(sprintf("[%s] ✔ 完成物种: %s (genes=%d)", end_ts, sp, gene_n))
    res
}

## ─── 主执行 ─────────────────────────────────────────────────────────────────
meta_all <- fread(meta_file)
setnames(meta_all, c("Run","Tissue.1"), c("Sample","Tissue"))
meta_all[, Tissue := gsub("\\s+","_", trimws(Tissue))]

species_list <- unique(meta_all$ScientificName)
message(sprintf("★ 总共需要处理物种数: %d", length(species_list)))

# 并行统计 --------------------------------------------------------------------
res_list <- bplapply(species_list, function(sp){
    isoform_stats(sp, meta_all[ScientificName == sp])
}, BPPARAM = MulticoreParam(workers = workers))

# 合并有效结果 ----------------------------------------------------------------
var_all_species <- rbindlist(Filter(Negate(is.null), res_list))

fwrite(var_all_species,
       file.path(out_root, "gene_isoform_variance_allSpecies.tsv"),
       sep = "\t", quote = FALSE)
saveRDS(var_all_species,"var_all_species.RDS")
## ———— 绘图 1：多物种密度曲线 ————————————
out_fig1 <- file.path(out_root, "IsoVar_density_allSpecies.pdf")
message("[PLOT] 绘制密度曲线 …")

p1 <- ggplot(var_all_species,
             aes(x = log10(IsoVar + eps), colour = Species)) +
  geom_density(linewidth = 1) +
  labs(title = "跨物种转录本多样性方差分布",
       x = "log10(Variance + ε)", y = "Density") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(out_fig1, p1, width = 7, height = 4, dpi = 300)

## ———— 绘图 2：Count + Density 双 Y 轴 —————————
message("[PLOT] 绘制柱形/密度双轴图 …")

df_plot <- var_all_species[, .(logVar = log10(IsoVar + eps))]

## 关键修正：剔除缺失或非有限值
df_plot <- df_plot[is.finite(logVar)]

n_bins      <- 60
density_max <- max(density(df_plot$logVar)$y)
hist_max    <- max(hist(df_plot$logVar, plot = FALSE, breaks = n_bins)$counts)
ratio       <- hist_max / density_max

out_fig2 <- file.path(out_root, "IsoVar_hist_density_allSpecies.pdf")
p2 <- ggplot(df_plot, aes(x = logVar)) +
  geom_histogram(aes(y = ..count..), bins = n_bins,
                 fill = "steelblue", alpha = 0.5) +
  geom_density(aes(y = ..density.. * ratio),
               color = "firebrick", linewidth = 1) +
  scale_y_continuous(
    name     = "基因数 (Count)",
    sec.axis = sec_axis(~ . / ratio, name = "密度 (Density)")
  ) +
  labs(x = "log10(Variance + ε)",
       title = "跨物种转录本多样性方差：Count 与 Density") +
  theme_bw(base_size = 12) +
  theme(axis.title.y.left  = element_text(color = "steelblue"),
        axis.title.y.right = element_text(color = "firebrick"))
ggsave(out_fig2, p2, width = 7, height = 4, dpi = 300)


message(sprintf("★ 全流程完成: 结果与图像已保存至 %s", out_root))

