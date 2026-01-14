##############################################################################
# 计算所有物种各基因跨组织 AS 计数方差，并汇总到一张密度曲线图
# —— 既不会因“缺失列”报错，也能自动识别每个 data.frame 中实际存在的组织列
##############################################################################

library(data.table)
library(ggplot2)

## —————————— 参数 ——————————
eps        <- 1e-5   # 避免 log10(0)
plot_bins  <- 60     # 若改直方图时的柱数
out_fig    <- "All_species_AS_variance_density.pdf"

## —————————— 核心函数：给定一个物种的列表 (组织→data.frame) 返回方差表 ——————————
calc_var_one_species <- function(sp_list) {
  # 1. 把同一物种所有组织 melt 成长表；仅保留非 NA 的事件（即存在于该组织）
  long_tbl <- rbindlist(
    lapply(sp_list, function(df) {
      dt  <- as.data.table(df)
      # 仅 melt 除 Event_ID 以外的真实列即可，na.rm=TRUE 自动去掉 NA
      melt(dt, id.vars = "Event_ID",
           variable.name = "Tissue", value.name = "val",
           na.rm = TRUE)[, GeneID := tstrsplit(Event_ID, "_", keep = 2)][,
                                                                         .N, by = .(GeneID, Tissue)]  # 每基因每组织 AS 计数
    }),
    use.names = TRUE
  )
  
  # 2. 宽表：GeneID × Tissue
  wide <- dcast(long_tbl, GeneID ~ Tissue, value.var = "N", fill = 0)
  tissue_cols <- setdiff(names(wide), "GeneID")
  
  # 3. 方差（只要 ≥2 个组织有值即可正常计算；否则为 0）
  wide[, Variance := ifelse(rowSums(.SD) > 0,
                            apply(.SD, 1, var), 0),
       .SDcols = tissue_cols]
  wide[, .(GeneID, Variance)]
}

## —————————— 主循环：遍历所有物种 ——————————
var_all_species <- rbindlist(
  lapply(names(all_species_tissue_data_withoutavg_1_8), function(sp) {
    message("Processing species: ", sp)
    res <- calc_var_one_species(all_species_tissue_data_withoutavg_1_8[[sp]])
    res[, Species := sp]
    res
  }),
  use.names = TRUE
)

## —————————— 绘图：多物种密度曲线 ——————————
p <- ggplot(var_all_species,
            aes(x = log10(Variance + eps), colour = Species)) +
  geom_density(linewidth = 1) +
  labs(title = "跨物种 AS 计数方差分布",
       x = "log10(Variance + ε)", y = "Density") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(out_fig, p, width = 7, height = 4, dpi = 300)
print(p)




#-----------------------------------------------
 
library(data.table)
library(ggplot2)

## —————— 参数设置 ——————
eps       <- 1e-5            # 避免 log10(0)
n_bins    <- 60              # 直方图柱数
out_file  <- "AS_variance_dual_axis.pdf"

## —————— 准备数据 ——————
df <- copy(var_all_species)   # 假设已包含 Species, GeneID, Variance
df[, logVar := log10(Variance + eps)]
df <- df[!is.na(logVar)]      # 去除缺失值

## —————— 计算缩放比率 ——————
# 直方图（不绘图）获取 counts
h <- hist(df$logVar, plot = FALSE, breaks = n_bins)
max_count   <- max(h$counts)

# 密度估算，去除 NA 后安全计算
dens <- density(df$logVar)
max_density <- max(dens$y)

ratio <- max_count / max_density

## —————— 绘图 —— 双 Y 轴 ——————
p <- ggplot(df, aes(x = logVar)) +
  # 左轴：Count
  geom_histogram(aes(y = ..count..),
                 bins = n_bins,
                 fill = "steelblue", alpha = 0.5) +
  # 右轴：Density
  geom_density(aes(y = ..density.. * ratio),
               color = "firebrick", size = 1) +
  # 双 Y 轴
  scale_y_continuous(
    name    = "基因数 (Count)",
    sec.axis = sec_axis(~ . / ratio,
                        name = "密度 (Density)")
  ) +
  labs(
    x     = "log10(Variance + ε)",
    title = "跨物种 AS 方差分布 —— Count 与 Density 双 Y 轴"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "steelblue"),
    axis.title.y.right = element_text(color = "firebrick")
  )

ggsave(out_file, p, width = 7, height = 4, dpi = 300)
print(p)
