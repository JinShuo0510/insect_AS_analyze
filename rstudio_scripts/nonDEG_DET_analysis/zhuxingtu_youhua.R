###############################################################################
# 网格堆积柱形图：在“Stage1 × Stage2 / Tissue1 × Tissue2”网格中，
# 每个单元用统一宽度的纵向堆积柱形条替代气泡，直观呈现 DEG / nonDEG 数量。
# ────────────────────────────────────────────────────────────────────────────
# 输出文件：
#   • p_age_gridbar.pdf
#   • p_tissue_gridbar.pdf
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ─────────────────────────── 1. 公共设置 ─────────────────────────────
my_colors <- c(
  "DTx_nonDEG" = "#2A5783",
  "DETx_DEG"   = "#7FB3D5"
)

# ─────────────────────────── 2. 读入或直接复制数据 ───────────────────
## 2.1 发育阶段
df_age <- tibble::tribble(
  ~comparison,                  ~DETx_DEG, ~DTx_nonDEG,
  "Adult_vs_Cell",                2628,       386,
  "Adult_vs_Egg",                 7041,       648,
  "Adult_vs_Larva",               3587,       276,
  "Adult_vs_Larva_and_Pupa",       158,        44,
  "Adult_vs_Pupa",                3172,       329,
  "Cell_vs_Egg",                  2536,      1281,
  "Cell_vs_Larva",                2415,       563,
  "Cell_vs_Larva_and_Pupa",        342,        48,
  "Cell_vs_Pupa",                 2822,       546,
  "Egg_vs_Larva",                 8452,       877,
  "Egg_vs_Larva_and_Pupa",         261,        38,
  "Egg_vs_Pupa",                  4727,       311,
  "Larva_vs_Larva_and_Pupa",       219,        26,
  "Larva_vs_Pupa",                3058,       303,
  "Larva_and_Pupa_vs_Pupa",        129,        40
)

## 2.2 组织
df_tissue <- tibble::tribble(
  ~comparison,                         ~DETx_DEG, ~DTx_nonDEG,
  "Developmental_tissue_vs_Digestive",  6851,       675,
  "Developmental_tissue_vs_Fat_body",   5730,      1000,
  "Developmental_tissue_vs_Gland",      8260,      1522,
  "Developmental_tissue_vs_Head",       3374,       371,
  "Developmental_tissue_vs_Ovary",      3779,       305,
  "Developmental_tissue_vs_Testis",     8253,       471,
  "Developmental_tissue_vs_Thorax",     1753,       346,
  "Developmental_tissue_vs_Whole_body", 3155,       582,
  "Digestive_vs_Fat_body",              5030,       160,
  "Digestive_vs_Gland",                 4995,       158,
  "Digestive_vs_Head",                  6284,       287,
  "Digestive_vs_Ovary",                 8501,       951,
  "Digestive_vs_Testis",                8359,       236,
  "Digestive_vs_Thorax",                3961,       537,
  "Digestive_vs_Whole_body",            4897,       244,
  "Fat_body_vs_Gland",                  3702,       112,
  "Fat_body_vs_Head",                   5245,       416,
  "Fat_body_vs_Ovary",                  7501,      1540,
  "Fat_body_vs_Testis",                 7731,       321,
  "Fat_body_vs_Thorax",                 3849,       964,
  "Fat_body_vs_Whole_body",             4382,       491,
  "Gland_vs_Head",                      6736,       754,
  "Gland_vs_Ovary",                     8227,      1608,
  "Gland_vs_Testis",                    7977,       340,
  "Gland_vs_Thorax",                    4350,      1754,
  "Gland_vs_Whole_body",                5203,       608,
  "Head_vs_Ovary",                      4216,       392,
  "Head_vs_Testis",                     7381,       373,
  "Head_vs_Thorax",                     1580,       508,
  "Head_vs_Whole_body",                 3316,       338,
  "Ovary_vs_Testis",                    6716,       669,
  "Ovary_vs_Thorax",                    1831,       551,
  "Ovary_vs_Whole_body",                5307,       669,
  "Testis_vs_Thorax",                   6289,       651,
  "Testis_vs_Whole_body",               7829,       406,
  "Thorax_vs_Whole_body",               2414,       297
)

# ─────────────────────────── 3. 长格式转换 ───────────────────────────
df_age_long <- df_age %>% 
  pivot_longer(cols = c(DETx_DEG, DTx_nonDEG),
               names_to = "type", values_to = "count") %>% 
  mutate(type = factor(type, levels = c("DTx_nonDEG", "DETx_DEG")))

df_tissue_long <- df_tissue %>% 
  pivot_longer(cols = c(DETx_DEG, DTx_nonDEG),
               names_to = "type", values_to = "count") %>% 
  mutate(type = factor(type, levels = c("DTx_nonDEG", "DETx_DEG")))

# ───────────────────── 4. 网格堆积柱形图函数 ─────────────────────────
plot_gridbar <- function(
    df_long, sep_regex, x_lab, y_lab, file_out,
    width_out = 8, height_out = 6,
    bar_width  = 0.1,                 # ① 更瘦的柱
    expand_mult = 0.9,                 # ② 左右留白比例
    panel_gap = 0.8                    # ③ 面板间距，单位行
) {
  df_sep <- df_long %>% 
    tidyr::separate(comparison, c("Cat1", "Cat2"), sep = sep_regex)
  
  p <- ggplot(df_sep, aes(x = 1, y = count, fill = type)) +
    geom_col(width = bar_width) +
    scale_fill_manual(values = my_colors,
                      labels = c("nonDEG", "DEG")) +
    
    # 每格 1 根柱，行列网格
    facet_grid(Cat2 ~ Cat1, switch = "y") +
    
    # 给 x 轴左右各留 expand_mult × 总宽
    scale_x_continuous(expand = expansion(mult = expand_mult)) +
    scale_y_continuous(expand = c(0, 0)) +
    
    labs(x = x_lab, y = y_lab) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "top",
      strip.text.x = element_text(angle = 45, hjust = 0),
      strip.text.y.left = element_text(angle = 0),
      strip.background = element_rect(fill = NA, colour = NA),
      panel.spacing = unit(panel_gap, "lines")  # 面板间距
    )
  
  ggsave(file_out, p, width = width_out, height = height_out, dpi = 300)
  p
}

## 调用示例（发育阶段 / 组织）
p_age_gridbar <- plot_gridbar(df_age_long,    "_vs_", "阶段 A", "阶段 B",
                              "p_age_gridbar.pdf",
                              width_out = 8, height_out = 6,
                              bar_width = 0.35, expand_mult = 0.4)

p_tissue_gridbar <- plot_gridbar(df_tissue_long, "_vs_", "组织 A", "组织 B",
                                 "p_tissue_gridbar.pdf",
                                 width_out = 8, height_out = 8,
                                 bar_width = 0.35, expand_mult = 0.4)

print(p_age_gridbar)
print(p_tissue_gridbar)
