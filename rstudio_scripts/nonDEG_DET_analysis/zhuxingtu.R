# 加载必要的 R 包
library(dplyr)
library(tidyr)
library(ggplot2)

# 自定义配色
my_colors <- c(
  "DTx_nonDEG" = "#2A5783",
  "DETx_DEG"   = "#7FB3D5"
)

######################
# 发育阶段数据准备
######################
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
) %>%
  pivot_longer(
    cols = c(DETx_DEG, DTx_nonDEG),
    names_to  = "type",
    values_to = "count"
  ) %>%
  mutate(type = factor(type, levels = c("DTx_nonDEG", "DETx_DEG")))

######################
# 组织数据准备
######################
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
) %>%
  pivot_longer(
    cols = c(DETx_DEG, DTx_nonDEG),
    names_to  = "type",
    values_to = "count"
  ) %>%
  mutate(type = factor(type, levels = c("DTx_nonDEG", "DETx_DEG")))

######################
# 气泡图绘制
######################

## 发育阶段气泡图
df_age_bubble <- df_age %>%
  separate(comparison, into = c("Stage1", "Stage2"), sep = "_vs_")

p_age_bubble <- ggplot(df_age_bubble, aes(x = Stage1, y = Stage2)) +
  geom_point(aes(size = count, color = type), alpha = 0.7) +
  scale_size_continuous(range = c(3, 12), name = "DETx 数量") +
  scale_color_manual(values = my_colors, name = NULL, labels = c("nonDEG", "DEG")) +
  labs(x = "阶段 A", y = "阶段 B") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

## 组织气泡图
df_tissue_bubble <- df_tissue %>%
  separate(comparison, into = c("Tissue1", "Tissue2"), sep = "_vs_")

p_tissue_bubble <- ggplot(df_tissue_bubble, aes(x = Tissue1, y = Tissue2)) +
  geom_point(aes(size = count, color = type), alpha = 0.7) +
  scale_size_continuous(range = c(3, 12), name = "DETx 数量") +
  scale_color_manual(values = my_colors, name = NULL, labels = c("nonDEG", "DEG")) +
  labs(x = "组织 A", y = "组织 B") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# 绘图展示及保存
print(p_age_bubble)
print(p_tissue_bubble)
ggsave("p_age_bubble.pdf",    p_age_bubble,    width = 8, height = 6, dpi = 300)
ggsave("p_tissue_bubble.pdf", p_tissue_bubble, width = 8, height = 6, dpi = 300)
