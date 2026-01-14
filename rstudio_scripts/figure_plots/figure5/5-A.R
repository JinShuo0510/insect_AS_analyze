library(tidyr)
library(dplyr)
library(ggplot2)

df_sum <- final_df %>%
  mutate(Species = sapply(Species, create_abbreviation)) %>%
  group_by(Species, Tissue) %>%
  summarise(
    Detected_AS = sum(Event_count, na.rm = TRUE),
    DevAS       = sum(devAS_count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols      = c("Detected_AS", "DevAS"),
    names_to  = "Category",
    values_to = "Count"
  )

df_sum$Species <- factor(df_sum$Species,
                         levels = sort(unique(df_sum$Species)))

df_sum <- df_sum %>%
  mutate(Tissue = fct_recode(Tissue, "Wing" = "wing"))
# tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
df_sum$Tissue  <- factor(df_sum$Tissue,
                         levels = names(tissue_color_mapping))


p <- ggplot(df_sum,
            aes(x = Species,
                y = Count,
                fill = Tissue,
                alpha = Category)) +
  geom_col(width = 0.7) +
  facet_grid(~ Tissue,
             scales = "free_x",
             space  = "free_x",
             switch = "x") +
  scale_fill_manual(
    values = tissue_color_mapping,
    name   = "Tissue"                # 设置图例标题
  ) + 
  scale_alpha_manual(
    values = c("Detected_AS" = 0.4,
               "DevAS"       = 1),
    name   = "Category"
  ) +
  # —— 将 base_size 调大到 18 —— 
  theme_bw(base_size = 18) +
  theme(
    # 清空面板与分面背景
    panel.background    = element_blank(),
    panel.border        = element_blank(),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    strip.background    = element_blank(),
    strip.placement     = "outside",
    strip.text.x.top    = element_blank(),
    # 底部分面标签
    strip.text.x.bottom = element_text(
      size = 16,    # 分面标签字体
      vjust = 1
    ),
    # 仅保留 x、y 轴线
    axis.line.x         = element_line(),
    axis.line.y         = element_line(),
    axis.ticks          = element_line(),
    # 坐标轴标题
    axis.title.x        = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y        = element_text(size = 20, margin = margin(r = 10)),
    # 坐标轴刻度标签
    axis.text.x         = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 14
    ),
    axis.text.y         = element_text(size = 14),
    # 图例文字（仅 alpha 图例）
    legend.title        = element_text(size = 16),
    legend.text         = element_text(size = 14),
    panel.spacing       = unit(0.2, "lines")
  ) +
  labs(
    x = "Species",
    y = "Number of Microexons"
  )

print(p)
ggsave("5-A.pdf",p,width=13, height=6,dpi=300)