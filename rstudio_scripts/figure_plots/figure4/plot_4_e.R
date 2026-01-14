library(ggplot2)
library(dplyr)
library(ggpattern)
library(scales)  # alpha() 函数来自该包
library(tidyr)   # complete() 函数用于补全数据

#-------------------------
# 1. 读取原始数据
#-------------------------
df <- read.csv("all_tissues_enriched_in_foreground_vs_background.csv_annotated_results.csv")
df <- subset(df, region == "downstream")
names(df)[names(df) == "regulation"] <- "direction"

#-------------------------
# 2. 构建组织-颜色映射
#-------------------------
all_tissues <- unique(df$Tissue)

if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
}

# 手动映射已知组织
tissue_color_mapping["Bombyx_mori_as_filter_devas_uniq"] = "gray50"
tissue_color_mapping["Bombyx_mori_exons.sorted_filter_as"] = "black"

# 候选颜色
color_pool <- c(
  "#E41A1C", "#377EB8", "#984EA3", "#8DA0CB", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#A6D854", "#FFD92F", "#E5C494",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)

# 找出尚未指定颜色的组织
new_tissues <- setdiff(all_tissues, names(tissue_color_mapping))

# 分配新组织颜色
if (length(new_tissues) > 0) {
  available_colors <- setdiff(color_pool, tissue_color_mapping)
  new_mapping <- setNames(available_colors[1:length(new_tissues)], new_tissues)
  tissue_color_mapping <- c(tissue_color_mapping, new_mapping)
}

#-------------------------
# 3. 数据汇总
#-------------------------
df_summarized <- df %>%
  filter(significant..FDR...0.05. == "True") %>%
  mutate(category = if_else(
    Matched_RBP_Name == "No_Match" & Matched_Motif_ID == "No_Match",
    "Unknown", "Known"
  )) %>%
  group_by(Tissue, direction, category) %>%
  summarise(hexamer_count = n(), .groups = "drop")


# 补全所有组合，缺失补0
df_summarized <- df_summarized %>%
  complete(Tissue = all_tissues, direction = c("up", "down"), category = c("Known", "Unknown"),
           fill = list(hexamer_count = 0))

# 为每条记录匹配颜色：Known 用原色，Unknown 用 alpha 降低透明度
df_summarized <- df_summarized %>%
  mutate(
    tissue_color = tissue_color_mapping[Tissue],
    plot_fill = if_else(category == "Known",
                        tissue_color,
                        alpha(tissue_color, 0.5))
  )

#-------------------------
# 4. 分拆 up/down数据及添加 x 轴索引
#-------------------------
df_up <- df_summarized %>% filter(direction == "up")
df_down <- df_summarized %>% filter(direction == "down")

# 明确因子顺序
df_summarized <- df_summarized %>%
  mutate(Tissue_idx = as.numeric(factor(Tissue, levels = all_tissues)))
df_up <- df_up %>%
  mutate(Tissue_idx = as.numeric(factor(Tissue, levels = all_tissues)))
df_down <- df_down %>%
  mutate(Tissue_idx = as.numeric(factor(Tissue, levels = all_tissues)))

# 计算每个 Tissue 与 direction 的总计数，用于添加数字标注
df_labels <- df_summarized %>%
  group_by(Tissue, direction) %>%
  summarise(total = sum(hexamer_count), .groups = "drop") %>%
  complete(Tissue = all_tissues, direction = c("up", "down"), fill = list(total = 0)) %>%
  mutate(Tissue_idx = as.numeric(factor(Tissue, levels = all_tissues)),
         x_offset = if_else(direction == "up", Tissue_idx - 0.2, Tissue_idx + 0.2))

#-------------------------
# 5. 绘图并设置字体大小
#-------------------------
p <- ggplot() +
  # up 柱子
  geom_col(
    data = df_up,
    aes(x = Tissue_idx - 0.2, y = hexamer_count, fill = plot_fill),
    position = "stack",
    width = 0.35
  ) +
  # down 柱子，添加斜线图案
  geom_col_pattern(
    data = df_down,
    aes(x = Tissue_idx + 0.2, y = hexamer_count, fill = plot_fill),
    position = "stack",
    width = 0.35,
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.05,
    pattern_spacing = 0.05,
    pattern_colour = "black",
    pattern_fill = NA
  ) +
  # 添加柱子顶部数字标注
  geom_text(
    data = df_labels,
    aes(x = x_offset, y = total, label = total),
    vjust = -0.5,
    size = 5  # 调整标注字体大小，可根据需要调整
  ) +
  scale_x_continuous(
    breaks = seq_along(all_tissues),
    labels = all_tissues
  ) +
  scale_fill_identity() +
  labs(
    x = "Tissue",
    y = "Hexamer Count",
    fill = "Category"
  ) +
  theme_bw(base_size = 12) +  # 提高整体基准字号
  theme(
    axis.title.x = element_text(size = 22),   # 横轴标题字号
    axis.title.y = element_text(size = 22),   # 纵轴标题字号
    axis.text.x  = element_text(size = 18),     # 横轴刻度文字字号
    axis.text.y  = element_text(size = 18)      # 纵轴刻度文字字号
  )

print(p)


ggsave("4-downstream.pdf",p,width=8,height=5,dpi=300)
