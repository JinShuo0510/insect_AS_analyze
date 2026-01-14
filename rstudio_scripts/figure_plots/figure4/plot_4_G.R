# 加载库
library(ggplot2)
library(dplyr)
library(readr)
library(ggbreak)
library(tidyr)
library(purrr)

# --- 配置 ---
input_csv_file <- "aggregated_motif_TCTCTC_positions.csv"
target_tissue <- "Whole_body" # <--- 设置你要绘制的组织
target_motif <- "AACTTT" # 你的目标基序 (用于标题)
# 文件名
motif_length <- 6 # 基序长度 (hexamer)
output_plot_file <- paste0("motif_", target_motif, "_stacked_proportion_covered_", target_tissue, ".pdf")

# 从用户提供的信息中获取目标组织的外显子总数
tissue_exon_counts <- list(
  "Fat_body" = 366,
  "Gland" = 482,
  "Head" = 1243,
  "Whole_body" = 1792,
  "Ovary" = 177,
  "Testis" = 698
)
# 检查目标组织是否存在于计数列表中
if (!target_tissue %in% names(tissue_exon_counts)) {
  stop("错误: 目标组织 '", target_tissue, "' 的外显子总数未知。请在 'tissue_exon_counts' 列表中添加该信息。")
}
total_exons_in_target_tissue <- tissue_exon_counts[[target_tissue]]
print(paste("信息: 使用", target_tissue, "组织的外显子总数:", total_exons_in_target_tissue))

# --- 数据加载和处理 ---
# 检查文件是否存在
if (!file.exists(input_csv_file)) {
  stop("错误: 输入的 CSV 文件 '", input_csv_file, "' 未找到。")
}
# 使用 read_csv 加载数据
tryCatch({
  all_data <- read_csv(input_csv_file, show_col_types = FALSE)
}, error = function(e) {
  stop("错误: 无法读取 CSV 文件。请检查文件格式和内容。错误信息: ", e$message)
})

# 检查必要的列是否存在
required_cols <- c("Tissue", "Region", "RelativePositionSpliceSite", "Direction", "Count")
missing_cols <- setdiff(required_cols, colnames(all_data))
if (length(missing_cols) > 0) {
  stop("错误: 输入的 CSV 文件缺少必要的列: ", paste(missing_cols, collapse=", "))
}

# 筛选目标组织的数据
tissue_data_filtered <- all_data %>%
  filter(Tissue == target_tissue) %>%
  filter(!is.na(RelativePositionSpliceSite)) # 确保相对位置有效

# 处理 'Direction' 列，将 NA 或其他值归为 "Other"
tissue_data_processed <- tissue_data_filtered %>%
  mutate(
    Direction_cat = as.character(Direction),
    Direction_cat = case_when(
      Direction_cat == "up" ~ "up",
      Direction_cat == "down" ~ "down",
      is.na(Direction_cat) ~ "Other",
      Direction_cat == "" ~ "Other",
      TRUE ~ "Other"
    ),
    Direction_cat = factor(Direction_cat, levels = c("up", "down", "Other"))
  )

# 定义区域分割
# upstream: -100 到 -1
# exon: 0 到 +50 
# downstream: +51 到 +100
upstream_range <- c(-100, -1)
exon_range <- c(0, 50)
downstream_range <- c(51, 150)

# 展开数据，使每个基序覆盖其所有位置
expanded_data <- tissue_data_processed %>%
  # 为每个基序生成其覆盖的位置序列 (从起始位置到起始位置 + 长度 - 1)
  mutate(CoveredPosition = map(RelativePositionSpliceSite, ~ seq(from = .x, to = .x + motif_length - 1, by = 1))) %>%
  # 展开数据框，使得每个覆盖的位置都有一行
  unnest(CoveredPosition)

# 根据位置添加区域分类
expanded_data <- expanded_data %>%
  mutate(
    Region_cat = case_when(
      CoveredPosition >= upstream_range[1] & CoveredPosition <= upstream_range[2] ~ "Upstream Intron",
      CoveredPosition >= exon_range[1] & CoveredPosition <= exon_range[2] ~ "Exon",
      CoveredPosition >= downstream_range[1] & CoveredPosition <= downstream_range[2] ~ "Downstream Intron",
      TRUE ~ "Other"
    ),
    Region_cat = factor(Region_cat, levels = c("Upstream Intron", "Exon", "Downstream Intron", "Other"))
  ) %>%
  filter(Region_cat != "Other") # 过滤掉不在这三个主要区域的数据

# 检查展开和过滤后是否有数据
if (nrow(expanded_data) == 0) {
  stop("错误: 在指定的区域范围内，没有基序覆盖任何位置。无法生成图形。")
} else {
  print(paste("信息：展开数据后，在区域范围内有", nrow(expanded_data), "条基序覆盖位置记录。"))
}

# --- 计算比例数据 ---
# 计算每个覆盖位置、每个区域和方向分类的基序贡献计数
motif_counts_per_covered_pos <- expanded_data %>%
  group_by(CoveredPosition, Region_cat, Direction_cat) %>%
  summarise(MotifContributionCount = sum(Count), .groups = 'drop')

# 计算比例
proportion_data <- motif_counts_per_covered_pos %>%
  mutate(
    TotalExonsInTissue = total_exons_in_target_tissue,
    Proportion = MotifContributionCount / TotalExonsInTissue
  ) %>%
  # 确保所有可能的X轴位置都有数据（即使计数为0）
  complete(
    CoveredPosition = c(
      seq(upstream_range[1], upstream_range[2], by = 1),
      seq(exon_range[1], exon_range[2], by = 1),
      seq(downstream_range[1], downstream_range[2], by = 1)
    ),
    Region_cat,
    Direction_cat,
    fill = list(MotifContributionCount = 0, Proportion = 0, TotalExonsInTissue = total_exons_in_target_tissue)
  )

# 检查计算后是否有数据
if (nrow(proportion_data) == 0) {
  stop("错误：计算比例后没有数据。请检查过滤和计算步骤。")
}

# --- 绘图 ---

# 使用柔和的颜色
direction_colors <- c(
  "up" = "#6495ED", # 柔和的蓝色
  "down" = "#F08080", # 柔和的红色
  "Other" = "#B0C4DE" # 柔和的灰蓝色
)

p <- ggplot(proportion_data, aes(x = CoveredPosition, y = Proportion, fill = Direction_cat)) +
  geom_col(position = "stack", width = 1) + # width=1 使柱子紧邻
  scale_fill_manual(values = direction_colors, name = "Direction") +
  # 添加区域标签
  annotate("text", x = mean(upstream_range), y = Inf, label = "Upstream Intron", 
           hjust = 0.5, vjust = -0.5, size = 4, fontface = "bold", alpha = 0.8) +
  annotate("text", x = mean(exon_range), y = Inf, label = "Exon", 
           hjust = 0.5, vjust = -0.5, size = 4, fontface = "bold", alpha = 0.8) +
  annotate("text", x = mean(downstream_range), y = Inf, label = "Downstream Intron", 
           hjust = 0.5, vjust = -0.5, size = 4, fontface = "bold", alpha = 0.8) +
  # 设置标签和标题
  labs(
    x = "Position Relative to Splice Sites", 
    y = "Proportion of Total Exons in Tissue") +
  # 应用主题
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 15),  # 增大图例标题字体大小
    legend.text = element_text(size = 13),   # 增大图例项目字体大小
    axis.title.x = element_text(size = 15),  # 增大x轴标题字体大小
    axis.title.y = element_text(size = 15),  # 增大y轴标题字体大小
    axis.text.x = element_text(size = 15),   # 增大x轴刻度标签字体大小
    axis.text.y = element_text(size = 15),   # 增大y轴刻度标签字体大小
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# 打印图形到屏幕
print(p)


# 保存图形
ggsave(output_plot_file, plot = p, width = 8, height = 5, dpi = 300)
print(paste("图形已保存到:", output_plot_file))