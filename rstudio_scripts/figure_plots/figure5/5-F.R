# 初始化空列表存储结果
head_data_list <- list()

# 遍历每个物种
for (species in names(find_devas_results_1_8_full_deleteNULL)) {
  species_data <- find_devas_results_1_8_full_deleteNULL[[species]]
  
  # 检查是否存在Head组织且Event_ID列非空
  if ("Head" %in% names(species_data)) {
    head_data <- species_data[["Head"]]
    
    # 确保Event_ID列存在且非空
    if ("Event_ID" %in% names(head_data) && length(head_data$Event_ID) > 0) {
      event_ids <- head_data$Event_ID
      head_df <- data.frame(
        Species = species,
        Event_ID = event_ids,
        stringsAsFactors = FALSE
      )
      head_data_list[[species]] <- head_df
    }
  }
}

# 合并所有有效数据框（跳过空数据）
final_head_table <- do.call(rbind, head_data_list)

# 输出结果
head(final_head_table)


# 辅助函数：提取外显子信息并计算微外显子和长度是否能被3整除
extract_exon_info <- function(event_id) {
  # 分割Event_ID获取相关信息
  parts <- strsplit(event_id, "_")[[1]]
  
  # 获取方向后的两个数值（外显子起始和结束位置）
  direction_index <- which(parts %in% c("+", "-"))
  if (length(direction_index) == 0) {
    return(data.frame(
      exon_start = NA,
      exon_end = NA,
      is_microexon = NA,
      length_divisible_by_3 = NA
    ))
  }
  
  exon_start <- as.numeric(parts[direction_index + 1])
  exon_end <- as.numeric(parts[direction_index + 2])
  exon_length <- abs(exon_end - exon_start) + 1
  
  # 计算是否为微外显子和长度是否能被3整除
  is_microexon <- exon_length < 28
  length_divisible_by_3 <- exon_length %% 3 == 0
  
  return(data.frame(
    exon_start = exon_start,
    exon_end = exon_end,
    is_microexon = is_microexon,
    length_divisible_by_3 = length_divisible_by_3
  ))
}

# 提取外显子信息并合并到原数据框
exon_info_list <- lapply(final_head_table$Event_ID, extract_exon_info)
exon_info_df <- do.call(rbind, exon_info_list)

# 合并到原数据框
final_head_table <- cbind(final_head_table, exon_info_df)

# 查看结果
head(final_head_table)

library(dplyr)
library(ggplot2)

# 假设你的数据为 df
# df 的字段包括：
#   Species, is_microexon, length_divisible_by_3, 以及其他列

# 1. 按物种和微外显子/大外显子进行分组，汇总总数与可被3整除的数目
df_summary <- final_head_table %>%
  group_by(Species, is_microexon) %>%
  summarise(
    total_exons = n(),
    total_3N = sum(length_divisible_by_3),
    .groups = "drop"
  ) %>%
  # 2. 计算比例和置信区间
  rowwise() %>%
  mutate(
    proportion_3N = total_3N / total_exons,
    ci = list(binom.test(total_3N, total_exons, conf.level = 0.95)$conf.int),
    ci_lower = ci[[1]],
    ci_upper = ci[[2]]
  ) %>%
  select(-ci)  # 删除临时列表列

# 假设要去除的物种列表
# species_to_remove <- c("Bicyclus_anynana", "Drosophila_suzukii", "Acyrthosiphon_pisum")  # 替换为实际需要去除的物种名称
# 方法1：直接过滤
# df_summary <- df_summary %>%
#   filter(!Species %in% species_to_remove)


# 2. 为方便绘图，将比例转换为百分比，并建立用于形状映射的变量 Type
combined_df <- df_summary %>%
  mutate(
    pct = proportion_3N * 100,
    lower = ci_lower * 100,
    upper = ci_upper * 100,
    Type = ifelse(is_microexon, "Microexons", "Macroexons")
  )

# 3. 定义颜色映射：TRUE 为 "#8DA0CB"，FALSE 为 "#ADD8E6"
my_colors <- c("TRUE" = "#66CC33",  
               "FALSE" = "#FF9900")

p <- ggplot(combined_df, aes(x = Species, y = pct/100, color = as.factor(is_microexon), shape = Type)) +
  geom_hline(yintercept = seq(0, 0.65, 0.1), color = "gray90") +
  geom_vline(xintercept = 1:length(unique(combined_df$Species)), color = "gray90") +
  
  geom_errorbar(aes(ymin = lower/100, ymax = upper/100),
                position = position_dodge(width = 0.5),
                width = 0.3,
                size = 1) +
  
  geom_point(position = position_dodge(width = 0.5),
             size = 3.5) +
  
  # 隐藏 shape 图例，保留颜色图例
  scale_shape_manual(values = c("Microexons" = 16, "Macroexons" = 16), guide = "none") +
  scale_color_manual(values = my_colors) +
  
  scale_y_continuous(limits = c(0, 0.65), 
                     breaks = seq(0, 0.65, 0.1),
                     expand = c(0.02, 0)) +
  
  labs(x = "",
       y = "Percentage of exons divisible by 3",
       color = "Exon Type") +
  
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(0.8, 0.9),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(colour = "black")
  )

print(p)
ggsave("Fig5F.pdf", plot = p, width = 8, height = 5, dpi = 300)
