# 载入必要的包
library(dplyr)
library(tidyr)
library(ggplot2)

##############################################################################
# 1. 假设merged_df已经加载到R环境中
##############################################################################

# 指定需要过滤direction为NA的组织
tissues_to_filter <- c("Fat_body", "Gland", "Head", "Ovary", "Testis", "Whole_body")

# 数据预处理：仅去除指定组织中direction为NA的行
df <- merged_df %>%
  filter(!(tissue %in% tissues_to_filter & is.na(direction)))

# 将NA方向标记为"neutral"（如果需要保留NA值）
df$direction <- ifelse(is.na(df$direction), "neutral", df$direction)

# # 仅保留direction为up或down的数据
# df <- df %>%
#   filter(direction %in% c("up", "down"))

##############################################################################
# 2. 计算每个组织中不同direction下的统计量
##############################################################################

df_summary <- df %>%
  group_by(tissue, direction) %>%
  summarize(
    n = n(),
    mean_score = mean(Average_PhastCons_Score, na.rm = TRUE),
    sd_score = sd(Average_PhastCons_Score, na.rm = TRUE),
    se_score = sd_score / sqrt(n),
    conf.low = mean_score - 1.96 * se_score,  # 95% 置信区间下界
    conf.high = mean_score + 1.96 * se_score, # 95% 置信区间上界
    .groups = "drop"
  )

##############################################################################
# 3. 组织配色映射
##############################################################################

# 创建或加载自定义颜色映射
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
}
  # 创建初始颜色映射，可根据需要修改默认颜色
  tissue_color_mapping["Bombyx_mori_as_filter_devas_uniq"] = "gray50"
  tissue_color_mapping["Bombyx_mori_exons.sorted_filter_as"] = "black"


# 获取当前所有组织
current_tissues <- unique(df_summary$tissue)
# 找出映射中尚未出现的组织
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

# 定义颜色候选池
color_pool <- c(
  "#E41A1C", "#377EB8", "#984EA3", "#8DA0CB", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#A6D854", "#FFD92F", "#E5C494",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)

# 为其他组织分配颜色
if(length(new_tissues) > 0) {
  available_colors <- setdiff(color_pool, tissue_color_mapping)
  new_mapping <- setNames(available_colors[1:length(new_tissues)], new_tissues)
  tissue_color_mapping <- c(tissue_color_mapping, new_mapping)
}

# 保存组织-颜色映射
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

##############################################################################
# 4. 计算显著性差异
##############################################################################

# 将 df_summary 数据转换为宽格式，使每个组织同时拥有 up 和 down 两组数据
df_wide <- df_summary %>%
  pivot_wider(names_from = direction, 
              values_from = c(n, mean_score, sd_score, se_score, conf.low, conf.high))

# 对于同时存在 up 与 down 数据的组织，计算统计显著性
df_sig <- df %>%
  filter(direction %in% c("up", "down")) %>%
  group_by(tissue) %>%
  do({
    up_data <- filter(., direction == "up")$Average_PhastCons_Score
    down_data <- filter(., direction == "down")$Average_PhastCons_Score
    
    # 只有当两组都有足够的数据时才进行检验
    if(length(up_data) > 1 && length(down_data) > 1) {
      test_result <- t.test(up_data, down_data)
      data.frame(
        p_value = test_result$p.value,
        signif = case_when(
          test_result$p.value < 0.001 ~ "***",
          test_result$p.value < 0.01 ~ "**",
          test_result$p.value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    } else {
      data.frame(p_value = NA, signif = NA)
    }
  }) %>%
  ungroup()

# 合并显著性结果到df_wide
df_wide <- left_join(df_wide, df_sig, by = "tissue")

# 确定显著性标注的 y 坐标
df_wide <- df_wide %>%
  mutate(
    y_pos = pmax(conf.high_up, conf.high_down, na.rm = TRUE) + 0.05
  )

##############################################################################
# 5. 定义特定组织的显示顺序（若需要）
##############################################################################

# 定义特定组织的顺序（可根据需要调整）
priority_tissues <- c("Bombyx_mori_exons.sorted_filter_as", "Bombyx_mori_as_filter_devas_uniq")
# 获取其他组织的原始顺序（排除已定义的优先组织）
other_tissues <- setdiff(unique(df_summary$tissue), priority_tissues)
# 组合成最终的组织顺序
final_tissue_order <- c(priority_tissues, other_tissues)
# 将tissue列转换为因子并设置顺序
df_summary$tissue <- factor(df_summary$tissue, levels = final_tissue_order)
df_wide$tissue <- factor(df_wide$tissue, levels = final_tissue_order)

##############################################################################
# 6. 绘图
##############################################################################

# 设置位置偏移参数
pd <- position_dodge(width = 0.5)

# 定义方向形状映射
direction_shapes <- c("up" = 16, "down" = 17, "neutral" = 20)

# 绘制图形
p <- ggplot(df_summary, aes(x = tissue, y = mean_score, color = tissue, shape = direction)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, size = 0.6, position = pd) +
  scale_color_manual(values = tissue_color_mapping) +
  scale_shape_manual(values = direction_shapes) +
  labs(
    x = "Tissue",
    y = "Mean PhastCons Score with 95% CI",
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(data = df_wide, 
            aes(x = tissue, y = y_pos, label = signif),
            color = "black", size = 6, position = position_dodge(width = 0), vjust = 0,
            inherit.aes = FALSE)

p <- p + theme(
  axis.title.x = element_text(size = 18),  # x轴标题大小
  axis.title.y = element_text(size = 18),  # y轴标题大小
  axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # x轴刻度标签大小（保持角度）
  axis.text.y = element_text(size = 14)    # y轴刻度标签大小
)

# 显示图形
print(p)

# 保存图形到文件
ggsave("PhastCons_Score_UpDown_Plot.pdf", plot = p, width = 6, height = 6)
