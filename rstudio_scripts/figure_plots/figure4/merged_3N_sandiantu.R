# 加载必要的R包
library(dplyr)
library(tidyr)
library(ggplot2)

#####################
# 第一步：读取与整理数据
#####################

# 假设结果数据已加载到工作环境中的 merged_df 对象中，
# 原始文件列依次为：tissue, chr, start, end, gene, divisible, direction
df_raw <- merged_df

# 给列命名
# colnames(df_raw) <- c("tissue", "chr", "start", "end", "gene", "divisible", "direction")

# 仅保留关心的列：组织(tissue)、方向(direction)、能否被3整除(divisible)
df <- df_raw %>%
  select(tissue, direction, divisible)

# 去除指定组织中 direction 为 NA 的数据行
tissues_to_filter <- c("Fat_body", "Gland", "Head", "Ovary", "Testis", "Whole_body")
df <- df %>%
  filter(!(tissue %in% tissues_to_filter & is.na(direction)))

# 转换 divisible 列为数值型（若尚未转换）
df$divisible <- as.numeric(df$divisible)

#####################
# 第二步：分组统计，计算平均值和二项分布95%置信区间
#####################

# 将NA方向标记为"neutral"
df$direction <- ifelse(is.na(df$direction), "neutral", df$direction)

df_summary <- df %>%
  group_by(tissue, direction) %>%
  # summarize(
  #   n = n(),
  #   x = sum(divisible),
  #   p = x / n,  # 1的比例作为平均值
  #   conf.low = prop.test(x, n)$conf.int[1],
  #   conf.high = prop.test(x, n)$conf.int[2],
  #   .groups = "drop"
  # )
  summarize(
    n = n(),
    x = sum(divisible),
    p = x / n,
    # 使用 Wilson 方法计算95%置信区间
    conf.low = binom.confint(x, n, conf.level = 0.95, methods = "wilson")$lower,
    conf.high = binom.confint(x, n, conf.level = 0.95, methods = "wilson")$upper,
    .groups = "drop"
  )

#####################
# 第三步：组织配色映射
#####################

# 创建自定义颜色映射
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
}

# 手动映射已知组织
tissue_color_mapping["Bombyx_mori_as_filter_devas"] = "gray50"
tissue_color_mapping["Bombyx_mori_exons_filter_as"] = "black"

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

#####################
# 第四步：计算显著性差异并绘图
#####################

# 将 df_summary 数据转换为宽格式，使每个组织同时拥有 up 和 down 两组数据
df_wide <- df_summary %>%
  pivot_wider(names_from = direction, 
              values_from = c(n, x, p, conf.low, conf.high))

# 对于同时存在 up 与 down 数据的组织，计算比例差异的显著性
df_sig <- df_wide %>%
  filter(!is.na(n_up) & !is.na(n_down)) %>%  # 保证两个方向都有数据
  rowwise() %>%
  mutate(
    # 使用 prop.test 比较 up 与 down 两组，获取 p 值
    p_value = prop.test(c(x_up, x_down), c(n_up, n_down))$p.value,
    # 根据 p 值设定显著性标记
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup() %>%
  # 确定显著性标注的 y 坐标：取 up 与 down 组中较高的置信区间值，再加上适当偏移
  mutate(
    y_pos = pmax(conf.high_up, conf.high_down, na.rm = TRUE) + 0.05
  )

# 设置位置偏移参数
pd <- position_dodge(width = 0.5)

# 定义方向形状映射
direction_shapes <- c("up" = 16, "down" = 17, "neutral" = 20)

#定义特定组织的顺序
priority_tissues <- c("Bombyx_mori_exons_filter_as", "Bombyx_mori_as_filter_devas")
# 获取其他组织的原始顺序（排除已定义的优先组织）
other_tissues <- setdiff(unique(df_summary$tissue), priority_tissues)
# 组合成最终的组织顺序
final_tissue_order <- c(priority_tissues, other_tissues)
# 将tissue列转换为因子并设置顺序
df_summary$tissue <- factor(df_summary$tissue, levels = final_tissue_order)
df_sig$tissue <- factor(df_sig$tissue, levels = final_tissue_order)
# 绘图代码保持不变（使用调整顺序后的数据）
p <- ggplot(df_summary, aes(x = tissue, y = p, color = tissue, shape = direction)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, size = 0.6, position = pd) +
  scale_color_manual(values = tissue_color_mapping) +
  scale_shape_manual(values = direction_shapes) +
  labs(
    x = "Tissue",
    y = "Proportion (divisible by 3)",
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(data = df_sig, 
            aes(x = tissue, y = y_pos, label = signif),
            color = "black", size = 6, position = pd, vjust = 0,
            inherit.aes = FALSE)
p <- p + theme(
  axis.title.x = element_text(size = 18),  # x轴标题大小
  axis.title.y = element_text(size = 18),  # y轴标题大小
  axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # x轴刻度标签大小（保持角度）
  axis.text.y = element_text(size = 14)    # y轴刻度标签大小
)

# p + coord_cartesian(ylim = c(0, 上限值))
# 显示图形
# 保存图形到文件
ggsave("3N_output_plot_up_down.pdf", plot = p, width = 6, height = 6, dpi=300)

print(p)

