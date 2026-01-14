# 加载必要的 R 包
library(dplyr)
library(binom)
library(ggplot2)

# ============================
# 1. 数据处理部分：统计每个Tissue中总数以及direction为"up"的数量
# ============================
# 对microexon数据进行统计
micro_by_tissue <- micro_exon_results %>%
  # 过滤掉不需要的组织
  filter(Tissue != "Head_and_Thorax") %>%
  filter(Tissue != "Leg") %>%
  filter(Tissue != "Neuron") %>%
  filter(Tissue != "Excretory") %>%
  filter(Tissue != "Hemolymph") %>%
  filter(Tissue != "Abdomen") %>%
  filter(Tissue != "Bacteriocyte") %>%
  filter(Tissue != "Thorax") %>%
  

  group_by(Tissue) %>%
  summarise(
    micro_total = n(),                               # 每个组织的总样本数
    micro_up    = sum(direction == "up", na.rm = TRUE) # direction为"up"的数量
  ) %>%
  ungroup()

  combined_df <- bind_rows(Bombyx_mori_event_fits_all_exon, .id = "Tissue")
# 对macroexon数据进行统计
macro_by_tissue <- Fout_time_all_species_results  %>%
  filter(Tissue != "Head_and_Thorax") %>%
  filter(Tissue != "Leg") %>%
  filter(Tissue != "Neuron") %>%
  filter(Tissue != "Excretory") %>%
  group_by(Tissue) %>%
  summarise(
    macro_total = n(),                               # 每个组织的总样本数
    macro_up    = sum(direction == "up", na.rm = TRUE) # direction为"up"的数量
  ) %>%
  ungroup()

# 将两部分数据按 Tissue 合并
tissue_df <- inner_join(micro_by_tissue, macro_by_tissue, by = "Tissue")

# 计算每个组织中的"up"比例、95%置信区间（Wilson法）以及Fisher检验的p值
tissue_df <- tissue_df %>%
  rowwise() %>%
  mutate(
    # 计算"up"比例
    micro_prop = micro_up / micro_total,
    macro_prop = macro_up / macro_total,
    # 计算95%置信区间，使用 Wilson 方法（这里对每一行取第1个值）
    micro_low  = binom.confint(micro_up, micro_total, conf.level = 0.95, methods = "wilson")$lower[1],
    micro_high = binom.confint(micro_up, micro_total, conf.level = 0.95, methods = "wilson")$upper[1],
    macro_low  = binom.confint(macro_up, macro_total, conf.level = 0.95, methods = "wilson")$lower[1],
    macro_high = binom.confint(macro_up, macro_total, conf.level = 0.95, methods = "wilson")$upper[1],
    # 构造2×2列联表并进行Fisher精确检验（行：microexon vs. macroexon；列："up"与非"up"）
    fisher_p   = fisher.test(
      matrix(c(micro_up, micro_total - micro_up,
               macro_up, macro_total - macro_up),
             nrow = 2, byrow = TRUE)
    )$p.value
  ) %>%
  ungroup() %>%
  select(Tissue, 
         micro_total, micro_up, micro_prop, micro_low, micro_high,
         macro_total, macro_up, macro_prop, macro_low, macro_high, fisher_p)

# ============================
# 2. 数据转换为长格式以便绘图
# ============================
micro_plot <- tissue_df %>%
  select(Tissue, fisher_p, p = micro_prop, low = micro_low, high = micro_high) %>%
  mutate(ExonType = "Microexon")

macro_plot <- tissue_df %>%
  select(Tissue, fisher_p, p = macro_prop, low = macro_low, high = macro_high) %>%
  mutate(ExonType = "Macroexon")

plot_df <- bind_rows(micro_plot, macro_plot)

# ============================
# 3. 绘图部分：构建图形
# ============================
# 定义位置调整参数，控制点和误差线的水平位置
pd <- position_dodge(width = 0.6)

# 为每个组织添加显著性标记，根据 Fisher 检验的 p 值
significance_data <- tissue_df %>%
  mutate(
    stars = case_when(
      fisher_p < 0.001 ~ "***",
      fisher_p < 0.01  ~ "**",
      fisher_p < 0.05  ~ "*",
      TRUE             ~ "ns"
    ),
    star_y = 0.8  # 固定标记显示位置
  )

# （请确保事先定义好tissue_color_mapping，例如：
#  tissue_color_mapping <- c("Whole_body" = "black", "某组织1" = "red", "某组织2" = "blue", ...)
#  或根据你的实际组织设置颜色。）

# 绘制图形
p <- ggplot(plot_df, aes(x = Tissue, y = p, color = Tissue, fill = Tissue, shape = ExonType)) +
  # 绘制误差线
  geom_errorbar(aes(ymin = low, ymax = high), 
                width = 0.4,    # 误差线宽度
                size = 1.0,     # 误差线粗细
                position = pd) +
  # 绘制数据点
  geom_point(size = 4, position = pd) +
  # 添加显著性标记（星号），固定显示在 y = 0.5
  geom_text(data = significance_data, 
            aes(x = Tissue, y = star_y, label = stars),
            inherit.aes = FALSE, size = 6, color = "black") +
  # 设置点的形状：Microexon 为实心点，Macroexon 为空心点
  scale_shape_manual(values = c("Microexon" = 16, "Macroexon" = 1), name = NULL) +
  scale_color_manual(values = tissue_color_mapping) +
  scale_fill_manual(values = tissue_color_mapping) +
  scale_y_continuous(limits = c(0, 1)) +  # 调整y轴范围
  labs(x = "", y = "Proportion of 'up'") +
  # 调整图例显示，仅保留形状图例，去除颜色与填充图例
  guides(color = "none", fill = "none", shape = guide_legend(override.aes = list(color = "black"))) +
  theme_bw(base_size = 15) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

# 显示图形并保存为 PDF
print(p)
ggsave("Fig5C_by_Tissue_all.pdf", plot = p, width = 8, height = 6, dpi = 300)
