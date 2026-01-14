# 加载必要的 R 包
library(dplyr)
library(tidyr)
library(binom)
library(ggplot2)

# ============================
# 1~3. 数据处理部分保持不变，删除Head_and_Thorax组织
# ============================

# 微外显子数据按组织汇总 SE_count 与 SE_devAS_count
micro_by_tissue <- final_df %>%
  group_by(Tissue) %>%
  summarise(micro_total = sum(SE_count, na.rm = TRUE),
            micro_devAS = sum(SE_devAS_count, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Tissue != "Head_and_Thorax") # 过滤掉Head_and_Thorax

# 常规（macro）外显子数据按组织汇总 SE_count 与 SE_devAS_count
macro_by_tissue <- all_events_and_devAS_count %>%
  group_by(Tissue) %>%
  summarise(macro_total = sum(SE_count, na.rm = TRUE),
            macro_devAS = sum(SE_devAS_count, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Tissue != "Head_and_Thorax") # 过滤掉Head_and_Thorax

# 将两部分数据按 Tissue 合并
tissue_df <- inner_join(micro_by_tissue, macro_by_tissue, by = "Tissue")

# 计算每个组织中 devAS 比例、95% CI（Wilson 法）及 Fisher 检验 p 值
tissue_df <- tissue_df %>%
  rowwise() %>%
  mutate(
    # 计算 devAS 比例
    micro_prop = micro_devAS / micro_total,
    macro_prop = macro_devAS / macro_total,
    # 计算 95% 置信区间（Wilson 方法）
    micro_low  = binom.confint(micro_devAS, micro_total, conf.level = 0.95, methods = "wilson")$lower[1],
    micro_high = binom.confint(micro_devAS, micro_total, conf.level = 0.95, methods = "wilson")$upper[1],
    macro_low  = binom.confint(macro_devAS, macro_total, conf.level = 0.95, methods = "wilson")$lower[1],
    macro_high = binom.confint(macro_devAS, macro_total, conf.level = 0.95, methods = "wilson")$upper[1],
    # 构造 2×2 列联表，行：微外显子与常规外显子，列：devAS 与 non-devAS
    fisher_p = fisher.test(
      matrix(c(micro_devAS, micro_total - micro_devAS,
               macro_devAS, macro_total - macro_devAS),
             nrow = 2, byrow = TRUE)
    )$p.value
  ) %>%
  ungroup() %>%
  select(Tissue, micro_total, micro_devAS, micro_prop, micro_low, micro_high,
         macro_total, macro_devAS, macro_prop, macro_low, macro_high, fisher_p)

# 数据转换为长格式以便绘图
micro_plot <- tissue_df %>%
  select(Tissue, fisher_p, p = micro_prop, low = micro_low, high = micro_high) %>%
  mutate(ExonType = "Microexon")

macro_plot <- tissue_df %>%
  select(Tissue, fisher_p, p = macro_prop, low = macro_low, high = macro_high) %>%
  mutate(ExonType = "Macroexon")

plot_df <- bind_rows(micro_plot, macro_plot)

# ============================
# 4. 新的绘图代码：增强误差线显示，删除样本量标签
# ============================

# 定义位置调整参数
pd <- position_dodge(width = 0.6)

# 为每个组织添加显著性标记
significance_data <- tissue_df %>%
  mutate(
    # 根据p值确定显著性星号
    stars = case_when(
      fisher_p < 0.001 ~ "***",
      fisher_p < 0.01 ~ "**",
      fisher_p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    # 设置星号固定位置y=0.5
    star_y = 0.5
  )

# 绘制主图
p <- ggplot(plot_df, aes(x = Tissue, y = p, color = Tissue, fill = Tissue, shape = ExonType)) +
  # 首先绘制误差线，增加宽度和视觉效果
  geom_errorbar(aes(ymin = low, ymax = high), 
                width = 0.4,      # 增加误差线的宽度
                size = 1.0,      # 增加误差线的粗细
                position = pd) +
  # 然后绘制点，使其位于误差线上方
  geom_point(size = 4, position = pd) +
  # 添加显著性标记（星号）固定在y=0.5位置
  geom_text(data = significance_data, 
            aes(x = Tissue, y = star_y, label = stars),
            inherit.aes = FALSE, size = 6, color = "black") +
  # 设置形状 - 微外显子为实心点，宏外显子为空心点
  scale_shape_manual(values = c("Microexon" = 16, "Macroexon" = 1),
                     name = NULL) + # 移除图例标题
  scale_color_manual(values = tissue_color_mapping) +
  scale_fill_manual(values = tissue_color_mapping) +
  scale_y_continuous(limits = c(0, 0.6)) + # 调整y轴范围，去掉负值部分
  labs(x = "", y = "Proportion of devAS") +
  # 只保留形状图例，移除颜色图例
  guides(color = "none", fill = "none", shape = guide_legend(override.aes = list(color = "black"))) +
  theme_bw(base_size = 15) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1), # 保留x轴标签倾斜
    axis.text.y = element_text(size = 14),
    legend.position = "top", # 图例置顶
    legend.title = element_blank(), # 确保图例没有标题
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1) # 加粗边框
  )

# 显示图形并保存为 PDF
print(p)
ggsave("Fig5B_by_Tissue_improved.pdf", plot = p, width = 8, height = 6, dpi = 300)