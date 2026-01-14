# 首先过滤掉不需要的物种
species_to_remove <- c("Acyrthosiphon_pisum", "Bicyclus_anynana", 
                       "Drosophila_simulans", "Neodiprion_lecontei", 
                       "Drosophila_suzukii")

# 过滤数据
head_micro_data <- head_micro_data %>%
  filter(!Species %in% species_to_remove)
head_data <- head_data %>%
  filter(!Species %in% species_to_remove)

# 计算统计数据
up_stats <- head_data %>%
  group_by(Species) %>%
  summarise(
    total_events = n(),
    up_events = sum(direction == "down", na.rm = TRUE),
    up_percentage = up_events / total_events * 100
  ) %>%
  arrange(desc(up_percentage))

up_micro_stats <- head_micro_data %>%
  group_by(Species) %>%
  summarise(
    total_events = n(),
    up_events = sum(direction == "down", na.rm = TRUE),
    up_percentage = up_events / total_events * 100
  ) %>%
  arrange(desc(up_percentage))

# 计算置信区间
micro_ci <- calc_ci_up(up_micro_stats) %>% mutate(Type = "Microexons")
macro_ci <- calc_ci_up(up_stats) %>% mutate(Type = "Macroexons")

# 合并数据
combined_df <- bind_rows(micro_ci, macro_ci)
combined_df$Species <- factor(combined_df$Species,
                              levels = unique(combined_df$Species))

# 优化后的绘图代码
p <- ggplot(combined_df, aes(x = Species, y = pct/100, shape = Type)) + # 将百分比转换为小数
  # 添加网格线
  geom_hline(yintercept = seq(0, 0.65, 0.1), color = "gray90") +
  geom_vline(xintercept = 1:length(unique(combined_df$Species)), color = "gray90") +
  
  # 设置误差线
  geom_errorbar(aes(ymin = lower/100, ymax = upper/100),  # 将百分比转换为小数
                position = position_dodge(width = 0.5),
                width = 0.3,  # 增加横线宽度
                size = 1,    # 增加线条粗细
                color = "#8DA0CB") +
  
  # 设置数据点
  geom_point(position = position_dodge(width = 0.5),
             size = 3.5,     # 增加点的大小
             color = "#8DA0CB") +
  
  scale_shape_manual(values = c("Microexons" = 16, "Macroexons" = 1)) +
  scale_y_continuous(limits = c(0, 0.65), 
                     breaks = seq(0, 0.65, 0.1),
                     expand = c(0.02, 0)) +
  
  labs(x = "",
       y = "Percentage of 'down' exons",
       shape = "") +
  
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
ggsave("Fig5D.pdf", plot = p, width = 8, height = 5, dpi = 300)

