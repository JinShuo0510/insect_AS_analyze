library(ggplot2)
library(dplyr)


filtered_df <- results_df[results_df$Num_Samples <= 50, ]

averaged_df <- filtered_df %>%
  group_by(Num_Samples) %>%
  summarise(Mean_Num_Events = mean(Num_Events, na.rm = TRUE))


# 按 Num_Samples 分组，去除每个组内的极端值
cleaned_avg_df <- filtered_df %>%
  group_by(Num_Samples) %>%
  mutate(
    Q1 = quantile(Num_Events, 0.25, na.rm = TRUE),  # 第一四分位数
    Q3 = quantile(Num_Events, 0.75, na.rm = TRUE),  # 第三四分位数
    IQR = Q3 - Q1,  # 四分位距
    lower_bound = Q1 - 1.5 * IQR,  # 下界
    upper_bound = Q3 + 1.5 * IQR   # 上界
  ) %>%
  filter(Num_Events >= lower_bound & Num_Events <= upper_bound) %>%  # 去除极端值
  summarise(Mean_Num_Events = mean(Num_Events, na.rm = TRUE)) %>%  # 计算平均值
  ungroup()


# 查看清理后的数据
head(cleaned_avg_df)

filtered_df$Log_Mean_Num_Events <- log(filtered_df$Num_Events + 1)  # 加 1 避免 log(0)


p <- ggplot(filtered_df, aes(x = Num_Samples, y = Log_Mean_Num_Events)) +
  geom_point(color = "blue", alpha = 0.6) +  # 绘制散点图
  geom_smooth(method = "loess", span = 0.75, color = "purple", se = TRUE) +  # 添加LOESS拟合曲线
  labs(title = "LOESS Fit: Num_Samples (0-50) vs Log(Mean_Num_Events)",
       x = "Number of Samples",
       y = "Log(Mean_Num_Events)") +
  theme_minimal()



print(p)
