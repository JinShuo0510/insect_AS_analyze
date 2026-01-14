library(dplyr)
library(tidyr)


Bor_wholebody_data <- all_species_tissue_data_withoutavg$Aquatica_lateralis$Whole_body


# 筛选出除第一列外，所有列值都不为NA的行
filtered_data <- Bor_wholebody_data %>%
  filter(if_all(-1, ~ !is.na(.)))

# 查看筛选后的数据
str(filtered_data)
print(colnames(filtered_data))

#转换为长数据格式
df_long <- filtered_data %>%
  pivot_longer(
    cols = -Event_ID,
    names_to = "stage_sample",
    values_to = "psi_string"
  ) %>%
  separate(stage_sample, into = c("stage", "sample_type"), sep = "-") %>%
  select(-sample_type)

# 分割 PSI 值字符串并转换为数值
df_long <- df_long %>%
  mutate(psi = strsplit(psi_string, ",")) %>%
  unnest(psi) %>%
  mutate(psi = as.numeric(psi))

# 1. Kruskal-Wallis 检验
kw_results <- df_long %>%
  group_by(Event_ID) %>%
  summarise(
    kw_pvalue = kruskal.test(psi ~ stage)$p.value
  )


# 2. 多重检验校正 (例如，使用 BH 方法)
kw_results$p_value_adj <- p.adjust(kw_results$kw_pvalue, method = "BH")


# 3. 筛选显著事件 (例如，校正后 P 值 < 0.05)
significant_events <- kw_results$Event_ID[kw_results$p_value_adj < 0.05]


# # 4. 两两比较 (Wilcoxon 秩和检验)
# pairwise_results <- df_long %>%
#   filter(Event_ID %in% significant_events) %>%
#   group_by(Event_ID) %>%
#   do(
#     pairwise.wilcox.test(.$psi, .$stage, p.adjust.method = "BH") %>%
#       broom::tidy() # 使用 broom::tidy() 函数将结果转换为数据框
#   )

# 检查 significant_events 是否为空
if (length(significant_events) == 0) {
  print("没有显著事件，跳过两两比较。")
  pairwise_results <- data.frame()  # 返回空数据框
} else {
  # 执行两两比较 (Wilcoxon 秩和检验)
  pairwise_results <- df_long %>%
    filter(Event_ID %in% significant_events) %>%
    group_by(Event_ID) %>%
    do({
      # 捕获错误，避免程序中断
      result <- tryCatch({
        test <- pairwise.wilcox.test(.$psi, .$stage, p.adjust.method = "BH")
        # 手动提取结果
        data.frame(
          group1 = rownames(test$p.value),
          group2 = colnames(test$p.value),
          p.value = as.vector(test$p.value)
        )
      }, error = function(e) {
        # 如果出错，返回空数据框
        data.frame(group1 = NA, group2 = NA, p.value = NA)
      })
      result
    }) %>%
    filter(!is.na(group1))  # 过滤掉出错的行
}

# # 计算 dPSI
# dpsi_results <- df_long %>%
#   filter(Event_ID %in% significant_events) %>%
#   group_by(Event_ID) %>%
#   summarise(
#     dPSI = max(psi, na.rm = TRUE) - min(psi, na.rm = TRUE)
#   )

# 计算 dPSI
dpsi_results <- df_long %>%
  filter(Event_ID %in% significant_events) %>%
  group_by(Event_ID) %>%
  summarise(
    dPSI = if (all(is.na(psi))) {
      NA  # 如果 psi 列中所有值都是 NA，返回 NA
    } else {
      max(psi, na.rm = TRUE) - min(psi, na.rm = TRUE)  # 否则计算 dPSI
    }
  )


# 结合统计检验结果和 dPSI 值
final_results <- kw_results %>%
  filter(Event_ID %in% significant_events) %>%
  left_join(dpsi_results, by = "Event_ID")

# 筛选 DevAS 事件 (例如，校正后 P 值 < 0.05 且 dPSI > 0.2)
devAS_events <- final_results %>%
  filter(p_value_adj < 0.05, dPSI > 0.2)
