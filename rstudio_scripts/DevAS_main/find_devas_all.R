library(dplyr)
library(tidyr)
library(purrr)

# 定义一个函数来执行你的分析
analyze_tissue_data <- function(data) {
  # 筛选出除第一列外，所有列值都不为NA的行
  filtered_data <- data %>%
    filter(if_all(-1, ~ !is.na(.)))
  
  # 转换为长数据格式
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
  
  # 计算每个时期的样本数量
  stage_sample_counts <- df_long %>%
    group_by(stage) %>%
    summarise(n_samples = n_distinct(Event_ID))
  
  # 过滤样本数量不足的时期
  valid_stages <- stage_sample_counts %>%
    filter(n_samples >= 2) %>%
    pull(stage)
  
  # 记录被过滤的时期
  excluded_stages <- stage_sample_counts %>%
    filter(n_samples < 2) %>%
    pull(stage)
  
  if (length(excluded_stages) > 0) {
    message("The following stages were excluded due to insufficient samples: ", paste(excluded_stages, collapse = ", "))
  }
  
  # 仅保留有效时期的数据
  df_long_filtered <- df_long %>%
    filter(stage %in% valid_stages)
  
  # 检查剩余时期的数量
  num_valid_stages <- length(valid_stages)
  
  if (num_valid_stages == 1) {
    message("Only one valid stage found after filtering. Skipping statistical tests.")
    return(NULL)
  }
  
  # 检查 stage 列中是否至少有两个不同的值
  if (length(unique(df_long_filtered$stage)) < 2) {
    message("Only one stage found in the data. Skipping statistical tests.")
    return(NULL)
  }
  
  if (num_valid_stages == 2) {
    # 只有两个有效时期，使用 Wilcoxon 秩和检验
    kw_results <- df_long_filtered %>%
      group_by(Event_ID) %>%
      summarise(
        kw_pvalue = wilcox.test(psi ~ stage)$p.value
      )
  } else {
    # 大于两个有效时期，使用 Kruskal-Wallis 检验
    kw_results <- df_long_filtered %>%
      group_by(Event_ID) %>%
      summarise(
        kw_pvalue = kruskal.test(psi ~ stage)$p.value
      )
  }
  
  # 多重检验校正
  kw_results$p_value_adj <- p.adjust(kw_results$kw_pvalue, method = "BH")
  
  # 筛选显著事件
  significant_events <- kw_results$Event_ID[kw_results$p_value_adj < 0.05]
  
  # 检查 significant_events 是否包含 NA
  if (all(is.na(significant_events))) {
    message("No significant events found after multiple testing correction.")
    return(NULL)
  }
  
  # 如果 significant_events 包含 NA，则移除 NA
  significant_events <- significant_events[!is.na(significant_events)]
  
  # 如果显著事件为空，直接返回空
  if (length(significant_events) == 0) {
    message("No significant events found after multiple testing correction.")
    return(NULL)
  }
  
  # 两两比较（仅适用于大于两个时期）
  if (num_valid_stages > 2) {
    pairwise_results <- df_long_filtered %>%
      filter(Event_ID %in% significant_events) %>%
      group_by(Event_ID) %>%
      do(
        pairwise.wilcox.test(.$psi, .$stage, p.adjust.method = "BH") %>%
          broom::tidy()
      )
  } else {
    pairwise_results <- NULL
  }
  
  # 计算 dPSI
  dpsi_results <- df_long_filtered %>%
    filter(Event_ID %in% significant_events) %>%
    group_by(Event_ID) %>%
    summarise(
      dPSI = max(psi, na.rm = TRUE) - min(psi, na.rm = TRUE)
    )
  
  # 结合统计检验结果和 dPSI 值
  final_results <- kw_results %>%
    filter(Event_ID %in% significant_events) %>%
    left_join(dpsi_results, by = "Event_ID")
  
  # 筛选 DevAS 事件
  devAS_events <- final_results %>%
    filter(p_value_adj < 0.05, dPSI > 0.2)
  
  return(devAS_events)
}


# 遍历所有物种和组织
results <- all_species_tissue_data_withoutavg %>%
  map(~ map(.x, analyze_tissue_data))

# 查看结果
str(results)
