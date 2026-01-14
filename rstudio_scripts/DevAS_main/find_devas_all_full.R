library(dplyr)
library(tidyr)
library(broom)

analyze_tissue_data <- function(data) {
  # 转换为长数据格式
  df_long <- data %>%
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
  
  # 检查每个 Event_ID 在不同时期的有效样本数量
  event_stage_counts <- df_long %>%
    group_by(Event_ID, stage) %>%
    summarise(n_samples = sum(!is.na(psi)), .groups = 'drop') %>%
    filter(n_samples > 0)  # 保留至少有一个有效样本的时期
  
  # 筛选出至少有两个有效时期的 Event_ID
  valid_events <- event_stage_counts %>%
    group_by(Event_ID) %>%
    summarise(n_stages = n_distinct(stage), .groups = 'drop') %>%
    filter(n_stages >= 2)  # 至少需要两个有效时期
  
  # 如果没有任何 Event_ID 满足条件，返回 NULL
  if (nrow(valid_events) == 0) {
    message("No events found with at least two valid stages.")
    return(NULL)
  }
  
  # 仅保留有效 Event_ID 的数据
  df_long_filtered <- df_long %>%
    filter(Event_ID %in% valid_events$Event_ID)
  
  # 检查 stage 列中是否至少有两个不同的值
  if (length(unique(df_long_filtered$stage)) < 2) {
    message("Only one stage found in the data. Skipping statistical tests.")
    return(NULL)
  }
  
  # 执行统计检验
  kw_results <- df_long_filtered %>%
    group_by(Event_ID) %>%
    summarise(
      kw_pvalue = tryCatch({
        if (length(unique(stage)) == 2) {
          # 只有两个有效时期，使用 Wilcoxon 秩和检验
          wilcox.test(psi ~ stage, exact = FALSE)$p.value
        } else {
          # 大于两个有效时期，使用 Kruskal-Wallis 检验
          kruskal.test(psi ~ stage)$p.value
        }
      }, error = function(e) NA_real_)  # 如果出错，返回 NA
    )
  
  # 多重检验校正（手动处理 NA 值）
  kw_results <- kw_results %>%
    mutate(p_value_adj = ifelse(
      is.na(kw_pvalue),
      NA_real_,
      p.adjust(kw_pvalue, method = "BH")
    ))
  
  # 筛选显著事件
  significant_events <- kw_results$Event_ID[kw_results$p_value_adj < 0.05 & !is.na(kw_results$p_value_adj)]
  
  # 检查 significant_events 是否包含 NA
  if (all(is.na(significant_events))) {
    message("No significant events found after multiple testing correction.")
    return(NULL)
  }
  
  # 如果显著事件为空，直接返回空
  if (length(significant_events) == 0) {
    message("No significant events found after multiple testing correction.")
    return(NULL)
  }
  
  # 处理两两比较或秩和检验结果
  pairwise_results <- df_long_filtered %>%
    filter(Event_ID %in% significant_events) %>%
    group_by(Event_ID) %>%
    do({
      stages <- unique(.$stage)
      if (length(stages) == 2) {
        # 只有两个时期，使用 Wilcoxon 秩和检验
        result <- tryCatch({
          wilcox.test(psi ~ stage, data = ., exact = FALSE) %>%
            broom::tidy() %>%
            select(p.value) %>%
            mutate(group1 = stages[1], group2 = stages[2])  # 添加组别信息
        }, error = function(e) data.frame(group1 = character(), group2 = character(), p.value = numeric()))
      } else {
        # 大于两个时期，使用 pairwise.wilcox.test
        result <- tryCatch({
          pairwise.wilcox.test(.$psi, .$stage, p.adjust.method = "BH", exact = FALSE) %>%
            broom::tidy()
        }, error = function(e) data.frame(group1 = character(), group2 = character(), p.value = numeric()))
      }
      if (nrow(result) == 0) {
        data.frame(group1 = character(), group2 = character(), p.value = numeric())  # 确保返回包含 p.value 列的空数据框
      } else {
        result
      }
    }) %>%
    filter(nrow(.) > 0) %>%  # 移除空数据框
    { if ("p.value" %in% colnames(.)) arrange(., Event_ID, p.value) else . }  # 如果 p.value 列存在，按 p.value 排序
  
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
  
  # 如果进行了两两比较或秩和检验，将结果整理为宽格式并合并到最终结果中
  if (!is.null(pairwise_results)) {
    pairwise_wide <- pairwise_results %>%
      group_by(Event_ID) %>%
      mutate(pairwise_id = row_number()) %>%  # 为每对比较生成唯一ID
      pivot_wider(
        names_from = pairwise_id,
        values_from = c(group1, group2, p.value),
        names_glue = "pairwise_{pairwise_id}_{.value}"
      )
    
    final_results <- final_results %>%
      left_join(pairwise_wide, by = "Event_ID")
  }
  
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
