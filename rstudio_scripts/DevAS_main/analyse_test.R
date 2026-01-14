library(dplyr)

# 假设 all_species_tissue_data_withoutavg 是一个包含多个物种数据的列表
# 每个物种的数据是一个 tibble，包含多个时期的 PSI 值

# 初始化一个空的列表来存储结果
results <- list()

# 遍历每个物种
for (species_name in names(all_species_tissue_data_withoutavg)) {
  species_data <- all_species_tissue_data_withoutavg[[species_name]]
  
  # 初始化一个空的列表来存储当前物种的结果
  species_results <- list()
  
  # 遍历每个时期
  for (tissue_name in names(species_data)) {
    tissue_data <- species_data[[tissue_name]]
    
    # 统计可变剪切事件的数量（行数）
    num_events <- nrow(tissue_data)
    
    # 统计样本数量（仅计算第一行的 PSI 值数量）
    num_samples <- tissue_data %>%
      slice(1) %>%  # 只取第一行
      select(-Event_ID) %>%  # 排除 Event_ID 列
      mutate(across(everything(), ~ ifelse(is.na(.x), 0, length(unlist(strsplit(.x, ",")))))) %>%
      max()  # 取最大值
    
    # 将结果存储在 species_results 中
    species_results[[tissue_name]] <- list(
      num_events = num_events,
      num_samples = num_samples
    )
  }
  
  # 将当前物种的结果存储在 results 中
  results[[species_name]] <- species_results
}

# 打印结果
print(results)
