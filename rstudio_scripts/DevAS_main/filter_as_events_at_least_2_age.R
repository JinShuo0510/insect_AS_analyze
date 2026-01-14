library(dplyr)
library(purrr)
library(furrr)  # 用于并行计算

# 启用多核并行计算
plan(multisession, workers = availableCores() - 1)  # 使用所有可用核心减一

# 定义一个函数来筛选至少两个时期有值的可变剪切事件，并保留所有列
filter_events <- function(tissue_data) {
  # 获取所有时期列的名称（排除 Event_ID）
  period_cols <- setdiff(names(tissue_data), "Event_ID")
  
  # 如果时期列少于两列，跳过筛选
  if (length(period_cols) < 2) {
    return(tissue_data)  # 返回原始数据
  }
  
  # 计算每个事件在多少时期有值
  num_non_na <- rowSums(!is.na(tissue_data[, period_cols]))
  
  # 筛选至少两个时期有值的事件，并保留所有列
  tissue_data[num_non_na >= 2, ]
}


# 遍历每个物种和组织，应用筛选函数（并行化处理）
#这段代码是不正确的，其筛选的是物种下的组织个数。
#filtered_data <- future_map(all_species_tissue_data_withoutavg, ~ map(.x, filter_events))
filtered_data <- future_map(all_species_tissue_data_withoutavg_1_8, function(species_data) {
  map(species_data, filter_events)  # 对每个组织的数据应用筛选函数
})

# 输出筛选后的结果
print(filtered_data)