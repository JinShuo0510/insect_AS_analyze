library(dplyr)
library(tidyr)
library(tibble)

# # 假设 all_species_tissue_data_withoutavg 是一个包含多个物种数据的列表
# # 每个物种的数据是一个 tibble，包含多个时期的 PSI 值
# 
# # 初始化一个空的列表来存储结果
# results <- list()
# 
# # 遍历每个物种
# for (species_name in names(filt)) {
#   species_data <- all_species_tissue_data_withoutavg[[species_name]]
# 
#   # 初始化一个空的列表来存储当前物种的结果
#   species_results <- list()
# 
#   # 遍历每个时期
#   for (tissue_name in names(species_data)) {
#     tissue_data <- species_data[[tissue_name]]
# 
#     # 统计可变剪切事件的数量（行数）
#     num_events <- nrow(tissue_data)
# 
#     # 统计样本数量（每行PSI值的个数的和，选最长）
#     # 使用向量化操作计算每行的 PSI 值数量
#     max_samples <- tissue_data %>%
#       select(-Event_ID) %>%  # 排除 Event_ID 列
#       mutate(across(everything(), ~ sapply(.x, function(x) ifelse(is.na(x), 0, length(unlist(strsplit(x, ","))))))) %>%
#       rowwise() %>%
#       mutate(max_samples = max(c_across(everything()))) %>%
#       ungroup() %>%
#       pull(max_samples) %>%
#       max()
# 
#     # 将结果存储在 species_results 中
#     species_results[[tissue_name]] <- list(
#       num_events = num_events,
#       max_samples = max_samples
#     )
#   }
# 
#   # 将当前物种的结果存储在 results 中
#   results[[species_name]] <- species_results
# }
# 
# # 打印结果
# print(results)
# 
# 
# results_df <- results %>%
#   # 将嵌套列表转换为长格式
#   enframe(name = "Species", value = "Tissue") %>%
#   # 展开 Tissue 列中的嵌套列表
#   unnest_longer(Tissue) %>%
#   # 将 Tissue 列中的嵌套列表进一步展开
#   unnest_wider(Tissue) %>%
#   # 重命名列
#   rename(
#     Tissue = Tissue_id,  # 如果 Tissue_id 是默认列名
#     Num_Events = num_events,
#     Num_Samples = num_samples
#   ) %>%
#   # 调整列顺序，将 Tissue 列放到第二列
#   select(Species, Tissue, Num_Events, Num_Samples)


# num_samples_file <- read.csv("sample_statistics_by_tissue.csv")  # 替换为实际文件路径
# 
# results_df <- results_df %>%
#   left_join(num_samples_file, by = c("Species", "Tissue"))
# 
# print(results_df)



