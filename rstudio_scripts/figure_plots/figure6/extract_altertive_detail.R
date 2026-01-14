library(dplyr)
library(tidyr)

# 步骤 1：从 candidate_details_alternative 中提取所有唯一的物种
all_species <- candidate_details_alternative %>%
  distinct(species) %>%
  pull(species)

# 初始化一个空的数据框来存储结果
result_table <- data.frame()

# 步骤 2 和 3：对每个物种执行操作
for (species in all_species) {
  # 提取包含当前物种的 ortholog_group
  ortholog_groups <- candidate_details_alternative %>%
    filter(species == !!species) %>%  # 使用 !! 确保 species 被正确解析
    pull(ortholog_group)
  
  # 过滤 results$match_results 数据并保留所有列
  filtered_SE_matches <- results$match_results %>%
    filter(ortholog_group %in% ortholog_groups & species == !!species) 
    # 确保每个基因只保留一条唯一记录
    # group_by(gene_id)
    # distinct(gene_id, .keep_all = TRUE)
  # slice(1) %>%
  # ungroup()
  
  # 将结果添加到主数据框
  result_table <- bind_rows(result_table, filtered_SE_matches)
}

# 结果按物种和基因ID排序
result_table_alternative_detail <- result_table %>%
  arrange(species, gene_id)

# 查看结果
print(result_table_alternative_detail)