library(dplyr)
library(tidyr)
# 步骤 1：从 candidate_details 中提取所有唯一的物种
all_species <- candidate_details %>%
  separate_rows(present_species_str, sep = ", ") %>%
  distinct(present_species_str) %>%
  pull(present_species_str)
# 初始化一个空的列表来存储每个物种的 gene_id
species_gene_ids <- list()
# 步骤 2 和 3：对每个物种执行操作
for (species in all_species) {
  # 提取包含当前物种的 ortholog_group
  ortholog_groups <- candidate_details %>%
    filter(grepl(species, present_species_str)) %>%
    pull(ortholog_group)
  # 过滤 orthogroup_SE_matches 数据
  filtered_SE_matches <- results$match_results %>%
    filter(ortholog_group %in% ortholog_groups & species == species)
  # 提取 gene_id 并去重
  gene_ids <- filtered_SE_matches %>%
    pull(gene_id) %>%
    unique()
  # 将结果存储到列表中
  species_gene_ids[[species]] <- gene_ids
}
# 将结果整理为一个表格（使用基础 R 替代 enframe）
result_table <- data.frame(
  species = rep(names(species_gene_ids), sapply(species_gene_ids, length)),
  gene_id = unlist(species_gene_ids),
  stringsAsFactors = FALSE
)
# 查看结果
print(result_table)
