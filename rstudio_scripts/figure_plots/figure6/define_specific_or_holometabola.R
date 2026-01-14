library(dplyr)
library(tidyr)
library(stringr) # 用于字符串处理

# 定义目标物种
target_species <- c("Bombyx_mori", "Drosophila_mojavensis", 
                    "Helicoverpa_armigera", "Tribolium_castaneum", "Aedes_aegypti")

# --- 步骤 1：处理所有数据并添加外显子类型 ---
# (从上一步复制代码，确保results$match_results中的PSI列被保留)
# 更严格的判断逻辑
result_table_intermediate <- results$match_results %>%
  # 首先按ortholog_group分组
  group_by(ortholog_group) %>%
  mutate(
    # 计算目标物种数量
    target_species_count = n_distinct(species[species %in% target_species]),
    # 计算总物种数量
    total_species_count = n_distinct(species),
    # 标记当前记录是否为目标物种
    is_target_species = species %in% target_species,
    # 判断外显子类型
    exon_type = case_when(
      # Species-specific条件：
      # 1. 当前记录是目标物种
      # 2. 该ortholog_group只在一个目标物种中出现
      # 3. (可选)可以添加总物种数量的限制
      is_target_species & 
        target_species_count == 1 ~ "Species-specific",
      
      # Holometabola-specific条件：
      # 1. 在两个或更多目标物种中出现
      target_species_count >= 2 ~ "Holometabola-specific",
      
      # 其他情况
      TRUE ~ "Other"
    )
  ) %>%
  ungroup()

# --- 步骤 2：识别PSI列并计算最高表达组织 ---

# 识别所有PSI列
psi_cols <- names(result_table_intermediate)[grepl("_PSI$", names(result_table_intermediate))]

# 创建一个唯一的行标识符，以防后续操作打乱顺序
result_table_intermediate <- result_table_intermediate %>% mutate(row_id = row_number())

# 处理PSI值：长格式化 -> 提取组织 -> 计算平均值 -> 找最大值
psi_analysis <- result_table_intermediate %>%
  select(row_id, all_of(psi_cols)) %>%
  # 将PSI数据转为长格式
  pivot_longer(
    cols = all_of(psi_cols),
    names_to = "psi_column",
    values_to = "psi_value",
    values_drop_na = TRUE # 忽略NA值
  ) %>%
  # 提取组织名称 (移除 Age- 前缀和 _PSI 后缀)
  mutate(
    Tissue = str_remove(psi_column, "^[^-]+-"), # 移除第一个'-'之前的部分和'-'本身
    Tissue = str_remove(Tissue, "_PSI$")       # 移除_PSI后缀
  ) %>%
  # 按原始行（外显子）和组织分组，计算平均PSI（处理不同时期的同组织）
  group_by(row_id, Tissue) %>%
  summarise(avg_psi = mean(psi_value, na.rm = TRUE), .groups = 'drop') %>%
  # 对每个原始行（外显子），找到平均PSI最高的组织
  group_by(row_id) %>%
  # 使用slice_max处理可能的平局，with_ties=FALSE只选一个（通常是第一个）
  slice_max(order_by = avg_psi, n = 1, with_ties = FALSE) %>% 
  # 如果需要报告所有平局的组织，可以用下面这行代替slice_max
  # summarise(max_psi_tissue = paste(sort(Tissue), collapse = "; "), max_avg_psi = first(avg_psi)) %>%
  ungroup() %>%
  # 选择需要的列进行合并
  select(row_id, max_psi_tissue = Tissue, max_avg_psi = avg_psi)

# --- 步骤 3：合并结果并整理 ---

# 将最高PSI组织信息合并回主表
final_result <- result_table_intermediate %>%
  left_join(psi_analysis, by = "row_id") %>%
  # 选择并排序列
  select(
    species,
    gene_id,
    ortholog_group,
    transcript_id,
    SE_chr,
    SE_strand,
    SE_exon_start,
    SE_exon_end,
    match_type,
    exon_type,
    max_psi_tissue, # 新增列：最高PSI组织
    max_avg_psi,    # 新增列：对应的最高平均PSI值
    target_species_count,
    is_target_species,
    # 可以按需保留或移除 row_id
    -row_id 
  ) %>%
  arrange(species, ortholog_group, gene_id)

# 显示结果的一部分
print(head(final_result))

# 显示目标物种的分类统计 (包含新列的概要)
summary_stats_psi <- final_result %>%
  filter(is_target_species) %>%
  group_by(species, exon_type, max_psi_tissue) %>%
  summarise(
    count = n(),
    mean_max_psi = mean(max_avg_psi, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(species, exon_type, desc(count))

print(summary_stats_psi)

# 保存详细结果
# write.csv(final_result, "annotated_exons_with_max_psi_tissue.csv", row.names = FALSE)
