########################################################################
# B. 标记 Species-specific / Holometabola-specific + 所有有 PSI 值的组织
########################################################################
library(stringr)
library(dplyr)
library(tidyr)

### 1. 标记外显子类型 -----------------------------------------------------
target_species <- c("Bombyx_mori", "Drosophila_mojavensis", 
                    "Helicoverpa_armigera", "Tribolium_castaneum", 
                    "Aedes_aegypti")

result_table_intermediate <- result_table_alternative_detail %>% 
  group_by(ortholog_group) %>% 
  mutate(
    target_species_count = n_distinct(species[species %in% target_species]),
    is_target_species    = species %in% target_species,
    exon_type = case_when(
      is_target_species & target_species_count == 1 ~ "Species-specific",
      target_species_count >= 2                     ~ "Holometabola-specific",
      TRUE                                          ~ "Other"
    )
  ) %>% 
  ungroup() %>% 
  mutate(row_id = row_number())                    # 稳定行索引

### 2. 标记每条外显子在有 PSI 值的所有组织 -----------------------------------
psi_cols <- names(result_table_intermediate)[str_detect(names(result_table_intermediate), "_PSI$")]

psi_analysis <- result_table_intermediate %>% 
  select(row_id, all_of(psi_cols)) %>% 
  pivot_longer(
    cols = all_of(psi_cols),
    names_to  = "psi_column",
    values_to = "psi_value",
    values_drop_na = TRUE
  ) %>% 
  mutate(
    Tissue = str_remove(psi_column, "^[^-]+-"),   # 去掉前缀 “Age-” 或 “Stage-”
    Tissue = str_remove(Tissue, "_PSI$")          # 去掉后缀 “_PSI”
  ) %>% 
  distinct(row_id, Tissue) %>%                   # 每个 row_id–Tissue 组合只保留一次
  mutate(
    max_psi_tissue = Tissue,                     # 重用列名保持接口一致
    max_avg_psi    = 1                            # PSI 存在则记为 1
  ) %>% 
  select(row_id, max_psi_tissue, max_avg_psi)

### 3. 合并并输出 ---------------------------------------------------------
final_result <- result_table_intermediate %>% 
  left_join(psi_analysis, by = "row_id") %>% 
  select(
    species, gene_id, ortholog_group, transcript_id,
    SE_chr, SE_strand, SE_exon_start, SE_exon_end,
    match_type, exon_type, max_psi_tissue, max_avg_psi,
    target_species_count, is_target_species
  ) %>% 
  arrange(species, ortholog_group, gene_id)

cat("#B  完整结果行数 = ", nrow(final_result), "\n")

### 4. 统计摘要（可选） ----------------------------------------------------
summary_stats_psi_alternative <- final_result %>% 
  filter(is_target_species) %>% 
  group_by(species, exon_type, max_psi_tissue, .drop = FALSE) %>% 
  summarise(
    count        = n(),
    mean_max_psi = mean(max_avg_psi, na.rm = TRUE),
    .groups      = "drop"
  ) %>% 
  arrange(species, exon_type, desc(count))

print(summary_stats_psi_alternative)

### 5. 若需保存 ------------------------------------------------------------
# write.csv(final_result, "annotated_exons_with_max_psi_tissue.csv", row.names = FALSE)
