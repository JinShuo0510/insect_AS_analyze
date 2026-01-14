library(dplyr)

# 假设你的数据结构为 all_species_tissue_data_withoutavg
all_species_tissue_data_withoutavg_1_8 <- lapply(all_species_tissue_data_withoutavg_1_8, function(species) {
  lapply(species, function(tissue) {
    tissue %>%
      filter(rowSums(!is.na(select(., -Event_ID))) > 0)  # 向量化操作，检查每行是否有非 NA 值
  })
})

saveRDS(all_species_tissue_data_withoutavg_1_8, file = "all_species_tissue_data_1_8.RDS")