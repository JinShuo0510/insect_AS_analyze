# 筛选出组织数大于5的物种名称
species_with_more_than_5_tissues <- names(all_species_tissue_data_withoutavg_1_8)[
  sapply(all_species_tissue_data_withoutavg_1_8, function(species) length(species) > 5)
]

# 输出满足条件的物种名称
print(species_with_more_than_5_tissues)

# 统计满足条件的物种个数
count <- length(species_with_more_than_5_tissues)
print(count)
