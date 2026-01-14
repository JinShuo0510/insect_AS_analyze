# 加载dplyr包
library(dplyr)
# 统计每个ortholog_group下非重复的species个数
species_count <- SE_matches %>%
  group_by(ortholog_group) %>%
  summarise(unique_species_count = n_distinct(species))
# 查看结果
print(species_count)