# 删除为NULL的组织
result <- lapply(find_devas_results_1_8_full, function(species) {
  species[!sapply(species, is.null)]
})

# 删除所有组织为NULL的物种
result <- result[sapply(result, length) > 0]

# 查看结果
print(result)

saveRDS(result,file="find_devas_results_1_8_full_deleteNULL.RDS")