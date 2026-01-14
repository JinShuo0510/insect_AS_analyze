# -------------------- 1. 遍历 results_by_orthogroup 筛选 DevAS 事件 --------------------
results_by_orthogroup_devas <- list()

# 遍历每个 orthogroup
for (og_id in names(results_by_orthogroup)) {
  og_df <- results_by_orthogroup[[og_id]]
  
  # 如果当前 orthogroup 没有数据，跳过
  if (is.null(og_df)) next
  
  # 得到该 orthogroup 包含的物种
  species_list <- unique(og_df$Species)
  filtered_species_df <- list()
  
  # 针对每个物种
  for (sp in species_list) {
    # 检查 DevAS 变量中是否有该物种的数据
    if (!(sp %in% names(find_devas_results_1_8_full_deleteNULL))) {
      message(sprintf("警告：DevAS数据中未找到物种 %s 的数据", sp))
      next
    }
    
    # 对该物种所有组织的 DevAS 事件提取基因ID
    gene_ids_devas <- c()
    sp_devas <- find_devas_results_1_8_full_deleteNULL[[sp]]
    for (tissue in names(sp_devas)) {
      df_tissue <- sp_devas[[tissue]]
      # Event_ID 是以下划线分隔，第二个元素是 geneID
      gene_ids <- sapply(df_tissue$Event_ID, function(x) {
        parts <- strsplit(x, "_")[[1]]
        if(length(parts) >= 2) parts[2] else NA
      })
      gene_ids_devas <- unique(c(gene_ids_devas, gene_ids))
    }
    
    # 筛选当前 orthogroup 中该物种的记录
    sp_rows <- og_df[og_df$Species == sp, ]
    sp_filtered <- sp_rows[sp_rows$GeneID %in% gene_ids_devas, ]
    
    if (nrow(sp_filtered) > 0) {
      filtered_species_df[[sp]] <- sp_filtered
    }
  }
  
  # 如果该 orthogroup 下有至少一个物种通过筛选，则保存该 orthogroup 的结果
  if (length(filtered_species_df) > 0) {
    results_by_orthogroup_devas[[og_id]] <- dplyr::bind_rows(filtered_species_df)
  } else {
    results_by_orthogroup_devas[[og_id]] <- NULL
  }
}

# -------------------- 2. 合并所有筛选后的 orthogroup 结果 --------------------
combined_results_devas <- dplyr::bind_rows(results_by_orthogroup_devas, .id = "Orthogroup")

# 查看结果
print(combined_results_devas)

# -------------------- 3. 保存结果 --------------------
saveRDS(combined_results_devas, file = "combined_results_devas_byOrthogroup.RDS")
saveRDS(results_by_orthogroup_devas, file = "results_devas_byOrthogroup.RDS")
