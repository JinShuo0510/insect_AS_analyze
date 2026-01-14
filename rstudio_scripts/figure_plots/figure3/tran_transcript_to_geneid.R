# -------------------- 1. 读取 Orthogroup 文件 --------------------
og <- read.delim("filtered_orthogroups_full.txt", header = TRUE, stringsAsFactors = FALSE)

# -------------------- 2. 提取物种名称 --------------------
species_names <- sub("\\.protein$", "", colnames(og)[-1])
print(species_names)

# -------------------- 3. 提取各物种对应的 GeneID --------------------
gene_ids_list <- list()
for(i in seq_along(species_names)) {
  sp <- species_names[i]
  # 对应列下标为 i+1（因为第一列是 Orthogroup）
  col_values <- og[[i + 1]]
  # 对每个值拆分逗号，去除多余空白，并合并为一个向量
  gene_ids <- unlist(lapply(col_values, function(x) {
    # 去除首尾空白
    x <- trimws(x)
    # 如果非空则拆分，否则返回 NULL
    if(nzchar(x)) {
      ids <- unlist(strsplit(x, ",\\s*"))
      trimws(ids)
    } else {
      NULL
    }
  }))
  # 去重并剔除空字符
  gene_ids <- unique(gene_ids[gene_ids != ""])
  gene_ids_list[[sp]] <- gene_ids
}

# 检查结果
str(gene_ids_list)

# -------------------- 4. 筛选新数据结构中的SE事件 --------------------
# 准备结果列表，保持原始数据结构
results <- list()

for(sp in species_names) {
  # 检查 all_species_tissue_data_1 中是否包含该物种数据
  if(sp %in% names(all_species_tissue_data_1)) {
    # 为该物种创建结果子列表
    results[[sp]] <- list()
    
    # 获取该物种下所有组织
    tissues <- names(all_species_tissue_data_1[[sp]])
    
    for(tissue in tissues) {
      # 获取该组织下的事件数据框
      tissue_df <- all_species_tissue_data_1[[sp]][[tissue]]
      
      # 从Event_ID列中提取GeneID
      event_genes <- sapply(tissue_df$Event_ID, function(x) {
        # 假设格式为"SE_GeneID_其他信息"
        parts <- strsplit(x, "_")[[1]]
        if(length(parts) >= 2) {
          return(parts[2])  # 返回基因ID
        } else {
          return(NA)
        }
      })
      
      # 过滤符合条件的事件
      matching_genes <- gene_ids_list[[sp]]
      matched_rows <- event_genes %in% matching_genes
      
      if(any(matched_rows)) {
        # 保持原始数据结构
        results[[sp]][[tissue]] <- tissue_df[matched_rows, ]
      } else {
        # 如果没有匹配的行，仍然创建空的数据框保持结构一致性
        results[[sp]][[tissue]] <- tissue_df[0, ]
      }
    }
  } else {
    message(sprintf("警告：数据中未找到物种 %s 的信息", sp))
  }
}

# -------------------- 5. 输出结果 --------------------
# 保留原始输出结构
saveRDS(results, file = "filtered_species_tissue_data.RDS")