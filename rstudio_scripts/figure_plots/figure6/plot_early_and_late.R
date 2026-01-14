# 定义需要处理的物种
species_list <- c("Bombyx_mori", 
                  "Drosophila_mojavensis", "Helicoverpa_armigera", "Tribolium_castaneum")

# 定义组织和阶段
tissues <- c("Head", "Whole_body")
stages <- c("early", "late")

# 创建一个空列表来存储所有基因文件的数据
all_gene_data <- list()

# 读取所有物种的文件
for (species in species_list) {
  for (tissue in tissues) {
    for (stage in stages) {
      file_name <- paste0(species, "_", tissue, "_", stage, "_genes.csv")
      file_path <- file_name  # 如需指定路径，请修改此行
      
      # 尝试读取文件，如果文件不存在则跳过
      tryCatch({
        file_data <- read.csv(file_path, stringsAsFactors = FALSE)
        
        # 将数据添加到列表中
        key <- paste0(species, "_", tissue, "_", stage)
        all_gene_data[[key]] <- file_data
        cat("Successfully loaded:", file_name, "\n")
      }, error = function(e) {
        cat("Warning: Could not find or read file:", file_name, "\n")
      })
    }
  }
}

# 创建函数检查result_table中的基因是否在早期/晚期基因文件中
check_genes_multi_species <- function(result_table, all_gene_data, species_list) {
  # 创建结果数据框
  result_df <- data.frame(
    row_name = rownames(result_table),
    species = result_table$species,
    gene_id = result_table$gene_id,
    Head_early = FALSE,
    Head_late = FALSE,
    Whole_body_early = FALSE,
    Whole_body_late = FALSE,
    stringsAsFactors = FALSE
  )
  
  # 遍历result_table中的每一行
  for (i in 1:nrow(result_df)) {
    species_name <- result_df$species[i]
    gene_id <- result_df$gene_id[i]
    
    # 检查该基因是否在目标物种列表中
    if (species_name %in% species_list) {
      # 检查该基因在各个文件中的存在情况
      for (tissue in tissues) {
        for (stage in stages) {
          key <- paste0(species_name, "_", tissue, "_", stage)
          
          # 如果有该物种-组织-阶段的数据
          if (key %in% names(all_gene_data)) {
            gene_file <- all_gene_data[[key]]
            
            # 检查基因ID是否在Gene列中
            if (gene_id %in% gene_file$Gene) {
              col_name <- paste0(tissue, "_", stage)
              result_df[i, col_name] <- TRUE
            }
          }
        }
      }
    }
  }
  
  return(result_df)
}

# 执行检查
gene_results <- check_genes_multi_species(result_table, all_gene_data, species_list)

# 添加一个总结列，表示基因是否在任何早期或晚期文件中出现
gene_results$any_early_or_late <- gene_results$Head_early | 
  gene_results$Head_late | 
  gene_results$Whole_body_early | 
  gene_results$Whole_body_late

# 筛选出至少在一个早期或晚期文件中出现的基因
filtered_genes <- gene_results[gene_results$any_early_or_late, ]

# 查看结果
print(head(filtered_genes))
cat("Total genes identified:", nrow(filtered_genes), "\n")

# 添加早期和晚期总结列
gene_results$any_early <- gene_results$Head_early | gene_results$Whole_body_early
gene_results$any_late <- gene_results$Head_late | gene_results$Whole_body_late

# 再次筛选出至少在一个早期或晚期文件中出现的基因
filtered_genes <- gene_results[gene_results$any_early_or_late, ]

# 保存结果到CSV文件
write.csv(filtered_genes, "multi_species_stage_genes.csv", row.names = FALSE)

# 创建一个更简洁的结果表
concise_results <- data.frame(
  row_name = filtered_genes$row_name,
  species = filtered_genes$species,
  gene_id = filtered_genes$gene_id,
  any_early = filtered_genes$any_early,
  any_late = filtered_genes$any_late,
  stringsAsFactors = FALSE
)

# 保存简洁版结果
write.csv(concise_results, "multi_species_stage_genes_concise.csv", row.names = FALSE)

# 输出结果的摘要
cat("\n====== 结果摘要 ======\n")
for (species in species_list) {
  # 筛选当前物种的基因
  species_genes <- concise_results[concise_results$species == species, ]
  
  early_count <- sum(species_genes$any_early, na.rm = TRUE)
  late_count <- sum(species_genes$any_late, na.rm = TRUE)
  
  cat(species, ":\n")
  cat("  早期基因数量:", early_count, "\n")
  cat("  晚期基因数量:", late_count, "\n")
  cat("  总计:", early_count + late_count, "\n\n")
}

# 按物种统计基因数量
species_count <- table(filtered_genes$species)
cat("\n各物种基因数量统计:\n")
print(species_count)

# 创建详细报告，包含每个组织-阶段的统计信息
detailed_summary <- data.frame(
  Species = character(),
  Head_early = numeric(),
  Head_late = numeric(),
  Whole_body_early = numeric(),
  Whole_body_late = numeric(),
  Total = numeric(),
  stringsAsFactors = FALSE
)

for (species in species_list) {
  species_genes <- filtered_genes[filtered_genes$species == species, ]
  
  if (nrow(species_genes) > 0) {
    head_early_count <- sum(species_genes$Head_early, na.rm = TRUE)
    head_late_count <- sum(species_genes$Head_late, na.rm = TRUE)
    body_early_count <- sum(species_genes$Whole_body_early, na.rm = TRUE)
    body_late_count <- sum(species_genes$Whole_body_late, na.rm = TRUE)
    total_count <- nrow(species_genes)
    
    detailed_summary <- rbind(detailed_summary, data.frame(
      Species = species,
      Head_early = head_early_count,
      Head_late = head_late_count,
      Whole_body_early = body_early_count,
      Whole_body_late = body_late_count,
      Total = total_count,
      stringsAsFactors = FALSE
    ))
  }
}

# 输出详细报告
cat("\n====== 详细统计 ======\n")
print(detailed_summary)
write.csv(detailed_summary, "species_detailed_summary.csv", row.names = FALSE)