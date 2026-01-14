# 整理早期和晚期基因ID - 修复组织名提取问题
library(dplyr)

# 读取所有组织的早期和晚期基因文件
consolidate_gene_ids <- function() {
  # 获取工作目录下所有CSV文件名
  all_files <- list.files(pattern = ".*_(early|late)_genes\\.csv$")
  
  # 创建空的结果数据框
  all_genes <- data.frame(
    Gene = character(0),
    Trend = character(0),
    Tissue = character(0),
    stringsAsFactors = FALSE
  )
  
  # 处理每个文件
  for (file in all_files) {
    # 正确从文件名中提取组织名和趋势
    if (grepl("early_genes.csv$", file)) {
      trend <- "Early"
      # 提取组织名 (移除 "_early_genes.csv")
      tissue <- sub("_early_genes\\.csv$", "", file)
    } else if (grepl("late_genes.csv$", file)) {
      trend <- "Late"
      # 提取组织名 (移除 "_late_genes.csv")
      tissue <- sub("_late_genes\\.csv$", "", file)
    } else {
      # 跳过不匹配的文件
      next
    }
    
    message(paste("处理文件:", file, "组织:", tissue, "趋势:", trend))
    
    # 读取文件
    tryCatch({
      gene_data <- read.csv(file)
      
      # 提取基因ID
      if ("Gene" %in% colnames(gene_data)) {
        # 创建简化的数据框
        simplified <- data.frame(
          Gene = gene_data$Gene,
          Trend = trend,
          Tissue = tissue,
          stringsAsFactors = FALSE
        )
        
        # 添加到结果中
        all_genes <- rbind(all_genes, simplified)
      } else {
        message(paste("文件", file, "没有Gene列"))
      }
    }, error = function(e) {
      message(paste("处理文件", file, "时出错:", e$message))
    })
  }
  
  # 处理综合文件中的数据（如果存在）
  if (file.exists("early_genes_all_tissues.csv")) {
    early_all <- read.csv("early_genes_all_tissues.csv")
    if ("Gene" %in% colnames(early_all) && "Tissue" %in% colnames(early_all)) {
      # 确保从原始文件中保留完整的组织名
      simplified <- data.frame(
        Gene = early_all$Gene,
        Trend = "Early",
        Tissue = early_all$Tissue,
        stringsAsFactors = FALSE
      )
      all_genes <- rbind(all_genes, simplified)
    }
  }
  
  if (file.exists("late_genes_all_tissues.csv")) {
    late_all <- read.csv("late_genes_all_tissues.csv")
    if ("Gene" %in% colnames(late_all) && "Tissue" %in% colnames(late_all)) {
      simplified <- data.frame(
        Gene = late_all$Gene,
        Trend = "Late",
        Tissue = late_all$Tissue,
        stringsAsFactors = FALSE
      )
      all_genes <- rbind(all_genes, simplified)
    }
  }
  
  # 去除可能的重复项
  all_genes <- distinct(all_genes)
  
  # 保存结果
  write.csv(all_genes, "all_early_late_genes_summary.csv", row.names = FALSE)
  
  # 创建简化版本，仅包含唯一基因ID和它们的趋势
  unique_genes <- all_genes %>%
    group_by(Gene, Trend) %>%
    summarise(Tissues = paste(unique(Tissue), collapse = ","), .groups = "drop")
  
  write.csv(unique_genes, "unique_early_late_genes.csv", row.names = FALSE)
  
  # 统计信息
  stats <- all_genes %>%
    group_by(Tissue, Trend) %>%
    summarise(GeneCount = n_distinct(Gene), .groups = "drop")
  
  write.csv(stats, "early_late_gene_counts_by_tissue.csv", row.names = FALSE)
  
  return(all_genes)
}

# 执行统计
gene_summary <- consolidate_gene_ids()