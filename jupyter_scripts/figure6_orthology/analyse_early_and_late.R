# 加载所需包
library(Mfuzz)
library(edgeR)
library(Biobase)

# 定义物种列表
species_list <- c("Acyrthosiphon_pisum", "Aedes_aegypti", "Apis_mellifera", 
                  "Blattella_germanica", "Bombyx_mori", "Drosophila_mojavensis", 
                  "Gryllus_bimaculatus", "Helicoverpa_armigera", "Tribolium_castaneum")

# 遍历每个物种
for (species in species_list) {
  input_file <- paste0(species, "_GeneExpression_GroupedData_origin.tsv")
  output_prefix <- species
  
  # 检查文件是否存在
  if (!file.exists(input_file)) {
    message(paste("文件", input_file, "不存在，跳过物种", species))
    next
  }
  
  message(paste("开始处理物种:", species))
  
  # 1. 读入数据
  data <- read.table(input_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  # 2. 定义函数，将逗号分隔的字符串转换为数值向量的平均值
  convert_to_numeric <- function(x) {
    if (is.numeric(x)) return(x)
    vals <- as.numeric(unlist(strsplit(as.character(x), split = ",")))
    return(mean(vals, na.rm = TRUE))
  }
  
  # 针对数据框中的每一列进行转换
  data_numeric <- data
  for(i in 1:ncol(data_numeric)) {
    data_numeric[, i] <- sapply(data_numeric[, i], convert_to_numeric)
  }
  
  # 3. 过滤低表达基因（总表达量 > 10）
  data_filtered <- data_numeric[rowSums(data_numeric, na.rm = TRUE) > 10, ]
  
  # 4. 对数据进行 log2 转换
  data_log <- log2(data_filtered + 1)
  
  # 5. 提取组织和阶段信息
  tissue_names <- sapply(colnames(data_log), function(x) {
    parts <- strsplit(x, "__")[[1]]
    return(tail(parts, n=1))
  })
  
  # 新增代码：将Egg和Pupa重命名为Whole_body
  tissue_names[tissue_names %in% c("Egg", "Pupa")] <- "Whole_body"
  
  stage_names <- sapply(colnames(data_log), function(x) {
    parts <- strsplit(x, "__")[[1]]
    return(parts[1])
  })
  
  unique_tissues <- unique(tissue_names)
  print(unique_tissues)
  
  # 定义期望的阶段顺序
  desired_order <- c("Egg", "Larva", "Pupa", "Adult")
  
  # 6. 对每个组织进行聚类分析
  default_cluster_num <- 4
  cluster_results <- list()
  
  # 创建存储早期和晚期基因的列表
  early_genes_by_tissue <- list()
  late_genes_by_tissue <- list()
  all_early_genes <- list()
  all_late_genes <- list()
  
  # 提前创建空的摘要数据框
  summary_data <- data.frame(
    Tissue = character(),
    Cluster = integer(),
    Trend = character(),
    GeneCount = integer(),
    stringsAsFactors = FALSE
  )
  
  for (tissue in unique_tissues) {
    # 找出该组织对应的列索引
    tissue_cols <- which(tissue_names == tissue)
    n_samples <- length(tissue_cols)
    
    if (n_samples < 2) {
      message(paste("跳过组织", tissue, "：样本数小于2"))
      next
    }
    
    # 提取并排序阶段信息
    tissue_stage <- stage_names[tissue_cols]
    tissue_stage_factor <- factor(tissue_stage, levels = desired_order)
    ord <- order(tissue_stage_factor, na.last = TRUE)
    tissue_cols <- tissue_cols[ord]
    ordered_stages <- tissue_stage[ord]
    
    # 确定聚类数
    cluster_num <- ifelse(n_samples < default_cluster_num, n_samples, default_cluster_num)
    message(paste("处理组织", tissue, "：样本数 =", n_samples, "，聚类数 =", cluster_num))
    
    # 提取组织特异的表达数据
    data_tissue <- data_log[, tissue_cols, drop = FALSE]
    
    # 数据清洗：移除存在 NA、NaN 或 Inf 的行
    data_tissue <- data_tissue[apply(data_tissue, 1, function(x) all(is.finite(x))), ]
    
    # 去除方差为0的行
    row_variances <- apply(data_tissue, 1, var)
    data_tissue <- data_tissue[row_variances > 0, ]
    
    if(nrow(data_tissue) == 0) {
      message(paste("组织", tissue, "在清洗后无有效基因，跳过"))
      next
    }
    
    # 创建 ExpressionSet 对象并标准化
    eset <- new("ExpressionSet", exprs = as.matrix(data_tissue))
    eset_scaled <- standardise(eset)
    
    # 执行软聚类
    cl <- mfuzz(eset_scaled, c = cluster_num, m = 1.25)
    cluster_results[[tissue]] <- cl
    
    # 绘制聚类中心图
    output_file <- paste0(output_prefix, "_mfuzz_clusters_", tissue, ".png")
    png(output_file, width = 1200, height = 800, res = 150)
    matplot(t(cl$centers), type = "l", lty = 1,
            col = rainbow(nrow(cl$centers)), lwd = 2,
            xlab = "Developmental Stage (ordered: Egg, Larva, Pupa, Adult)",
            ylab = "Standardised Expression",
            main = paste("Cluster Centers for tissue:", tissue))
    legend("topright", legend = rownames(cl$centers), col = rainbow(nrow(cl$centers)), lty = 1, lwd = 2)
    dev.off()
    
    # 初始化组织特异早期和晚期基因列表
    early_genes_tissue <- list()
    late_genes_tissue <- list()
    
    # 提取高归属度基因并识别早期/晚期基因
    membership <- cl$membership
    
    # 对每个聚类进行处理
    for (i in 1:cluster_num) {
      cluster_center <- cl$centers[i, ]
      
      # 计算趋势方向
      if (length(cluster_center) >= 2) {
        # 简单趋势判断：比较首末点
        trend <- cluster_center[length(cluster_center)] - cluster_center[1]
        is_early <- trend < 0
        is_late <- trend > 0
        
        # 更复杂的趋势判断：线性回归斜率
        if (length(cluster_center) > 2) {
          x <- 1:length(cluster_center)
          model <- lm(cluster_center ~ x)
          slope <- coef(model)[2]
          is_early <- slope < -0.1
          is_late <- slope > 0.1
        }
      } else {
        is_early <- FALSE
        is_late <- FALSE
      }
      
      # 设定归属度阈值并筛选高归属度基因
      membership_threshold <- 0.7
      high_membership_genes <- rownames(membership)[membership[, i] >= membership_threshold]
      
      # 如果基因太少，降低阈值
      if (length(high_membership_genes) < 10 && membership_threshold > 0.5) {
        membership_threshold <- 0.5
        high_membership_genes <- rownames(membership)[membership[, i] >= membership_threshold]
      }
      
      if (length(high_membership_genes) > 0) {
        # 创建基因数据框
        gene_data <- data.frame(
          Gene = high_membership_genes,
          Membership = membership[high_membership_genes, i],
          Cluster = i,
          Trend = ifelse(is_early, "Early", ifelse(is_late, "Late", "Neither")),
          Tissue = tissue,
          stringsAsFactors = FALSE
        )
        
        # 添加原始表达数据
        gene_expression <- data_tissue[high_membership_genes, , drop = FALSE]
        gene_data <- cbind(gene_data, gene_expression)
        
        # 分别处理早期和晚期基因
        if (is_early) {
          early_genes <- gene_data[gene_data$Trend == "Early", , drop = FALSE]
          if (nrow(early_genes) > 0) {
            early_genes_tissue[[paste0("Cluster_", i)]] <- early_genes
            all_early_genes[[paste0(tissue, "_Cluster_", i)]] <- early_genes
          }
        }
        
        if (is_late) {
          late_genes <- gene_data[gene_data$Trend == "Late", , drop = FALSE]
          if (nrow(late_genes) > 0) {
            late_genes_tissue[[paste0("Cluster_", i)]] <- late_genes
            all_late_genes[[paste0(tissue, "_Cluster_", i)]] <- late_genes
          }
        }
        
        # 更新趋势统计
        early_count <- sum(gene_data$Trend == "Early")
        late_count <- sum(gene_data$Trend == "Late")
        neither_count <- sum(gene_data$Trend == "Neither")
        
        if (early_count > 0) {
          summary_data <- rbind(summary_data, data.frame(
            Tissue = tissue, Cluster = i, 
            Trend = "Early", GeneCount = early_count,
            stringsAsFactors = FALSE
          ))
        }
        if (late_count > 0) {
          summary_data <- rbind(summary_data, data.frame(
            Tissue = tissue, Cluster = i, 
            Trend = "Late", GeneCount = late_count,
            stringsAsFactors = FALSE
          ))
        }
        if (neither_count > 0) {
          summary_data <- rbind(summary_data, data.frame(
            Tissue = tissue, Cluster = i, 
            Trend = "Neither", GeneCount = neither_count,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # 绘制带趋势标记的聚类图
    output_file_trend <- paste0(output_prefix, "_mfuzz_clusters_trend_", tissue, ".png")
    png(output_file_trend, width = 1200, height = 800, res = 150)
    matplot(t(cl$centers), type = "l", lty = 1,
            col = rainbow(nrow(cl$centers)), lwd = 2,
            xlab = "Developmental Stage (ordered: Egg, Larva, Pupa, Adult)",
            ylab = "Standardised Expression",
            main = paste("Cluster Centers for tissue:", tissue))
    
    # 添加趋势标记
    for (i in 1:cluster_num) {
      cluster_center <- cl$centers[i, ]
      if (length(cluster_center) >= 2) {
        trend <- cluster_center[length(cluster_center)] - cluster_center[1]
        is_early <- trend < 0
        is_late <- trend > 0
        if (is_early) {
          text(1, cluster_center[1], "Early", col = rainbow(nrow(cl$centers))[i], pos = 3)
        } else if (is_late) {
          text(length(cluster_center), cluster_center[length(cluster_center)], 
               "Late", col = rainbow(nrow(cl$centers))[i], pos = 3)
        }
      }
    }
    
    legend("topright", legend = rownames(cl$centers), col = rainbow(nrow(cl$centers)), lty = 1, lwd = 2)
    dev.off()
    
    # 合并该组织的所有早期和晚期基因
    if (length(early_genes_tissue) > 0) {
      all_early_tissue <- do.call(rbind, early_genes_tissue)
      early_genes_by_tissue[[tissue]] <- all_early_tissue
      write.csv(all_early_tissue, paste0(output_prefix, "_", tissue, "_early_genes.csv"), row.names = FALSE)
      message(paste("已保存组织", tissue, "的早期基因至", paste0(output_prefix, "_", tissue, "_early_genes.csv")))
    }
    
    if (length(late_genes_tissue) > 0) {
      all_late_tissue <- do.call(rbind, late_genes_tissue)
      late_genes_by_tissue[[tissue]] <- all_late_tissue
      write.csv(all_late_tissue, paste0(output_prefix, "_", tissue, "_late_genes.csv"), row.names = FALSE)
      message(paste("已保存组织", tissue, "的晚期基因至", paste0(output_prefix, "_", tissue, "_late_genes.csv")))
    }
  }
  
  # 保存趋势摘要
  write.csv(summary_data, paste0(output_prefix, "_gene_cluster_trends_summary.csv"), row.names = FALSE)
  
  # 标准化并合并所有组织的早期/晚期基因
  process_and_save_gene_list <- function(gene_list, output_filename) {
    if (length(gene_list) == 0) {
      message(paste("没有找到", output_filename, "的基因"))
      return(NULL)
    }
    
    # 获取所有列名
    all_colnames <- unique(unlist(lapply(gene_list, colnames)))
    
    # 标准化所有数据框
    standardized_list <- lapply(gene_list, function(df) {
      missing_cols <- setdiff(all_colnames, colnames(df))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {
          df[[col]] <- NA
        }
      }
      df <- df[, all_colnames]
      return(df)
    })
    
    # 合并数据框并保存
    combined_data <- do.call(rbind, standardized_list)
    write.csv(combined_data, output_filename, row.names = FALSE)
    message(paste("已保存", output_filename))
    return(combined_data)
  }
  
  # 保存合并的早期/晚期基因列表
  process_and_save_gene_list(all_early_genes, paste0(output_prefix, "_early_genes_all_tissues.csv"))
  process_and_save_gene_list(all_late_genes, paste0(output_prefix, "_late_genes_all_tissues.csv"))
  
  message(paste("物种", species, "分析完成！所有结果已保存。"))
}

