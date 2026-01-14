# 安装并加载必要的包
if (!require("data.table")) install.packages("data.table")
library(data.table)

process_alternative_splicing_optimized <- function(expression_file, as_data) {
  # 开始计时，用于总体进度报告
  start_time <- Sys.time()
  cat("开始处理可变剪切数据分析...\n")
  
  # 读取基因表达数据并转换为data.table
  cat("读取基因表达文件:", expression_file, "...\n")
  gene_expr <- fread(expression_file)
  cat("成功加载基因表达数据，共", nrow(gene_expr), "行。\n")
  
  # 预先添加结果列
  gene_expr[, `:=`(AS_count = 0L, dPSI_count = 0L)]
  
  # 转换as_data为data.table并创建索引
  cat("正在预处理可变剪切数据并创建索引...\n")
  as_data_dt <- list()
  tissue_count <- length(names(as_data))
  
  for (i in seq_along(as_data)) {
    tissue_name <- names(as_data)[i]
    cat(sprintf("[%d/%d] 处理组织 '%s'...\n", i, tissue_count, tissue_name))
    
    if (is.data.frame(as_data[[tissue_name]])) {
      dt <- as.data.table(as_data[[tissue_name]])
      cat("  - 提取的事件数量:", nrow(dt), "\n")
      
      # 提取基因ID部分用于匹配
      cat("  - 解析事件ID并创建索引...\n")
      dt[, gene_id := gsub("^[^_]+_([^_]+)_.*$", "\\1", Event_ID)]
      setkey(dt, gene_id)
      as_data_dt[[tissue_name]] <- dt
      cat("  - 索引创建完成!\n")
    }
  }
  cat("预处理完成，共处理", length(as_data_dt), "个组织。\n")
  
  # 设置进度显示间隔和计数器
  total_rows <- nrow(gene_expr)
  report_interval <- max(1, round(total_rows / 20)) # 报告大约20次进度
  last_report_time <- start_time
  
  cat(sprintf("\n开始逐行处理基因表达数据 (共%d行)...\n", total_rows))
  
  # 批量处理表达数据
  for (i in 1:total_rows) {
    # 显示进度更新（每5%或至少每10秒显示一次）
    current_time <- Sys.time()
    if (i %% report_interval == 0 || i == total_rows || 
        as.numeric(difftime(current_time, last_report_time, units="secs")) >= 10) {
      
      percent_done <- round(i / total_rows * 100, 1)
      elapsed <- round(as.numeric(difftime(current_time, start_time, units="mins")), 2)
      
      if (i < total_rows) {
        est_total <- elapsed / (i / total_rows)
        est_remaining <- est_total - elapsed
        cat(sprintf("[%s] 处理进度: %d/%d (%.1f%%) - 已用时间: %.2f分钟 - 预计剩余: %.2f分钟\n", 
                    format(current_time, "%H:%M:%S"), i, total_rows, percent_done, elapsed, est_remaining))
      } else {
        cat(sprintf("[%s] 处理完成: %d/%d (100%%) - 总用时: %.2f分钟\n", 
                    format(current_time, "%H:%M:%S"), i, total_rows, elapsed))
      }
      
      last_report_time <- current_time
    }
    
    # 提取当前行数据
    tissue <- gene_expr$tissue[i]
    gene_id <- gene_expr$GeneID[i]
    prev_main <- gene_expr$prev_main_stage[i]
    prev_sub <- gene_expr$prev_sub_stage[i]
    current_main <- gene_expr$main_stage[i]
    current_sub <- gene_expr$sub_stage[i]
    
    # 清除基因ID中可能的前缀（如"+"）
    clean_gene_id <- gsub("^\\+", "", gene_id)
    
    # 构建前一个和当前阶段的列名
    prev_column <- paste0(prev_main, "-", prev_sub, "-", tissue)
    current_column <- paste0(current_main, "-", current_sub, "-", tissue)
    
    # 检查组织是否存在
    if (!tissue %in% names(as_data_dt)) {
      if (i %% report_interval == 0) cat(sprintf("  - 跳过行 %d: 未找到组织 '%s'\n", i, tissue))
      next
    }
    
    tissue_dt <- as_data_dt[[tissue]]
    
    # 检查列名是否存在
    if (!(prev_column %in% names(tissue_dt)) || !(current_column %in% names(tissue_dt))) {
      if (i %% report_interval == 0) cat(sprintf("  - 跳过行 %d: 未找到列 '%s' 或 '%s'\n", i, prev_column, current_column))
      next
    }
    
    # 使用索引高效查找与基因相关的事件
    gene_events <- tissue_dt[clean_gene_id]
    
    if (nrow(gene_events) == 0) {
      if (i %% report_interval == 0) cat(sprintf("  - 行 %d: 基因 '%s' 无相关AS事件\n", i, clean_gene_id))
      next
    }
    
    # 记录AS事件总数
    as_count <- nrow(gene_events)
    gene_expr$AS_count[i] <- as_count
    
    # 使用向量化操作计算dPSI并统计
    prev_psi <- gene_events[[prev_column]]
    current_psi <- gene_events[[current_column]]
    
    # 同时存在有效PSI值的行
    valid_rows <- !is.na(prev_psi) & !is.na(current_psi)
    dpsi_values <- abs(current_psi[valid_rows] - prev_psi[valid_rows])
    
    # 计算dPSI > 0.2的数量
    dpsi_count <- sum(dpsi_values > 0.2, na.rm = TRUE)
    gene_expr$dPSI_count[i] <- dpsi_count
    
    if (i %% report_interval == 0) {
      cat(sprintf("  - 行 %d: 基因 '%s' 共有 %d 个AS事件, %d 个dPSI>0.2\n", 
                  i, clean_gene_id, as_count, dpsi_count))
    }
  }
  
  # 报告总结
  end_time <- Sys.time()
  total_time <- round(as.numeric(difftime(end_time, start_time, units="mins")), 2)
  cat(sprintf("\n处理完成! 总用时: %.2f分钟\n", total_time))
  cat(sprintf("总基因数: %d, 有AS事件的基因数: %d, 有显著dPSI(>0.2)的基因数: %d\n",
              nrow(gene_expr), 
              sum(gene_expr$AS_count > 0), 
              sum(gene_expr$dPSI_count > 0)))
  
  return(gene_expr)
}


# 使用示例
result <- process_alternative_splicing_optimized("gene_expression_fold_changes.csv", all_species_tissue_data_1[["Bombyx_mori"]])
