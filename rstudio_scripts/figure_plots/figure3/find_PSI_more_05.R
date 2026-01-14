library(data.table)

# 读取第一个文件
gene_tissue_data <- fread("all_early_late_genes_summary.csv")
setDT(gene_tissue_data)

# 创建简单的进度显示函数
progress_tracker <- function(total, current, label = "Processing") {
  if(current %% max(1, floor(total/20)) == 0 || current == total) {
    percent <- round(current/total * 100)
    cat(sprintf("\r%s: %d/%d (%d%%)", label, current, total, percent))
    if(current == total) cat("\n")
    flush.console()
  }
}

# 处理每个组织和状态组合的函数
process_tissue_status <- function(tissue, status) {
  cat(sprintf("\nProcessing tissue: %s, status: %s\n", tissue, status))
  
  # 获取该组织下特定状态的基因ID
  gene_ids <- gene_tissue_data[Tissue == tissue & Trend == status, Gene]
  cat(sprintf("Found %d genes with status %s in %s\n", length(gene_ids), status, tissue))
  
  # 如果没有符合条件的基因，返回零结果
  if (length(gene_ids) == 0) {
    return(data.table(
      Tissue = tissue,
      Status = status,
      TotalEvents = 0,
      EventsWithHighDPSI = 0
    ))
  }
  
  # 检查R变量中是否存在该组织的数据
  if (!tissue %in% names(all_species_tissue_data_1_8[["Bombyx_mori"]])) {
    cat(sprintf("Tissue %s not found in data\n", tissue))
    return(data.table(
      Tissue = tissue,
      Status = status,
      TotalEvents = 0,
      EventsWithHighDPSI = 0
    ))
  }
  
  # 使用data.table提高效率
  tissue_data <- as.data.table(all_species_tissue_data_1_8[["Bombyx_mori"]][[tissue]])
  cat(sprintf("Loaded tissue data with %d events\n", nrow(tissue_data)))
  
  # 获取所有PSI列
  psi_columns <- grep("^(Larva|Adult|Pupa|Egg|Cell)-", names(tissue_data), value = TRUE)
  
  # 将Event_ID列提取出来以加快处理速度
  event_ids <- tissue_data$Event_ID
  
  # 分批处理基因ID，每次处理50个
  batch_size <- 500
  num_batches <- ceiling(length(gene_ids) / batch_size)
  
  # 初始化追踪匹配的事件
  matched_events_idx <- integer(0)
  
  # 分批处理基因ID
  cat(sprintf("Processing genes in %d batches\n", num_batches))
  for (batch in 1:num_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, length(gene_ids))
    current_genes <- gene_ids[start_idx:end_idx]
    
    # 针对每个基因ID单独匹配，然后合并结果
    batch_matches <- integer(0)
    for (i in 1:length(current_genes)) {
      gene_id <- current_genes[i]
      gene_pattern <- paste0("_", gene_id, "_")
      current_matches <- which(grepl(gene_pattern, event_ids, fixed = TRUE))
      batch_matches <- c(batch_matches, current_matches)
    }
    
    # 合并唯一匹配的索引
    matched_events_idx <- unique(c(matched_events_idx, batch_matches))
    progress_tracker(num_batches, batch, "Processing gene batches")
  }
  cat("\n")
  
  # 统计匹配事件
  total_events <- length(matched_events_idx)
  cat(sprintf("Found %d matching events\n", total_events))
  
  # 如果没有匹配的事件，返回零结果
  if (total_events == 0) {
    return(data.table(
      Tissue = tissue,
      Status = status,
      TotalEvents = 0,
      EventsWithHighDPSI = 0
    ))
  }
  
  # 处理找到的事件
  high_dpsi_events <- 0
  
  # 创建更小的子集来进行处理
  matched_data <- tissue_data[matched_events_idx]
  
  cat("Calculating dPSI values...\n")
  # 为每个匹配的事件计算最大dPSI
  for (i in 1:total_events) {
    max_dpsi <- 0
    
    # 遍历所有PSI列计算最大dPSI
    for (psi_col in psi_columns) {
      psi_value <- matched_data[[psi_col]][i]
      if (!is.na(psi_value)) {
        # 解析逗号分隔的PSI值
        psi_values <- as.numeric(unlist(strsplit(psi_value, ",")))
        if (length(psi_values) > 1) {
          # 计算最大PSI差异
          current_dpsi <- max(abs(diff(range(psi_values, na.rm = TRUE))), na.rm = TRUE)
          if (current_dpsi > max_dpsi) {
            max_dpsi <- current_dpsi
          }
        }
      }
    }
    
    # 统计dPSI > 0.5的事件
    if (max_dpsi > 0.5) {
      high_dpsi_events <- high_dpsi_events + 1
    }
    
    progress_tracker(total_events, i, "Calculating dPSI")
  }
  cat("\n")
  
  cat(sprintf("Results for %s (%s): %d total events, %d with high dPSI\n", 
              tissue, status, total_events, high_dpsi_events))
  
  return(data.table(
    Tissue = tissue,
    Status = status,
    TotalEvents = total_events,
    EventsWithHighDPSI = high_dpsi_events
  ))
}

# 获取唯一组织列表和状态列表
tissues <- unique(gene_tissue_data$Tissue)
statuses <- unique(gene_tissue_data$Trend)
cat(sprintf("Processing %d tissues with %d statuses\n", length(tissues), length(statuses)))

# 创建所有组织和状态的组合
tissue_status_combinations <- expand.grid(
  Tissue = tissues,
  Status = statuses,
  stringsAsFactors = FALSE
)
setDT(tissue_status_combinations)

# 总共需要处理的组合数
total_combinations <- nrow(tissue_status_combinations)
cat(sprintf("Total combinations to process: %d\n", total_combinations))

# 初始化结果列表
results <- list()
result_count <- 0

# 处理每个组织-状态组合
for (i in 1:total_combinations) {
  tissue <- tissue_status_combinations$Tissue[i]
  status <- tissue_status_combinations$Status[i]
  
  cat(sprintf("\n[%d/%d] Processing combination: %s - %s\n", 
              i, total_combinations, tissue, status))
  
  # 处理当前组合
  combination_result <- process_tissue_status(tissue, status)
  result_count <- result_count + 1
  results[[result_count]] <- combination_result
  
  # 显示当前累积结果
  current_results <- rbindlist(results[1:result_count])
  cat("\nCurrent results:\n")
  print(current_results)
  
  # 每处理5个组合或所有组合都处理完后保存一次中间结果
  if (i %% 5 == 0 || i == total_combinations) {
    intermediate_file <- paste0("tissue_splicing_results_", i, "_of_", total_combinations, ".csv")
    fwrite(current_results, intermediate_file)
    cat(sprintf("Intermediate results saved to %s\n", intermediate_file))
  }
}

# 合并所有结果
final_result <- rbindlist(results)

# 排序结果，先按组织，再按状态
setorder(final_result, Tissue, Status)

# 输出最终结果
fwrite(final_result, "tissue_splicing_events_summary_by_status.csv")
cat("\nFinal results:\n")
print(final_result)

# 创建摘要表
summary_table <- dcast(final_result, Tissue ~ Status, 
                       value.var = c("TotalEvents", "EventsWithHighDPSI"))
fwrite(summary_table, "tissue_splicing_summary_by_status_wide.csv")
cat("\nSummary table:\n")
print(summary_table)