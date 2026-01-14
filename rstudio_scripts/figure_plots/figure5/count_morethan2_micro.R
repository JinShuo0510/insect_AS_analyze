# 初始化结果数据框
result_df <- data.frame(
  Species = character(),
  Tissue = character(),
  Age = character(),
  Event_count = integer(),
  SE_count = integer(),
  MXE_count = integer(),
  RI_count = integer(),
  A3SS_count = integer(),
  A5SS_count = integer(),
  devAS_count = integer(),
  SE_devAS_count = integer(),
  MXE_devAS_count = integer(),
  RI_devAS_count = integer(),
  A3SS_devAS_count = integer(),
  A5SS_devAS_count = integer(),
  stringsAsFactors = FALSE
)

# 辅助函数：检查外显子长度是否小于28nt
check_exon_length <- function(event_id) {
  # 分割Event_ID获取相关信息
  parts <- strsplit(event_id, "_")[[1]]
  # 获取方向后的两个数值（外显子起始和结束位置）
  direction_index <- which(parts %in% c("+", "-"))
  if (length(direction_index) == 0) return(FALSE) # 如果没有方向信息，跳过
  
  exon_start <- as.numeric(parts[direction_index + 1])
  exon_end <- as.numeric(parts[direction_index + 2])
  exon_length <- abs(exon_end - exon_start) + 1 # 计算外显子长度
  
  return(exon_length < 28)
}

# 处理每个物种和组织
for (species in names(find_devas_results_1_8_full_deleteNULL)) {
  for (tissue in names(find_devas_results_1_8_full_deleteNULL[[species]])) {
    tryCatch({
      # 获取DevAS数据
      devAS_data <- find_devas_results_1_8_full_deleteNULL[[species]][[tissue]]
      
      # 如果没有DevAS数据，跳过此组织
      if (nrow(devAS_data) == 0) next
      
      # 获取filtered_list中对应的数据
      filtered_data <- filtered_list[[species]][[tissue]]
      
      # 对filtered_data按Event_ID去重
      filtered_data <- filtered_data[!duplicated(filtered_data$Event_ID), ]
      
      # 对devAS_data按Event_ID去重
      devAS_data <- devAS_data[!duplicated(devAS_data$Event_ID), ]
      
      # 提取时期信息
      col_names <- colnames(filtered_data)[-1]
      ages <- unique(sapply(strsplit(col_names, "-"), `[`, 1))
      
      # 获取所有事件的类型和ID
      all_event_ids <- filtered_data$Event_ID
      all_event_types <- sapply(strsplit(all_event_ids, "_"), `[`, 1)
      
      # 筛选外显子长度小于28nt的事件
      short_exon_events <- sapply(all_event_ids, check_exon_length)
      all_event_ids <- all_event_ids[short_exon_events]
      all_event_types <- all_event_types[short_exon_events]
      
      # 获取在filtered_list中存在的DevAS事件的ID和类型
      devAS_event_ids <- intersect(devAS_data$Event_ID, all_event_ids)
      devAS_event_types <- sapply(strsplit(devAS_event_ids, "_"), `[`, 1)
      
      # 对每个时期进行统计
      for (age in ages) {
        # 获取当前时期的列
        age_cols <- grep(paste0("^", age, "-"), col_names, value = TRUE)
        
        # 检查当前时期的列是否有非NA值（所有事件）
        all_non_na_rows <- rowSums(!is.na(filtered_data[, age_cols, drop = FALSE])) > 0
        
        # 进一步筛选外显子长度小于28nt的事件
        all_non_na_rows <- all_non_na_rows & (filtered_data$Event_ID %in% all_event_ids)
        
        # 预设所有事件类型计数为0
        all_event_counts <- c(SE = 0, MXE = 0, RI = 0, A3SS = 0, A5SS = 0)
        
        # 统计所有事件类型的数量
        temp_counts <- table(all_event_types[all_non_na_rows[short_exon_events]])
        all_event_counts[names(temp_counts)] <- temp_counts
        
        # 检查当前时期的列是否有非NA值（DevAS事件）
        devAS_rows <- filtered_data$Event_ID %in% devAS_event_ids
        devAS_non_na_rows <- all_non_na_rows & devAS_rows
        
        # 预设DevAS事件类型计数为0
        devAS_counts <- c(SE = 0, MXE = 0, RI = 0, A3SS = 0, A5SS = 0)
        
        # 统计DevAS事件类型的数量
        temp_devAS_counts <- table(all_event_types[devAS_non_na_rows[short_exon_events]])
        devAS_counts[names(temp_devAS_counts)] <- temp_devAS_counts
        
        # 检查并警告异常情况
        for (event_type in names(devAS_counts)) {
          if (devAS_counts[event_type] > all_event_counts[event_type]) {
            warning(paste("Warning: In", species, tissue, age, "- DevAS", event_type, "count (", devAS_counts[event_type], ") is greater than all event", event_type, "count (", all_event_counts[event_type], ")"))
          }
        }
        
        # 添加日志信息
        cat(sprintf("Species: %s, Tissue: %s, Age: %s\n", species, tissue, age))
        cat(sprintf("Total events: %d, DevAS events: %d\n", sum(all_non_na_rows), sum(devAS_non_na_rows)))
        cat(sprintf("Event counts: %s\n", paste(all_event_counts, collapse = ", ")))
        cat(sprintf("DevAS counts: %s\n", paste(devAS_counts, collapse = ", ")))
        cat("\n")
        
        # 将结果添加到数据框中
        result_df <- rbind(result_df, data.frame(
          Species = species,
          Tissue = tissue,
          Age = age,
          Event_count = sum(all_non_na_rows),
          SE_count = all_event_counts["SE"],
          MXE_count = all_event_counts["MXE"],
          RI_count = all_event_counts["RI"],
          A3SS_count = all_event_counts["A3SS"],
          A5SS_count = all_event_counts["A5SS"],
          devAS_count = sum(devAS_non_na_rows),
          SE_devAS_count = devAS_counts["SE"],
          MXE_devAS_count = devAS_counts["MXE"],
          RI_devAS_count = devAS_counts["RI"],
          A3SS_devAS_count = devAS_counts["A3SS"],
          A5SS_devAS_count = devAS_counts["A5SS"],
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      cat("Error occurred for species:", species, "and tissue:", tissue, "\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Traceback:", conditionCall(e), "\n")
    })
  }
}

# 计算百分比
final_df <- result_df
final_df$DevAS_percent <- (final_df$devAS_count / final_df$Event_count) * 100
final_df$SE_DevAS_percent <- (final_df$SE_devAS_count / final_df$SE_count) * 100
final_df$MXE_DevAS_percent <- (final_df$MXE_devAS_count / final_df$MXE_count) * 100
final_df$RI_DevAS_percent <- (final_df$RI_devAS_count / final_df$RI_count) * 100
final_df$A3SS_DevAS_percent <- (final_df$A3SS_devAS_count / final_df$A3SS_count) * 100
final_df$A5SS_DevAS_percent <- (final_df$A5SS_devAS_count / final_df$A5SS_count) * 100

# 打印最终结果
print(final_df)

saveRDS(final_df, "all_events_and_devAS_microexon_counts.RDS")

