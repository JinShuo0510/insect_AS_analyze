library(dplyr)
library(tidyr)
library(purrr)

# 定义常量
meta_cols <- c("GeneID","geneSymbol","chr","strand",
               "exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")

# 获取所有物种名称
species_names <- names(all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]])

# 定义方向判断函数
judge_direction <- function(diff_vec, threshold = 0.05) {
  # 如果输入为NA或空，直接返回"none"
  if(length(diff_vec) == 0 || all(is.na(diff_vec))) {
    return("none")
  }
  
  # 过滤小幅变化
  sign_vec <- ifelse(abs(diff_vec) < threshold, 0, sign(diff_vec))
  
  # 如果所有差值均为非负且至少有正值，则整体上调
  if(all(sign_vec >= 0) && any(sign_vec == 1)) {
    return("up")
  }
  # 如果所有差值均为非正且至少有负值，则整体下调
  if(all(sign_vec <= 0) && any(sign_vec == -1)) {
    return("down")
  }
  
  # 如果既有正又有负，则判断第一次出现的方向
  pos_inds <- which(sign_vec == 1)
  neg_inds <- which(sign_vec == -1)
  if(length(pos_inds) > 0 && length(neg_inds) > 0) {
    if(min(pos_inds) < min(neg_inds)) {
      return("up-down")
    } else {
      return("down-up")
    }
  }
  
  return("none")
}

# 创建一个空的结果列表，用于存储每个物种处理后的数据框
result_list <- list()

# 循环处理每个物种
for (species in species_names) {
  tryCatch({
    # 读取当前物种的数据
    df_raw <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]][[species]]
    
    # 如果数据为空，则跳过该物种
    if (is.null(df_raw) || nrow(df_raw) == 0) {
      message(paste0("Skipping empty data for species: ", species))
      next
    }
    
    # 1. 转换为长表格
    df_long <- df_raw %>%
      # 将除了 meta_cols 以外的列都 pivot_longer
      pivot_longer(
        cols      = -all_of(meta_cols),
        names_to  = "Age_Tissue",  # 原列名，如 "Pupa-Pupa", "Larva-Whole_body" 等
        values_to = "PSI_string"   # 其中存的是如 "1.0,1.0,NA,..." 之类的字符串
      ) 
    
    # 检查是否能够通过"-"分割
    if(all(grepl("-", df_long$Age_Tissue))) {
      df_long <- df_long %>%
        # 分割 Age_Tissue => Age 和 Tissue
        separate(
          col    = "Age_Tissue",
          into   = c("Age", "Tissue"),
          sep    = "-",
          remove = FALSE
        ) %>%
        # 将 Tissue == "Egg" 或 Tissue == "Pupa" 的行改为 "Whole_body"
        mutate(
          Tissue = if_else(Tissue %in% c("Egg","Pupa"), "Whole_body", Tissue)
        )
    } else {
      # 如果无法分割，创建默认值
      message(paste0("Cannot split Age_Tissue for species: ", species, ". Using defaults."))
      df_long <- df_long %>%
        mutate(
          Age = "Unknown",
          Tissue = "Unknown"
        )
    }
    
    # 创建唯一事件ID
    df_long <- df_long %>%
      mutate(event_id = paste(GeneID, chr, strand, exonStart_0base, exonEnd,
                              upstreamES, upstreamEE, downstreamES, downstreamEE,
                              sep = "__"))
    
    # 2. 数据预处理
    df_long_parsed <- df_long %>%
      # 只保留感兴趣的阶段
      filter(Age %in% c("Egg", "Larva", "Nymph", "Pupa", "Adult")) %>%
      # 将幼虫期名称统一，Nymph 也视为 Larva
      mutate(
        Age = if_else(Age %in% c("Larva", "Nymph"), "Larva", Age)
      )
    
    # 检查是否有数据
    if(nrow(df_long_parsed) == 0) {
      message(paste0("No data with target ages for species: ", species))
      next
    }
    
    df_long_parsed <- df_long_parsed %>%
      # 将逗号分隔的 PSI_string 拆分为多行，处理警告
      separate_rows(PSI_string, sep = ",") %>%
      # 安全地转换为数值型
      mutate(PSI_value = suppressWarnings(as.numeric(PSI_string))) %>%
      # 过滤掉 NA 值
      filter(!is.na(PSI_value)) %>%
      # 对同一事件、组织、阶段取平均值
      group_by(event_id, Tissue, Age) %>%
      summarise(PSI_avg = mean(PSI_value, na.rm = TRUE), .groups = "drop")
    
    # 如果处理后没有数据，跳过该物种
    if (nrow(df_long_parsed) == 0) {
      message(paste0("No valid data after parsing for species: ", species))
      next
    }
    
    # 3. 将 Age 定义为因子，保证顺序为 Egg, Larva, Pupa, Adult
    df_long_parsed <- df_long_parsed %>%
      mutate(Age = factor(Age, levels = c("Egg", "Larva", "Pupa", "Adult")))
    
    # 4. 转换为宽表，每个事件在同一组织下得到对应阶段数据
    df_wide <- df_long_parsed %>%
      pivot_wider(names_from = Age, values_from = PSI_avg, values_fill = list(PSI_avg = NA))
    
    # 检查哪些年龄阶段列实际存在
    age_columns <- intersect(c("Egg", "Larva", "Pupa", "Adult"), names(df_wide))
    
    # 如果少于2个阶段，则跳过该物种
    if(length(age_columns) < 2) {
      message(paste0("Only ", length(age_columns), " age stages available for species: ", species))
      next
    }
    
    # 5. 计算非缺失阶段数
    df_wide <- df_wide %>%
      mutate(non_missing = rowSums(!is.na(select(., all_of(age_columns)))))
    
    # 仅保留非缺失阶段数达到要求的事件(至少有2个阶段，或者尽可能多的阶段)
    min_required_stages <- min(3, length(age_columns))
    df_wide <- df_wide %>%
      filter(non_missing >= min_required_stages)
    
    # 如果筛选后没有数据，跳过该物种
    if (nrow(df_wide) == 0) {
      message(paste0("No valid data after filtering for species: ", species))
      next
    }
    
    # 6. 计算方向
    df_wide <- df_wide %>%
      rowwise() %>%
      mutate(
        # 按照顺序提取存在的各阶段数值
        stage_values = list(unlist(across(all_of(age_columns)))),
        # 仅保留非 NA 数值
        stage_values = list(stage_values[!is.na(stage_values)]),
        # 计算连续差值：如果 stage_values 长度大于 1
        diffs = list(if(length(stage_values) > 1) diff(stage_values) else numeric(0)),
        direction = judge_direction(unlist(diffs), threshold = 0.05)
      ) %>%
      ungroup()
    
    # 清理中间变量
    df_wide <- df_wide %>%
      select(-stage_values, -diffs, -non_missing)
    
    # 添加物种信息
    df_wide <- df_wide %>%
      mutate(Species = species)
    
    # 存储结果
    result_list[[species]] <- df_wide
    
    # 输出进度信息
    message(paste0("Processed ", species, ": found ", nrow(df_wide), " events with ", 
                   paste(age_columns, collapse=", "), " stages"))
    
  }, error = function(e) {
    message(paste0("Error processing species ", species, ": ", e$message))
  })
}

# 合并所有结果
all_species_results <- bind_rows(result_list)
saveRDS(all_species_results,"Fout_time_all_species_results.RDS")

# 可以查看总结果
print(paste0("Total processed events across all species: ", nrow(all_species_results)))

# 或者按物种查看事件数量
summary_by_species <- all_species_results %>%
  group_by(Species) %>%
  summarise(event_count = n()) %>%
  arrange(desc(event_count))

print(summary_by_species)

# 也可以查看方向分布
direction_summary <- all_species_results %>%
  group_by(Species, direction) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Species, desc(count))

print(direction_summary)

saveRDS(direction_summary,"Fout_time_dircetion_summary.RDS")