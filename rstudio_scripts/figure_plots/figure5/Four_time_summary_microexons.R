# 加载必要的库
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
  if(length(diff_vec) == 0 || all(is.na(diff_vec))) {
    return("none")
  }
  
  sign_vec <- ifelse(abs(diff_vec) < threshold, 0, sign(diff_vec))
  
  if(all(sign_vec >= 0) && any(sign_vec == 1)) {
    return("up")
  }
  if(all(sign_vec <= 0) && any(sign_vec == -1)) {
    return("down")
  }
  
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

# 创建空的结果列表
result_list <- list()

# 循环处理每个物种
for (species in species_names) {
  tryCatch({
    # 读取当前物种的数据
    df_raw <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]][[species]]
    
    if (is.null(df_raw) || nrow(df_raw) == 0) {
      message(paste0("[WARN] Skipping empty data for species: ", species))
      next
    }
    
    # 1. 转换为长表格
    df_long <- df_raw %>%
      pivot_longer(
        cols = -all_of(meta_cols),
        names_to = "Age_Tissue",
        values_to = "PSI_string"
      ) 
    
    # 分割Age和Tissue
    if(all(grepl("-", df_long$Age_Tissue))) {
      df_long <- df_long %>%
        separate(
          col = "Age_Tissue",
          into = c("Age", "Tissue"),
          sep = "-",
          remove = FALSE
        ) %>%
        mutate(
          Tissue = if_else(Tissue %in% c("Egg","Pupa"), "Whole_body", Tissue)
        )
    } else {
      message(paste0("[WARN] Cannot split Age_Tissue for species: ", species, ". Using defaults."))
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
      filter(Age %in% c("Egg", "Larva", "Nymph", "Pupa", "Adult")) %>%
      mutate(
        Age = if_else(Age %in% c("Larva", "Nymph"), "Larva", Age)
      )
    
    if(nrow(df_long_parsed) == 0) {
      message(paste0("[WARN] No data with target ages for species: ", species))
      next
    }
    
    df_long_parsed <- df_long_parsed %>%
      separate_rows(PSI_string, sep = ",") %>%
      mutate(PSI_value = suppressWarnings(as.numeric(PSI_string))) %>%
      filter(!is.na(PSI_value)) %>%
      group_by(event_id, Tissue, Age) %>%
      summarise(PSI_avg = mean(PSI_value, na.rm = TRUE), .groups = "drop")
    
    if (nrow(df_long_parsed) == 0) {
      message(paste0("[WARN] No valid data after parsing for species: ", species))
      next
    }
    
    # 3. 将Age定义为因子
    df_long_parsed <- df_long_parsed %>%
      mutate(Age = factor(Age, levels = c("Egg", "Larva", "Pupa", "Adult")))
    
    # 4. 转换为宽表
    df_wide <- df_long_parsed %>%
      pivot_wider(names_from = Age, values_from = PSI_avg, values_fill = list(PSI_avg = NA))
    
    # 检查哪些年龄阶段列实际存在
    age_columns <- intersect(c("Egg", "Larva", "Pupa", "Adult"), names(df_wide))
    
    # 新增：统计组织数量和各组织事件数
    tissue_counts <- df_wide %>% 
      group_by(Tissue) %>% 
      summarise(total_events = n(), .groups = "drop")
    
    message(paste0("[INFO] Species ", species, " has ", nrow(tissue_counts), 
                   " tissues: ", paste(tissue_counts$Tissue, collapse = ", ")))
    
    # 打印每个组织的事件数
    for(i in 1:nrow(tissue_counts)) {
      message(paste0("[INFO] -- Tissue ", tissue_counts$Tissue[i], 
                     " has ", tissue_counts$total_events[i], " total events"))
    }
    
    if(length(age_columns) < 2) {
      message(paste0("[WARN] Only ", length(age_columns), 
                     " age stages available for species: ", species))
      next
    }
    
    # 5. 计算非缺失阶段数
    df_wide <- df_wide %>%
      mutate(non_missing = rowSums(!is.na(select(., all_of(age_columns)))))
    
    min_required_stages <- min(2, length(age_columns))
    df_wide <- df_wide %>%
      filter(non_missing >= min_required_stages)
    
    # 新增：统计过滤后各组织的事件数
    filtered_tissue_counts <- df_wide %>% 
      group_by(Tissue) %>% 
      summarise(filtered_events = n(), .groups = "drop")
    
    message(paste0("[INFO] After filtering (>= ", min_required_stages, 
                   " stages), species ", species, " retains:"))
    
    for(i in 1:nrow(filtered_tissue_counts)) {
      message(paste0("[INFO] -- Tissue ", filtered_tissue_counts$Tissue[i], 
                     ": ", filtered_tissue_counts$filtered_events[i], 
                     " events (", round(filtered_tissue_counts$filtered_events[i]/tissue_counts$total_events[i]*100, 1),
                     "% of original)"))
    }
    
    if (nrow(df_wide) == 0) {
      message(paste0("[WARN] No valid data after filtering for species: ", species))
      next
    }
    
    # 6. 计算方向
    df_wide <- df_wide %>%
      rowwise() %>%
      mutate(
        stage_values = list(unlist(across(all_of(age_columns)))),
        stage_values = list(stage_values[!is.na(stage_values)]),
        diffs = list(if(length(stage_values) > 1) diff(stage_values) else numeric(0)),
        direction = judge_direction(unlist(diffs), threshold = 0.05)
      ) %>%
      ungroup()
    
    # 添加物种信息
    df_wide <- df_wide %>%
      mutate(Species = species)
    
    # 存储结果
    result_list[[species]] <- df_wide
    
  }, error = function(e) {
    message(paste0("[ERROR] Processing species ", species, ": ", e$message))
  })
}

# 合并所有结果
all_species_results <- bind_rows(result_list)

# 最终统计输出
message("\n[SUMMARY] Final statistics:")
message(paste0("Total processed species: ", length(species_names)))
message(paste0("Species with valid data: ", length(unique(all_species_results$Species))))

# 按物种和组织统计事件数
final_summary <- all_species_results %>%
  group_by(Species, Tissue) %>%
  summarise(events = n(), .groups = "drop")

message("\n[SUMMARY] Events per species and tissue:")
print(final_summary, n = Inf)

# 保存结果
saveRDS(all_species_results, "Fout_time_all_species_results.RDS")


# 加载必要包（若未安装需先执行 install.packages("dplyr")）
library(dplyr)
# 筛选微外显子事件
micro_exon_results <- all_species_results %>%
  filter(sapply(event_id, function(x) {
    parts <- strsplit(x, "__")[[1]]
    if(length(parts) >= 7 && parts[3] %in% c("+", "-")) {
      exon_length <- abs(as.numeric(parts[7]) - as.numeric(parts[6])) + 1
      return(exon_length < 28)
    }
    FALSE
  }))
# 验证结果
dim(micro_exon_results)  # 显示筛选后数据维度
head(micro_exon_results) # 查看前6条记录
saveRDS(micro_exon_results, "Fout_time_all_species_results_micro.RDS")


