# 载入必要的包
if(!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

# 1. 函数：使用向量化操作计算PSI字符串的平均值
calculate_psi_mean_vectorized <- function(psi_string) {
  # 检查是否为NA
  if(is.na(psi_string)) return(NA_real_)
  
  # 分割字符串并转换为数值
  values <- tryCatch({
    # 使用strsplit分割字符串
    split_values <- strsplit(psi_string, ",", fixed = TRUE)[[1]]
    
    # 过滤掉"NA"字符串
    numeric_values <- as.numeric(split_values[split_values != "NA"])
    
    # 如果没有有效值，返回NA
    if(length(numeric_values) == 0) return(NA_real_)
    
    # 返回平均值
    mean(numeric_values)
  }, error = function(e) {
    # 如果处理过程中出错，返回NA
    NA_real_
  })
  
  return(values)
}

# 2. 函数：获取所有目标物种中唯一的PSI条件名称
get_all_psi_conditions <- function(se_events_list, target_species) {
  all_conditions <- character(0)
  for (species in target_species) {
    if (species %in% names(se_events_list)) {
      species_se <- se_events_list[[species]]
      # 假设PSI列从第11列开始
      if (ncol(species_se) > 10) {
        psi_cols <- names(species_se)[11:ncol(species_se)]
        all_conditions <- c(all_conditions, psi_cols)
      }
    }
  }
  unique_conditions <- unique(all_conditions)
  cat("  Found", length(unique_conditions), "unique PSI conditions across target species.\n")
  return(unique_conditions)
}


# 3. 函数：从SE事件数据中提取PSI信息，仅针对目标物种，使用通用条件名
extract_psi_info_optimized <- function(se_events_list, target_species) {
  # 创建一个列表来存储目标物种的PSI数据
  psi_data_list <- list()
  
  # 处理每个目标物种
  cat("Extracting PSI data for target species:", paste(target_species, collapse=", "), "\n")
  for(species in target_species) {
    # 检查该物种是否存在于se_events_list中
    if (species %in% names(se_events_list)) {
      cat("  Processing PSI for species:", species, "\n")
      
      # 转换为data.table
      species_se <- as.data.table(se_events_list[[species]])
      
      # 找出PSI列（除去前10列）
      psi_cols <- character(0)
      if (ncol(species_se) > 10) {
        psi_cols <- names(species_se)[11:ncol(species_se)]
      }
      
      # 如果没有PSI列，创建一个空的占位符并跳过计算
      if(length(psi_cols) == 0) {
        cat("    No PSI columns found for species:", species, "\n")
        psi_dt <- data.table(
          species = species,
          GeneID = species_se$GeneID,
          exonStart_0base = species_se$exonStart_0base,
          exonEnd = species_se$exonEnd,
          upstreamES = species_se$upstreamES,
          upstreamEE = species_se$upstreamEE,
          downstreamES = species_se$downstreamES,
          downstreamEE = species_se$downstreamEE
        )
        psi_data_list[[species]] <- psi_dt
        next
      }
      
      # 创建一个新的数据表，包含基因ID和坐标
      psi_dt <- data.table(
        species = species,
        GeneID = species_se$GeneID,
        exonStart_0base = species_se$exonStart_0base,
        exonEnd = species_se$exonEnd,
        upstreamES = species_se$upstreamES,
        upstreamEE = species_se$upstreamEE,
        downstreamES = species_se$downstreamES,
        downstreamEE = species_se$downstreamEE
      )
      
      # 使用lapply高效计算每个PSI列的平均值
      # 列名直接使用原始的条件名（例如 "Pupa-Pupa"）
      for(col in psi_cols) {
        set(psi_dt, j = col, 
            value = vapply(species_se[[col]], calculate_psi_mean_vectorized, numeric(1)))
      }
      
      # 添加到列表
      psi_data_list[[species]] <- psi_dt
      cat("    Processed", length(psi_cols), "PSI columns for", nrow(psi_dt), "events\n")
      
    } else {
      cat("  Species", species, "not found in SE events list. Skipping PSI calculation.\n")
      # 可以选择添加一个空的占位符，以便后续合并时不会出错
      psi_data_list[[species]] <- NULL 
    }
  }
  
  return(psi_data_list)
}

# 4. 函数：更高效地读取和处理原始正交组数据
read_ortholog_groups_alt <- function(input_file) {
  # 读取原始数据
  data <- fread(input_file, header = TRUE, sep = "\t")
  
  # 初始化结果data.table
  result_dt <- data.table(
    ortholog_group = character(),
    species_count = integer(),
    species = character(),
    gene_id = character(),
    transcript_id = character(),
    start = integer(),
    end = integer()
  )
  
  # 预分配内存以提高性能
  total_entries <- sum(sapply(data$ortholog_group, function(x) length(unlist(strsplit(x, ",")))))
  result_dt <- data.table(
    ortholog_group = character(total_entries),
    species_count = integer(total_entries),
    species = character(total_entries),
    gene_id = character(total_entries),
    transcript_id = character(total_entries),
    start = integer(total_entries),
    end = integer(total_entries)
  )
  
  # 使用计数器跟踪当前行
  counter <- 0
  
  # 处理每一行，避免多次rbindlist
  for (i in 1:nrow(data)) {
    # 生成唯一OG ID
    og_id <- paste0("OG", i)
    species_count <- data$species_count[i]
    og_string <- data$ortholog_group[i]
    
    # 分割ortholog_group字符串
    entries <- unlist(strsplit(og_string, ","))
    
    # 提取每个条目的信息
    for (entry in entries) {
      # 使用正则表达式提取信息
      parts <- strsplit(entry, "::", fixed = TRUE)[[1]]
      if (length(parts) == 3) {
        species <- parts[1]
        gene_id <- parts[2]
        
        # 分割最后一部分以获取transcript_id和坐标
        transcript_parts <- strsplit(parts[3], "_pos", fixed = TRUE)[[1]]
        if (length(transcript_parts) == 2) {
          transcript_id <- transcript_parts[1]
          
          # 分割坐标
          coords <- strsplit(transcript_parts[2], "-", fixed = TRUE)[[1]]
          if (length(coords) == 2) {
            start <- as.integer(coords[1])
            end <- as.integer(coords[2])
            
            # 增加计数器
            counter <- counter + 1
            
            # 直接设置值，而不是使用rbindlist
            set(result_dt, counter, "ortholog_group", og_id)
            set(result_dt, counter, "species_count", species_count)
            set(result_dt, counter, "species", species)
            set(result_dt, counter, "gene_id", gene_id)
            set(result_dt, counter, "transcript_id", transcript_id)
            set(result_dt, counter, "start", start)
            set(result_dt, counter, "end", end)
          }
        }
      }
    }
    
    # 定期报告进度
    if (i %% 1000 == 0 || i == nrow(data)) {
      cat(sprintf("Processed %d/%d rows\n", i, nrow(data)))
    }
  }
  
  # 如果计数器小于预分配的大小，需要裁剪
  if(counter < total_entries) {
    result_dt <- result_dt[1:counter]
  }
  
  return(result_dt)
}

# 5. 函数：使用data.table优化匹配正交组和SE事件，并添加通用PSI信息
optimize_with_data_table_and_psi <- function(all_groups_dt, se_events_list, psi_data_list, all_psi_conditions, tolerance = 10) {
  # 创建结果列表
  results_list <- list()
  
  # 只处理在正交组数据和SE事件列表中都存在的物种
  valid_species <- intersect(unique(all_groups_dt$species), names(se_events_list))
  
  # 对每个有效物种单独处理
  for(species in valid_species) {
    cat("Processing species:", species, "\n")
    
    # 提取当前物种的正交组数据
    species_ogs <- all_groups_dt[species == species]
    
    # 提取当前物种的SE事件并转换为data.table
    species_se <- as.data.table(se_events_list[[species]])
    
    # 获取当前物种的PSI数据 (如果存在于 psi_data_list 中)
    species_psi <- if(species %in% names(psi_data_list)) psi_data_list[[species]] else NULL
    
    # 确保坐标是数值型
    species_ogs[, start := as.numeric(start)]
    species_ogs[, end := as.numeric(end)]
    species_se[, exonStart_0base := as.numeric(exonStart_0base)]
    species_se[, exonEnd := as.numeric(exonEnd)]
    species_se[, upstreamES := as.numeric(upstreamES)]
    species_se[, upstreamEE := as.numeric(upstreamEE)]
    species_se[, downstreamES := as.numeric(downstreamES)]
    species_se[, downstreamEE := as.numeric(downstreamEE)]
    
    # 基因ID过滤（如果适用）
    if(all(c("gene_id", "GeneID") %in% c(names(species_ogs), names(species_se)))) {
      common_genes <- intersect(species_ogs$gene_id, species_se$GeneID)
      if(length(common_genes) > 0) {
        species_ogs_filtered <- species_ogs[gene_id %in% common_genes]
        species_se_filtered <- species_se[GeneID %in% common_genes]
        
        # 只有当过滤后的数据集不为空时才使用过滤后的数据
        if(nrow(species_ogs_filtered) > 0 && nrow(species_se_filtered) > 0) {
          species_ogs <- species_ogs_filtered
          species_se <- species_se_filtered
        }
      }
    }
    
    # 准备结果数据表
    matches_dt <- data.table()
    
    # 批量处理以提高效率，同时减轻内存负担
    batch_size <- 500  # 可以根据物种数据量调整
    total_batches <- ceiling(nrow(species_ogs) / batch_size)
    
    for(batch_idx in 1:total_batches) {
      start_row <- (batch_idx - 1) * batch_size + 1
      end_row <- min(batch_idx * batch_size, nrow(species_ogs))
      
      batch_ogs <- species_ogs[start_row:end_row]
      cat(sprintf("  Processing batch %d/%d (rows %d-%d of %d)\n", 
                  batch_idx, total_batches, start_row, end_row, nrow(species_ogs)))
      
      # 预先计算所有可能的匹配
      og_matches <- data.table(
        og_idx = integer(),
        se_idx = integer(),
        match_type = character()
      )
      
      # 为每个SE事件创建索引范围
      se_indices <- data.table(
        se_idx = 1:nrow(species_se),
        # exon区域
        exon_start_min = species_se$exonStart_0base - tolerance,
        exon_start_max = species_se$exonStart_0base + tolerance,
        exon_end_min = species_se$exonEnd - tolerance,
        exon_end_max = species_se$exonEnd + tolerance,
        # upstream区域
        up_start_min = species_se$upstreamES - tolerance,
        up_start_max = species_se$upstreamES + tolerance,
        up_end_min = species_se$upstreamEE - tolerance,
        up_end_max = species_se$upstreamEE + tolerance,
        # downstream区域
        down_start_min = species_se$downstreamES - tolerance,
        down_start_max = species_se$downstreamES + tolerance,
        down_end_min = species_se$downstreamEE - tolerance,
        down_end_max = species_se$downstreamEE + tolerance
      )
      
      # 使用向量化操作批量查找匹配
      for(i in 1:nrow(batch_ogs)) {
        og_row <- batch_ogs[i]
        og_start <- og_row$start
        og_end <- og_row$end
        
        exon_matches <- which(
          og_start >= se_indices$exon_start_min & 
            og_start <= se_indices$exon_start_max &
            og_end >= se_indices$exon_end_min & 
            og_end <= se_indices$exon_end_max
        )
        if(length(exon_matches) > 0) {
          og_matches <- rbindlist(list(og_matches, data.table(og_idx = i, se_idx = exon_matches, match_type = "match")))
        }
        
        upstream_matches <- which(
          og_start >= se_indices$up_start_min & 
            og_start <= se_indices$up_start_max &
            og_end >= se_indices$up_end_min & 
            og_end <= se_indices$up_end_max
        )
        if(length(upstream_matches) > 0) {
          og_matches <- rbindlist(list(og_matches, data.table(og_idx = i, se_idx = upstream_matches, match_type = "upstream")))
        }
        
        downstream_matches <- which(
          og_start >= se_indices$down_start_min & 
            og_start <= se_indices$down_start_max &
            og_end >= se_indices$down_end_min & 
            og_end <= se_indices$down_end_max
        )
        if(length(downstream_matches) > 0) {
          og_matches <- rbindlist(list(og_matches, data.table(og_idx = i, se_idx = downstream_matches, match_type = "downstream")))
        }
      }
      
      # 如果找到匹配，生成结果
      if(nrow(og_matches) > 0) {
        # 为每个匹配创建结果行
        batch_results <- list() # Use a list to collect results efficiently
        
        for(i in 1:nrow(og_matches)) {
          og_idx <- og_matches$og_idx[i]
          se_idx <- og_matches$se_idx[i]
          match_type <- og_matches$match_type[i]
          
          og_row <- batch_ogs[og_idx]
          se_row <- species_se[se_idx]
          
          # 创建基本匹配结果 (as a list first)
          match_row_list <- list(
            ortholog_group = og_row$ortholog_group,
            species_count = og_row$species_count,
            species = species,
            gene_id = og_row$gene_id,
            transcript_id = og_row$transcript_id,
            start = og_row$start,
            end = og_row$end,
            match_type = match_type,
            SE_gene_id = se_row$GeneID,
            SE_chr = se_row$chr,
            SE_strand = se_row$strand,
            SE_exon_start = se_row$exonStart_0base,
            SE_exon_end = se_row$exonEnd,
            SE_upstream_start = se_row$upstreamES,
            SE_upstream_end = se_row$upstreamEE,
            SE_downstream_start = se_row$downstreamES,
            SE_downstream_end = se_row$downstreamEE
          )
          
          # 查找对应的PSI信息 (仅当species_psi存在时)
          if (!is.null(species_psi)) {
            psi_key <- species_psi[GeneID == se_row$GeneID & 
                                     exonStart_0base == se_row$exonStart_0base & 
                                     exonEnd == se_row$exonEnd]
            
            # 如果找到PSI信息，添加到结果
            if(nrow(psi_key) > 0) {
              # 获取当前物种的PSI列
              species_psi_cols <- setdiff(names(psi_key), 
                                          c("species", "GeneID", "exonStart_0base", "exonEnd", 
                                            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"))
              
              # 遍历所有全局PSI条件
              for(condition in all_psi_conditions) {
                # 生成通用PSI列名
                psi_col_name <- paste0(condition, "_PSI")
                # 如果当前物种有这个条件的数据，则添加值，否则为NA
                if(condition %in% species_psi_cols) {
                  match_row_list[[psi_col_name]] <- psi_key[[condition]][1]
                } else {
                  match_row_list[[psi_col_name]] <- NA_real_
                }
              }
            } else {
              # 如果没找到PSI信息，为所有通用PSI列添加NA
              for(condition in all_psi_conditions) {
                psi_col_name <- paste0(condition, "_PSI")
                match_row_list[[psi_col_name]] <- NA_real_
              }
            }
          } else {
            # 如果该物种没有PSI数据，为所有通用PSI列添加NA
            for(condition in all_psi_conditions) {
              psi_col_name <- paste0(condition, "_PSI")
              match_row_list[[psi_col_name]] <- NA_real_
            }
          }
          
          # 添加到批次结果列表
          batch_results[[length(batch_results) + 1]] <- match_row_list
        }
        
        # 将批次结果列表高效合并为data.table
        batch_dt <- rbindlist(batch_results, use.names = TRUE, fill = TRUE)
        
        # 添加到总结果
        matches_dt <- rbindlist(list(matches_dt, batch_dt), use.names = TRUE, fill = TRUE)
      }
      
      # 清理内存
      rm(batch_ogs, og_matches, se_indices, batch_results)
      gc()
    }
    
    # 保存当前物种的结果
    if(nrow(matches_dt) > 0) {
      results_list[[species]] <- as.data.frame(matches_dt)
    }
  }
  
  # 合并所有物种的结果
  if(length(results_list) > 0) {
    final_results <- rbindlist(lapply(results_list, as.data.table), use.names = TRUE, fill = TRUE)
    return(as.data.frame(final_results))
  } else {
    return(data.frame())
  }
}

# 6. 函数：统计和分析匹配结果
analyze_matches <- function(match_results) {
  if(is.null(match_results) || nrow(match_results) == 0) {
    cat("No matches found.\n")
    return(NULL)
  }
  
  # 确保使用data.table处理
  match_dt <- as.data.table(match_results)
  
  # 基本统计
  cat("Total matches found:", nrow(match_dt), "\n")
  
  # 按匹配类型统计
  match_type_counts <- match_dt[, .N, by = match_type]
  cat("\nMatches by type:\n")
  print(match_type_counts)
  
  # 按物种统计
  species_counts <- match_dt[, .N, by = species]
  cat("\nMatches by species:\n")
  print(species_counts)
  
  # 按物种计数统计
  species_count_stats <- match_dt[, .N, by = species_count]
  cat("\nMatches by species count:\n")
  print(species_count_stats)
  
  # 按正交组统计
  og_counts <- match_dt[, .N, by = ortholog_group]
  setorder(og_counts, -N)
  cat("\nTop 10 ortholog groups with most matches:\n")
  print(head(og_counts, 10))
  
  # PSI值的基本统计（如果存在）
  # 查找通用PSI列名
  psi_cols <- names(match_dt)[endsWith(names(match_dt), "_PSI")]
  if(length(psi_cols) > 0) {
    cat("\nPSI value statistics:\n")
    # 使用data.table高效计算PSI统计
    psi_stats_list <- list()
    for(col in psi_cols) {
      stats <- match_dt[, .(
        mean = mean(get(col), na.rm = TRUE),
        median = median(get(col), na.rm = TRUE),
        min = min(get(col), na.rm = TRUE),
        max = max(get(col), na.rm = TRUE),
        na_count = sum(is.na(get(col)))
      )]
      cat("  ", col, ":\n")
      cat("    Mean:", stats$mean, "\n")
      cat("    Median:", stats$median, "\n")
      cat("    Min:", stats$min, "\n")
      cat("    Max:", stats$max, "\n")
      cat("    NA count:", stats$na_count, "\n")
      psi_stats_list[[col]] <- stats
    }
    psi_stats <- rbindlist(psi_stats_list, idcol = "condition")
  } else {
    psi_stats <- NULL
  }
  
  # 检查是否有同一正交组中同一物种有多个基因的情况
  og_species_genes <- match_dt[, .(gene_count = uniqueN(gene_id)), by = .(ortholog_group, species)]
  multi_gene_groups <- og_species_genes[gene_count > 1]
  
  if(nrow(multi_gene_groups) > 0) {
    cat("\nWARNING: Found ortholog groups with multiple genes from the same species!\n")
    cat("This should not happen in 1:1 ortholog groups.\n")
    cat("Number of such cases:", nrow(multi_gene_groups), "\n")
    cat("Sample of problematic groups:\n")
    print(head(multi_gene_groups, 10))
  } else {
    cat("\nAll ortholog groups maintain 1:1 relationship (no multiple genes from same species).\n")
  }
  
  # 返回分析结果
  return(list(
    match_type_counts = match_type_counts,
    species_counts = species_counts,
    species_count_stats = species_count_stats,
    og_counts = og_counts,
    multi_gene_groups = multi_gene_groups,
    psi_stats = psi_stats
  ))
}

# 7. 函数：保存结果
save_results <- function(match_results, file_prefix = "orthogroup_SE_matches") {
  if(is.null(match_results) || nrow(match_results) == 0) {
    cat("No matches to save.\n")
    return(NULL)
  }
  
  # 保存为CSV
  csv_file <- paste0(file_prefix, ".csv")
  fwrite(match_results, csv_file)  # 使用data.table的fwrite代替write.csv
  cat("Results saved to", csv_file, "\n")
  
  # 保存为RDS（R数据格式，更高效读取）
  rds_file <- paste0(file_prefix, ".rds")
  saveRDS(match_results, rds_file)
  cat("Results saved to", rds_file, "\n")
  
  # 如果安装了arrow包，也可以保存为parquet格式（高效存储和读取）
  if(require(arrow)) {
    parquet_file <- paste0(file_prefix, ".parquet")
    write_parquet(match_results, parquet_file)
    cat("Results saved to", parquet_file, "\n")
  }
  
  return(TRUE)
}

# 8. 主执行函数：整合所有步骤，包括通用PSI信息
run_full_analysis_with_psi <- function(ortholog_file, se_events_list, tolerance = 10, output_prefix = NULL) {
  # 设置默认输出前缀
  if(is.null(output_prefix)) {
    output_prefix <- "orthogroup_SE_matches_with_psi"
  }
  
  # 开始计时
  start_time <- Sys.time()
  cat("Analysis started at:", format(start_time), "\n")
  
  # 步骤1：读取和处理正交组数据
  cat("Step 1: Reading and processing ortholog groups data...\n")
  all_groups_dt <- read_ortholog_groups_alt(ortholog_file)
  cat("  Processed", nrow(all_groups_dt), "entries from", length(unique(all_groups_dt$ortholog_group)), "ortholog groups\n")
  
  # 获取目标物种列表 (来自正交组数据)
  target_species <- unique(all_groups_dt$species)
  cat("  Target species:", paste(target_species, collapse=", "), "\n")
  
  # 步骤2：获取所有PSI条件名称
  cat("\nStep 2: Identifying all unique PSI conditions...\n")
  all_psi_conditions <- get_all_psi_conditions(se_events_list, target_species)
  
  # 步骤3：提取PSI信息 (仅针对目标物种，使用通用条件名)
  cat("\nStep 3: Extracting PSI information for target species...\n")
  psi_data_list <- extract_psi_info_optimized(se_events_list, target_species)
  cat("  Extracted PSI data for", length(psi_data_list), "species\n")
  
  # 步骤4：匹配正交组与SE事件，并添加通用PSI信息
  cat("\nStep 4: Matching ortholog groups to SE events and adding PSI information...\n")
  match_results <- optimize_with_data_table_and_psi(all_groups_dt, se_events_list, psi_data_list, all_psi_conditions, tolerance)
  
  # 步骤5：分析匹配结果
  cat("\nStep 5: Analyzing match results...\n")
  analysis <- analyze_matches(match_results)
  
  # 步骤6：保存结果
  cat("\nStep 6: Saving results...\n")
  save_results(match_results, output_prefix)
  
  # 结束计时并显示总用时
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  cat("\nAnalysis completed in:", format(elapsed), "\n")
  
  # 返回结果和分析
  return(list(
    match_results = match_results,
    analysis = analysis,
    runtime = elapsed
  ))
}

# 9. 使用示例
# 替换为您的实际文件路径和数据
ortholog_file <- "one_to_one_orthologs.txt"
se_events_list <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]]

# 执行完整分析，包括PSI信息
results <- run_full_analysis_with_psi(ortholog_file, se_events_list)

# 查看匹配结果
head(results$match_results)

# 查看PSI统计
results$analysis$psi_stats


matched_rows <- subset(results$match_results, match_type == "match")

saveRDS(results$match_results,"analyse_result.RDS")
saveRDS(matched_rows,"analyse_match_result.RDS")