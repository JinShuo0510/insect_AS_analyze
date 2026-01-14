# -------------------- 1. 读取 Orthogroup 文件 --------------------
og <- read.delim("important_species_orthogroups_at_lest_7_species.txt", 
                 header = TRUE, stringsAsFactors = FALSE)

# -------------------- 2. 提取物种名称 --------------------
# 除去第一列"Orthogroup"，其他列名去除后缀".protein"
species_names <- sub("\\.protein$", "", colnames(og)[-1])
print(species_names)

# -------------------- 3. 预加载各物种的 GTF 映射（transcript -> gene） --------------------
mapping_list <- list()
for (sp in species_names) {
  gtf_file <- paste0(sp, ".exon.gtf")
  if (!file.exists(gtf_file)) {
    message(sprintf("警告：文件 %s 不存在！", gtf_file))
    next
  }
  gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#", 
                    sep = "\t", stringsAsFactors = FALSE)
  gtf$transcript_id <- sub(".*transcript_id\\s+([^;]+);.*", "\\1", gtf$V9)
  gtf$gene_id       <- sub(".*gene_id\\s+([^;]+);.*", "\\1", gtf$V9)
  mapping <- unique(gtf[, c("transcript_id", "gene_id")])
  mapping_list[[sp]] <- mapping
}

# -------------------- 4. 按 Orthogroup 生成结果变量 --------------------
results_by_orthogroup <- list()

for (i in 1:nrow(og)) {
  og_id <- og[i, "Orthogroup"]
  orthogroup_results <- list()  # 用于存储该 orthogroup 下各物种的筛选结果
  
  # 针对每个物种
  for (sp in species_names) {
    # 这里列名假设为物种名，如若列名带后缀可调整为 paste0(sp, ".protein")
    col_name <- sp  
    transcript_str <- og[i, col_name]
    
    # 检查是否为 NA 或空字符串
    if (length(transcript_str) == 0 || is.na(transcript_str)) next
    transcript_str <- trimws(transcript_str)
    if (!nzchar(transcript_str)) next
    
    # 拆分逗号得到 transcript IDs
    transcripts <- unique(unlist(strsplit(transcript_str, ",\\s*")))
    
    # 利用预加载的 mapping 将 transcript 映射为 gene IDs
    if (!is.null(mapping_list[[sp]])) {
      mapping <- mapping_list[[sp]]
      gene_ids_mapped <- unique(mapping$gene_id[mapping$transcript_id %in% transcripts])
      
      # 检查 all_event_wide_lists_withoutavg_nondeNA_1_8 中是否包含该物种数据（假定在 "SE" 层级中）
      if (sp %in% names(all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]])) {
        sp_df <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]][[sp]]
        subset_df <- sp_df[sp_df$GeneID %in% gene_ids_mapped, ]
        if (nrow(subset_df) > 0) {
          # 添加物种信息列（确保非合并结果中每条数据都能知道所属物种）
          subset_df$Species <- sp  
          subset_df <- subset_df[, c("Species", setdiff(names(subset_df), "Species"))]
          orthogroup_results[[sp]] <- subset_df
        }
      } else {
        message(sprintf("警告：R变量中未找到物种 %s 的数据", sp))
      }
    }
  }
  
  # 如果该 orthogroup 下有至少一个物种的结果，则合并存入列表
  if (length(orthogroup_results) > 0) {
    results_by_orthogroup[[og_id]] <- dplyr::bind_rows(orthogroup_results)
  } else {
    results_by_orthogroup[[og_id]] <- NULL
  }
}

# -------------------- 5. 查看结果 --------------------
# 此时 results_by_orthogroup 列表中每个数据框都包含 "Species" 列，
# 你可以通过以下方式检查某个 orthogroup 的结果
str(results_by_orthogroup[["OG0000969"]])

# -------------------- 6. 合并并保存结果 --------------------
combined_results_by_og <- dplyr::bind_rows(results_by_orthogroup, .id = "Orthogroup")
print(combined_results_by_og)

saveRDS(combined_results_by_og, file = "combined_results_important_at7_byOrthogroup.RDS")
saveRDS(results_by_orthogroup, file = "results_important_at7_byOrthogroup.RDS")
