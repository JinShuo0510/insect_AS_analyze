
# -------------------- 1. 读取 Orthogroup 文件 --------------------
# 假设文件路径为 "orthogroup_file.txt"
og <- read.delim("important_species_orthogroups_at_lest_7_species.txt", header = TRUE, stringsAsFactors = FALSE)

# -------------------- 2. 提取物种名称 --------------------
# 除去第一列"Orthogroup"，其他列名去除后缀".protein"
species_names <- sub("\\.protein$", "", colnames(og)[-1])
# 查看提取的物种名称
print(species_names)

# -------------------- 3. 提取各物种对应的 GeneID --------------------
# 对于每个物种（即每个列），拆分逗号分隔的多个基因ID，并整理成列表
gene_ids_list <- list()
for(i in seq_along(species_names)) {
  sp <- species_names[i]
  # 对应列下标为 i+1（因为第一列是 Orthogroup）
  col_values <- og[[i + 1]]
  
  # 对每个值拆分逗号，去除多余空白，并合并为一个向量
  gene_ids <- unlist(lapply(col_values, function(x) {
    # 去除首尾空白
    x <- trimws(x)
    # 如果非空则拆分，否则返回 NA（可选）
    if(nzchar(x)) {
      ids <- unlist(strsplit(x, ",\\s*"))
      trimws(ids)
    } else {
      NULL
    }
  }))
  
  # 去重并剔除空字符
  gene_ids <- unique(gene_ids[gene_ids != ""])
  gene_ids_list[[sp]] <- gene_ids
}
# 检查结果
str(gene_ids_list)

#########################################
# 2. 利用 GTF 文件建立 transcript->gene 映射 #
#########################################
# 准备结果列表，用于保存筛选后的 R 变量数据
results <- list()

for(sp in species_names) {
  # 构建当前物种的 GTF 文件名，格式为 "{species}.exon.gtf"
  gtf_file <- paste0(sp, ".exon.gtf")
  if(!file.exists(gtf_file)) {
    message(sprintf("警告：文件 %s 不存在！", gtf_file))
    next
  }
  
  # 读取 GTF 文件（过滤注释行，以 tab 分隔，字段无表头）
  gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
  # GTF 文件标准字段：V1=seqname, V2=source, V3=feature, V4=start, V5=end,
  # V6=score, V7=strand, V8=frame, V9=attributes
  
  # 从属性列中提取 transcript_id 与 gene_id
  # 此处使用正则表达式分别提取属性中的 transcript_id 和 gene_id
  gtf$transcript_id <- sub(".*transcript_id\\s+([^;]+);.*", "\\1", gtf$V9)
  gtf$gene_id       <- sub(".*gene_id\\s+([^;]+);.*", "\\1", gtf$V9)
  
  # 取唯一映射
  mapping <- unique(gtf[, c("transcript_id", "gene_id")])
  
  # 获取 orthogroup 文件中该物种的转录本 ID 列表
  transcripts <- gene_ids_list[[sp]]
  
  # 根据 transcript_id 映射得到对应的 gene_id（有可能有部分转录本在 GTF 中未找到）
  gene_ids_mapped <- mapping$gene_id[mapping$transcript_id %in% transcripts]
  gene_ids_mapped <- unique(gene_ids_mapped)
  
  ###############################
  # 3. 利用 gene_id 筛选 R 变量数据 #
  ###############################
  # 检查 all_event_wide_lists_withoutavg_nondeNA_1_8 中是否包含该物种数据
  if(sp %in% names(all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]])) {
    sp_df <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]][[sp]]
    # 假设 sp_df 中基因的标识在列 "GeneID"
    subset_df <- sp_df[sp_df$GeneID %in% gene_ids_mapped, ]
    results[[sp]] <- subset_df
  } else {
    message(sprintf("警告：R变量中未找到物种 %s 的数据", sp))
  }
}

# -------------------- 5. 输出结果 --------------------
# 方案1：结果为列表，每个物种对应一个 data.frame
results


library(dplyr)

combined_results <- bind_rows(results, .id = "Species")

# 查看合并后的结果
print(combined_results)

saveRDS(combined_results,file = "combined_results_important_at7.RDS")
saveRDS(results,file = "results_important_at7.RDS")

