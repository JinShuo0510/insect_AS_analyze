# 加载必要的包
library(stringr)
library(progress)  # 用于显示进度条
library(data.table)

# 函数：格式化时间差
format_time <- function(start_time) {
  time_diff <- difftime(Sys.time(), start_time, units = "secs")
  sprintf("%.2f seconds", as.numeric(time_diff))
}

# 记录开始时间
total_start_time <- Sys.time()
cat(sprintf("\n[%s] 开始处理数据...\n", format(total_start_time, "%Y-%m-%d %H:%M:%S")))

# 初始化空列表存储结果
head_data_list <- list()

# 获取物种总数用于进度显示
total_species <- length(names(all_species_tissue_data_1_8))
cat(sprintf("[%s] 共发现 %d 个物种需要处理\n", 
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
            total_species))

# 创建进度条
pb <- progress_bar$new(
  format = "物种处理进度 [:bar] :percent eta: :eta",
  total = total_species,
  clear = FALSE,
  width = 60
)

# 遍历每个物种
for (species in names(all_species_tissue_data_1_8)) {
  species_data <- all_species_tissue_data_1_8[[species]]
  if ("Head" %in% names(species_data)) {
    head_data <- species_data[["Head"]]
    if ("Event_ID" %in% names(head_data) && length(head_data$Event_ID) > 0) {
      event_ids <- head_data$Event_ID
      head_df <- data.frame(
        Species = species,
        Event_ID = event_ids,
        stringsAsFactors = FALSE
      )
      head_data_list[[species]] <- head_df
    }
  }
  pb$tick()
}

# 合并所有有效数据框
merge_start_time <- Sys.time()
cat(sprintf("\n[%s] 开始合并数据框...\n", format(merge_start_time, "%Y-%m-%d %H:%M:%S")))

final_head_table <- do.call(rbind, head_data_list)
total_rows <- nrow(final_head_table)

cat(sprintf("[%s] 数据框合并完成，共 %d 行数据\n", 
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
            total_rows))

# 向量化提取外显子信息的函数
extract_exon_info_vectorized <- function(event_ids, chunk_size = 100000) {
  cat(sprintf("[%s] 开始处理外显子信息，总数据量: %d\n", 
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
              length(event_ids)))
  
  # 计算需要处理的块数
  n_chunks <- ceiling(length(event_ids) / chunk_size)
  
  # 创建进度条
  pb <- progress_bar$new(
    format = "处理外显子信息 [:bar] :percent eta: :eta",
    total = n_chunks,
    clear = FALSE,
    width = 60
  )
  
  # 初始化结果列表
  result_list <- vector("list", n_chunks)
  
  # 分块处理数据
  for (i in 1:n_chunks) {
    start_idx <- ((i-1) * chunk_size) + 1
    end_idx <- min(i * chunk_size, length(event_ids))
    chunk_ids <- event_ids[start_idx:end_idx]
    
    # 处理当前块
    parts_list <- str_split(chunk_ids, "_")
    direction_positions <- sapply(parts_list, function(parts) {
      which(parts %in% c("+", "-"))[1]
    })
    
    n <- length(chunk_ids)
    exon_start <- rep(NA_real_, n)
    exon_end <- rep(NA_real_, n)
    is_microexon <- rep(NA, n)
    
    valid_idx <- !is.na(direction_positions)
    if (any(valid_idx)) {
      valid_parts <- parts_list[valid_idx]
      valid_dir_pos <- direction_positions[valid_idx]
      
      exon_start[valid_idx] <- as.numeric(sapply(seq_along(valid_parts), function(i) {
        valid_parts[[i]][valid_dir_pos[i] + 1]
      }))
      exon_end[valid_idx] <- as.numeric(sapply(seq_along(valid_parts), function(i) {
        valid_parts[[i]][valid_dir_pos[i] + 2]
      }))
      
      exon_length <- abs(exon_end[valid_idx] - exon_start[valid_idx]) + 1
      is_microexon[valid_idx] <- exon_length < 28
    }
    
    result_list[[i]] <- data.frame(
      exon_start = exon_start,
      exon_end = exon_end,
      is_microexon = is_microexon
    )
    
    pb$tick()
  }
  
  # 合并所有结果
  cat(sprintf("[%s] 合并处理结果...\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  final_result <- do.call(rbind, result_list)
  return(final_result)
}

# 提取外显子信息
exon_start_time <- Sys.time()
cat(sprintf("\n[%s] 开始提取外显子信息...\n", format(exon_start_time, "%Y-%m-%d %H:%M:%S")))

# 使用较大的chunk_size来处理数据
exon_info_df <- extract_exon_info_vectorized(final_head_table$Event_ID, chunk_size = 100000)

# 合并到原数据框
merge_exon_time <- Sys.time()
cat(sprintf("\n[%s] 开始合并外显子信息到原数据框...\n", format(merge_exon_time, "%Y-%m-%d %H:%M:%S")))

final_head_table <- cbind(final_head_table, exon_info_df)

# 输出总处理时间
total_end_time <- Sys.time()
total_time <- difftime(total_end_time, total_start_time, units = "mins")
cat(sprintf("\n[%s] 处理完成！总耗时: %.2f 分钟\n", 
            format(total_end_time, "%Y-%m-%d %H:%M:%S"), 
            as.numeric(total_time)))

# 输出结果统计
cat(sprintf("\n结果统计：\n"))
cat(sprintf("总行数: %d\n", nrow(final_head_table)))
cat(sprintf("微外显子数量: %d\n", sum(final_head_table$is_microexon, na.rm = TRUE)))
cat(sprintf("有效外显子数量: %d\n", sum(!is.na(final_head_table$exon_start))))

# 查看结果样本
cat("\n结果预览：\n")
print(head(final_head_table))




# 使用向量化操作处理数据
library(stringr)

# 筛选 Bombyx mori 的数据
bombyx_data <- final_head_table[final_head_table$Species == "Bombyx_mori", ]

cat(sprintf("[%s] 开始处理 Bombyx mori 的数据，共 %d 条记录...\n", 
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
            nrow(bombyx_data)))

# 向量化解析 Event_ID
event_ids <- bombyx_data$Event_ID
parts_list <- str_split(event_ids, "_")

# 提取 GeneID (第二个元素)
gene_ids <- sapply(parts_list, function(x) x[2])

# 提取染色体信息
chrs <- sapply(parts_list, function(x) {
  na_index <- which(x == "NA")
  if(length(na_index) > 0) {
    # 获取目标元素并处理第一个元素
    element1 <- substr(x[na_index + 1], 4, nchar(x[na_index + 1]))  # 去除前3个字符
    element2 <- x[na_index + 2]
    paste(element1, element2, sep = "_")  # 用下划线连接两个元素
  } else {
    NA
  }
})

# 创建bed格式数据框
bed_data <- data.frame(
  chr = chrs,
  start = bombyx_data$exon_start,
  end = bombyx_data$exon_end,
  gene_id = gene_ids,
  stringsAsFactors = FALSE
)

# 输出到BED文件
output_file <- "bombyx_mori_exons.bed"
cat(sprintf("[%s] 正在写入BED文件: %s\n", 
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
            output_file))

write.table(bed_data, 
            file = output_file, 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)

# 输出统计信息
cat(sprintf("\n处理完成！\n"))
cat(sprintf("总共处理了 %d 条 Bombyx mori 的记录\n", nrow(bed_data)))
cat(sprintf("结果已保存到: %s\n", output_file))

# 预览输出的前几行
cat("\n文件内容预览：\n")
system(sprintf("head -n 5 %s", output_file))



