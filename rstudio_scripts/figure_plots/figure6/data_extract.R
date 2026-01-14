# 读取数据
# 假设数据保存在名为"orthogroups.txt"的文件中
data <- read.table("one_to_one_orthologs.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 创建结果列表来存储每个ortholog group的信息
result_list <- list()

# 遍历每一行
for(i in 1:nrow(data)) {
  # 获取物种数量和ortholog group信息
  species_count <- data$species_count[i]
  orthogroup_info <- data$ortholog_group[i]
  
  # 分割ortholog group信息，获取每个物种的信息
  species_entries <- unlist(strsplit(orthogroup_info, ","))
  
  # 创建一个data frame来存储当前ortholog group的信息
  orthogroup_data <- data.frame(
    species = character(),
    gene_id = character(),
    transcript_id = character(),
    start = numeric(),
    end = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 解析每个物种条目
  for(entry in species_entries) {
    # 使用正则表达式提取物种、基因ID、转录本ID和坐标信息
    pattern <- "(.+)::(.+)::(.+)_pos(\\d+)-(\\d+)"
    matches <- regexec(pattern, entry)
    extracted <- regmatches(entry, matches)[[1]]
    
    if(length(extracted) == 6) {
      species <- extracted[2]
      gene_id <- extracted[3]
      transcript_id <- extracted[4]
      start <- as.numeric(extracted[5])
      end <- as.numeric(extracted[6])
      
      # 添加到当前ortholog group的data frame
      orthogroup_data <- rbind(orthogroup_data, data.frame(
        species = species,
        gene_id = gene_id,
        transcript_id = transcript_id,
        start = start,
        end = end,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 创建一个包含物种数量和详细信息的列表
  group_info <- list(
    species_count = species_count,
    species_details = orthogroup_data
  )
  
  # 将当前ortholog group的信息添加到结果列表
  group_id <- paste0("OG", i)
  result_list[[group_id]] <- group_info
}

# 查看结果
# 打印每个ortholog group的信息
for(og_id in names(result_list)) {
  cat("Ortholog Group:", og_id, "\n")
  cat("Species Count:", result_list[[og_id]]$species_count, "\n")
  cat("Species Details:\n")
  print(result_list[[og_id]]$species_details)
  cat("\n")
}

# 示例：访问特定ortholog group的信息
# 例如，查看第一个ortholog group的所有物种
# print(result_list[["OG1"]]$species_details)

# 保存结果为RDS文件
saveRDS(result_list, "ortholog_groups_analysis.rds")

# 可选：导出为CSV文件
# 创建一个包含所有信息的data frame
all_groups_df <- data.frame(
  ortholog_group = character(),
  species_count = numeric(),
  species = character(),
  gene_id = character(),
  transcript_id = character(),
  start = numeric(),
  end = numeric(),
  stringsAsFactors = FALSE
)

for(og_id in names(result_list)) {
  group_df <- result_list[[og_id]]$species_details
  group_df$ortholog_group <- og_id
  group_df$species_count <- result_list[[og_id]]$species_count
  
  all_groups_df <- rbind(all_groups_df, group_df[, c("ortholog_group", "species_count", 
                                                     "species", "gene_id", "transcript_id", 
                                                     "start", "end")])
}

# 保存为CSV
write.csv(all_groups_df, "ortholog_groups_details.csv", row.names = FALSE)
saveRDS(all_groups_df, "ortholog_groups_details.RDS")