# 定义事件类型信息
event_types_info <- list(
  SE = list(
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "exonStart_0base", "exonEnd",
                        "upstreamES", "upstreamEE",
                        "downstreamES", "downstreamEE"),
    coord_cols = c("exonStart_0base", "exonEnd")
  ),
  A3SS = list(
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "longExonStart_0base", "longExonEnd",
                        "shortES", "shortEE",
                        "flankingES", "flankingEE"),
    coord_cols = c("longExonStart_0base", "longExonEnd")
  ),
  A5SS = list(
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "longExonStart_0base", "longExonEnd",
                        "shortES", "shortEE",
                        "flankingES", "flankingEE"),
    coord_cols = c("longExonStart_0base", "longExonEnd")
  ),
  RI = list(
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "riExonStart_0base", "riExonEnd",
                        "upstreamES", "upstreamEE",
                        "downstreamES", "downstreamEE"),
    coord_cols = c("riExonStart_0base", "riExonEnd")
  ),
  MXE = list(
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "1stExonStart_0base", "1stExonEnd",
                        "2ndExonStart_0base", "2ndExonEnd",
                        "upstreamES", "upstreamEE",
                        "downstreamES"),
    coord_cols = c("1stExonStart_0base", "1stExonEnd",
                   "2ndExonStart_0base", "2ndExonEnd")
  )
)

# 假设 four_four_all_event_wide_lists 包含所有物种的数据
# 获取所有物种的名称
all_species <- unique(unlist(lapply(four_all_event_wide_lists, names)))

# 初始化一个列表，用于存储所有物种的数据
all_species_data <- list()

# 遍历所有物种
for (species_name in all_species) {
  # 初始化一个列表，用于存储当前物种的所有事件数据
  species_data <- list()
  
  # 遍历所有事件类型
  for (event_type in names(four_all_event_wide_lists)) {
    # 检查当前事件类型中是否存在当前物种的数据
    if (species_name %in% names(four_all_event_wide_lists[[event_type]])) {
      # 提取当前物种的数据框
      species_data[[event_type]] <- four_all_event_wide_lists[[event_type]][[species_name]]
    }
  }
  
  # 将当前物种的数据存储到 all_species_data 中
  all_species_data[[species_name]] <- species_data
  
  # 保存当前物种的数据到文件
  # saveRDS(species_data, file = paste0(species_name, "_data.rds"))
}

# 定义生成 Event_ID 的函数
generate_event_id <- function(event_type, row, event_info) {
  gene_id <- row[["GeneID"]]
  coord_values <- sapply(event_info$coord_cols, function(col) row[[col]])
  coord_str <- paste(coord_values, collapse = "_")
  event_id <- paste(event_type, gene_id, coord_str, sep = "_")
  return(event_id)
}

# 初始化一个列表，用于存储所有物种的 Event_ID 数据
all_species_event_data_with_id <- list()

# 遍历每个物种
for (species_name in all_species) {
  # 初始化一个列表，用于存储当前物种的 Event_ID 数据
  species_event_data_with_id <- list()
  
  # 遍历每个事件类型
  for (event_type in names(event_types_info)) {
    # 获取当前事件类型的数据框
    event_data <- all_species_data[[species_name]][[event_type]]
    
    # 检查 event_data 是否为空
    if (is.null(event_data) || nrow(event_data) == 0) {
      message(paste("Skipping", event_type, "for species", species_name, "due to empty data."))
      next  # 如果为空，跳过当前事件类型
    }
    
    # 获取当前事件类型的规则
    event_info <- event_types_info[[event_type]]
    
    # 为每一行生成 Event_ID
    event_data$Event_ID <- apply(event_data, 1, function(row) {
      generate_event_id(event_type, row, event_info)
    })
    
    # 提取 PSI 值列（假设 PSI 值列的名称以 "PSI_" 开头）
    psi_cols <- grep("-", names(event_data), value = TRUE)
    
    # 保留 Event_ID 和 PSI 值列
    species_event_data_with_id[[event_type]] <- event_data[, c("Event_ID", psi_cols)]
  }
  
  # 将当前物种的 Event_ID 数据存储到 all_species_event_data_with_id 中
  all_species_event_data_with_id[[species_name]] <- species_event_data_with_id
}

# 合并所有物种的 Event_ID 数据
all_species_combined_data <- lapply(all_species, function(species_name) {
  species_event_data_with_id <- all_species_event_data_with_id[[species_name]]
  all_cols <- unique(unlist(lapply(species_event_data_with_id, colnames)))
  
  for (event_type in names(species_event_data_with_id)) {
    missing_cols <- setdiff(all_cols, colnames(species_event_data_with_id[[event_type]]))
    species_event_data_with_id[[event_type]][missing_cols] <- NA
  }
  
  combined_data <- do.call(rbind, species_event_data_with_id)
  return(combined_data)
})

# 命名列表
names(all_species_combined_data) <- all_species

# 查看合并后的表格
print(all_species_combined_data)

# 解析列名并提取 Tissue
extract_tissue <- function(colname) {
  parts <- unlist(strsplit(colname, "-"))
  tissue <- parts[length(parts)]  # 最后一个部分是 Tissue
  
  # 将 Egg 和 Pupa 归类到 Whole_body
  if (tissue %in% c("Egg", "Pupa")) {
    tissue <- "Whole_body"
  }
  
  return(tissue)
}

# 初始化一个列表，用于存储所有物种的 Tissue 数据
all_species_tissue_data <- list()

# 遍历每个物种
for (species_name in all_species) {
  combined_data <- all_species_combined_data[[species_name]]
  
  # 获取所有列名（不包括 Event_ID）
  all_cols <- colnames(combined_data)[-1]  # 排除 Event_ID 列
  
  # 按 Tissue 分组
  tissue_groups <- split(all_cols, sapply(all_cols, extract_tissue))
  
  # 创建 [Tissue][dataframe] 结构
  tissue_data <- lapply(names(tissue_groups), function(tissue) {
    cols_to_keep <- c("Event_ID", tissue_groups[[tissue]])
    combined_data[, cols_to_keep, drop = FALSE]
  })
  
  # 命名列表
  names(tissue_data) <- names(tissue_groups)
  
  # 将当前物种的 Tissue 数据存储到 all_species_tissue_data 中
  all_species_tissue_data[[species_name]] <- tissue_data
  
  # 保存当前物种的 Tissue 数据到文件
  # saveRDS(tissue_data, file = paste0(species_name, "_tissue_data.rds"))
}

# 保存所有物种的 Tissue 数据到文件
saveRDS(all_species_tissue_data, file = "all_species_tissue_data.rds")

# 查看结果
print(all_species_tissue_data)
