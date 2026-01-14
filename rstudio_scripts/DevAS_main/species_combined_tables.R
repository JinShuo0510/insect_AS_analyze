# # 假设家蚕的物种名称为 "Bombyx_mori"
# species_name <- "Bombyx_mori"
# 
# # 初始化一个列表，用于存储家蚕的所有事件数据
# bombyx_mori_data <- list()
# 
# # 遍历所有事件类型
# for (event_type in names(all_event_wide_lists)) {
#   # 检查当前事件类型中是否存在家蚕的数据
#   if (species_name %in% names(all_event_wide_lists[[event_type]])) {
#     # 提取家蚕的数据框
#     bombyx_mori_data[[event_type]] <- all_event_wide_lists[[event_type]][[species_name]]
#   }
# }
# 
# # 查看提取的家蚕数据
# print(bombyx_mori_data)
# 
# saveRDS(bombyx_mori_data, file="bombyx_mori_data.rds")
# 
# 
# event_types_info <- list(
#   SE = list(
#     identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
#                         "exonStart_0base", "exonEnd",
#                         "upstreamES", "upstreamEE",
#                         "downstreamES", "downstreamEE"),
#     coord_cols = c("exonStart_0base", "exonEnd")
#   ),
#   A3SS = list(
#     identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
#                         "longExonStart_0base", "longExonEnd",
#                         "shortES", "shortEE",
#                         "flankingES", "flankingEE"),
#     coord_cols = c("longExonStart_0base", "longExonEnd")
#   ),
#   A5SS = list(
#     identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
#                         "longExonStart_0base", "longExonEnd",
#                         "shortES", "shortEE",
#                         "flankingES", "flankingEE"),
#     coord_cols = c("longExonStart_0base", "longExonEnd")
#   ),
#   RI = list(
#     identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
#                         "riExonStart_0base", "riExonEnd",
#                         "upstreamES", "upstreamEE",
#                         "downstreamES", "downstreamEE"),
#     coord_cols = c("riExonStart_0base", "riExonEnd")
#   ),
#   MXE = list(
#     identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
#                         "1stExonStart_0base", "1stExonEnd",
#                         "2ndExonStart_0base", "2ndExonEnd",
#                         "upstreamES", "upstreamEE",
#                         "downstreamES"),
#     coord_cols = c("1stExonStart_0base", "1stExonEnd",
#                    "2ndExonStart_0base", "2ndExonEnd")
#   )
# )
# 
# 
# # 定义生成 Event_ID 的函数
# generate_event_id <- function(event_type, row, event_info) {
#   gene_id <- row[["GeneID"]]
#   coord_values <- sapply(event_info$coord_cols, function(col) row[[col]])
#   coord_str <- paste(coord_values, collapse = "_")
#   event_id <- paste(event_type, gene_id, coord_str, sep = "_")
#   return(event_id)
# }
# 
# # 初始化一个列表，用于存储添加了 Event_ID 的数据框
# event_data_with_id <- list()
# 
# # 遍历每个事件类型
# for (event_type in names(event_types_info)) {
#   # 获取当前事件类型的数据框
#   event_data <- all_event_wide_lists[[event_type]][["Bombyx_mori"]]
# 
#   # 获取当前事件类型的规则
#   event_info <- event_types_info[[event_type]]
# 
#   # 为每一行生成 Event_ID
#   event_data$Event_ID <- apply(event_data, 1, function(row) {
#     generate_event_id(event_type, row, event_info)
#   })
# 
#   # 提取 PSI 值列（假设 PSI 值列的名称以 "PSI_" 开头）
#   psi_cols <- grep("-", names(event_data), value = TRUE)
# 
#   # 保留 Event_ID 和 PSI 值列
#   event_data_with_id[[event_type]] <- event_data[, c("Event_ID", psi_cols)]
# }
# 
# all_cols <- unique(unlist(lapply(event_data_with_id, colnames)))
# 
# for (event_type in names(event_data_with_id)) {
#   missing_cols <- setdiff(all_cols, colnames(event_data_with_id[[event_type]]))
#   event_data_with_id[[event_type]][missing_cols] <- NA
# }
# # 合并所有事件类型的数据框
# combined_data <- do.call(rbind, event_data_with_id)
# 
# # 查看合并后的表格
# print(combined_data)


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

# 查看结果
print(tissue_data)


