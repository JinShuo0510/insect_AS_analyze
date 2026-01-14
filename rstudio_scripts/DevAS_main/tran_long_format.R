# 1. 定义事件类型及其对应的信息
event_types_info <- list(
  SE = list(
    df_name = "SE_events_detail",
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "exonStart_0base", "exonEnd",
                        "upstreamES", "upstreamEE",
                        "downstreamES", "downstreamEE"),
    coord_cols = c("exonStart_0base", "exonEnd")
  ),
  A3SS = list(
    df_name = "A3SS_events_detail",
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "longExonStart_0base", "longExonEnd",
                        "shortES", "shortEE",
                        "flankingES", "flankingEE"),
    coord_cols = c("longExonStart_0base", "longExonEnd")
  ),
  A5SS = list(
    df_name = "A5SS_events_detail",
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "longExonStart_0base", "longExonEnd",
                        "shortES", "shortEE",
                        "flankingES", "flankingEE"),
    coord_cols = c("longExonStart_0base", "longExonEnd")
  ),
  RI = list(
    df_name = "RI_events_detail",
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "riExonStart_0base", "riExonEnd",
                        "upstreamES", "upstreamEE",
                        "downstreamES", "downstreamEE"),
    coord_cols = c("riExonStart_0base", "riExonEnd")
  ),
  MXE = list(
    df_name = "MXE_events_detail",
    identifier_cols = c("GeneID", "geneSymbol", "chr", "strand",
                        "1stExonStart_0base", "1stExonEnd",
                        "2ndExonStart_0base", "2ndExonEnd",
                        "upstreamES", "upstreamEE",
                        "downstreamES"),
    coord_cols = c("1stExonStart_0base", "1stExonEnd",
                   "2ndExonStart_0base", "2ndExonEnd")
  )
)

# 2. 加载必要的包
library(dplyr)
library(tidyr)

# 3. 假设 all_event_wide_lists 已经按照之前的步骤生成
# 结构为 all_event_wide_lists[[event_type]][[species]]

# 4. 为每个事件类型和物种添加 eventID
# 初始化一个新列表，用于存储包含 eventID 的数据
all_event_wide_with_id <- all_event_wide_lists

# 遍历每种事件类型
for(event in names(all_event_wide_with_id)) {
  
  # 获取当前事件类型的坐标列
  coord_cols <- event_types_info[[event]]$coord_cols
  
  # 遍历每个物种
  for(species in names(all_event_wide_with_id[[event]])) {
    
    # 获取当前物种的数据框
    df <- all_event_wide_with_id[[event]][[species]]
    
    # 检查是否存在所有坐标列
    missing_coords <- setdiff(coord_cols, colnames(df))
    if(length(missing_coords) > 0){
      warning(paste("事件类型", event, "的物种", species, "缺少坐标列:", paste(missing_coords, collapse=", ")))
      next
    }
    
    # 生成 eventID
    if(event != "MXE") {
      df <- df %>%
        mutate(eventID = paste(GeneID, event, 
                               !!sym(coord_cols[1]), !!sym(coord_cols[2]),
                               sep = "_"))
    } else {
      df <- df %>%
        mutate(eventID = paste(GeneID, event, 
                               !!sym(coord_cols[1]), !!sym(coord_cols[2]),
                               !!sym(coord_cols[3]), !!sym(coord_cols[4]),
                               sep = "_"))
    }
    
    # 将 eventID 放在第一列
    df <- df %>%
      select(eventID, everything())
    
    # 更新列表
    all_event_wide_with_id[[event]][[species]] <- df
  }
}

# 5. 为每个物种和事件类型生成综合表格
# 初始化一个列表，用于存储每个物种的各事件类型的综合表格
species_event_combined_list <- list()

# 获取所有物种名称（跨事件类型）
all_species <- unique(unlist(lapply(all_event_wide_with_id, names)))

# 遍历每个物种
for(species in all_species) {
  
  # 初始化一个子列表，用于存储不同事件类型的数据框
  species_event_combined_list[[species]] <- list()
  
  # 遍历每种事件类型
  for(event in names(all_event_wide_with_id)) {
    
    # 检查当前事件类型是否包含该物种
    if(species %in% names(all_event_wide_with_id[[event]])) {
      
      # 获取当前事件类型和物种的数据框
      df <- all_event_wide_with_id[[event]][[species]]
      
      # 提取 eventID 和 PSI 列
      eventID_col <- "eventID"
      psi_cols <- setdiff(colnames(df), eventID_col)
      
      # 从 PSI 列名中提取 Tissue 信息
      tissues <- sapply(psi_cols, function(col_name) {
        # 假设 Age_Tissue 格式为 "Age_AgeX_Age_xifenX_Tissue"
        parts <- unlist(strsplit(col_name, "-"))
        if(length(parts) >= 1){
          return(tail(parts, 1))  # 取最后一个部分作为 Tissue
        } else {
          return(NA)
        }
      })
      
      # 创建 PSI 列与 Tissue 的映射
      psi_tissue_map <- data.frame(
        psi_col = psi_cols,
        Tissue = tissues,
        stringsAsFactors = FALSE
      )
      
      # 获取唯一的 Tissue 列表，保持出现顺序
      unique_tissues <- unique(psi_tissue_map$Tissue)
      
      # 按照 Tissue 排列 PSI 列
      ordered_psi_cols <- sapply(unique_tissues, function(tissue) {
        cols <- psi_tissue_map$psi_col[psi_tissue_map$Tissue == tissue]
        return(cols)
      })
      
      # 展平列表
      ordered_psi_cols <- unlist(ordered_psi_cols, use.names = FALSE)
      
      # 确保所有 PSI 列都存在于数据框中
      ordered_psi_cols <- intersect(ordered_psi_cols, colnames(df))
      
      # 按照 Tissue 排列 PSI 列，并将 eventID 放在第一列
      df <- df %>%
        select(eventID, all_of(ordered_psi_cols))
      
      # 将处理后的数据框添加到物种的子列表中
      species_event_combined_list[[species]][[event]] <- df
    }
  }
}

# 6. （可选）将生成的表格导出为 CSV 文件
# 定义输出主目录
output_dir <- "species_combined_tables"

# 创建输出主目录（如果不存在）
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 遍历每个物种
for(species in names(species_event_combined_list)) {
  
  # 定义物种的输出子目录
  species_dir <- file.path(output_dir, gsub(" ", "_", species))
  
  # 创建物种的输出子目录（如果不存在）
  if(!dir.exists(species_dir)) {
    dir.create(species_dir)
  }
  
  # 遍历每个事件类型
  for(event in names(species_event_combined_list[[species]])) {
    
    # 获取当前事件类型的数据框
    df <- species_event_combined_list[[species]][[event]]
    
    # 定义文件名，格式为 "species_event.csv"
    file_name <- paste0("combined_", event, ".csv")
    
    # 定义文件路径
    file_path <- file.path(species_dir, file_name)
    
    # 写入 CSV 文件
    write.csv(df, file = file_path, row.names = FALSE)
  }
}

cat("所有综合表格已导出至 'species_combined_tables' 目录，每个物种包含各自的事件类型子目录。\n")
