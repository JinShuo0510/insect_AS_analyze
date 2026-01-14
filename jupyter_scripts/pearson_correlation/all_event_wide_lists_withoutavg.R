SE_events_detail <- readRDS("SE_events_detail_withoutavg_1_8.RDS")
A3SS_events_detail <- readRDS("A3SS_events_detail_withoutavg_1_8.RDS")
A5SS_events_detail <- readRDS("A5SS_events_detail_withoutavg_1_8.RDS")
RI_events_detail <- readRDS("RI_events_detail_withoutavg_1_8.RDS")
MXE_events_detail <- readRDS("MXE_events_detail_withoutavg_1_8.RDS")

# 1. 定义事件类型及其对应的数据框和唯一标识列
event_types_info <- list(
  SE = list(
    df_name = "SE_events_detail",
    identifier_cols = c(
      "GeneID", "geneSymbol", "chr", "strand",
      "exonStart_0base", "exonEnd",
      "upstreamES", "upstreamEE",
      "downstreamES", "downstreamEE"
    )
  ),
  A3SS = list(
    df_name = "A3SS_events_detail",
    identifier_cols = c(
      "GeneID", "geneSymbol", "chr", "strand",
      "longExonStart_0base", "longExonEnd",
      "shortES", "shortEE",
      "flankingES", "flankingEE"
    )
  ),
  A5SS = list(
    df_name = "A5SS_events_detail",
    identifier_cols = c(
      "GeneID", "geneSymbol", "chr", "strand",
      "longExonStart_0base", "longExonEnd",
      "shortES", "shortEE",
      "flankingES", "flankingEE"
    )
  ),
  RI = list(
    df_name = "RI_events_detail",
    identifier_cols = c(
      "GeneID", "geneSymbol", "chr", "strand",
      "riExonStart_0base", "riExonEnd",
      "upstreamES", "upstreamEE",
      "downstreamES", "downstreamEE"
    )
  ),
  MXE = list(
    df_name = "MXE_events_detail",
    identifier_cols = c(
      "GeneID", "geneSymbol", "chr", "strand",
      "1stExonStart_0base", "1stExonEnd",
      "2ndExonStart_0base", "2ndExonEnd",
      "upstreamES", "upstreamEE",
      "downstreamES"
    )
  )
)

# 2. 加载必要的包
library(dplyr)
library(tidyr)

# 3. 初始化一个空列表，用于存储所有事件类型的结果
all_event_wide_lists <- list()

# 4. 遍历每种事件类型，处理数据
for (event in names(event_types_info)) {
  
  # 获取当前事件类型的信息
  event_info <- event_types_info[[event]]
  df_name <- event_info$df_name
  identifier_cols <- event_info$identifier_cols
  
  # 检查数据框是否存在
  if(!exists(df_name)){
    warning(paste("数据框", df_name, "不存在。跳过事件类型", event, "。"))
    next
  }
  
  # 获取数据框
  event_df <- get(df_name)
  
  # 创建 Age_Tissue 列
  event_df <- event_df %>%
	  mutate(Age_Tissue = paste(Age, Age_Xifen, Tissue, sep = "__"))

  print(event_df)
  # 定义需要选择的列，包括唯一标识列、Age_Tissue、avg_psi 和 Species
#  required_cols <- c(identifier_cols, "Age_Tissue", "avg_psi", "Species")
  required_cols <- c(identifier_cols, "Age_Tissue", "IncLevel1", "Species")
  
  # 检查所有必要的列是否存在
  missing_cols <- setdiff(required_cols, colnames(event_df))
  
  if(length(missing_cols) > 0){
    warning(paste("在事件类型", event, "中缺少以下列:", paste(missing_cols, collapse = ", "), "。跳过该事件类型。"))
    next
  }
  
  # 筛选必要的列并去除重复行
  print(sum(duplicated(event_df)))

  event_data <- event_df %>%
    select(all_of(required_cols)) %>%
    distinct()
   
  # 按 Species 分割数据，使用 split() 函数
  species_split <- split(event_data, event_data$Species)
  
  # 对每个物种的数据进行宽格式转换
  species_wide_list <- lapply(names(species_split), function(species) {
    df <- species_split[[species]]
    
    # 转换为宽格式
    df_wide <- df %>%
      select(-Species) %>%
      pivot_wider(
        names_from = Age_Tissue,
#        values_from = avg_psi
	values_from = IncLevel1
      )
    
    return(df_wide)
  })
  
  # 为列表元素命名为对应的物种名称
  names(species_wide_list) <- names(species_split)
  
  # 将当前事件类型的结果添加到主列表中
  all_event_wide_lists[[event]] <- species_wide_list
}


saveRDS(all_event_wide_lists, file = "all_event_wide_lists_withoutavg_1_8.RDS")
