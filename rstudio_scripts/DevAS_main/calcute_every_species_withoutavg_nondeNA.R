library(dplyr)
library(readr)

file_names <- c("SE_events_details.csv", "A5SS_events_details.csv", "A3SS_events_details.csv", "MXE_events_details.csv", "RI_events_details.csv")

# 计算 NA 值比例
calculate_na_ratio <- function(psi_string) {
  if (is.na(psi_string) || psi_string == "") {
    return(1)
  }
  psi_values <- as.numeric(unlist(strsplit(psi_string, ",")))
  na_ratio <- sum(is.na(psi_values)) / length(psi_values)
  return(na_ratio)
}

process_file <- function(file_name) {
  if (!file.exists(file_name)) {
    stop(paste("文件不存在:", file_name))
  }
  
  data <- read_csv(file_name)
  
  # 计算 NA 值比例（直接对 IncLevel1 列计算）
  data$na_ratio <- sapply(data$IncLevel1, calculate_na_ratio)
  
  # 如果需要过滤数据，可以取消注释以下行
  # filtered_data <- data[data$na_ratio <= 0.5, ]
  filtered_data <- data
  
  # 移除 ID 列（如果存在）
  if ("ID" %in% colnames(filtered_data)) {
    df <- filtered_data %>% dplyr::select(-ID)
  } else {
    df <- filtered_data
  }
  
  # 保存为 RDS 文件
  event_type <- sub("_.*", "", file_name)
  saveRDS(df, file = paste0(event_type, "_events_detail_withoutavg_nondeNA.RDS"))
  
  return(df)
}

# 处理文件
for (file_name in file_names) {
  tryCatch({
    df <- process_file(file_name)
    print(df)
  }, error = function(e) {
    message("处理文件时出错: ", file_name)
    message("错误信息: ", e$message)
  })
}

