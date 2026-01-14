# 加载必要的库
library(dplyr)

# 定义一个函数用于提取 SE 事件的基因 ID
extract_se_geneIDs <- function(species_data) {
  se_geneIDs <- c()
  
  for (tissue in names(species_data)) {
    tissue_data <- species_data[[tissue]]
    
    if (!is.null(tissue_data)) {
      event_ids <- tissue_data$Event_ID
      se_event_ids <- event_ids[grep("^SE_", event_ids)]
      geneIDs <- gsub("^SE_(.+?)_.+", "\\1", se_event_ids)
      se_geneIDs <- c(se_geneIDs, geneIDs)
    }
  }
  
  unique(se_geneIDs)
}

# 提取 Apis_mellifera 的 SE 事件基因 ID
apis_geneIDs <- extract_se_geneIDs(find_devas_results_1_8_full[["Apis_mellifera"]])

# 提取 Bombyx_mori 的 SE 事件基因 ID
bombyx_geneIDs <- extract_se_geneIDs(find_devas_results_1_8_full[["Bombyx_mori"]])

# 打印结果
cat("Apis_mellifera SE events - Unique GeneIDs:\n")
print(apis_geneIDs)
cat("Total unique GeneIDs for Apis_mellifera:", length(apis_geneIDs), "\n\n")

cat("Bombyx_mori SE events - Unique GeneIDs:\n")
print(bombyx_geneIDs)
cat("Total unique GeneIDs for Bombyx_mori:", length(bombyx_geneIDs), "\n\n")

# 合并两个物种的基因 ID
combined_geneIDs <- unique(c(apis_geneIDs, bombyx_geneIDs))

cat("Combined unique GeneIDs for both species:\n")
print(combined_geneIDs)
cat("Total combined unique GeneIDs:", length(combined_geneIDs), "\n")


# 加载必要的库
library(dplyr)

# 1. 读取同源基因表
orthogroups <- read.table("converted_orthogroups.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# 2. 筛选基因 ID 并保留三列信息
# 筛选 Apis_mellifera 的基因 ID
apis_in_table <- orthogroups %>% filter(Apis_mellifera %in% apis_geneIDs)

# 筛选 Bombyx_mori 的基因 ID
bombyx_in_table <- orthogroups %>% filter(Bombyx_mori %in% bombyx_geneIDs)

# 3. 取两个物种的并集
combined_table <- bind_rows(apis_in_table, bombyx_in_table) %>%
  distinct(Orthogroup, .keep_all = TRUE)

# 4. 展示筛选后的完整结果
cat("Filtered Orthogroups with corresponding gene IDs:\n")
print(combined_table)

# 将结果写入文件（可选）
write.table(combined_table, file = "filtered_orthogroups_full.txt", sep = "\t", row.names = FALSE, quote = FALSE)

