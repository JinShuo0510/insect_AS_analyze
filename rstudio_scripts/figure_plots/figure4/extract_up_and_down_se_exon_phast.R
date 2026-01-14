# 加载必要的包
library(dplyr)
library(tidyr)
library(stringr)


df_raw <- read.table("phastcons_results.txt", header = TRUE, stringsAsFactors = FALSE)
# df_raw <- read.table("phastcons_results_phyloP.txt", header = TRUE, stringsAsFactors = FALSE)

# 将非数值字符串替换为 NA，然后强制转换为 numeric
df_raw$Average_PhastCons_Score <- gsub("^NA$|^N/A$", "", df_raw$Average_PhastCons_Score) 
df_raw$Average_PhastCons_Score <- as.numeric(df_raw$Average_PhastCons_Score)


# 假设你的数据框名为 df_raw
df_raw <- df_raw %>%
  extract(
    File,
    into = c("gene", "region", "chr", "start", "end"),
    regex = "^([^_]+)_([^_]+)_((?:[^_]+_[^_]+))_([0-9]+)[-_]([0-9]+)(?:_[0-9]+)?\\.maf$"
  )
names(df_raw)[1] <- "tissue"


# 使用字符串列名引用 + 解除分组
extract_events <- function(tissue_data, tissue_name) {
  tissue_data %>%
    dplyr::ungroup() %>%  # 关键步骤：解除分组
    dplyr::transmute(
      tissue = tissue_name,
      chr = .data[["chr"]],
      start = .data[["exonStart_0base"]],
      end = .data[["exonEnd"]],
      gene = .data[["GeneID"]],
      direction = .data[["direction"]]
    ) %>%
  dplyr::filter(.data[["direction"]] %in% c("up", "down"))
}

# 执行提取（确保输出仅6列）
final_result <- Bombyx_mori_event_fits_all_exon  %>% 
  imap_dfr(extract_events) %>% 
  dplyr::select(tissue, chr, gene, direction)

# 验证结构
str(final_result)



library(dplyr)

final_result$chr <- gsub("^chr", "", final_result$chr)

merged_df <- merge(df_raw, final_result, 
                   by = c("tissue", "chr", "gene"),
                   all.x = TRUE)  # 保留第一个表的所有行

print(merged_df)

unique(merged_df$direction)