library(dplyr)
library(purrr)

df_raw <- read.table("merged_3N_output.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(df_raw) <- c("tissue","chr","start","end","gene","divisible")


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
final_result <- Bombyx_mori_event_fits_all  %>% 
  imap_dfr(extract_events) %>% 
  dplyr::select(tissue, chr, start, end, gene, direction)

# 验证结构
str(final_result)



library(dplyr)

final_result$chr <- gsub("^chr", "", final_result$chr)

merged_df <- merge(df_raw, final_result, 
                   by = c("tissue", "chr", "start", "end", "gene"),
                   all.x = TRUE)  # 保留第一个表的所有行

