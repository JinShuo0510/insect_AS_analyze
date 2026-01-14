result_merge_all <- lapply(result_merge_all, function(df) {
  df %>% 
    mutate(divable = ifelse((exonEnd - exonStart_0base) %% 3 == 0, 1, 0))
})

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
      direction = .data[["direction"]],
      divisible = .data[["divable"]]
    ) %>%
    dplyr::filter(.data[["direction"]] %in% c("up", "down"))
}

# 执行提取（确保输出仅6列）
final_result <- result_merge_all %>% 
  imap_dfr(extract_events) %>% 
  dplyr::select(tissue, chr, start, end, gene, direction, divisible)

