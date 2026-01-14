library(dplyr)
library(stringr)
library(tidyr)

# 获取所有组织名称
tissues <- names(Bombyx_mori_event_fits_all)

# 对每个组织进行遍历，并返回对应的 merged_DevAS_all，保持列表结构
result_merge_all <- lapply(tissues, function(tissue) {
  
  # 检查当前组织在 Bombyx_mori_event_fits_all 中是否存在有效数据
  current_fits <- Bombyx_mori_event_fits_all[[tissue]]
  if (is.null(current_fits)) {
    message(paste("组织", tissue, "在 Bombyx_mori_event_fits_all 中不存在，跳过该组织。"))
    return(NULL)
  }
  
  # 检查当前组织在 find_devas_results_1_8_full_deleteNULL 中是否存在有效数据
  current_devas <- find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][[tissue]]
  if (is.null(current_devas)) {
    message(paste("组织", tissue, "在 find_devas_results_1_8_full_deleteNULL 中不存在，跳过该组织。"))
    return(NULL)
  }

  # 1. 提取当前组织的 devAS_events_every
  devAS_events_every <- current_fits %>%
    filter(dPSI > 0.2, adj_p < 0.05) %>%
    select(GeneID, chr, strand, exonStart_0base, exonEnd,
           upstreamES, upstreamEE, downstreamES, downstreamEE,
           dPSI, direction, raw_p, adj_p)

  # # 2. 从 find_devas 数据中提取 SE 类型事件，并整理变量
  meta_cols <- c("GeneID", "geneSymbol", "chr", "strand",
                 "exonStart_0base", "exonEnd",
                 "upstreamES", "upstreamEE",
                 "downstreamES", "downstreamEE")

  filtered_events <- current_devas %>%
    filter(str_detect(Event_ID, "^SE")) %>%
    separate(
      col = Event_ID,
      into = c("Event_type", "GeneID", "geneSymbol", "chr_part1", "chr_part2", "strand",
               "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
               "downstreamES", "downstreamEE"),
      sep = "_",
      remove = FALSE
    ) %>%
    mutate(
      chr = paste(chr_part1, chr_part2, sep = "_"),
      .after = geneSymbol
    ) %>%
    select(-chr_part1, -chr_part2) %>%
    select(Event_type, all_of(meta_cols)) %>%
    select(-1)

  # # 3. 将部分变量转换为数值型，确保后续匹配时类型一致
  filtered_events_fixed <- filtered_events %>%
    mutate(
      exonStart_0base = as.numeric(exonStart_0base),
      exonEnd = as.numeric(exonEnd),
      upstreamES = as.numeric(upstreamES),
      upstreamEE = as.numeric(upstreamEE),
      downstreamES = as.numeric(downstreamES),
      downstreamEE = as.numeric(downstreamEE)
    )

  # # 4. 利用 semi_join 进行事件匹配
  matched_events <- current_fits %>%
    semi_join(filtered_events_fixed,
              by = c("GeneID", "chr", "strand",
                     "exonStart_0base", "exonEnd",
                     "upstreamES", "upstreamEE",
                     "downstreamES", "downstreamEE")) %>%
    select(GeneID, chr, strand, exonStart_0base, exonEnd,
           upstreamES, upstreamEE, downstreamES, downstreamEE,
           dPSI, direction, raw_p, adj_p)

  # 5. 合并 devAS_events_every 与 matched_events，并去除重复记录
  merged_DevAS_all <- bind_rows(devAS_events_every, matched_events) %>% distinct()
  # merged_DevAS_all <- devAS_events_every %>% distinct()
  # 返回结果，保持原始变量结构（单个数据框）
  return(merged_DevAS_all)
})

# 为列表元素命名，保持与组织名称一致
names(result_merge_all) <- tissues
result_merge_all <- discard(result_merge_all, is.null)

# 保存结果（列表结构）
saveRDS(result_merge_all, "Bombyx_mori_event_fits_all_merge.RDS")

# 打印每个组织的结果结构
lapply(names(result_merge_all), function(tissue) {
  cat("组织：", tissue, "\n")
  print(result_merge_all[[tissue]])
  cat("\n")
})
