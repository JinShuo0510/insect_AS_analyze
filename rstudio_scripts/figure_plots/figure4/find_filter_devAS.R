# 创建一个函数来处理每个物种和组织的匹配
match_events_across_datasets <- function(find_devas_results, all_species_results) {
  # 创建一个空列表存储结果
  matched_results <- list()
  
  # 遍历每个物种
  for(species_name in names(find_devas_results)) {
    matched_results[[species_name]] <- list()
    
    # 遍历该物种的每个组织
    for(tissue_name in names(find_devas_results[[species_name]])) {
      # 筛选SE事件并提取必要信息
      filtered_events <- find_devas_results[[species_name]][[tissue_name]] %>%
        filter(str_detect(Event_ID, "^SE")) %>%
        separate(
          col = Event_ID,
          into = c("Event_type", "GeneID", "geneSymbol", "chr_part1", "chr_part2", "strand",
                   "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
          sep = "_",
          remove = FALSE
        ) %>%
        mutate(
          chr = if_else(is.na(chr_part2), chr_part1, paste(chr_part1, chr_part2, sep = "_")),
          .after = geneSymbol
        ) %>%
        select(-chr_part1, -chr_part2) %>%
        select(Event_ID, GeneID, geneSymbol, chr, strand, exonStart_0base, exonEnd, 
               upstreamES, upstreamEE, downstreamES, downstreamEE) %>%
        mutate(
          exonStart_0base = as.numeric(exonStart_0base),
          exonEnd = as.numeric(exonEnd),
          upstreamES = as.numeric(upstreamES),
          upstreamEE = as.numeric(upstreamEE),
          downstreamES = as.numeric(downstreamES),
          downstreamEE = as.numeric(downstreamEE)
        )
      
      # 创建事件标识符，与Fout_time_all_species_results格式匹配
      filtered_events <- filtered_events %>%
        mutate(event_id = paste(GeneID, chr, strand, exonStart_0base, exonEnd, 
                                upstreamES, upstreamEE, downstreamES, downstreamEE, sep = "__"))
      
      # 在all_species_results中查找匹配项
      species_results <- all_species_results %>%
        filter(Species == species_name, Tissue == tissue_name)
      
      # 将事件ID拆分为组成部分
      matched_events <- filtered_events %>%
        inner_join(species_results, by = "event_id") %>%
        select(GeneID, geneSymbol, chr, strand, exonStart_0base, exonEnd, 
               upstreamES, upstreamEE, downstreamES, downstreamEE, 
               direction, Larva, Pupa, Adult, Egg)
      
      matched_results[[species_name]][[tissue_name]] <- matched_events
    }
  }
  
  return(matched_results)
}

# 使用函数
matched_events_all <- match_events_across_datasets(
  find_devas_results_1_8_full_deleteNULL, 
  Fout_time_all_species_results
)

# 为匹配结果创建方向统计摘要
generate_direction_summary <- function(matched_events_all) {
  # 创建一个空的tibble用于存储所有结果
  all_summaries <- tibble(
    Species = character(),
    Tissue = character(),
    direction = character(),
    count = integer()
  )
  
  # 遍历每个物种和组织
  for (species_name in names(matched_events_all)) {
    for (tissue_name in names(matched_events_all[[species_name]])) {
      # 获取当前物种和组织的匹配事件
      current_data <- matched_events_all[[species_name]][[tissue_name]]
      
      if (nrow(current_data) > 0) {
        # 生成方向摘要
        direction_summary <- current_data %>%
          group_by(direction) %>%
          summarise(count = n(), .groups = "drop") %>%
          mutate(
            Species = species_name,
            Tissue = tissue_name
          ) %>%
          select(Species, Tissue, direction, count)
        
        # 合并到总结果中
        all_summaries <- bind_rows(all_summaries, direction_summary)
      }
    }
  }
  
  # 整理并排序结果
  all_summaries <- all_summaries %>%
    arrange(Species, Tissue, desc(count))
  
  return(all_summaries)
}

# 使用函数生成摘要
direction_summary <- generate_direction_summary(matched_events_all)

# 查看总体摘要统计
species_tissue_summary <- direction_summary %>%
  group_by(Species, Tissue) %>%
  summarise(
    total_events = sum(count),
    up_events = sum(count[direction == "up"], na.rm = TRUE),
    down_events = sum(count[direction == "down"], na.rm = TRUE),
    up_down_events = sum(count[direction == "up-down"], na.rm = TRUE),
    down_up_events = sum(count[direction == "up-down"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_events))

# 查看跨物种的方向摘要
species_direction_summary <- direction_summary %>%
  group_by(Species, direction) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(Species, desc(count))