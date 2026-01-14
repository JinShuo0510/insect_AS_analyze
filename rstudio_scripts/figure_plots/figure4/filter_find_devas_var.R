library(dplyr)
library(stringr)

meta_cols <- c("GeneID","geneSymbol","chr","strand",
               "exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")

# 假设你的 tibble 变量名为 df
filtered_events <- find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Head"]] %>%
  filter(str_detect(Event_ID, "^SE")) %>%
  separate(
    col = Event_ID,
    into = c("Event_type", "GeneID", "geneSymbol", "chr_part1", "chr_part2", "strand",
             "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
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


# 显示提取的结果
print(filtered_events)

library(dplyr)


filtered_events_fixed <- filtered_events %>%
  mutate(
    exonStart_0base = as.numeric(exonStart_0base),
    exonEnd = as.numeric(exonEnd),
    upstreamES = as.numeric(upstreamES),
    upstreamEE = as.numeric(upstreamEE),
    downstreamES = as.numeric(downstreamES),
    downstreamEE = as.numeric(downstreamEE)
  )

matched_events <- results[["Head"]] %>%
  semi_join(filtered_events_fixed, 
            by = c("GeneID", "chr", "strand", "exonStart_0base", "exonEnd", 
                   "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")) %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE,
         downstreamES, downstreamEE, dPSI, direction, raw_p, adj_p)

print(matched_events)


merged_DevAS_table <- bind_rows(devAS_events, matched_events)

merged_DevAS_table <- merged_DevAS_table %>% distinct()


print(merged_DevAS_table)


saveRDS(merged_DevAS_table,file="merged_DevAS_table.RDS")