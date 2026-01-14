library(stringr)

devAS_events_every <- Bombyx_mori_event_fits_all[["Testis"]] %>%
  filter(dPSI > 0.2, adj_p < 0.05) %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE,
         downstreamES, downstreamEE, dPSI, direction, raw_p, adj_p)


#find devAS
meta_cols <- c("GeneID","geneSymbol","chr","strand",
               "exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")

# 假设你的 tibble 变量名为 df
filtered_events <- find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Testis"]] %>%
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

matched_events <- Bombyx_mori_event_fits_all[["Testis"]] %>%
  semi_join(filtered_events_fixed, 
            by = c("GeneID", "chr", "strand", "exonStart_0base", "exonEnd", 
                   "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")) %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE,
         downstreamES, downstreamEE, dPSI, direction, raw_p, adj_p)

# print(matched_events)


merged_DevAS_all <- bind_rows(devAS_events_every, matched_events)

merged_DevAS_all <- merged_DevAS_all %>% distinct()




devAS_events_find <- matched_events %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE,
         downstreamES, downstreamEE, dPSI, direction, raw_p, adj_p)


result <- devAS_events_every %>%
       filter(direction %in% c("up", "down", "down-up", "up-down")) %>%
       group_by(direction) %>% 
       summarise(count = n())

result_find <- devAS_events_find %>%
  filter(direction %in% c("up", "down", "down-up", "up-down")) %>%
  group_by(direction) %>% 
  summarise(count = n())

result_merge <- merged_DevAS_all %>%
  filter(direction %in% c("up", "down", "down-up", "up-down")) %>%
  group_by(direction) %>% 
  summarise(count = n())


print(result)
print(result_find)
print(result_merge)