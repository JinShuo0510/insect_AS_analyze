# 加载必要的 R 包
library(dplyr)
library(readr)
library(tidyr)
library(purrr)

# —— 1. 解析 DevAS 的 Event_ID 并提取 SE 事件 —— 
devAS_parsed <- find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Head"]] %>%
  as_tibble() %>%
  mutate(Event_ID = as.character(Event_ID)) %>%
  extract(
    Event_ID,
    into = c("event_type","gene_id","chr","strand",
             "exon_start","exon_end",
             "upstream_start","upstream_end",
             "downstream_start","downstream_end"),
    regex = "^([^_]+)_([^_]+)_NA_chr([A-Za-z0-9_.]+)_([+-])_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)$",
    remove = FALSE
  ) %>%
  mutate_at(vars(exon_start:downstream_end), as.integer)

SE_events <- devAS_parsed %>%
  filter(event_type == "SE") %>%
  select(chr, exon_start, exon_end,
         upstream_start, upstream_end,
         downstream_start, downstream_end)

# —— 2. 读取两份 DTx 文件并合并 —— 
devDTx <- read_tsv(
  "Head/devDTx_nonDEG_with_exons.tsv",
  col_types = cols(
    id         = col_character(),
    logFC      = col_double(),
    FDR        = col_double(),
    gene_id    = col_character(),
    exon_chr   = col_character(),
    exon_start = col_integer(),
    exon_end   = col_integer()
  )
)

nonDTx <- read_tsv(
  "Head/DTx_nonDEG_with_exons.tsv",
  col_types = cols(
    id         = col_character(),
    logFC      = col_double(),
    FDR        = col_double(),
    gene_id    = col_character(),
    exon_chr   = col_character(),
    exon_start = col_integer(),
    exon_end   = col_integer()
  )
)

mergedDTx <- bind_rows(devDTx, nonDTx)

# —— 3. 匹配方式1：exon + upstream + downstream 全部匹配 ±10 bp —— 
mode1 <- SE_events %>%
  pmap_dfr(function(chr,
                    exon_start, exon_end,
                    upstream_start, upstream_end,
                    downstream_start, downstream_end) {
    df_chr <- mergedDTx %>% filter(exon_chr == chr)
    ids_exon <- df_chr %>%
      filter(abs(exon_start - !!exon_start) <= 10,
             abs(exon_end   - !!exon_end)   <= 10) %>%
      pull(id)
    ids_up   <- df_chr %>%
      filter(abs(exon_start - !!upstream_start) <= 10,
             abs(exon_end   - !!upstream_end)   <= 10) %>%
      pull(id)
    ids_down <- df_chr %>%
      filter(abs(exon_start - !!downstream_start) <= 10,
             abs(exon_end   - !!downstream_end)   <= 10) %>%
      pull(id)
    common   <- Reduce(intersect, list(ids_exon, ids_up, ids_down))
    if (length(common) > 0) tibble(id = common, is_mode1 = TRUE)
    else tibble(id = character(0), is_mode1 = logical(0))
  }) %>%
  distinct()

# —— 4. 匹配方式2：仅 upstream + downstream 匹配 ±10 bp —— 
mode2 <- SE_events %>%
  pmap_dfr(function(chr,
                    exon_start, exon_end,
                    upstream_start, upstream_end,
                    downstream_start, downstream_end) {
    df_chr <- mergedDTx %>% filter(exon_chr == chr)
    ids_up   <- df_chr %>%
      filter(abs(exon_start - !!upstream_start) <= 10,
             abs(exon_end   - !!upstream_end)   <= 10) %>%
      pull(id)
    ids_down <- df_chr %>%
      filter(abs(exon_start - !!downstream_start) <= 10,
             abs(exon_end   - !!downstream_end)   <= 10) %>%
      pull(id)
    common   <- intersect(ids_up, ids_down)
    if (length(common) > 0) tibble(id = common, is_mode2 = TRUE)
    else tibble(id = character(0), is_mode2 = logical(0))
  }) %>%
  distinct()

# —— 5. 合并匹配结果并添加 match_type —— 
merged_matched <- mergedDTx %>%
  left_join(mode1, by = "id") %>%
  left_join(mode2, by = "id") %>%
  mutate(
    is_mode1   = replace_na(is_mode1, FALSE),
    is_mode2   = replace_na(is_mode2, FALSE),
    match_type = case_when(
      is_mode1 ~ "mode1",
      is_mode2 ~ "mode2",
      TRUE     ~ "none"
    )
  )

# —— 6. 输出最终表格 —— 
write_tsv(merged_matched, "merged_DTx_all_with_match_head.tsv")
