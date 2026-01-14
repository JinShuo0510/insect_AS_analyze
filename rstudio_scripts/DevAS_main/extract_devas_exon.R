
# 加载所需包
library(dplyr)
library(stringr)
library(purrr)
library(readr)


bombyx_mori_devas <- find_devas_results_1_8_deleteNULL['Bombyx_mori']

# 定义最终输出的变量名称（包含 Event_type 和处理后的 chr）
identifier_cols <- c(
  "Event_type", "GeneID", "geneSymbol", "chr", "strand",
  "exonStart_0base", "exonEnd",
  "upstreamES", "upstreamEE",
  "downstreamES", "downstreamEE"
)

bombyx_mori_devas_new <- lapply(bombyx_mori_devas$Bombyx_mori, function(df) {
  df %>%
    # 筛选 Event_ID 以 "SE" 开头的行
    filter(str_detect(Event_ID, "^SE")) %>%
    # 按 "_" 拆分 Event_ID，生成列表列 split_ID
    mutate(split_ID = str_split(Event_ID, "_")) %>%
    # 提取各字段，同时对染色体字段进行处理：
    # 将第4和第5项合并，并去除合并后字符串开头的 "chr"
    mutate(
      Event_type      = map_chr(split_ID, ~ .x[1]),
      GeneID          = map_chr(split_ID, ~ .x[2]),
      geneSymbol      = map_chr(split_ID, ~ .x[3]),
      chr             = map_chr(split_ID, ~ str_remove(paste0(.x[4], "_", .x[5]), "^chr")),
      strand          = map_chr(split_ID, ~ .x[6]),
      exonStart_0base = map_chr(split_ID, ~ .x[7]),
      exonEnd         = map_chr(split_ID, ~ .x[8]),
      upstreamES      = map_chr(split_ID, ~ .x[9]),
      upstreamEE      = map_chr(split_ID, ~ .x[10]),
      downstreamES    = map_chr(split_ID, ~ .x[11]),
      downstreamEE    = map_chr(split_ID, ~ .x[12])
    ) %>%
    # 调整输出列顺序：原始 Event_ID 与所有新提取的变量
    select(Event_ID, all_of(identifier_cols))
})

bombyx_mori_devas_new$Developmental_tissue <- NULL


# 使用新列表的名称作为组织标识
tissue_names <- names(bombyx_mori_devas_new)

# 合并各组织数据，并增加新列 tissue 标识组织名称
combined_df <- map2_dfr(bombyx_mori_devas_new, tissue_names, ~ mutate(.x, tissue = .y))

# 保存为 CSV 文件
write_csv(combined_df, "bombyx_mori_devas_new.csv")
