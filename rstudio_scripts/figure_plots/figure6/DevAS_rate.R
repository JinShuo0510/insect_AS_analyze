library(dplyr)
library(dplyr)
library(stringr)

# 1. 从 find_devas_results_1_8_full_deleteNULL 中提取每个物种所有 DevAS 的 SE_gene_id


# 修正 devAS_df：同时抽出 SE_gene_id
devAS_df <- bind_rows(
  lapply(names(find_devas_results_1_8_full_deleteNULL), function(sp) {
    bind_rows(
      lapply(find_devas_results_1_8_full_deleteNULL[[sp]], function(tbl) {
        tbl %>%
          filter(str_starts(Event_ID, "SE_")) %>%
          # 正则提取 SE_gene_id（Event_ID 第二段）
          mutate(
            SE_gene_id = str_match(Event_ID, "^SE_([^_]+)_")[,2]
          ) %>%
          select(SE_gene_id)
      })
    ) %>%
      distinct() %>%
      mutate(species = sp)
  })
) %>%
  distinct() %>%
  mutate(devAS_flag = TRUE)

# 1. 定义需要处理的阶段
ages <- c("Cell","Egg","Nymph","Larva","Pupa","Adult")


library(dplyr)
library(data.table)

## ---------- 1. 读入 CDS 区间 ----------
# 按列名读取；BED 为 0-based，end 为 open，所以稍后要做 +1 ↔ −1 的统一
cds_dt <- fread(
  "combined_cds_output.bed",
  col.names = c("chrom","cds_start0","cds_end0","tx_id","score","strand","species")
)[ , .(
  species,
  chrom_std = sub("^chr", "", chrom),
  strand,
  cds_start = cds_start0 + 1,   # 1-based closed
  cds_end   = cds_end0
)] |> unique() |> as.data.table()
setkey(cds_dt, species, chrom_std, strand, cds_start, cds_end)



# 2. 确保 result_table_detail 中已有 SE_gene_id 列，无需额外转换

# 3. 按 species + SE_gene_id 左连接，并打标
result_table_detail <- result_table_detail %>%
  left_join(devAS_df, by = c("species", "SE_gene_id")) %>%
  mutate(isdevAS = if_else(is.na(devAS_flag), "no", "yes")) %>%
  select(-devAS_flag)


# 2. 使用 rowwise + c_across 对每个阶段判断：全局是 DevAS 且该阶段任意 PSI 非 NA
result_table_detail <- result_table_detail %>%
  rowwise() %>%
  mutate(
    Cell_is_devas  = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Cell-.*_PSI$")))),  "yes", "no"),
    Egg_is_devas   = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Egg-.*_PSI$")))),   "yes", "no"),
    Nymph_is_devas = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Nymph-.*_PSI$")))), "yes", "no"),
    Larva_is_devas = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Larva-.*_PSI$")))), "yes", "no"),
    Pupa_is_devas  = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Pupa-.*_PSI$")))),  "yes", "no"),
    Adult_is_devas = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Adult-.*_PSI$")))), "yes", "no")
  ) %>%
  ungroup()


result_table_detail <- result_table_detail %>%
  mutate(
    # 计算外显子长度：根据 match_type 选不同的 start/end
    exon_length = case_when(
      match_type == "match"     ~ as.integer(SE_exon_end)       - as.integer(SE_exon_start)       + 1,
      match_type == "upstream"  ~ as.integer(SE_upstream_end)   - as.integer(SE_upstream_start)   + 1,
      match_type == "downstream"~ as.integer(SE_downstream_end) - as.integer(SE_downstream_start) + 1,
      TRUE                       ~ NA_integer_
    ),
    # 判断是否能被3整除
    is_divide_by3 = if_else(!is.na(exon_length) & exon_length %% 3 == 0, "yes", "no")
  ) %>%
  select(-exon_length)  # 如不需要保留中间列，可删除


# ---------------- 1. 找到正确的染色体列 ----------------
result_dt <- as.data.table(result_table_detail)

# 自动识别列名（先尝试 SE_chr，再尝试 chrom）
chr_col <- intersect(c("SE_chr", "chrom"), names(result_dt))[1]
if (is.na(chr_col)) stop("未检测到染色体列 SE_chr 或 chrom")

# ---------------- 2. 统一坐标区间 ----------------
result_dt[, `:=`(
  exon_start = fifelse(match_type == "match",      as.integer(SE_exon_start),
                       fifelse(match_type == "upstream",  as.integer(SE_upstream_start),
                               as.integer(SE_downstream_start))),
  exon_end   = fifelse(match_type == "match",      as.integer(SE_exon_end),
                       fifelse(match_type == "upstream",  as.integer(SE_upstream_end),
                               as.integer(SE_downstream_end)))
)]

# 去掉 chr 前缀，生成 chrom_std
result_dt[, chrom_std := sub("^chr", "", get(chr_col))]

# strand 列同样动态处理
strand_col <- intersect(c("SE_strand","strand"), names(result_dt))[1]
setnames(result_dt, strand_col, "strand", skip_absent = TRUE)

result_dt[, row_id := .I]
setkey(result_dt, species, chrom_std, strand, exon_start, exon_end)

# ---------------- 3. 读取并标准化 CDS ----------------

# ---------------- 4. 判定重叠 ----------------
olap <- foverlaps(
  result_dt, cds_dt,
  by.x = c("species","chrom_std","strand","exon_start","exon_end"),
  by.y = c("species","chrom_std","strand","cds_start","cds_end"),
  type = "any", nomatch = 0L
)

result_dt[ , is_cds_overleap := ifelse(.I %in% unique(olap$row_id), "yes", "no")]

# ---------------- 5. 写回 ----------------
result_table_detail <- result_dt[
  , !c("row_id","chrom_std","exon_start","exon_end"), with = FALSE
] |> as_tibble()



result_dt <- as.data.table(result_table_detail)

for (age in c("Cell","Egg","Nymph","Larva","Pupa","Adult")) {
  cols <- grep(paste0("^", age, "-.*_PSI$"), names(result_dt), value = TRUE)
  result_dt[, paste0(age, "_exist") := ifelse(rowSums(!is.na(.SD)) > 0, "yes", "no"), .SDcols = cols]
}

result_table_detail <- as_tibble(result_dt)


# 删除所有以 “_PSI” 结尾的列
result_table_detail <- result_table_detail %>%
  # select(-matches("_PSI$")) %>%
  select(-matches("exon_chr$"))




result_table_alternative_detail <- result_table_alternative_detail %>%
  left_join(devAS_df, by = c("species", "SE_gene_id")) %>%
  mutate(isdevAS = if_else(is.na(devAS_flag), "no", "yes")) %>%
  select(-devAS_flag)

result_table_alternative_detail <- result_table_alternative_detail %>%
  rowwise() %>%
  mutate(
    Cell_is_devas  = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Cell-.*_PSI$")))),  "yes", "no"),
    Egg_is_devas   = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Egg-.*_PSI$")))),   "yes", "no"),
    Nymph_is_devas = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Nymph-.*_PSI$")))), "yes", "no"),
    Larva_is_devas = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Larva-.*_PSI$")))), "yes", "no"),
    Pupa_is_devas  = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Pupa-.*_PSI$")))),  "yes", "no"),
    Adult_is_devas = if_else(isdevAS == "yes" & any(!is.na(c_across(matches("^Adult-.*_PSI$")))), "yes", "no")
  ) %>%
  ungroup()

result_table_alternative_detail <- result_table_alternative_detail %>%
  mutate(
    # 计算外显子长度：根据 match_type 选不同的 start/end
    exon_length = case_when(
      match_type == "match"     ~ as.integer(SE_exon_end)       - as.integer(SE_exon_start)       + 1,
      match_type == "upstream"  ~ as.integer(SE_upstream_end)   - as.integer(SE_upstream_start)   + 1,
      match_type == "downstream"~ as.integer(SE_downstream_end) - as.integer(SE_downstream_start) + 1,
      TRUE                       ~ NA_integer_
    ),
    # 判断是否能被3整除
    is_divide_by3 = if_else(!is.na(exon_length) & exon_length %% 3 == 0, "yes", "no")
  ) %>%
  select(-exon_length)  # 如不需要保留中间列，可删除


# ---------------- 1. 找到正确的染色体列 ----------------
result_dt <- as.data.table(result_table_alternative_detail)

# 自动识别列名（先尝试 SE_chr，再尝试 chrom）
chr_col <- intersect(c("SE_chr", "chrom"), names(result_dt))[1]
if (is.na(chr_col)) stop("未检测到染色体列 SE_chr 或 chrom")

# ---------------- 2. 统一坐标区间 ----------------
result_dt[, `:=`(
  exon_start = fifelse(match_type == "match",      as.integer(SE_exon_start),
                       fifelse(match_type == "upstream",  as.integer(SE_upstream_start),
                               as.integer(SE_downstream_start))),
  exon_end   = fifelse(match_type == "match",      as.integer(SE_exon_end),
                       fifelse(match_type == "upstream",  as.integer(SE_upstream_end),
                               as.integer(SE_downstream_end)))
)]

# 去掉 chr 前缀，生成 chrom_std
result_dt[, chrom_std := sub("^chr", "", get(chr_col))]

# strand 列同样动态处理
strand_col <- intersect(c("SE_strand","strand"), names(result_dt))[1]
setnames(result_dt, strand_col, "strand", skip_absent = TRUE)

result_dt[, row_id := .I]
setkey(result_dt, species, chrom_std, strand, exon_start, exon_end)

# ---------------- 3. 读取并标准化 CDS ----------------

# ---------------- 4. 判定重叠 ----------------
olap <- foverlaps(
  result_dt, cds_dt,
  by.x = c("species","chrom_std","strand","exon_start","exon_end"),
  by.y = c("species","chrom_std","strand","cds_start","cds_end"),
  type = "any", nomatch = 0L
)

result_dt[ , is_cds_overleap := ifelse(.I %in% unique(olap$row_id), "yes", "no")]

# ---------------- 5. 写回 ----------------
result_table_alternative_detail <- result_dt[
  , !c("row_id","chrom_std","exon_start","exon_end"), with = FALSE
] |> as_tibble()


result_dt <- as.data.table(result_table_alternative_detail)

for (age in c("Cell","Egg","Nymph","Larva","Pupa","Adult")) {
  cols <- grep(paste0("^", age, "-.*_PSI$"), names(result_dt), value = TRUE)
  result_dt[, paste0(age, "_exist") := ifelse(rowSums(!is.na(.SD)) > 0, "yes", "no"), .SDcols = cols]
}

result_table_alternative_detail <- as_tibble(result_dt)


# 删除所有以 “_PSI” 结尾的列
result_table_alternative_detail <- result_table_alternative_detail %>%
  # select(-matches("_PSI$")) %>%
  select(-matches("exon_chr$"))

saveRDS(result_table_detail,"result_table_detail_result.RDS")
saveRDS(result_table_alternative_detail,"result_table_alternative_detail_result.RDS")