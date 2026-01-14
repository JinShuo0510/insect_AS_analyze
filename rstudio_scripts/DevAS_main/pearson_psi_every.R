

data_list <- all_event_wide_lists_withoutavg_1_8[["SE"]][["Bombyx_mori"]][1:1000, ]



########################################
## 1.  先挑出元信息列 + PSI列
########################################
# 你的数据前面有 10 列左右是基因位点相关信息(见 str 输出)，
# 然后每个 {时期}-{组织} 为一列。这里仅做示例，假设前 10 列为元信息:
meta_cols <- c("GeneID","geneSymbol","chr","strand",
               "exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")

# 如果实际的元信息列数量和名称有差异，请自行修正
psi_cols <- setdiff(colnames(data_list), meta_cols)

########################################
## 2.  将宽表转换为长表，并解析 PSI 值
########################################
# 思路：pivot_longer 将多个 "{时期}-{组织}" 列转为：行里存 {时期}-{组织} 以及 原始字符串PSI
data_long <- data_list %>%
  select(all_of(meta_cols), all_of(psi_cols)) %>%
  pivot_longer(
    cols      = all_of(psi_cols),         # 要从宽转长的列
    names_to  = "Stage_Tissue",           # 原本列名会存到这
    values_to = "PSI_string"              # 原本的PSI字符串
  )

# 现在每一行只有一个 PSI_string，但它可能有多个重复值，例如"0.962,0.968,1.0"
# 可以把这些字符串分割成数值，得到多行或取平均
# 下面示例是：先拆成多行，再对同一GeneID的重复做平均（如果你想保留每个重复，可不平均）
data_long_parsed <- data_long %>%
  # 把 PSI_string 用逗号分隔拆开到多行
  separate_rows(PSI_string, sep = ",") %>%
  # 把字符转换为数值(过滤真正的NA或空字符)
  mutate(PSI_value = as.numeric(PSI_string)) %>%
  # 可能会有一些NA，去掉
  filter(!is.na(PSI_value))



########################################
## 0. 加载依赖包
########################################
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)

########################################
## 1. 假设你已有 data_long_parsed:
##    它至少包含以下列：
##    - GeneID
##    - Stage_Tissue (如 "Larva-Head","Adult-Leg" 等)
##    - PSI_value (数值)
##
##    先从 Stage_Tissue 分割出 Stage / Tissue
########################################

data_long_parsed2 <- data_long_parsed %>%
  separate(
    col    = Stage_Tissue,
    into   = c("Stage","Stage-xifen", "Tissue"),
    sep    = "__",
    remove = FALSE  # 保留原列
  )


data_long_parsed2 <- data_long_parsed2 %>%
  mutate(
    Stage = case_when(
      Stage == "Cell" ~ "Egg",
      TRUE                        ~ Stage
    )
  ) %>%
  mutate(
    Tissue = case_when(
      Tissue %in% c("Egg", "Pupa") ~ "Whole_body",  # Egg/Pupa -> Whole_body
      TRUE                        ~ Tissue          # 其他情况保持不变
    )
  )

# 如果你希望给 Stage 设置固定的因子顺序(使画图时 x 轴合乎发育先后),
# 例如: Egg -> Larva -> Pupa -> Adult
stage_levels <- c( "Egg", "Larva", "Pupa", "Adult") 


data_long_parsed2 <- data_long_parsed2 %>%
  mutate(
    Stage = factor(Stage, levels = stage_levels),  # 大时期排序
    Combined_Stage = paste(Stage, `Stage-xifen`, sep = "-")  # 联合字段
  ) %>%
  arrange(Stage, `Stage-xifen`)  # 同一大时期内按细分时期排序

########################################
## 2. 定义函数：对某个 Tissue (组内数据) 找其“最早Stage”，
##    并对其他 Stage 做 bootstrap 相关
########################################
get_earliest_stage_cor_bootstrap <- function(df_grouped, n_boot = 1000, seed = 123) {
  set.seed(seed)
  
  #（a）检查数据是否为空
  df_valid <- df_grouped %>%
    filter(!is.na(Stage), !is.na(`Stage-xifen`))
  if (nrow(df_valid) == 0) {
    return(
      tibble(
        Stage      = NA,
        `Stage-xifen` = NA,
        Tissue     = unique(df_grouped$Tissue),
        cor_median = NA_real_,
        cor_lower  = NA_real_,
        cor_upper  = NA_real_,
        n_genes    = NA_integer_
      )
    )
  }
  
  #（b）按 {Stage, Stage-xifen} 找到最早的时期
  earliest_stage <- df_valid %>%
    arrange(Stage, `Stage-xifen`) %>%  # 先按大时期，再按细分时期排序
    slice(1) %>%
    select(Stage, `Stage-xifen`) %>%
    mutate(Combined_Stage = paste(Stage, `Stage-xifen`, sep = "-"))  # 联合时期
  
  print(earliest_stage)
  
  #（c）提取最早时期的 PSI 数据
  df_earliest <- df_valid %>%
    filter(Stage == earliest_stage$Stage[1] & `Stage-xifen` == earliest_stage$`Stage-xifen`[1]) %>%
    select(GeneID, PSI_earliest = PSI_value)
  
  print(df_earliest)
  
  #（d）将其他时期的数据与最早时期合并
  df_compare <- df_valid %>%
    filter(!(Stage == earliest_stage$Stage[1] & `Stage-xifen` == earliest_stage$`Stage-xifen`[1])) %>%
    left_join(df_earliest, by = "GeneID") %>%
    filter(!is.na(PSI_earliest))  # 去掉参考 PSI 为空的行
  
  print(df_compare)
  
  #（e）如果没有其他时期可比较，则只返回最早时期
  if (nrow(df_compare) == 0) {
    return(
      tibble(
        Stage      = earliest_stage$Stage[1],
        `Stage-xifen` = earliest_stage$`Stage-xifen`[1],
        Tissue     = unique(df_grouped$Tissue),
        cor_median = 1,
        cor_lower  = 1,
        cor_upper  = 1,
        n_genes    = nrow(df_earliest)
      )
    )
  }
  
  #（f）对每个细分时期做 bootstrap 相关
  cor_boot <- df_compare %>%
    group_by(Stage, `Stage-xifen`) %>%
    group_modify(~ {
      df_stage <- .x
      cor_vec <- numeric(n_boot)
      for (i in seq_len(n_boot)) {
        idx <- sample(seq_len(nrow(df_stage)), size = nrow(df_stage), replace = TRUE)
        sample_df <- df_stage[idx, ]
        if (length(unique(sample_df$PSI_value)) < 2 || length(unique(sample_df$PSI_earliest)) < 2) {
          cor_vec[i] <- NA_real_
        } else {
          cor_vec[i] <- cor(sample_df$PSI_value, sample_df$PSI_earliest, use = "pairwise.complete.obs")
        }
      }
      tibble(
        Tissue      = unique(df_stage$Tissue),
        cor_median  = median(cor_vec, na.rm = TRUE),
        cor_lower   = quantile(cor_vec, probs = 0.025, na.rm = TRUE),
        cor_upper   = quantile(cor_vec, probs = 0.975, na.rm = TRUE),
        n_genes     = nrow(df_stage)
      )
    })
  print(cor_boot)
  
  #（g）把最早时期也加回来
  cor_boot_earliest <- tibble(
    Stage      = earliest_stage$Stage[1],
    `Stage-xifen` = earliest_stage$`Stage-xifen`[1],
    Tissue     = unique(df_compare$Tissue),
    cor_median = 1,
    cor_lower  = 1,
    cor_upper  = 1,
    n_genes    = nrow(df_earliest)
  )
  
  result <- bind_rows(cor_boot_earliest, cor_boot) %>%
    arrange(Stage, `Stage-xifen`)
  
  print(result)
  return(result)
}

########################################
## 3. 对每个 Tissue 分组，做相关分析
########################################
cor_all <- data_long_parsed2 %>%
  group_by(Tissue) %>%
  nest() %>%
  mutate(
    cor_table = map(
      data,
      ~ get_earliest_stage_cor_bootstrap(.x, n_boot = 1000, seed = 123)
    )
  ) %>%
  select(-data) %>%
  unnest(cols = cor_table)


# cor_all 结构类似：
#  Tissue | Stage | cor_median | cor_lower | cor_upper | n_genes
#  Head   | Egg   | 1          | 1         | 1         | NA
#  Head   | Larva | 0.935      | 0.921     | 0.949     | 32000
#  ...    | ...   | ...        | ...       | ...       | ...


# 先数一下每个Tissue的阶段有多少个
tissue_stage_count <- cor_all %>%
  group_by(Tissue) %>%
  summarise(n_stage = n_distinct(Stage))

# 找到至少有3个不同Stage的组织
valid_tissue <- tissue_stage_count %>%
  filter(n_stage >= 2) %>%
  pull(Tissue)

# 筛选 cor_all
cor_all <- cor_all %>%
  filter(Tissue %in% valid_tissue)


########################################
## 4. 画图：不同 Tissue 多条线
##    x=Stage, y=相关系数, 彩色区分 Tissue。
##    带置信区间带 (geom_ribbon)
########################################
cor_all <- cor_all %>%
  mutate(Combined_Stage = paste(Stage, `Stage-xifen`, sep = "-"))

cor_all <- cor_all %>%
  mutate(
    Combined_Stage = factor(
      Combined_Stage,
      levels = unique(data_long_parsed2$Combined_Stage)  # 确保顺序一致
    )
  )



custom_order <- c(
  "Egg-0_days", "Egg-20-22_hours", "Egg-24_hours", "Egg-1_day", "Egg-2_days",  "Egg-36_hours", "Egg-3_days", "Egg-72_hours", "Egg-120_hours","Egg-168_hours",  "Egg-192_hours", "Egg-216_hours", 
  "Larva-0_hours", "Larva-1_day_of_1_instar", "Larva-48_hours", "Larva-96_hours" , "Larva-1_instar", "Larva-2_instar", "Larva-0_day_of_3_instar", "Larva-24_hours_of_4_instar", "Larva-1_day_of_4_instar", "Larva-2_day_of_4_instar", "Larva-3_day_of_4_instar", "Larva-4_day_of_4_instar", "Larva-4_instar" , "Larva-10_day_of_4_instar", "Larva-11_day_of_4_instar", "Larva-12_day_of_4_instar", "Larva-0_day_of_5_instar", "Larva-1_day_of_5_instar", "Larva-2_day_of_5_instar", "Larva-3_day_of_5_instar", "Larva-4_day_of_5_instar", "Larva-5_instar","Larva-6_day_of_5_instar", "Larva-7_day_of_5_instar", "Larva-8_day_of_5_instar" , "Larva-9_day_of_5_instar", "Larva-10_day_of_5_instar", "Larva-12_day_of_5_instar", "Larva-13_day_of_5_instar", "Larva-14_day_of_5_instar", "Larva-15_day_of_5_instar", "Larva-wandering_stage", "Larva-36h_wandering_stage","Larva-52h_wandering_stage", 
  "Pupa-stage_5", 
  "Adult-0h","Adult-1_day", "Adult-2_days", "Adult-4_days","Adult-5_days"
)
# 从数据中过滤掉不需要的时期
data_long_parsed2 <- data_long_parsed2 %>%
  filter(Combined_Stage %in% custom_order)

cor_all <- cor_all %>%
  filter(Combined_Stage %in% custom_order)
data_long_parsed2 <- data_long_parsed2 %>%
  mutate(
    Combined_Stage = factor(
      Combined_Stage,
      levels = custom_order  # 手动指定顺序
    )
  )

cor_all <- cor_all %>%
  mutate(
    Combined_Stage = factor(
      Combined_Stage,
      levels = custom_order  # 保持顺序一致
    )
  )

# 设置置信区间范围阈值
threshold <- 0.5

# 剔除置信区间范围过大的点
cor_all_filtered <- cor_all %>%
  filter((cor_upper - cor_lower) <= threshold | is.na(cor_upper - cor_lower))

library(zoo)

# 对每个 Tissue 分组，进行插值平滑处理
cor_all_smoothed <- cor_all_filtered %>%
  group_by(Tissue) %>%
  mutate(
    cor_median_smooth = na.approx(cor_median, na.rm = FALSE),  # 平滑中位数
    cor_lower_smooth = na.approx(cor_lower, na.rm = FALSE),    # 平滑下界
    cor_upper_smooth = na.approx(cor_upper, na.rm = FALSE)     # 平滑上界
  ) %>%
  ungroup()




p <- # 用平滑后的数据进行绘图
  ggplot(cor_all_smoothed, aes(x = Combined_Stage, y = cor_median_smooth, group = Tissue, color = Tissue)) +
  geom_ribbon(aes(ymin = cor_lower_smooth, ymax = cor_upper_smooth, fill = Tissue), alpha = 0.2, color = NA) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 14) +
  labs(
    title = "PSI correlation with filtered and smoothed confidence intervals",
    x = "Stage",
    y = "Pearson correlation (median ± 95% CI)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

ggsave("correlation_plot.png", p, width = 10, height = 6, dpi = 300)

# 
# library(tidyr)
# library(pheatmap)
# 
# ##################################################
# ## 2. 热图 - 以 GeneID 为行，Stage 为列
# ##    (示例：筛选 Tissue = "Whole_body")
# ##################################################
# data_long_parsed2 <- data_long_parsed2[!is.na(data_long_parsed2$Stage), ]
# 
# 
# 
# p_grouped_box <- data_long_parsed2 %>%
#   ggplot(aes(x = Stage, y = PSI_value, fill = Tissue)) +
#   geom_violin(position=position_dodge(width=0.8), alpha=0.4) +
#   geom_boxplot(position = position_dodge(width=0.8), outlier.alpha = 0.3) +
#   theme_bw(base_size = 14) +
#   labs(
#     title = "Grouped Boxplot of PSI: different Tissues across Stages",
#     x = "Stage",
#     y = "PSI"
#   )
# print(p_grouped_box)
# 
# ggsave("grouped_boxplot.png", p_grouped_box, width = 10, height = 6, dpi = 300)
