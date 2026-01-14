library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(splines)

# 读取并整理数据（你已有前置代码）
# all_event_wide_lists_withoutavg_1_8 <- readRDS("all_event_wide_lists_withoutavg_1_8.RDS")
data_list <- all_event_wide_lists_withoutavg_1_8[["SE"]][["Bombyx_mori"]]
meta_cols <- c("GeneID","geneSymbol","chr","strand",
               "exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")
psi_cols <- setdiff(colnames(data_list), meta_cols)

custom_order <- c(
  "Egg-0_days", "Egg-20-22_hours", "Egg-24_hours", "Egg-1_day", "Egg-2_days", "Egg-36_hours", 
  "Egg-3_days", "Egg-72_hours", "Egg-120_hours", "Egg-168_hours", "Egg-192_hours", "Egg-216_hours", 
  "Larva-0_hours", "Larva-1_day_of_1_instar", "Larva-48_hours", "Larva-96_hours", "Larva-1_instar", 
  "Larva-2_instar", "Larva-0_day_of_3_instar", "Larva-24_hours_of_4_instar", "Larva-1_day_of_4_instar", 
  "Larva-2_day_of_4_instar", "Larva-3_day_of_4_instar", "Larva-4_day_of_4_instar", "Larva-4_instar", 
  "Larva-10_day_of_4_instar", "Larva-11_day_of_4_instar", "Larva-12_day_of_4_instar", "Larva-0_day_of_5_instar", 
  "Larva-1_day_of_5_instar", "Larva-2_day_of_5_instar", "Larva-3_day_of_5_instar", "Larva-4_day_of_5_instar", 
  "Larva-5_instar", "Larva-6_day_of_5_instar", "Larva-7_day_of_5_instar", "Larva-8_day_of_5_instar", 
  "Larva-9_day_of_5_instar", "Larva-10_day_of_5_instar", "Larva-12_day_of_5_instar", "Larva-13_day_of_5_instar", 
  "Larva-14_day_of_5_instar", "Larva-15_day_of_5_instar", "Larva-wandering_stage", "Larva-36h_wandering_stage",
  "Larva-52h_wandering_stage", "Pupa-stage_5", 
  "Adult-0h", "Adult-1_day", "Adult-2_days", "Adult-4_days", "Adult-5_days"
)


data_long_parsed2 <- data_list %>%
  select(all_of(meta_cols), all_of(psi_cols)) %>%
  pivot_longer(
    cols      = all_of(psi_cols),
    names_to  = "Stage_Tissue",
    values_to = "PSI_string"
  ) %>%
  separate_rows(PSI_string, sep = ",") %>%
  mutate(PSI_value = as.numeric(PSI_string)) %>%
  filter(!is.na(PSI_value)) %>%
  group_by(GeneID, geneSymbol, chr, strand,
           exonStart_0base, exonEnd,
           upstreamES, upstreamEE, downstreamES, downstreamEE,
           Stage_Tissue) %>%    # 将所有meta信息也纳入分组
  summarise(PSI_value = mean(PSI_value), .groups="drop") %>%
  separate(
    col  = Stage_Tissue,
    into = c("Stage", "Stage_xifen", "Tissue"),
    sep  = "__",
    remove = FALSE
  ) %>%
  mutate(
    Stage = if_else(Stage == "Cell", "Egg", Stage),
    Tissue = if_else(Tissue %in% c("Egg","Pupa"), "Whole_body", Tissue),
    Stage = factor(Stage, levels = c("Egg","Larva","Pupa","Adult")),
    Combined_Stage = paste(Stage, Stage_xifen, sep = "-")
  ) %>%
  filter(Combined_Stage %in% custom_order) %>%
  mutate(Combined_Stage = factor(Combined_Stage, levels = custom_order)) %>%
  arrange(Stage, Stage_xifen)

saveRDS(data_long_parsed2, file="Bombyx_mori_data_long_parsed2.RDS")


# 提取 Whole_body 作为基准数据
df_whole <- data_long_parsed2 %>%
  filter(Tissue == "Whole_body") %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, Combined_Stage, PSI_value) %>%
  rename(PSI_whole = PSI_value) 


# 构建阶段到时间（小时）的映射（根据前面生成的索引）
time_index_mapping <- c(
  "Egg-0_days" = 2,
  "Egg-20-22_hours" = 21,
  "Egg-24_hours" = 24,
  "Egg-1_day" = 24,
  "Egg-2_days" = 48,
  "Egg-36_hours" = 36,
  "Egg-3_days" = 72,
  "Egg-72_hours" = 72,
  "Egg-120_hours" = 120,
  "Egg-168_hours" = 168,
  "Egg-192_hours" = 192,
  "Egg-216_hours" = 216,
  "Larva-0_hours" = 216,
  "Larva-1_day_of_1_instar" = 240,
  "Larva-48_hours" = 264,
  "Larva-96_hours" = 312,
  "Larva-1_instar" = 312,
  "Larva-2_instar" = 360,
  "Larva-0_day_of_3_instar" = 408,
  "Larva-24_hours_of_4_instar" = 456,
  "Larva-1_day_of_4_instar" = 480,
  "Larva-2_day_of_4_instar" = 504,
  "Larva-3_day_of_4_instar" = 528,
  "Larva-4_day_of_4_instar" = 552,
  "Larva-4_instar" = 552,
  "Larva-10_day_of_4_instar" = 696,
  "Larva-11_day_of_4_instar" = 720,
  "Larva-12_day_of_4_instar" = 744,
  "Larva-0_day_of_5_instar" = 768,
  "Larva-1_day_of_5_instar" = 792,
  "Larva-2_day_of_5_instar" = 816,
  "Larva-3_day_of_5_instar" = 840,
  "Larva-4_day_of_5_instar" = 864,
  "Larva-5_instar" = 864,
  "Larva-6_day_of_5_instar" = 912,
  "Larva-7_day_of_5_instar" = 936,
  "Larva-8_day_of_5_instar" = 960,
  "Larva-9_day_of_5_instar" = 984,
  "Larva-10_day_of_5_instar" = 1008,
  "Larva-12_day_of_5_instar" = 1056,
  "Larva-13_day_of_5_instar" = 1080,
  "Larva-14_day_of_5_instar" = 1104,
  "Larva-15_day_of_5_instar" = 1128,
  "Larva-wandering_stage" = 1152,
  "Larva-36h_wandering_stage" = 1188,
  "Larva-52h_wandering_stage" = 1204,
  "Pupa-stage_5" = 1248,
  "Adult-0h" = 1248,
  "Adult-1_day" = 1272,
  "Adult-2_days" = 1296,
  "Adult-4_days" = 1344,
  "Adult-5_days" = 1368
)

# 将 Combined_Stage 转换为数值时间索引
df_whole <- df_whole %>%
  mutate(time_index = as.numeric(recode(as.character(Combined_Stage), !!!time_index_mapping)))

# ---------------------------
# 转换时间单位及构造自变量：以天为单位，并对天数取对数
# ---------------------------
df_whole <- df_whole %>%
  mutate(day = time_index / 24)  # 将小时转换为天

# ---------------------------
# 针对每个可变剪切事件分别拟合三次样条模型
# 使用对数（天）作为自变量，自由度设为4
# ---------------------------
event_fits <- df_whole %>%
  group_by(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(PSI_whole ~ bs(log(day), df = 4),
                            family = quasibinomial, data = .x)),
    newdata = map(data, ~ {
      new_day <- seq(min(.x$day), max(.x$day), length.out = 1000)
      data.frame(day = new_day)
    }),
    PSI_pred = map2(model, newdata, ~ predict(.x, newdata = .y, type = "response"))
  )




# 计算每个事件的 dPSI（预测 PSI 的最大值与最小值之差）
event_fits <- event_fits %>%
  mutate(
    dPSI = map_dbl(PSI_pred, ~ max(.x) - min(.x))
  )

df_whole <- df_whole %>%
  mutate(event_id = paste(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, sep = "__"))

event_fits <- event_fits %>%
  mutate(event_id = paste(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, sep = "__"))

# 计算每个事件的 devAS 变化方向
event_fits <- event_fits %>%
  mutate(
    # 计算连续时间点之间的 PSI 差值
    delta_PSI = map(PSI_pred, diff),
    
    # 统计所有正差值的总和（up）和负差值绝对值的总和（down）
    up = map_dbl(delta_PSI, ~ sum(.x[.x > 0], na.rm = TRUE)),
    down = map_dbl(delta_PSI, ~ sum(abs(.x[.x < 0]), na.rm = TRUE)),
    
    # 提取每个事件对应的时间向量
    time_vec = map(newdata, ~ .x$day),
    
    # 计算正变化的时间加权均值 (up_timing)
    up_timing = map2_dbl(delta_PSI, time_vec, ~ {
      pos_inds <- which(.x > 0)
      if(length(pos_inds) > 0) {
        sum(.x[pos_inds] * .y[pos_inds]) / sum(.x[pos_inds])
      } else {
        NA_real_
      }
    }),
    
    # 计算负变化的时间加权均值 (down_timing)
    down_timing = map2_dbl(delta_PSI, time_vec, ~ {
      neg_inds <- which(.x < 0)
      if(length(neg_inds) > 0) {
        sum(abs(.x[neg_inds]) * .y[neg_inds]) / sum(abs(.x[neg_inds]))
      } else {
        NA_real_
      }
    }),
    
    # 计算 up 与 down 总变化的比例
    ratio = up / (up + down),
    
    # 根据比例和时间加权均值判断变化方向：
    # 如果 ratio > 0.7，则主要上调 ("up")
    # 如果 ratio < 0.3，则主要下调 ("down")
    # 否则，根据 up_timing 和 down_timing 判断：若 up_timing < down_timing，则先上调后下调 ("up-down")，反之 ("down-up")
    direction = case_when(
      ratio > 0.7 ~ "up",
      ratio < 0.3 ~ "down",
      TRUE ~ if_else(up_timing < down_timing, "up-down", "down-up")
    )
  )




#!--------------------计算DevAS事件
library(broom)  # 用于整洁模型输出

# 假设 event_fits 已经存在，且每个事件的模型存储在 model 列中

# 计算DevAS，这里选取的是整个模型的所有值
event_fits <- event_fits %>%
  mutate(
    null_model = map(data, ~ glm(PSI_whole ~ 1, family = quasibinomial, data = .x)),
    
    raw_p = map2_dbl(null_model, model, ~ {
      pval <- NA_real_
      aov_out <- tryCatch(anova(.x, .y, test = "F"), error = function(e) NULL)
      
      # 如果 anova 返回有效结果，且包含至少2行（null 和 full 模型），才提取 p 值
      if (!is.null(aov_out) && nrow(aov_out) >= 2) {
        test_val <- aov_out$`Pr(>F)`[2]
        if (!is.null(test_val) && length(test_val) == 1 && !is.na(test_val)) {
          pval <- test_val
        }
      }
      
      return(pval)
    }),
    
    adj_p = p.adjust(raw_p, method = "BH")
  )


saveRDS(event_fits,"Bombyx_mori_Whole_body_event_fits.RDS")




# 3. 筛选出符合 devAS 标准的事件（例如 dPSI > 0.2 且 adj_p < 0.05）
devAS_events <- event_fits %>%
  filter(dPSI > 0.2, adj_p < 0.05) %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE,
         downstreamES, downstreamEE, dPSI, direction, raw_p, adj_p)

# 查看筛选结果
print(devAS_events)









# 示例：绘制某个特定事件（GeneID）的拟合结果
selected_event <- "Ck2a__chrNC_085111.1__+__9743077__9743299__9742889__9742968__9743432__9743716"  # 请替换为你感兴趣的事件ID
df_plot <- df_whole %>% filter(event_id == selected_event) 
pred_data <- tibble(day = seq(min(df_plot$day), max(df_plot$day), length.out = 1000),
                    PSI_pred = event_fits %>% filter(event_id == selected_event) %>% pull(PSI_pred) %>% .[[1]])

# 绘制拟合结果
p <- ggplot() +
  geom_point(data = df_plot, aes(x = day, y = PSI_whole)) +
  geom_line(data = pred_data, aes(x = day, y = PSI_pred), color = "blue") +
  labs(title = paste("事件", selected_event, "的基于对数(天)自变量的三次样条拟合"),
       x = "Time Index", y = "PSI")

print(p)





