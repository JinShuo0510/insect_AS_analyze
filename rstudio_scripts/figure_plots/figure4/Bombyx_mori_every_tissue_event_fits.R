library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(splines)
library(broom)
library(progress)

# 定义对指定组织（tissue）进行事件拟合与 DevAS 分析的函数
run_event_analysis <- function(tissue_name) {
  
  message("开始处理组织: ", tissue_name)
  
  # 筛选该组织的数据，并统一重命名 PSI 列为 PSI
  df_tissue <- data_long_parsed2 %>%
    filter(Tissue == tissue_name) %>%
    select(GeneID, chr, strand, exonStart_0base, exonEnd,
           upstreamES, upstreamEE, downstreamES, downstreamEE,
           Combined_Stage, PSI_value) %>%
    rename(PSI = PSI_value)
  
  # 阶段到时间（小时）的映射
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
  
  # 将 Combined_Stage 转换为数值型时间索引，并将小时转为天
  df_tissue <- df_tissue %>%
    mutate(time_index = as.numeric(recode(as.character(Combined_Stage), !!!time_index_mapping)),
           day = time_index / 24)
  
  # 针对每个可变剪切事件构建拟合模型
  event_fits <- df_tissue %>%
    group_by(GeneID, chr, strand, exonStart_0base, exonEnd,
             upstreamES, upstreamEE, downstreamES, downstreamEE) %>%
    nest() %>%
    mutate(
      model = map(data, ~ glm(PSI ~ bs(log(day), df = 4),
                              family = quasibinomial, data = .x)),
      newdata = map(data, ~ {
        new_day <- seq(min(.x$day), max(.x$day), length.out = 1000)
        data.frame(day = new_day)
      }),
      PSI_pred = map2(model, newdata, ~ predict(.x, newdata = .y, type = "response"))
    ) %>%
    mutate(
      dPSI = map_dbl(PSI_pred, ~ max(.x) - min(.x))
    )
  
  # 构造事件标识符（便于后续对比）
  df_tissue <- df_tissue %>%
    mutate(event_id = paste(GeneID, chr, strand, exonStart_0base, exonEnd,
                            upstreamES, upstreamEE, downstreamES, downstreamEE, sep = "__"))
  
  event_fits <- event_fits %>%
    mutate(event_id = paste(GeneID, chr, strand, exonStart_0base, exonEnd,
                            upstreamES, upstreamEE, downstreamES, downstreamEE, sep = "__")) %>%
    mutate(
      delta_PSI = map(PSI_pred, diff),
      up = map_dbl(delta_PSI, ~ sum(.x[.x > 0], na.rm = TRUE)),
      down = map_dbl(delta_PSI, ~ sum(abs(.x[.x < 0]), na.rm = TRUE)),
      time_vec = map(newdata, ~ .x$day),
      up_timing = map2_dbl(delta_PSI, time_vec, ~ {
        pos_inds <- which(.x > 0)
        if(length(pos_inds) > 0) {
          sum(.x[pos_inds] * .y[pos_inds]) / sum(.x[pos_inds])
        } else {
          NA_real_
        }
      }),
      down_timing = map2_dbl(delta_PSI, time_vec, ~ {
        neg_inds <- which(.x < 0)
        if(length(neg_inds) > 0) {
          sum(abs(.x[neg_inds]) * .y[neg_inds]) / sum(abs(.x[neg_inds]))
        } else {
          NA_real_
        }
      }),
      ratio = up / (up + down),
      direction = case_when(
        ratio > 0.7 ~ "up",
        ratio < 0.3 ~ "down",
        TRUE ~ if_else(up_timing < down_timing, "up-down", "down-up")
      )
    )
  
  # 计算 DevAS 事件：构造 null 模型并用 anova 比较全模型与 null 模型
  event_fits <- event_fits %>%
    mutate(
      null_model = map(data, ~ glm(PSI ~ 1, family = quasibinomial, data = .x)),
      raw_p = map2_dbl(null_model, model, ~ {
        pval <- NA_real_
        aov_out <- tryCatch(anova(.x, .y, test = "F"), error = function(e) NULL)
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
  
  message("完成处理组织: ", tissue_name)
  return(event_fits)
}

# 获取 data_long_parsed2 中所有 Tissue 的唯一值
tissues <- unique(data_long_parsed2$Tissue)
tissues <- tissues[tissues != "Whole_body"]

# 创建进度条
pb <- progress_bar$new(total = length(tissues),
                       format = "进度: [:bar] :current/:total (:percent)")

# 遍历每个 Tissue，运行分析并将结果保存到列表中
results <- map(tissues, function(tissue) {
  res <- run_event_analysis(tissue)
  pb$tick()  # 更新进度条
  return(res)
})
names(results) <- tissues

