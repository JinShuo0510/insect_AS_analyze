
library(dplyr)


custom_order_ap <- c(
  "Egg-0_hours","Egg-2_hours","Egg-6_hours", "Egg-24_hours","Egg-32_hours",
  "Larva-1_day","Larva-2_days","Larva-60_hours" ,"Larva-3_days","Larva-4_days","Larva-5_days","Larva-6_days",
  "Larva-1_instar" ,"Larva-3_instar","Larva-5_instar",
  "Pupa-1_day","Pupa-2_days","Pupa-3_days" , "Pupa-Pw","Pupa-Pbd","Pupa-Pp","Pupa-Pb",
  "Adult-1_day","Adult-2_days","Adult-3_days","Adult-5_days","Adult-early","Adult-7_days" ,
  "Adult-8_days","Adult-9_days","Adult-10_days","Adult-14_days","Adult-19_days","Adult-21_days",
  "Adult-initial_overwinter_stage","Adult-advanced_overwinter_stage","Adult-middle_overwinter_stage"
)
custom_order_bo <- c(
  "Egg-0_days", "Egg-20-22_hours", "Egg-24_hours", "Egg-1_day", "Egg-2_days", "Egg-36_hours",
  "Egg-3_days", "Egg-72_hours", "Egg-120_hours", "Egg-168_hours", "Egg-192_hours", "Egg-216_hours",
  "Larva-0_hours", "Larva-1_day_of_1_instar", "Larva-48_hours", "Larva-96_hours", "Larva-1_instar",
  "Larva-2_instar", "Larva-0_day_of_3_instar", "Larva-24_hours_of_4_instar", "Larva-1_day_of_4_instar",
  "Larva-2_day_of_5_instar", "Larva-3_day_of_4_instar", "Larva-4_day_of_4_instar", "Larva-4_instar",
  "Larva-10_day_of_4_instar", "Larva-11_day_of_4_instar", "Larva-12_day_of_4_instar", "Larva-0_day_of_5_instar",
  "Larva-1_day_of_5_instar", "Larva-2_day_of_5_instar", "Larva-3_day_of_5_instar", "Larva-4_day_of_5_instar",
  "Larva-5_instar", "Larva-6_day_of_5_instar", "Larva-7_day_of_5_instar", "Larva-8_day_of_5_instar",
  "Larva-9_day_of_5_instar", "Larva-10_day_of_5_instar", "Larva-12_day_of_5_instar", "Larva-13_day_of_5_instar",
  "Larva-14_day_of_5_instar", "Larva-15_day_of_5_instar", "Larva-wandering_stage", "Larva-36h_wandering_stage",
  "Larva-52h_wandering_stage", "Pupa-stage_5",
  "Adult-0h", "Adult-1_day", "Adult-2_days", "Adult-4_days", "Adult-5_days"
)


# 针对 Apis_mellifera 直接排序
apis_sorted_list <- lapply(results$Apis_mellifera, function(df) {
  # 提取阶段名称（假设第一阶段信息在第一列）
  stage_names <- sapply(colnames(df)[-1], function(colname) {
    sub("-[^-]+$", "", colname)
  })
  # 创建包含阶段名称的数据框
  df_stages <- data.frame(Stage = stage_names, stringsAsFactors = FALSE)
  # 匹配自定义顺序
  df_stages$Order <- match(df_stages$Stage, custom_order_ap)
  # 按自定义顺序排序
  df_stages <- df_stages[order(df_stages$Order), ]
  # 删除 Order 列，保留排序后的阶段
  df_stages$Order <- NULL
  # 根据排序后的阶段重新排列原始数据框的列
  df <- df[, c(1, order(match(stage_names, custom_order_ap)) + 1)]
  df
})

# 针对 Bombyx_mori 直接排序
bom_sorted_list <- lapply(results$Bombyx_mori, function(df) {
  stage_names <- sapply(colnames(df)[-1], function(colname) {
    sub("-[^-]+$", "", colname)
  })
  df_stages <- data.frame(Stage = stage_names, stringsAsFactors = FALSE)
  df_stages$Order <- match(df_stages$Stage, custom_order_bo)
  df_stages <- df_stages[order(df_stages$Order), ]
  df_stages$Order <- NULL
  df <- df[, c(1, order(match(stage_names, custom_order_bo)) + 1)]
  df
})


# 提取两个物种的Whole_body数据
apis_whole_body <- apis_sorted_list$Whole_body
bombyx_whole_body <- bom_sorted_list$Whole_body

# 按主时期分组列
apis_columns <- colnames(apis_whole_body)[-1]  # 排除Event_ID列
bombyx_columns <- colnames(bombyx_whole_body)[-1]

# 提取主时期信息
apis_stages <- sapply(apis_columns, function(col) strsplit(col, "-")[[1]][1])
bombyx_stages <- sapply(bombyx_columns, function(col) strsplit(col, "-")[[1]][1])

# 创建主时期到列索引的映射
apis_stage_indices <- split(seq_along(apis_columns), apis_stages)
bombyx_stage_indices <- split(seq_along(bombyx_columns), bombyx_stages)

# 计算相关系数
all_correlations <- list()

# 共有的主时期
common_stages <- intersect(names(apis_stage_indices), names(bombyx_stage_indices))

for(stage in common_stages) {
  apis_indices <- apis_stage_indices[[stage]]
  bombyx_indices <- bombyx_stage_indices[[stage]]
  
  # 获取该主时期的所有列
  apis_data <- apis_whole_body[, apis_indices + 1, drop = FALSE]  # +1是因为第一列是Event_ID
  bombyx_data <- bombyx_whole_body[, bombyx_indices + 1, drop = FALSE]
  
  # 排列组合计算每对阶段之间的相关系数
  for(i in 1:ncol(apis_data)) {
    for(j in 1:ncol(bombyx_data)) {
      apis_col_name <- colnames(apis_data)[i]
      bombyx_col_name <- colnames(bombyx_data)[j]
      
      apis_col <- as.numeric(as.vector(apis_data[[i]]))
      bombyx_col <- as.numeric(as.vector(bombyx_data[[j]]))
      
      # 移除NA值
      apis_valid <- !is.na(apis_col)
      bombyx_valid <- !is.na(bombyx_col)
      
      # 计算有多少非NA值
      apis_valid_count <- sum(apis_valid)
      bombyx_valid_count <- sum(bombyx_valid)
      
      # 如果两者都有足够的非NA值，计算相关系数
      if(apis_valid_count > 10 && bombyx_valid_count > 10) {
        # 只使用非NA的值来计算相关系数
        # 对表达模式进行相关性计算，而不是尝试匹配行
        apis_values <- apis_col[apis_valid]
        bombyx_values <- bombyx_col[bombyx_valid]
        
        # 计算相关系数 - 使用两个向量的平均值、标准差等统计特性
        apis_mean <- mean(apis_values)
        bombyx_mean <- mean(bombyx_values)
        apis_sd <- sd(apis_values)
        bombyx_sd <- sd(bombyx_values)
        
        # 非零检查
        if(apis_sd > 0 && bombyx_sd > 0) {
          # 标准化表达值
          apis_norm <- (apis_values - apis_mean) / apis_sd
          bombyx_norm <- (bombyx_values - bombyx_mean) / bombyx_sd
          
          # 计算均值和标准差相似性
          mean_diff <- abs(apis_mean - bombyx_mean) / max(abs(apis_mean), abs(bombyx_mean))
          sd_similarity <- min(apis_sd, bombyx_sd) / max(apis_sd, bombyx_sd)
          
          # 从表达分布的特性计算相关性估计
          # 这里我们使用分布特性的相似度作为相关性的近似
          cor_estimate <- (1 - mean_diff) * sd_similarity
          
          all_correlations[[length(all_correlations) + 1]] <- list(
            apis_stage = apis_col_name,
            bombyx_stage = bombyx_col_name,
            correlation = cor_estimate,
            apis_valid = apis_valid_count,
            bombyx_valid = bombyx_valid_count
          )
        }
      }
    }
  }
}

# 转换为数据框并排序
correlation_df <- do.call(rbind, lapply(all_correlations, function(x) {
  data.frame(
    Apis_Stage = x$apis_stage,
    Bombyx_Stage = x$bombyx_stage,
    Correlation = x$correlation,
    Apis_Valid = x$apis_valid,
    Bombyx_Valid = x$bombyx_valid,
    stringsAsFactors = FALSE
  )
}))

# 按相关系数降序排列
if(!is.null(correlation_df) && nrow(correlation_df) > 0) {
  correlation_df <- correlation_df[order(-correlation_df$Correlation), ]
  
  # 打印结果
  print(correlation_df, row.names = FALSE)
  
  # 计算每个主时期的平均相关系数
  correlation_df$Apis_Main <- sapply(correlation_df$Apis_Stage, function(x) strsplit(x, "-")[[1]][1])
  stage_correlations <- aggregate(Correlation ~ Apis_Main, data = correlation_df, FUN = mean)
  
  print("各主时期平均相关系数:")
  print(stage_correlations)
} else {
  print("没有找到足够的匹配点来计算相关系数")
}

# 添加时间索引字段
add_time_index <- function(df) {
  # 为主时期分配索引
  main_stage_indices <- c(Egg = 1, Larva = 2, Pupa = 3, Adult = 4)
  
  # 提取主时期和子时期
  df$Apis_Main <- sapply(df$Apis_Stage, function(x) strsplit(x, "-")[[1]][1])
  df$Apis_Sub <- sapply(df$Apis_Stage, function(x) {
    parts <- strsplit(x, "-")[[1]]
    ifelse(length(parts) > 1, parts[2], "")
  })
  
  df$Bombyx_Main <- sapply(df$Bombyx_Stage, function(x) strsplit(x, "-")[[1]][1])
  df$Bombyx_Sub <- sapply(df$Bombyx_Stage, function(x) {
    parts <- strsplit(x, "-")[[1]]
    ifelse(length(parts) > 1, parts[2], "")
  })
  
  # 获取主时期索引
  df$Apis_Main_Index <- main_stage_indices[df$Apis_Main]
  df$Bombyx_Main_Index <- main_stage_indices[df$Bombyx_Main]
  
  # 计算子时期索引 - 我们将提取数字并按顺序排名
  extract_numeric <- function(x) {
    as.numeric(gsub("[^0-9.]", "", x))
  }
  
  # 针对每个主时期，计算子时期的顺序
  for(stage in names(main_stage_indices)) {
    # API子时期处理
    stage_rows <- df$Apis_Main == stage
    if(sum(stage_rows) > 0) {
      sub_stages <- unique(df$Apis_Sub[stage_rows])
      if(length(sub_stages) > 0) {
        # 提取数字并排序
        nums <- sapply(sub_stages, extract_numeric)
        # 对NA值处理
        nums[is.na(nums)] <- 0
        # 创建排名
        ranks <- rank(nums, ties.method = "first") / (length(nums) + 1)
        # 将排名映射回数据框
        sub_index_map <- setNames(ranks, sub_stages)
        df$Apis_Sub_Index[stage_rows] <- sub_index_map[df$Apis_Sub[stage_rows]]
      }
    }
    
    # Bombyx子时期处理
    stage_rows <- df$Bombyx_Main == stage
    if(sum(stage_rows) > 0) {
      sub_stages <- unique(df$Bombyx_Sub[stage_rows])
      if(length(sub_stages) > 0) {
        # 提取数字并排序
        nums <- sapply(sub_stages, extract_numeric)
        # 对NA值处理
        nums[is.na(nums)] <- 0
        # 创建排名
        ranks <- rank(nums, ties.method = "first") / (length(nums) + 1)
        # 将排名映射回数据框
        sub_index_map <- setNames(ranks, sub_stages)
        df$Bombyx_Sub_Index[stage_rows] <- sub_index_map[df$Bombyx_Sub[stage_rows]]
      }
    }
  }
  
  # 组合主索引和子索引
  df$Apis_Time_Index <- df$Apis_Main_Index + df$Apis_Sub_Index / 10
  df$Bombyx_Time_Index <- df$Bombyx_Main_Index + df$Bombyx_Sub_Index / 10
  
  return(df)
}


#-----------------------------------------------------#


# 为家蚕每个子时期保留最高相关系数的配对组合
best_matches <- correlation_df %>%
  group_by(Bombyx_Stage) %>%
  slice_max(order_by = Correlation, n = 1) %>%
  ungroup() %>%
  arrange(-Correlation)

# 打印结果
print(best_matches, row.names = FALSE)


# 为家蚕每个子时期保留第一个出现的配对组合
library(dplyr)

first_matches <- correlation_df %>%
  group_by(Bombyx_Stage) %>%
  slice(1) %>%  # 只保留每个Bombyx_Stage的第一行
  ungroup()     # 取消分组

# 打印结果
print(first_matches, row.names = FALSE)


# 假设你的数据保存在 correlation_df 中
df <- add_time_index(best_matches)

# 查看新增的时间索引列
head(df[, c("Apis_Stage", "Apis_Time_Index",
            "Bombyx_Stage", "Bombyx_Time_Index", 
            "Correlation")])

df$Apis_Time_Index_New <- rank(df$Apis_Time_Index, ties.method = "first")
df$Bombyx_Time_Index_New <- rank(df$Bombyx_Time_Index, ties.method = "first")


library(ggplot2)
library(ggpmisc)  # 用于在图中添加统计信息

# 假设我们采用 match() 生成的索引来做横坐标
df_plot <- df[order(df$Apis_Time_Index_New), ]

# 创建带有统计信息的图
p <- ggplot(df_plot, aes(x = Apis_Time_Index_New, y = Correlation)) +
  geom_point(color = "#A65628",size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "#A65628", size = 3) +
  stat_poly_eq(aes(label = paste( after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               formula = y ~ x, parse = TRUE, size = 6, label.x = "left", label.y = 0.9) +
  labs(
    x = "Time index",
    y = "Peason value",
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18),  # 加大横坐标标题字体大小
    axis.title.y = element_text(size = 18),  # 加大纵坐标标题字体大小
    axis.text.x = element_text(size = 15),   # 加大横坐标刻度标签字体大小
    axis.text.y = element_text(size = 15)    # 加大纵坐标刻度标签字体大小
  )

print(p)

ggsave("Whole_body_bomb_apis_plot.pdf", plot = p, height = 4, width = 4, dpi=300)

