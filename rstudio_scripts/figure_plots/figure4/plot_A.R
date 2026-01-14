
library(tidyverse)
library(ggplot2)
library(RColorBrewer)  # For extended color palettes
library(viridis)  
library(scales)  # 确保载入 scales 包


# 生成所有组织按照时间索引的数据
# 提取 Whole_body 作为基准数据
df <- data_long_parsed2 %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, Tissue, Combined_Stage, PSI_value) 


# 构建阶段到时间（小时）的映射（根据前面生成的索引）
day_mapping <- c(
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
df <- df %>%
  mutate(day = as.numeric(recode(as.character(Combined_Stage), !!!day_mapping)))

# ---------------------------
# 转换时间单位及构造自变量：以天为单位，并对天数取对数
# ---------------------------
df <- df %>%
  mutate(day = day / 24)  # 将小时转换为天


# 提取所有组织中DevAS所对应的事件
up_devAS <- merged_DevAS_all %>%
  filter(direction == "up-down")
  
up_devAS_events <- df %>%
  semi_join(up_devAS, 
            by = c("GeneID", "chr", "strand", "exonStart_0base", "exonEnd", 
                   "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")) %>%
  select(GeneID, chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE,
         downstreamES, downstreamEE, Tissue, day, PSI_value)

# Group by Tissue and day, calculating mean PSI_value for each group
tissue_time_psi <- up_devAS_events %>%
  group_by(Tissue, day) %>%
  summarize(mean_PSI = mean(PSI_value, na.rm = TRUE), .groups = "drop")

# Check the summarized data
head(tissue_time_psi)

# Count unique tissues to ensure we have enough colors
n_tissues <- length(unique(tissue_time_psi$Tissue))
cat("Number of unique tissues:", n_tissues, "\n")

# 如果已有颜色映射文件，则读取，否则创建
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

# 获取当前数据中的组织
current_tissues <- unique(tissue_time_psi$Tissue)

# 判断哪些组织在已有映射中还不存在
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

# 颜色候选池
color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB","#4DAF4A","#FF7F00","#FFFF33","#A65628",
  "#F781BF","#999999","#66C2A5","#FC8D62","#A6D854","#FFD92F","#E5C494",
  "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"
)

# 更新颜色映射函数
update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  used_colors <- unname(existing_mapping)
  available_colors <- setdiff(color_pool, used_colors)
  if(length(new_tissues) > length(available_colors)) {
    # 如果颜色不够用，可以循环使用
    available_colors <- rep(available_colors, length.out = length(new_tissues))
  }
  new_mapping <- setNames(available_colors[1:length(new_tissues)], new_tissues)
  updated_mapping <- c(existing_mapping, new_mapping)
  return(updated_mapping)
}

# 更新并保存组织-颜色映射
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

# Group by Tissue and day, calculating mean PSI_value for each group
tissue_time_psi <- up_devAS_events %>%
  group_by(Tissue, day) %>%
  summarize(mean_PSI = mean(PSI_value, na.rm = TRUE), .groups = "drop")

# 筛选出每个组织包含数据点数 >= 3 的数据
tissue_time_psi_filtered <- tissue_time_psi %>%
  group_by(Tissue) %>%
  filter(n() >= 3) %>%
  ungroup()

# 更新组织列表及对应颜色，仅保留过滤后的组织
current_tissues <- unique(tissue_time_psi_filtered$Tissue)
tissue_colors <- tissue_color_mapping[current_tissues]

#组织分开绘制
# 将数据拆分为 Whole_body 和非 Whole_body 两部分
data_whole <- tissue_time_psi_filtered %>% filter(Tissue == "Testis")
data_other <- tissue_time_psi_filtered %>% filter(Tissue != "Testis")

# 获取原始的 Whole_body 颜色
whole_body_color <- tissue_color_mapping["Testis"]

# 创建最终颜色映射
final_colors <- tissue_color_mapping  # 复制原始颜色映射

# 将非 Whole_body 组织颜色设为半透明
for(tissue in names(final_colors)) {
  if(!is.na(tissue) && tissue != "Testis") {
    final_colors[tissue] <- adjustcolor(final_colors[tissue], alpha.f = 0.3)
  }
}

# 创建绘图
p <- ggplot() +
  # 绘制非Whole_body组织点
  geom_point(data = data_other, aes(x = day, y = mean_PSI, color = Tissue),
             size = 2, alpha = 0.3) +
  
  # 绘制非Whole_body组织线条
  geom_smooth(data = data_other, aes(x = day, y = mean_PSI, color = Tissue),
              method = "loess", span = 1, se = FALSE, linewidth = 2) +
  
  # 绘制Whole_body组织点
  geom_point(data = data_whole, aes(x = day, y = mean_PSI, color = Tissue),
             size = 2) +
  
  # 绘制Whole_body组织线条
  geom_smooth(data = data_whole, aes(x = day, y = mean_PSI, color = Tissue),
              method = "loess", span = 0.9, se = FALSE, linewidth = 3) +
  
  # 应用修改后的颜色映射
  scale_color_manual(values = final_colors) +
  
  labs(
    x = "Development",
    y = "PSI",
    color = "Tissue Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "down",  # 修正图例位置为"bottom"
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),  # 增大x轴标签字体大小
    axis.text.y = element_text(size = 12),  # 增大y轴标签字体大小
    axis.title.x = element_text(size = 14),  # 增大x轴标题字体大小
    axis.title.y = element_text(size = 14)   # 增大y轴标题字体大小
  )# 打印图形
print(p)


ggsave("Testis_down_plot.pdf", p, width = 4, height = 4, dpi=300)


down_devAS <- merged_DevAS_table %>%
  filter(direction == "down")

down_up_devAS <- merged_DevAS_table %>%
  filter(direction == "down-up")

up_down_devAS <- merged_DevAS_table %>%
  filter(direction == "up-down")
