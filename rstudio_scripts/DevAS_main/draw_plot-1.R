library(tidyr)
library(dplyr)
library(ggplot2)
library(colorRamp2)  # 用于颜色插值

# 修改缩写函数，只使用两个字段的首字母
create_abbreviation <- function(name) {
  parts <- strsplit(name, " ")[[1]]
  if (length(parts) >= 2) {
    paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))
  } else {
    substr(name, 1, 2)  # 如果只有一个单词，则取前两个字母
  }
}

# 假设你的数据框名为 all_events_and_devAS_count
df_long <- all_events_and_devAS_count %>%
  select(
    Species, Tissue, Age,
    MXE_count, SE_count, A3SS_count, A5SS_count, RI_count
  ) %>%
  # 重分类 Age
  mutate(Age = case_when(
    Age %in% c("Egg", "Larva", "Nymph", "Pupa", "Adult") ~ Age,
    TRUE ~ "Other"
  )) %>%
  # 简化物种名称
  mutate(Species = sapply(Species, create_abbreviation)) %>%
  # pivot_longer 将 5 种事件计数列转换为长表
  pivot_longer(
    cols = ends_with("_count"),
    names_to = "Event_type",
    values_to = "Count",
    names_pattern = "(.*)_count" # 去掉"_count"保留前缀做事件类型
  )

# 设置因子顺序
df_long$Age <- factor(df_long$Age, levels = c("Egg", "Larva", "Nymph", "Pupa", "Adult", "Other"))
df_long$Species <- factor(df_long$Species, levels = sort(unique(df_long$Species)))

# 计算每个 Age 包含的物种数目
age_species_count <- df_long %>%
  group_by(Age) %>%
  summarise(Species_count = n_distinct(Species)) %>%
  ungroup()

# 创建一个基础的自定义颜色集
color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB", "#4DAF4A",  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62",  "#A6D854", "#FFD92F", "#E5C494",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)


# 函数：更新颜色映射
update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  # 获取已使用的颜色
  used_colors <- unname(existing_mapping)
  
  # 找出未使用的颜色
  available_colors <- setdiff(color_pool, used_colors)
  
  # 为新组织分配颜色
  new_mapping <- setNames(
    available_colors[1:length(new_tissues)],
    new_tissues
  )
  
  # 合并现有映射和新映射
  updated_mapping <- c(existing_mapping, new_mapping)
  
  return(updated_mapping)
}

# 尝试加载现有的颜色映射
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

# 获取当前数据中的所有唯一组织
current_tissues <- unique(df_long$Tissue)

# 找出新的组织
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

# 更新颜色映射
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)

# 保存更新后的颜色映射
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

# 绘图
p <- ggplot(df_long, aes(x = Species, y = Count,
                         color = Tissue, shape = Event_type)) +
  geom_point(
    position = position_dodge(width = 0.7),
    size = 2
  ) +
  facet_grid(~ Age, scales = "free_x", space = "free_x") +
  scale_shape_manual(values = c(
    "MXE" = 21,
    "SE" = 19,
    "A3SS" = 24,
    "A5SS" = 25,
    "RI"  = 23
  )) +
  scale_color_manual(values = tissue_color_mapping) +
  scale_y_continuous(trans = "log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    strip.background = element_blank()
  ) +
  labs(
    x = "Species",
    y = "Number of Selected Events",
    color = "Tissue",
    shape = "Event Type"
  )

print(p)

ggsave("plotA-1.pdf", plot = p, width = 15, height = 6, units = "in", dpi = 300)

