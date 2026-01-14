library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)
library(smacof)
library(ggrepel)  # 用于智能标签放置
library(stringr)  # 添加这一行以使用str_extract和str_replace函数

data_list <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]][["Helicoverpa_armigera"]]
psi_data <- data_list[, grepl("Adult|Pupa|Larva|Egg|Cell", colnames(data_list))]
psi_data <- psi_data[rowMeans(is.na(psi_data)) < 0.2, ]

# 1. 确定每列的最大值数量
max_values <- sapply(psi_data, function(col) {
  max(sapply(strsplit(col, ","), length))
})

# 2. 动态拆分列
for (col_name in names(psi_data)) {
  # 生成动态列名
  new_col_names <- paste0(col_name, "_", seq_len(max_values[col_name]))
  # 拆分列
  psi_data <- psi_data %>%
    separate(
      col = col_name,
      into = new_col_names,
      sep = ",",
      fill = "right"  # 如果值不足，填充 NA
    )
}

# 查看结果
print(psi_data)
psi_data <- psi_data %>% mutate_all(as.numeric)

# 计算 Pearson 相关系数矩阵，忽略 NA
cor_mat <- cor(psi_data,
               use = "pairwise.complete.obs",
               method = "pearson")

# 检查哪些行全是NA
rows_all_na <- apply(cor_mat, 1, function(row) all(is.na(row)))
# 检查哪些列全是NA
cols_all_na <- apply(cor_mat, 2, function(col) all(is.na(col)))

cor_mat_clean <- cor_mat[!rows_all_na, !cols_all_na]
cor_mat_clean <-  na.omit(cor_mat_clean)
dist_mat <- as.dist(1 - cor_mat_clean)
mds_res <- cmdscale(dist_mat, k = 2)

#<!---------------------------plot-----------------------------------!>
# 在获得 MDS 结果后，我们需要处理点的标签
mds_df <- as.data.frame(mds_res)
mds_df$label <- rownames(mds_res)

# 解析标签，包括 "Cell" 阶段
mds_df <- mds_df %>%
  mutate(
    stage = str_extract(label, "^(Adult|Larva|Egg|Pupa|Cell)"),
    tissue = str_replace(str_extract(label, "-.*_"), "^-|_$", ""),  # 移除开头的连字符和结尾的下划线
    tissue = str_replace(tissue, "_$", ""),  # 移除最后一个下划线
    number = as.numeric(str_extract(label, "\\d+$"))
  )

# 创建生命阶段的因子，按指定顺序，包括 "Cell"
mds_df$stage <- factor(mds_df$stage, levels = c("Egg", "Cell", "Larva", "Pupa", "Adult"))

# 创建点大小的映射，包括 "Cell"
size_mapping <- c("Egg" = 2, "Cell" = 2.5, "Larva" = 3, "Pupa" = 4, "Adult" = 5)

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
current_tissues <- unique(mds_df$tissue)

# 找出新的组织
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

# 更新颜色映射
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)

# 保存更新后的颜色映射
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

# 绘图
p <- ggplot(mds_df, aes(x = V1, y = V2, color = tissue, size = stage)) +
  geom_point(alpha = 0.7) +
  scale_size_manual(values = size_mapping) +
  scale_color_manual(values = tissue_color_mapping) +  # 使用更新后的颜色映射
  geom_text_repel(aes(label = paste(stage, tissue, sep = "-")),
                  size = 3, box.padding = 0.5, point.padding = 0.5,
                  segment.color = "grey50", show.legend = FALSE) +
  labs(x = "MDS1", y = "MDS2",
       title = "MDS plot of Bombyx mori life stages and tissues",
       subtitle = "Based on 1 - Pearson's r distance",
       color = "Tissue", size = "Life Stage") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5),  # 减少图例边距
    legend.key.size = unit(0.6, "cm"),  # 减小图例键的大小
    legend.text = element_text(size = 6),  # 减小图例文本大小
    legend.title = element_text(size = 8, face = "bold"),  # 减小图例标题大小
    legend.text.align = 0
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),  # 减小颜色图例的点大小
    size = guide_legend(override.aes = list(color = "black"))
  )
  # xlim(-0.2, 0.2) +  # 限制x轴范围
  # ylim(-0.3, 0.3)    # 限制y轴范围

print(p)
ggsave("Helicoverpa_armigera_plot1-3.pdf", plot = p, width = 8, height = 6, units = "in", dpi = 300)