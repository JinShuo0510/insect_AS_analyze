library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)
library(smacof)
library(ggrepel)  # 用于智能标签放置
library(stringr)
library(readr)


# 假设 results 是包含各物种数据的列表
results <- results_important
# 对每个物种的数据进行处理、拆分并重命名列（列名格式： species-原列名_序号）
processed_list <- lapply(names(results), function(species) {
  df <- results[[species]]
  
  # 提取包含 Adult, Pupa, Larva, Egg, Cell 的列
  psi_data <- df %>% dplyr::select(matches("Adult|Pupa|Larva|Egg|Cell"))
  
  # 计算每列中以逗号拆分后最多的 token 数，用于后续拆分
  max_values <- sapply(psi_data, function(col) {
    max(sapply(strsplit(col, ","), length))
  })
  
  # 循环每一列，动态拆分并重命名，新列名格式为 "species-原列名_序号"
  for (col_name in names(psi_data)) {
    new_col_names <- paste0(species, "-", col_name, "_", seq_len(max_values[col_name]))
    psi_data <- psi_data %>% separate(
      col = col_name,
      into = new_col_names,
      sep = ",",
      fill = "right"  # 若值不足则填 NA
    )
  }
  
  # 将所有拆分后的数据转换为数值型
  psi_data <- psi_data %>% mutate(across(everything(), ~ parse_number(.)))
  
  # 对每一列，将非 NA 的值提到上方，其余 NA 补在后面（行顺序不再关心）
  psi_data <- as.data.frame(lapply(psi_data, function(col) {
    non_na <- col[!is.na(col)]
    na_part <- rep(NA, length(col) - length(non_na))
    c(non_na, na_part)
  }))
  
  return(psi_data)
})

# 各物种数据的行数可能不同，将每个数据框扩充为相同行数（不足的行补 NA）
max_rows <- max(sapply(processed_list, nrow))
processed_list <- lapply(processed_list, function(df) {
  if(nrow(df) < max_rows) {
    df[(nrow(df)+1):max_rows, ] <- NA
  }
  return(df)
})

# 合并所有物种数据（列拼接，行的顺序不必保持原始对应关系）
merged_data <- do.call(cbind, processed_list)
psi_data <- merged_data %>% mutate_all(as.numeric)

rows_all_na <- apply(psi_data, 1, function(row) all(is.na(row)))
cols_all_na <- apply(psi_data, 2, function(col) all(is.na(col)))
psi_data_clean <- psi_data[!rows_all_na, !cols_all_na]




# 设定NA比例阈值
threshold <- 0.99

# 对每一列计算NA比例，小于阈值的列保留
cols_to_keep <- apply(psi_data_clean, 2, function(col) mean(is.na(col)) < threshold)
psi_data_clean2 <- psi_data_clean[, cols_to_keep]

psi_data_logit <- psi_data_clean %>%
  mutate(across(everything(), ~ {
    x <- .
    # 避免精确 0 或 1
    x <- pmin(pmax(x, 1e-6), 1 - 1e-6)
    log(x / (1 - x))
  }))


#psi_data_clean <- na.omit(psi_data_clean)

 
# cor_mat_spearman <- cor(
#   psi_data_clean, 
#   use = "pairwise.complete.obs", 
#   method = "spearman"
# )

# cor_mat <- cor(psi_data_clean2[, nonzero_sd], use = "pairwise.complete.obs", method = "pearson")

# 计算 Pearson 相关系数矩阵，忽略 NA
cor_mat <- cor(psi_data_logit, use = "pairwise.complete.obs", method = "pearson")

# 清除全 NA 的行和列
rows_all_na <- apply(cor_mat, 1, function(row) all(is.na(row)))
cols_all_na <- apply(cor_mat, 2, function(col) all(is.na(col)))
cor_mat_clean <- cor_mat[!rows_all_na, !cols_all_na]
cor_mat_clean <- na.omit(cor_mat_clean)


#cor_mat_clean <- cor_mat

# 基于 1 - Pearson 相关系数构造距离矩阵，并进行 MDS 降维
dist_mat <- as.dist(1 - cor_mat_clean)
mds_res <- cmdscale(dist_mat, k = 2)

# 构造 MDS 数据框，并将 rownames 存入 label
mds_df <- as.data.frame(mds_res)
mds_df$label <- rownames(mds_res)

# 利用 tidyr::separate 根据分隔符正确分离 label
# 假设 label 格式为： "species.stage.tissue_number"，例如 "Bicyclus_anynana.Egg.Egg_1"
mds_df <- mds_df %>%
  separate(label, into = c("species", "stage", "tissue_num"), sep = "\\.", remove = FALSE) %>%
  separate(tissue_num, into = c("tissue", "number"), sep = "_", convert = TRUE)

# 将生命阶段转换为因子，指定顺序
mds_df$stage <- factor(mds_df$stage, levels = c("Egg", "Cell", "Larva", "Pupa", "Adult"))

# 定义各阶段的点大小映射
size_mapping <- c("Egg" = 2, "Cell" = 2.5, "Larva" = 3, "Pupa" = 4, "Adult" = 5)

# 加载或创建组织颜色映射（tissue_color_mapping），用于 scale_color_manual
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}
current_tissues <- unique(mds_df$tissue)
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))
update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  used_colors <- unname(existing_mapping)
  available_colors <- setdiff(color_pool, used_colors)
  new_mapping <- setNames(available_colors[1:length(new_tissues)], new_tissues)
  updated_mapping <- c(existing_mapping, new_mapping)
  return(updated_mapping)
}
color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628",
  "#F781BF","#999999","#66C2A5","#FC8D62","#A6D854","#FFD92F","#E5C494",
  "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"
)
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

# 定义物种形状映射，根据 mds_df 中提取的 species 信息
unique_species <- sort(unique(mds_df$species))
shapes_vec <- c(1, 16, 2, 17, 0, 15, 18, 5)  # ggplot2 中的 shape 数值
if(length(unique_species) > length(shapes_vec)) {
  stop("需要更多的形状定义!")
}
shape_mapping <- setNames(shapes_vec[seq_along(unique_species)], unique_species)

# 绘图，映射 color 为组织、size 为生命阶段、shape 为物种
p <- ggplot(mds_df, aes(x = V1, y = V2, color = tissue, size = stage, shape = species)) +
  geom_point(alpha = 0.7) +
  scale_size_manual(values = size_mapping) +
  scale_color_manual(values = tissue_color_mapping) +
  scale_shape_manual(values = shape_mapping) +
  # geom_text_repel(aes(label = paste(stage, tissue, sep = "-")),
  #                 size = 3, box.padding = 0.5, point.padding = 0.5,
  #                 segment.color = "grey50", show.legend = FALSE) +
  labs(x = "MDS1", y = "MDS2",
       title = "MDS图：不同生命阶段、组织及物种",
       subtitle = "基于 1 - Pearson 相关系数距离",
       color = "组织", size = "生命阶段", shape = "物种") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text.align = 0
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(color = "black"))
  ) 
#  xlim(-0.4, 0.3) +  # 限制x轴范围
#  ylim(-0.3, 0.3)    # 限制y轴范围
print(p)


ggsave("mds_le_important.pdf", plot = p, width = 12, height = 12)
