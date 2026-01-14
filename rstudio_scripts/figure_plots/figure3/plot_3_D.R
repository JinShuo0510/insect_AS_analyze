# 加载必要包
library(ggplot2)
library(ggpmisc)
library(dplyr)

#===================== 1. 读取数据 =====================#
# 假设汇总后的数据保存在 "summary_result.csv" 中，
# 数据中应包含：tissue、gene_count、total_dPSI_count 三个关键列
df <- read.csv("summary_result.csv", stringsAsFactors = FALSE)

#===================== 2. 读取或创建颜色映射 =====================#
# 检查/创建颜色映射文件
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

# 获取当前组织列表（排除 NA）
current_tissues <- na.omit(unique(df$tissue))

# 定义颜色候选池（保持原顺序）
color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB","#4DAF4A","#FF7F00","#FFFF33","#A65628",
  "#F781BF","#999999","#66C2A5","#FC8D62","#A6D854","#FFD92F","#E5C494",
  "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"
)

# 更新颜色映射函数（保持颜色池顺序）
update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  # 已用颜色（按颜色池顺序排列）
  used_colors <- intersect(color_pool, unname(existing_mapping))
  # 可用颜色（保持原顺序）
  available_colors <- setdiff(color_pool, used_colors)
  # 为新组织分配颜色（按颜色池顺序）
  new_mapping <- setNames(
    available_colors[seq_len(length(new_tissues))],
    new_tissues
  )
  # 合并时保留原有顺序
  updated_mapping <- c(existing_mapping, new_mapping)
  updated_mapping <- updated_mapping[order(match(names(updated_mapping), c(names(existing_mapping), new_tissues)))]
  return(updated_mapping)
}

# 更新颜色映射（将新组织加入映射中）
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

#===================== 3. 按组织生成图形并保存 =====================#
unique_tissues <- unique(df$tissue)
for (t in unique_tissues) {
  # 过滤当前组织的数据
  tissue_data <- df %>% filter(tissue == t)
  
  # 若无数据或未分配颜色，则跳过
  if (!t %in% names(tissue_color_mapping) || nrow(tissue_data) == 0) {
    next
  }
  
  # 获取当前组织的颜色
  current_color <- tissue_color_mapping[t]
  
  # 构建图形，采用提供的基本样式（除颜色外）
  p <- ggplot(tissue_data, aes(x = gene_count, y = total_dPSI_count)) +
    geom_point(color = current_color, size = 5) +
    geom_smooth(method = "lm", se = TRUE, color = current_color, size = 3) +
    stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
                 formula = y ~ x,
                 parse = TRUE,
                 size = 6,
                 label.x = "left",
                 label.y = 0.9) +
    labs(
      x = "Gene Count",
      y = "Total dPSI Count"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    )
  
  # 保存图形为 PDF 文件（每个组织一个图形）
  pdf_filename <- paste0("ScatterPlot_", t, ".pdf")
  pdf(pdf_filename, width = 4, height = 4)
  print(p)
  dev.off()
}
