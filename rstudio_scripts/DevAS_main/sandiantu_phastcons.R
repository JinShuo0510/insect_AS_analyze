# 载入必要的包
# install.packages("dplyr")    # 如果尚未安装，请取消注释此行
# install.packages("ggplot2")  # 如果尚未安装，请取消注释此行
library(dplyr)
library(ggplot2)

##############################################################################
# 1. 读取数据（假设文件名为 data.txt，包含三列：Tissue, File, Average_PhastCons_Score）
##############################################################################
df_raw <- read.table("phastcons_results.txt", header = TRUE, stringsAsFactors = FALSE)


# 将非数值字符串替换为 NA，然后强制转换为 numeric
df_raw$Average_PhastCons_Score <- gsub("^NA$|^N/A$", "", df_raw$Average_PhastCons_Score) 
df_raw$Average_PhastCons_Score <- as.numeric(df_raw$Average_PhastCons_Score)



# 只关注第一列 (Tissue) 和第三列 (Average_PhastCons_Score)
# df_raw 的列名请确保与实际文件一致，若不一致，请自行修改
df_summary <- df_raw %>%
  group_by(Tissue) %>%
  summarise(
    mean_score = mean(Average_PhastCons_Score, na.rm = TRUE),
    count      = n(),
    sd_score   = sd(Average_PhastCons_Score, na.rm = TRUE)
  ) %>%
  # 计算标准误并得到上下界：均值 ± 2×SE
  mutate(
    se_score = sd_score / sqrt(count),
    lower    = mean_score - 2 * se_score,
    upper    = mean_score + 2 * se_score
  ) %>%
  ungroup()

##############################################################################
# 2. 加载或更新组织颜色映射（tissue_color_mapping）
##############################################################################
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

current_tissues <- unique(df_summary$Tissue)
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  used_colors <- unname(existing_mapping)
  available_colors <- setdiff(color_pool, used_colors)
  new_mapping <- setNames(available_colors[1:length(new_tissues)], new_tissues)
  updated_mapping <- c(existing_mapping, new_mapping)
  return(updated_mapping)
}

color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB","#4DAF4A","#FF7F00","#FFFF33","#A65628",
  "#F781BF","#999999","#66C2A5","#FC8D62","#A6D854","#FFD92F","#E5C494",
  "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"
)

tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

##############################################################################
# 3. 绘制散点图
##############################################################################
p <- ggplot(df_summary, aes(x = Tissue, y = mean_score, color = Tissue)) +
  geom_point(size = 3) +
  # 误差棒，基于均值 ± 2×标准误
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.6) +
  # 使用手动映射，确保同一组织能获得一致的颜色
  scale_color_manual(values = tissue_color_mapping) +
  labs(
    x = "Tissue",
    y = "Mean PhastCons Score (± 2×SE)",
    title = "Mean PhastCons Score with 95% CI (±2×SE)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 将图形输出到窗口或保存到文件
print(p)


# 如果想将图像保存为 PDF 文件，可使用：
ggsave("PhastCons_Score_Plot.pdf", p, width = 6, height = 6)
