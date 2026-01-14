# 加载必要的R包
# install.packages("dplyr")   # 如果尚未安装，请取消注释此行
# install.packages("ggplot2") # 如果尚未安装，请取消注释此行
library(dplyr)
library(ggplot2)

#####################
# 第一步：读取与整理数据
#####################

# 假设你的结果文件为 "my_data.txt"
# 样例文件列示意：1=tissue, 2=chr, 3=start, 4=end, 5=gene, 6=divisible
df_raw <- read.table("merged_IDR_result.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# 为便于后续处理，给列命名
colnames(df_raw) <- c("tissue","gene","divisible")

# 仅保留关心的列：组织(tissue)和能否被3整除(divisible)，divisible=1或0
df <- df_raw %>%
  select(tissue, divisible)

# 转换 divisible 列为数值型（若已是数值可省略）
df$divisible <- as.numeric(df$divisible)

#####################
# 第二步：分组统计，计算平均值和二项分布95%置信区间
#####################

# 利用 prop.test 或 binom.test 计算二项分布置信区间，这里用 prop.test 做示例
# 对于每个组织，x 为 1 的个数，n 为总行数
df_summary <- df %>%
  group_by(tissue) %>%
  summarize(
    n = n(),
    x = sum(divisible),
    p = x / n,  # 1 的比例，作为平均值
    # 使用 prop.test 计算95%置信区间
    conf.low = prop.test(x, n)$conf.int[1],
    conf.high = prop.test(x, n)$conf.int[2]
  )

#####################
# 第三步：组织配色映射
#####################

# 如果已有 tissue_color_mapping.rds 文件，则读取，否则创建空映射
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

# 从 df_summary 中获取当前出现的组织
current_tissues <- unique(df_summary$tissue)
# 判断哪些组织在已有映射中还不存在
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

# 定义一个更新映射的函数
update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  used_colors <- unname(existing_mapping)
  available_colors <- setdiff(color_pool, used_colors)
  new_mapping <- setNames(available_colors[1:length(new_tissues)], new_tissues)
  updated_mapping <- c(existing_mapping, new_mapping)
  return(updated_mapping)
}

# 颜色候选池
color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB","#4DAF4A","#FF7F00","#FFFF33","#A65628",
  "#F781BF","#999999","#66C2A5","#FC8D62","#A6D854","#FFD92F","#E5C494",
  "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
  "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"
)

# 更新并保存组织-颜色映射
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

#####################
# 第四步：绘图
#####################

p <- ggplot(df_summary, aes(x = tissue, y = p, color = tissue)) +
  geom_point(size = 3) +
  # 误差棒，基于二项分布得到的conf.low和conf.high
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 0.6) +
  # 使用手动映射，确保同一组织能获得一致的颜色
  scale_color_manual(values = tissue_color_mapping) +
  labs(
    x = "Tissue",
    y = "Proportion (divisible by 3)",
    title = "Proportion of Divisible-by-3 Regions with 95% CI"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 将图形输出到窗口或保存到文件
print(p)
ggsave("IDR_output_plot.pdf", plot = p, width = 6, height = 6)
