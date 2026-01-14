# 安装并加载所需包
# install.packages(c("ggplot2", "binom", "colorspace"))
library(ggplot2)
library(binom)     # 用于计算二项分布置信区间
library(colorspace) # 用于颜色操作（调整亮度）

# 1. 准备数据 ----------------------------------------------------------
my_data <- data.frame(
  Tissue = c("Developmental_tissue", "Developmental_tissue",
             "Fat_body", "Fat_body",
             "Gland", "Gland",
             "Head", "Head",
             "Ovary", "Ovary",
             "Testis", "Testis",
             "Thorax", "Thorax",
             "Whole_body", "Whole_body"),
  Status = c("Early", "Late",
             "Early", "Late",
             "Early", "Late",
             "Early", "Late",
             "Early", "Late",
             "Early", "Late",
             "Early", "Late",
             "Early", "Late"),
  TotalEvents = c(122132, 49733,
                  119694, 102655,
                  543374, 73974,
                  148596, 136865,
                  50486, 79269,
                  204838, 130764,
                  23542, 40020,
                  194889, 271091),
  EventsWithHighDPSI = c(15651, 5214,
                         20192, 18302,
                         191422, 15717,
                         22792, 18745,
                         5910, 7460,
                         41371, 24812,
                         1737, 3350,
                         37048, 53534)
)

# 计算比例
my_data$Proportion <- my_data$EventsWithHighDPSI / my_data$TotalEvents

# 2. 计算95%置信区间 ---------------------------------------------------
for(i in 1:nrow(my_data)) {
  ci <- binom.confint(my_data$EventsWithHighDPSI[i], 
                      my_data$TotalEvents[i], 
                      conf.level = 0.95, 
                      methods = "wilson")  # Wilson方法适合大样本量
  
  my_data$lower[i] <- ci$lower
  my_data$upper[i] <- ci$upper
}

# 3. 设置颜色映射 ------------------------------------------------------
# 如果已有颜色映射文件，则读取，否则创建
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

# 获取当前数据中的组织
current_tissues <- unique(my_data$Tissue)

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

# 4. 创建与状态相关的颜色映射 ------------------------------------------
# 创建组合键：组织_状态
my_data$TissueStatus <- paste(my_data$Tissue, my_data$Status, sep = "_")

# 为每个组织创建轻（Early）深（Late）两种颜色
tissue_status_colors <- c()
for (tissue in current_tissues) {
  base_color <- tissue_color_mapping[tissue]
  # 使用colorspace包调整亮度：Early使用较亮版本，Late使用较暗版本
  light_color <- lighten(base_color, 0.3) # Early版本 - 较亮
  dark_color <- darken(base_color, 0.2)   # Late版本 - 较暗
  
  tissue_status_colors[paste(tissue, "Early", sep = "_")] <- light_color
  tissue_status_colors[paste(tissue, "Late", sep = "_")] <- dark_color
}

# 5. 绘制带有置信区间和自定义颜色的柱形图 ------------------------------
# 创建一个因子来控制图例顺序
my_data$TissueStatus <- factor(
  my_data$TissueStatus, 
  levels = paste(
    rep(current_tissues, each = 2), 
    c("Early", "Late"), 
    sep = "_"
  )
)

# 绘制图表
p <- ggplot(my_data, aes(x = Tissue, y = Proportion, fill = TissueStatus)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  # 添加误差线显示置信区间
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, 
                position = position_dodge(width = 0.9)) +
  # 使用自定义颜色
  scale_fill_manual(values = tissue_status_colors,
                    labels = function(x) sub(".*_", "", x)) + # 仅显示状态部分
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = "Tissue",
    y = "Proportion of Events with dPSI > 0.5",
    fill = "Status",
    title = "Proportion of High dPSI Events with 95% Confidence Intervals"
  ) +
  guides(fill = guide_legend(nrow = 1))

# 显示图形
print(p)

ggsave("3-B.pdf", plot = p, width = 8, height = 6, dpi = 300)
