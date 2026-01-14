# 加载需要的包
library(ggplot2)
library(tidyr)

# 构造示例数据（多个TissueTau点）
df <- data.frame(
  TissueTau = c("0.05", "0.15", "0.25", "0.45", 
                "0.55", "0.75", "0.85", "0.95"),
  DevAS     = c( 2200, 1400, 2000, 2700, 1800, 1600, 1200, 1000),
  AS        = c( 1300, 800,1200, 1000, 1500, 900,    600,  500),
  NonAS     = c(   550,  350,500,  450, 600, 400,    300,  250)
)

# 将宽表转换为长表
df_long <- pivot_longer(
  df, 
  cols = c("DevAS", "AS", "NonAS"), 
  names_to = "Category", 
  values_to = "NumberOfGenes"
)

# 强制设定 TissueTau 的顺序（从小到大）
df_long$TissueTau <- factor(df_long$TissueTau, 
                            levels = c("0.05", "0.15", "0.25", "0.35", "0.45", 
                                       "0.55", "0.65", "0.75", "0.85", "0.95"))

# 控制堆叠顺序：NonAS (底部) -> AS (中部) -> DevAS (顶部)
df_long$Category <- factor(df_long$Category, 
                           levels = c("DevAS","AS", "NonAS" ))

# 自定义配色
color_values <- c(
  "NonAS" = "#D9D9D9", # 浅蓝
  "AS"    = "#80B1D3", # 深蓝
  "DevAS" =  "#457B9D" # 粉/橙
)

# 绘图
p <- ggplot(df_long, aes(x = TissueTau, y = NumberOfGenes, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_values) +
  labs(x = "TissueTau", y = "Number of Genes") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 12),
    axis.title         = element_text(size = 14),
    legend.title       = element_blank()
  )
print(p)

ggsave("zhuxingtu.pdf", p, width = 8, height = 6, dpi = 300)
