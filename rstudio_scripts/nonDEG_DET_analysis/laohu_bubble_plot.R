# 1. 安装并加载包（若已安装可省略）
# install.packages("scatterpie")
library(scatterpie)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 2. 读取并整理数据（同原始逻辑）
df_pie <- tibble::tribble(
  ~comparison,           ~up,   ~down, ~total,
  "larva1_vs_adult",    2232,  2469,  4701,
  "larva2_vs_adult",    2498,  2902,  5400,
  "larva2_vs_larva1",    123,   205,   328,
  "larva3_vs_adult",    2285,  2182,  4467,
  "larva3_vs_larva1",    435,   302,   737,
  "larva3_vs_larva2",    177,   139,   316,
  "larva4_vs_adult",    3146,  3259,  6405,
  "larva4_vs_larva1",    846,   865,  1711,
  "larva4_vs_larva2",    260,   246,   506,
  "larva4_vs_larva3",    108,   122,   230,
  "pupa_vs_adult",      2203,  2487,  4690,
  "pupa_vs_larva1",     2236,  2451,  4687,
  "pupa_vs_larva2",     2551,  2450,  5001,
  "pupa_vs_larva3",     2362,  2699,  5061,
  "pupa_vs_larva4",     3233,  3242,  6475
) %>%
  separate(comparison, into = c("Stage1", "Stage2"), sep = "_vs_") %>%
  mutate(
    # 先按 development 顺序给阶段编码
    idx1 = as.numeric(factor(Stage1, levels = c("larva1","larva2","larva3","larva4","pupa"))),
    idx2 = as.numeric(factor(Stage2, levels = c("adult","larva1","larva2","larva3","larva4"))),
    # 半径映射逻辑保持不变
    radius = rescale(sqrt(total), to = c(0.1, 0.3))
  )

# 3. 绘图：交换 x 和 y
p_pie_ul <- ggplot() +
  # 3.1 用透明点生成大小图例：注意 aes(x=idx2, y=idx1)
  geom_point(
    data = df_pie,
    aes(x = idx2, y = idx1, size = total),
    color = "grey80",
    alpha = 0.3
  ) +
  # 3.2 画饼图气泡：aes(x=idx2, y=idx1, r=radius)
  geom_scatterpie(
    data = df_pie,
    aes(x = idx2, y = idx1, r = radius),
    cols = c("down", "up"),
    alpha = 0.8,
    color = "black",
    lwd = 0.5
  ) +
  # 4. 配色与图例
  scale_fill_manual(
    values = c("down" = "#2A5783", "up" = "#7FB3D5"),
    labels = c("Down", "Up"),
    name = "Direction"
  ) +
  scale_size_area(
    max_size = 8,
    name = "Count"
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
    size = guide_legend(order = 2, override.aes = list(alpha = 0.5, fill = "grey80"))
  ) +
  # 5. 坐标轴：x 对应原来的 Stage2，y 对应原来的 Stage1
  scale_x_continuous(
    breaks = 1:length(unique(df_pie$Stage2)),
    labels = unique(df_pie$Stage2),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = 1:length(unique(df_pie$Stage1)),
    labels = unique(df_pie$Stage1),
    expand = c(0.01, 0)
  ) +
  coord_fixed() +
  labs(x = "Stage B", y = "Stage A") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.box = "vertical"
  )

print(p_pie_ul)
ggsave("pie_plot.pdf",p_pie_ul,width = 6,height=5,dpi=300)