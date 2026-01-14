# 载入必要包
library(dplyr)
library(ggplot2)
library(ggrepel)

# 假设 merged_matched 已经在环境中
df <- merged_matched %>%
  mutate(
    negLogFDR = -log10(FDR),
    xval      = logFC,
    highlight = case_when(
      match_type == "mode1" ~ "mode1",
      match_type == "mode2" ~ "mode2",
      TRUE                  ~ "none"
    )
  )

# 选前 10 个高亮点用于标注
top10 <- df %>%
  filter(highlight %in% c("mode1","mode2")) %>%
  arrange(desc(negLogFDR)) %>%
  distinct(id, .keep_all = TRUE) %>%
  slice_head(n = 10)

# 绘制带图例的 Head 火山图
p_head <- ggplot(df, aes(x = xval, y = negLogFDR)) +
  geom_point(aes(fill = highlight),
             shape = 21, color = "white",
             size = 2.5, alpha = 0.8) +
  scale_fill_manual(
    name   = "匹配类型",
    values = c(
      none  = "grey80",
      mode1 = "#2A5783",  # 深蓝
      mode2 = "#7FB3D5"   # 浅蓝
    ),
    labels = c(
      none  = "无匹配",
      mode1 = "模式1",
      mode2 = "模式2"
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.4) +
  geom_vline(xintercept = c(-1, 1),    linetype = "dashed", size = 0.4) +
  geom_text_repel(
    data = top10, aes(label = id),
    size = 3, max.overlaps = 15, segment.size = 0.3
  ) +
  labs(
    title = "C. Head",
    x     = "log2(fold change)",
    y     = "-log10(FDR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(hjust = 0, vjust = 1, face = "bold"),
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "right"           # 显示图例并放到右侧
  )

print(p_head)

ggsave("match_devAS_Testis.pdf",p_head,,
       width=8, height=6, dpi=300)