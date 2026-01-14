# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

# ---------------------- #
# 处理第一个文件 (DevAS 外显子部分)
# ---------------------- #
data1 <- read.table("phastcons_bp_results_400bp_0415.txt", header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("N/A"))

data1 <- data1 %>%
  mutate(
    GeneID = sapply(strsplit(File_Name, "__"), `[`, 1),
    exonstart = as.numeric(sapply(strsplit(File_Name, "__"), `[`, 2)),
    exonend   = as.numeric(sapply(strsplit(File_Name, "__"), `[`, 3)),
    exon_length = exonend - exonstart,
    exon_type = ifelse(exon_length <= 28, "microexon", "macroexon"),
    dataset = "DevAS"  # 添加数据来源标识
  )

data_long1 <- data1 %>%
  pivot_longer(
    cols = starts_with("BP_"),
    names_to = "bp",
    values_to = "score",
    names_prefix = "BP_"
  ) %>%
  mutate(
    bp = as.numeric(bp),
    score = as.numeric(score)
  )

summary_data1 <- data_long1 %>%
  group_by(dataset, exon_type, bp) %>%
  summarise(
    n = sum(!is.na(score)),
    mean_score = mean(score, na.rm = TRUE),
    sd = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n > 0) %>%
  mutate(
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    lower = mean_score - ci,
    upper = mean_score + ci
  )

# ---------------------- #
# 处理第二个文件 (非 DevAS 外显子部分)
# ---------------------- #
data2 <- read.table("phastcons_bp_results_400bp_none.txt", header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("N/A"))

data2 <- data2 %>%
  mutate(
    GeneID = sapply(strsplit(File_Name, "__"), `[`, 1),
    exonstart = as.numeric(sapply(strsplit(File_Name, "__"), `[`, 2)),
    exonend   = as.numeric(sapply(strsplit(File_Name, "__"), `[`, 3)),
    exon_length = exonend - exonstart,
    exon_type = ifelse(exon_length <= 28, "microexon", "macroexon"),
    dataset = "nonDevAS"  # 添加数据来源标识
  )

data_long2 <- data2 %>%
  pivot_longer(
    cols = starts_with("BP_"),
    names_to = "bp",
    values_to = "score",
    names_prefix = "BP_"
  ) %>%
  mutate(
    bp = as.numeric(bp),
    score = as.numeric(score)
  )

summary_data2 <- data_long2 %>%
  group_by(dataset, exon_type, bp) %>%
  summarise(
    n = sum(!is.na(score)),
    mean_score = mean(score, na.rm = TRUE),
    sd = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n > 0) %>%
  mutate(
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    lower = mean_score - ci,
    upper = mean_score + ci
  )

# ---------------------- #
# 合并两份数据并添加 x 轴断层（代表 exon）
# ---------------------- #
gap <- 50  # 定义 gap 宽度

process_gap <- function(df) {
  df %>% 
    mutate(
      direction = if_else(bp < 0, "upstream", "downstream"),
      bp_plot = if_else(bp < 0, bp, bp + gap)
    )
}

summary_data1 <- process_gap(summary_data1)
summary_data2 <- process_gap(summary_data2)

summary_all <- bind_rows(summary_data1, summary_data2)

# ---------------------- #
# 绘图
# ---------------------- #
p <- ggplot(summary_all, 
            aes(x = bp_plot, y = mean_score, 
                group = interaction(dataset, exon_type, direction),
                color = dataset, linetype = exon_type)) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, size = 1) +
  scale_x_continuous(breaks = c(-200, -100, 0, gap + 1, gap + 100, gap + 200),
                     labels = c("-200", "-100", "0", "0", "100", "200")) +
  labs(
    x = "Postion",
    y = "Average Score",
    color = "Data Source",
    linetype = "Exon Type",
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(0.1, 0.7),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

print(p)
ggsave("Fig5E.pdf", plot = p, width = 9, height = 5, dpi = 300)
