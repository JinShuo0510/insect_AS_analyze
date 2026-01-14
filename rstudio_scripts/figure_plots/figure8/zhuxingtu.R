# 加载必要包
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 读取并清洗 metadata
meta <- read.csv(
  "2024_04_01_DATA_withGroupinfo_with_treatment_and_description_sort_final.csv",
  header      = TRUE,
  stringsAsFactors = FALSE,
  quote       = "\"",
  check.names = FALSE
)
meta_bm <- meta %>%
  filter(ScientificName == "Bombyx_mori") %>%
  mutate(
    Age    = gsub(" ", "_", trimws(Age)),
    Tissue = gsub(" ", "_", trimws(Tissue.1))
  ) %>%
  select(Run, Age, Tissue)

# 2. 读取 TPM 数据：先读表头，再跳过表头读数值
header_line <- readLines("Bombyx_mori_tpm_values.txt", n = 1)
runs <- strsplit(header_line, "\t")[[1]]
tpm_df <- read.table(
  "Bombyx_mori_tpm_values.txt",
  sep    = "\t",
  header = FALSE,
  skip   = 1,
  fill   = TRUE,
  stringsAsFactors = FALSE
)
colnames(tpm_df) <- c("transcript_id", runs)

# 3. 筛选目标转录本并转长格式
target_tx <- c(
  "NM_001173244.1",
  "XM_062674822.1"
)
tpm_long <- tpm_df %>%
  filter(transcript_id %in% target_tx) %>%
  pivot_longer(
    cols      = -transcript_id,
    names_to  = "Run",
    values_to = "TPM"
  ) %>%
  mutate(TPM = as.numeric(TPM))

# 4. 构建所有 transcript × Run 的完整组合，缺失 TPM 补 0
all_combos <- expand.grid(
  transcript_id = target_tx,
  Run           = unique(meta_bm$Run),
  stringsAsFactors = FALSE
)
expr_full <- all_combos %>%
  left_join(tpm_long, by = c("transcript_id", "Run")) %>%
  mutate(TPM = ifelse(is.na(TPM), 0, TPM)) %>%
  left_join(meta_bm, by = "Run")

# 5. 分别绘图并保存

# 1. 定义配色
academic_colors <- c("#E41A1C",  "#4DAF4A", "#377EB8", "#984EA3")
fill_colors <- rep(academic_colors,
                   length.out = length(unique(expr_full$transcript_id)))

# 2. 按发育时期 (Age) 绘制并保存
p_age <- ggplot(expr_full, aes(x = Age, y = TPM, fill = transcript_id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "各转录本在不同发育时期的 TPM",
    x     = "发育时期 (Age)",
    y     = "TPM",
    fill  = "转录本 ID"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
ggsave("CPR37_TPM_by_Age_academic_colors.pdf", p_age, width = 8, height = 5, dpi = 300)

# 3. 按组织 (Tissue) 绘制并保存
p_tissue <- ggplot(expr_full, aes(x = Tissue, y = TPM, fill = transcript_id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "各转录本在不同组织的 TPM",
    x     = "组织 (Tissue)",
    y     = "TPM",
    fill  = "转录本 ID"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
ggsave("CPR37_TPM_by_Tissue_academic_colors.pdf", p_tissue, width = 8, height = 5, dpi = 300)
