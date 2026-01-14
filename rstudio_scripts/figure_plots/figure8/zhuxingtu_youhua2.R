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
  "NM_001279379.1",
  "XM_062672840.1"
)
tpm_long <- tpm_df %>%
  filter(transcript_id %in% target_tx) %>%
  pivot_longer(
    cols      = -transcript_id,
    names_to  = "Run",
    values_to = "TPM"
  ) %>%
  mutate(TPM = as.numeric(TPM))

# 4. 构建所有 transcript × Run 的完整组合，缺失 TPM 补 0，并合并 metadata
all_combos <- expand.grid(
  transcript_id = target_tx,
  Run           = unique(meta_bm$Run),
  stringsAsFactors = FALSE
)
expr_full <- all_combos %>%
  left_join(tpm_long, by = c("transcript_id", "Run")) %>%
  mutate(TPM = ifelse(is.na(TPM), 0, TPM)) %>%
  left_join(meta_bm, by = "Run") %>%
  mutate(TPM_log10 = log10(TPM + 1))

# 5. 定义配色
academic_colors <- c("#E41A1C", "#4DAF4A", "#377EB8", "#984EA3")
fill_colors <- rep(academic_colors,
                   length.out = length(unique(expr_full$transcript_id)))

# 6. 将 Age 和 Tissue 明确设为因子（按出现顺序）
expr_full$Age    <- factor(expr_full$Age, levels = unique(expr_full$Age))
expr_full$Tissue <- factor(expr_full$Tissue, levels = unique(expr_full$Tissue))

# 7. 按发育时期 (Age) 绘制不重叠的柱状图并保存
p_age <- ggplot(expr_full, aes(x = Age, y = TPM_log10, fill = transcript_id)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width    = 0.8
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "isoform TPM in every tissue",
    x     = "Age",
    y     = "log10(TPM + 1)",
    fill  = "Transcript ID"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
ggsave(
  filename = "CYP367A1_log10TPM_by_Age_academic_colors_nododge_overlap.pdf",
  plot     = p_age,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# 8. 指定要筛选的组织列表
tissues_of_interest <- c("Brain", "Fatbody", "Midgut")

# 9. 对 expr_full 进行组织筛选，生成 expr_filtered
expr_filtered <- expr_full %>%
  filter(Tissue %in% tissues_of_interest) %>%
  mutate(Tissue = factor(Tissue, levels = tissues_of_interest))

# 10. 按组织 (Tissue) 绘制不重叠的柱状图并保存
p_tissue <- ggplot(expr_filtered, aes(x = Tissue, y = TPM_log10, fill = transcript_id)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width    = 0.8
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "isoform TPM in selected tissues",
    x     = "Tissue",
    y     = "log10(TPM + 1)",
    fill  = "Transcript ID"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )
ggsave(
  filename = "CYP367A1_log10TPM_by_SelectedTissues.pdf",
  plot     = p_tissue,
  width    = 8,
  height   = 5,
  dpi      = 300
)
