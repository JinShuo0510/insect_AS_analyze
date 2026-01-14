# 安装并加载必要的包
# install.packages(c("dplyr", "tidyr"))
library(dplyr)
library(tidyr)

# 假设你的原始数据在 results$match_results
df <- results$match_results

# 1. 数据预处理与长表转换
psi_cols <- grep("_PSI$", colnames(df), value = TRUE)
df_long <- df %>%
  select(ortholog_group, species, species_count, all_of(psi_cols)) %>%
  pivot_longer(
    cols = all_of(psi_cols),
    names_to = "condition",
    values_to = "psi_value",
    values_drop_na = TRUE
  )

# 2. 计算物种内 PSI 统计
MIN_VALID_OBS <- 2
psi_summary <- df_long %>%
  group_by(ortholog_group, species, species_count) %>%
  summarise(
    min_psi       = if (n() >= MIN_VALID_OBS) min(psi_value, na.rm = TRUE) else NA_real_,
    max_psi       = if (n() >= MIN_VALID_OBS) max(psi_value, na.rm = TRUE) else NA_real_,
    mean_psi      = if (n() >= MIN_VALID_OBS) mean(psi_value, na.rm = TRUE) else NA_real_,
    sd_psi        = if (n() >= MIN_VALID_OBS) sd(psi_value, na.rm = TRUE) else NA_real_,
    num_valid_psi = n(),
    .groups       = "drop"
  ) %>%
  filter(num_valid_psi >= MIN_VALID_OBS) %>%
  mutate(
    min_psi = ifelse(is.infinite(min_psi), NA_real_, min_psi),
    max_psi = ifelse(is.infinite(max_psi), NA_real_, max_psi)
  )

# 3. 设置阈值
CONSTITUTIVE_THRESHOLD    <- 0.90   # 组成型阈值
ALTERNATIVE_MAX_THRESHOLD <- 0     # 噪音剔除阈值（可设为0以只要有变动即纳入）
AMPLITUDE_THRESHOLD       <- 0.2   # 最小 dPSI 阈值

# 4. 分类剪接模式
classified_patterns <- psi_summary %>%
  mutate(
    splicing_pattern = case_when(
      min_psi >= CONSTITUTIVE_THRESHOLD                              ~ "Constitutive",
      (max_psi - min_psi) >= AMPLITUDE_THRESHOLD                     ~ "Alternative",
      TRUE                                                           ~ "Unclear"
    )
  )

# 5. 查看 Unclear 类别
unclear_events <- classified_patterns %>% filter(splicing_pattern == "Unclear")
cat("Unclear 事件数：", nrow(unclear_events), "\n")

# 6. 输出 Alternative 事件
alternative_events <- classified_patterns %>% filter(splicing_pattern == "Alternative")
cat("Alternative 事件数：", nrow(alternative_events), "\n")
print(alternative_events, n = Inf)

# 7. 跨物种模式比较
ancestral_ref_species <- c("Apis_mellifera", "Acyrthosiphon_pisum")
descendant_species_pool <- c(
  "Blattella_germanica", "Gryllus_bimaculatus",
  "Helicoverpa_armigera", "Bombyx_mori",
  "Tribolium_castaneum", "Drosophila_mojavensis",
  "Aedes_aegypti"
)

alternification_candidates <- classified_patterns %>%
  group_by(ortholog_group) %>%
  filter(
    sum(species %in% ancestral_ref_species) > 0,
    sum(species %in% descendant_species_pool) > 0,
    n_distinct(species) > 1
  ) %>%
  filter(
    all(splicing_pattern[species %in% ancestral_ref_species] == "Constitutive") &
      any(splicing_pattern[species %in% descendant_species_pool] == "Alternative")
  ) %>%
  ungroup() %>%
  distinct(ortholog_group)

cat("Found", nrow(alternification_candidates), "alternification candidates\n")
print(alternification_candidates, n = Inf)

# 8. 查看候选组详细信息
candidate_details_alternative <- classified_patterns %>%
  filter(ortholog_group %in% alternification_candidates$ortholog_group) %>%
  mutate(
    group_type = case_when(
      species %in% ancestral_ref_species    ~ "Ancestral_Ref",
      species %in% descendant_species_pool  ~ "Descendant_Pool",
      TRUE                                  ~ "Other"
    )
  ) %>%
  arrange(ortholog_group, group_type, species)

print(candidate_details_alternative, n = Inf)
