# 安装并加载必要的包
# install.packages("dplyr")
library(dplyr)

# 假设你的数据框叫 df
# df <- results$match_results # 使用你的实际变量名
 df <- analyse_result # 使用你的实际变量名

# --- 1. 定义物种组 (基于系统发育树, 同前) ---
ancestral_ref_species <- c("Apis_mellifera", "Acyrthosiphon_pisum")
descendant_species_pool <- c(
  "Blattella_germanica", "Gryllus_bimaculatus",
  "Helicoverpa_armigera", "Bombyx_mori",
  "Tribolium_castaneum",
  "Drosophila_mojavensis", "Aedes_aegypti"
)
all_defined_species <- c(ancestral_ref_species, descendant_species_pool)

# --- 2. 确定每个 ortholog_group 存在的物种列表 ---

# 获取每个 ortholog_group 实际存在的物种列表
# 注意：这里假设只要 df 中有记录，就代表该 group 在该 species 中存在
# 如果你的 df 可能包含没有实际外显子匹配但仍有记录的行，需要调整此逻辑
ortholog_presence <- df %>%
  # 确保我们只考虑定义内的物种，避免其他物种干扰（如果有的话）
  filter(species %in% all_defined_species) %>%
  # 对每个组，列出它存在于哪些物种中
  group_by(ortholog_group) %>%
  summarise(
    present_species = list(unique(species)),
    num_present_species = n_distinct(species), # 存在的物种数量
    .groups = 'drop'
  )

# --- 3. 筛选 "New Exon Birth" 候选 ---

# 核心逻辑:
# 1. 该 ortholog_group *不存在于任何* ancestral_ref_species 中。
# 2. 该 ortholog_group *至少存在于一个* descendant_species_pool 物种中。

new_exon_candidates <- ortholog_presence %>%
  filter(
    # 条件1: 检查是否 *所有* 祖先参考物种都不在 present_species 列表中
    # 使用 sapply 遍历祖先物种，检查它们是否都不在当前组的 present_species 里
    sapply(present_species, function(p_species) {
      all(!ancestral_ref_species %in% p_species)
    })
    &
      # 条件2: 检查是否 *至少一个* 衍生池物种在 present_species 列表中
      sapply(present_species, function(p_species) {
        any(descendant_species_pool %in% p_species)
      })
  ) %>%
  # 选取候选的 ortholog_group ID
  distinct(ortholog_group)

# --- 4. 输出结果 ---
cat("Species defined as Ancestral Reference:", paste(ancestral_ref_species, collapse=", "), "\n")
cat("Species defined as Descendant Pool:", paste(descendant_species_pool, collapse=", "), "\n\n")

cat("Potential 'New Exon Birth' Candidate Ortholog Groups (", nrow(new_exon_candidates), " found):\n", sep="")
if (nrow(new_exon_candidates) > 0) {
  print(new_exon_candidates)
  
  cat("\nDetails of Candidate Groups (showing species presence):\n")
  candidate_details <- ortholog_presence %>%
    filter(ortholog_group %in% new_exon_candidates$ortholog_group) %>%
    # 展开列表方便查看
    mutate(present_species_str = sapply(present_species, paste, collapse=", ")) %>%
    select(ortholog_group, num_present_species, present_species_str) %>%
    arrange(ortholog_group)
  print(candidate_details, n=Inf) # 打印所有行
  
  # 可选：进一步验证候选者确实只存在于衍生池物种中
  # verification <- candidate_details %>%
  #   left_join(ortholog_presence, by="ortholog_group") %>%
  #   mutate(
  #     only_in_descendants = sapply(present_species, function(p_species) {
  #       all(p_species %in% descendant_species_pool) & !any(p_species %in% ancestral_ref_species)
  #     })
  #   ) %>%
  #   filter(only_in_descendants) # 再次确认筛选逻辑
  # print("Verification check (should match candidate count):")
  # print(nrow(verification))
  
} else {
  cat("No candidate groups found matching the 'New Exon Birth' criteria.\n")
}