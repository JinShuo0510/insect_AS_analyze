# 加载必要的包
library(ggplot2)
library(dplyr)
library(binom)       # Wilson 置信区间
library(colorspace)  # lighten / darken

# 1. 计算每个物种的总 exon 数
species_totals <- summary_stats_psi_alternative %>%
  group_by(species) %>%
  summarise(
    species_total = sum(count, na.rm = TRUE),
    .groups      = "drop"
  )

# 2. 提取目标数据并计算比例
plot_data <- summary_stats_psi_alternative %>%
  filter(
    species        %in% c("Bombyx_mori",
                          "Drosophila_mojavensis",
                          "Helicoverpa_armigera",
                          "Tribolium_castaneum"),
    max_psi_tissue %in% c("Head", "Whole_body"),
    exon_type      %in% c("Holometabola-specific", "Species-specific")
  ) %>%
  left_join(species_totals, by = "species") %>%
  mutate(
    proportion = count / species_total
  )

# 3. 计算 Wilson 95% 置信区间
plot_data <- plot_data %>%
  rowwise() %>%
  mutate(
    ci        = list(binom.confint(count, species_total,
                                   conf.level = 0.95,
                                   methods    = "wilson")),
    lower     = ci$lower,
    upper     = ci$upper
  ) %>%
  ungroup() %>%
  select(-ci)

# 4. 计算合并检验的 P 值（使用所有物种的总分母一次性计算）
tot_species_total <- sum(species_totals$species_total)
test_results <- data.frame(
  Tissue     = character(),
  p_value    = numeric(),
  p_formatted= character(),
  stringsAsFactors = FALSE
)
for (tissue in c("Head", "Whole_body")) {
  d <- filter(plot_data, max_psi_tissue == tissue)
  total_h <- sum(d$count[d$exon_type == "Holometabola-specific"])
  total_s <- sum(d$count[d$exon_type == "Species-specific"])
  test    <- prop.test(
    x = c(total_h, total_s),
    n = c(tot_species_total, tot_species_total)
  )
  p       <- test$p.value
  p_fmt   <- if (p < 0.001) {
    f    <- format(p, scientific = TRUE, digits = 1)
    expo <- strsplit(f, "e")[[1]][2]
    paste0("P = ", substr(f, 1, 1), "×10^", expo)
  } else {
    paste0("P = ", round(p, 3))
  }
  test_results <- rbind(
    test_results,
    data.frame(Tissue = tissue,
               p_value = p,
               p_formatted = p_fmt,
               stringsAsFactors = FALSE)
  )
}

# 5. 美化物种名称与因子顺序
plot_data <- plot_data %>%
  mutate(
    Species = species %>%
      gsub("_", " ", .) %>%
      factor(levels = c("Bombyx mori",
                        "Drosophila mojavensis",
                        "Helicoverpa armigera",
                        "Tribolium castaneum"))
  )

# 6. 分组织循环绘图并保存
for (tissue in c("Head", "Whole_body")) {
  d          <- filter(plot_data, max_psi_tissue == tissue)
  base_color <- if (tissue == "Head") "#8DA0CB" else "#A65628"
  color_values <- c(
    "Holometabola-specific" = lighten(base_color, 0.3),
    "Species-specific"      = darken(base_color, 0.2)
  )
  p_text <- test_results$p_formatted[
    test_results$Tissue == tissue
  ]
  
  p <- ggplot(d, aes(x = Species, y = proportion, fill = exon_type)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8),
             width = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.2,
                  position = position_dodge(width = 0.8)) +
    geom_text(aes(label = count),
              position = position_dodge(width = 0.8),
              vjust = 1,    # 标注在柱子基部内侧
              size  = 3) +
    annotate("text", x = Inf, y = 0.6,
             label = p_text,
             hjust = 1.1, vjust = 1,
             size  = 4) +
    scale_y_continuous(limits = c(0, 0.6),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = color_values) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "top",
      legend.title     = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border     = element_rect(fill = NA, color = "black")
    ) +
    labs(
      x     = NULL,
      y     = "Percentage of exon",
      title = tissue
    )
  
  print(p)
  ggsave(
    filename = paste0("psi_alternative_type_", tissue, ".pdf"),
    plot     = p,
    width    = 6, height = 6, dpi = 300
  )
}
