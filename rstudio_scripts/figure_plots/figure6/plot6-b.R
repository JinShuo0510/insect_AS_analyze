# 加载必要的库
library(ggplot2)
library(dplyr)
library(colorspace)
library(binom)       # 二项式置信区间
library(stats)       # 包含prop.test函数

# 1. 读取并整理数据 ----------------------------------------------------
gene_data <- read.csv("csv/multi_species_stage_genes.csv", stringsAsFactors = FALSE)

#筛选新外显子的基因
filtered_data <- gene_data %>%
  semi_join(result_table_detail, by = c("species" = "species", "gene_id" = "gene_id"))
gene_data <- filtered_data

# 计算每个物种在每个组织和时期的基因计数
summary_data <- gene_data %>%
  filter(species %in% c("Acyrthosiphon_pisum", "Apis_mellifera", 
                        "Bombyx_mori", "Drosophila_mojavensis", 
                        "Helicoverpa_armigera", "Tribolium_castaneum")) %>%
  group_by(species) %>%
  summarize(
    Head_early_count = sum(Head_early, na.rm = TRUE),
    Head_late_count = sum(Head_late, na.rm = TRUE),
    Head_total = Head_early_count + Head_late_count,
    Whole_body_early_count = sum(Whole_body_early, na.rm = TRUE),
    Whole_body_late_count = sum(Whole_body_late, na.rm = TRUE),
    Whole_body_total = Whole_body_early_count + Whole_body_late_count
  )

# 创建长格式数据
plot_data <- data.frame()

for(i in 1:nrow(summary_data)) {
  species_name <- summary_data$species[i]
  
  # Head tissue
  if(summary_data$Head_total[i] > 0) {
    plot_data <- rbind(plot_data, data.frame(
      Species = species_name,
      Tissue = "Head",
      Stage = "early",  # 使用小写
      Count = summary_data$Head_early_count[i],
      Total = summary_data$Head_total[i],
      Proportion = summary_data$Head_early_count[i] / summary_data$Head_total[i]
    ))
    plot_data <- rbind(plot_data, data.frame(
      Species = species_name,
      Tissue = "Head",
      Stage = "late",   # 使用小写
      Count = summary_data$Head_late_count[i],
      Total = summary_data$Head_total[i],
      Proportion = summary_data$Head_late_count[i] / summary_data$Head_total[i]
    ))
  }
  
  # Whole body tissue
  if(summary_data$Whole_body_total[i] > 0) {
    plot_data <- rbind(plot_data, data.frame(
      Species = species_name,
      Tissue = "Whole_body",
      Stage = "early",  # 使用小写
      Count = summary_data$Whole_body_early_count[i],
      Total = summary_data$Whole_body_total[i],
      Proportion = summary_data$Whole_body_early_count[i] / summary_data$Whole_body_total[i]
    ))
    plot_data <- rbind(plot_data, data.frame(
      Species = species_name,
      Tissue = "Whole_body",
      Stage = "late",   # 使用小写
      Count = summary_data$Whole_body_late_count[i],
      Total = summary_data$Whole_body_total[i],
      Proportion = summary_data$Whole_body_late_count[i] / summary_data$Whole_body_total[i]
    ))
  }
}

# 替换为更友好的物种名称显示
plot_data$Species <- gsub("_", " ", plot_data$Species)

# 计算置信区间
plot_data$lower <- NA
plot_data$upper <- NA
for(i in 1:nrow(plot_data)) {
  if(plot_data$Total[i] > 0) {
    ci <- binom.confint(plot_data$Count[i],
                        plot_data$Total[i],
                        conf.level = 0.95,
                        methods = "wilson")
    plot_data$lower[i] <- ci$lower
    plot_data$upper[i] <- ci$upper
  }
}

# 设置物种顺序
plot_data$Species <- factor(plot_data$Species, 
                            levels = c("Acyrthosiphon pisum", "Apis mellifera", 
                                       "Bombyx mori", "Drosophila mojavensis",
                                       "Helicoverpa armigera", "Tribolium castaneum"))

# 2. 进行二项式检验 (跨物种合并) -------------------------------------
test_results <- data.frame(
  Tissue = character(),
  early_count = numeric(),
  late_count = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(tissue in c("Head", "Whole_body")) {
  # 筛选当前组织的数据
  tissue_data <- plot_data[plot_data$Tissue == tissue, ]
  
  # 合并所有物种的早期和晚期基因计数
  early_data <- tissue_data[tissue_data$Stage == "early", ]
  late_data <- tissue_data[tissue_data$Stage == "late", ]
  
  total_early <- sum(early_data$Count)
  total_late <- sum(late_data$Count)
  
  # 进行二项式检验
  test <- prop.test(c(total_early, total_late), 
                    c(total_early + total_late, total_early + total_late))
  
  # 保存结果
  test_results <- rbind(test_results, data.frame(
    Tissue = tissue,
    early_count = total_early,
    late_count = total_late,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  ))
}

# 格式化P值为科学计数法
format_p_scientific <- function(p) {
  if(p < 0.001) {
    formatted <- format(p, scientific = TRUE, digits = 1)
    exponent <- as.numeric(strsplit(formatted, "e")[[1]][2])
    return(paste0("P = ", substr(formatted, 1, 1), " × 10^", exponent))
  } else {
    return(paste0("P = ", round(p, 3)))
  }
}

test_results$p_formatted <- sapply(test_results$p_value, format_p_scientific)

# 3. 为每个组织创建图形 ----------------------------------------------
for (tissue in c("Head", "Whole_body")) {
  # 筛选数据
  tissue_data <- plot_data[plot_data$Tissue == tissue, ]
  
  # 设置基准颜色
  base_color <- if(tissue == "Head") "#8DA0CB" else "#A65628"
  
  # 创建颜色映射
  color_values <- c(
    "early" = lighten(base_color, 0.3),
    "late" = darken(base_color, 0.2)
  )
  
  # 获取当前组织的p值
  p_text <- test_results$p_formatted[test_results$Tissue == tissue]
  
  # 创建图形
  p <- ggplot(tissue_data, aes(x = Species, y = Proportion * 100, fill = Stage)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
              vjust = 1.5, color = "black", size = 3) +
    geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                  width = 0.2,
                  position = position_dodge(width = 0.9)) +
    # 添加P值注释
    annotate("text", x = Inf, y = Inf, 
             label = p_text, 
             parse = FALSE,
             hjust = 1.1, 
             vjust = 1.5, 
             size = 4, 
             fontface = "plain") +
    scale_fill_manual(values = color_values,
                      labels = c("early" = "Early", "late" = "Late")) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.margin = margin(10, 15, 10, 10)
    ) +
    labs(
      x = "",
      y = "Percentage of genes",
      title = paste(tissue, "Tissue"),
      fill = ""
    )
  
  # 显示和保存图形
  print(p)
  ggsave(paste0("gene_distribution_", tissue, ".pdf"), 
         p, width = 6, height = 6, dpi = 300)
}

# 4. 打印详细的统计结果 --------------------------------------------
for(tissue in c("Head", "Whole_body")) {
  tissue_data <- plot_data[plot_data$Tissue == tissue, ]
  early_count <- sum(tissue_data$Count[tissue_data$Stage == "early"])
  late_count <- sum(tissue_data$Count[tissue_data$Stage == "late"])
  total_count <- early_count + late_count
  
  cat("\n", tissue, "组织统计结果:\n")
  cat("早期基因数量:", early_count, "(", round(early_count/total_count*100, 1), "%)\n")
  cat("晚期基因数量:", late_count, "(", round(late_count/total_count*100, 1), "%)\n")
  cat("总基因数量:", total_count, "\n")
  cat("二项式检验p值:", format_p_scientific(test_results$p_value[test_results$Tissue == tissue]), "\n")
}