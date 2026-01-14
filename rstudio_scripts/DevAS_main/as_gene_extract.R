library(stringr)
library(tidyverse)

gene_ids_se <- all_event_wide_lists_withoutavg_nondeNA_1_8[["SE"]][["Bombyx_mori"]]$GeneID
gene_id_a3ss <- all_event_wide_lists_withoutavg_nondeNA_1_8[["A3SS"]][["Bombyx_mori"]]$GeneID
gene_id_a5ss <- all_event_wide_lists_withoutavg_nondeNA_1_8[["A5SS"]][["Bombyx_mori"]]$GeneID
gene_id_ri <- all_event_wide_lists_withoutavg_nondeNA_1_8[["RI"]][["Bombyx_mori"]]$GeneID
gene_id_mxe <- all_event_wide_lists_withoutavg_nondeNA_1_8[["MXE"]][["Bombyx_mori"]]$GeneID

gene_ids_as <- c(gene_ids_se, gene_id_a3ss, gene_id_a5ss, gene_id_ri, gene_id_mxe)


as_unique_gene_ids <- unique(gene_ids_as)



devas_gene_ids_Developmental_tissue <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Developmental_tissue"]]$Event_ID, "_", simplify = TRUE)[,2]
devas_gene_ids_Fat_body <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Fat_body"]]$Event_ID, "_", simplify = TRUE)[,2]
devas_gene_ids_Gland <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Gland"]]$Event_ID, "_", simplify = TRUE)[,2]
devas_gene_ids_Head <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Head"]]$Event_ID, "_", simplify = TRUE)[,2]
devas_gene_ids_Ovary <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Ovary"]]$Event_ID, "_", simplify = TRUE)[,2]

devas_gene_ids_Testis <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Testis"]]$Event_ID, "_", simplify = TRUE)[,2]
devas_gene_ids_Whole_body <- str_split(find_devas_results_1_8_full_deleteNULL[["Bombyx_mori"]][["Whole_body"]]$Event_ID, "_", simplify = TRUE)[,2]

devas_gene_ids <- c(devas_gene_ids_Developmental_tissue, devas_gene_ids_Fat_body, devas_gene_ids_Gland, devas_gene_ids_Head, devas_gene_ids_Ovary, devas_gene_ids_Testis, devas_gene_ids_Whole_body)

devas_unique_gene_ids <- unique(devas_gene_ids)

library(dplyr)

data <- read.delim("tau_output.txt", header = TRUE, stringsAsFactors = FALSE)

# 2. 使用 mutate() 添加新列
data <- data %>%
  mutate(
    AS = ifelse(GeneID %in% as_unique_gene_ids, "yes", "no"),
    DevAS = ifelse(GeneID %in% devas_unique_gene_ids, "yes", "no")
  )

library(ggplot2)

# 1. 为数据添加新分类变量：DevAS > AS > NonAS
data$Category <- ifelse(data$DevAS == "yes", "DevAS",
                        ifelse(data$AS == "yes", "AS", "NonAS"))

# 2. 根据 Tau 值分箱（0 到 1，每 0.1 一个区间，标签为区间中值）
data$Tau_bin <- cut(data$Tau, breaks = seq(0, 1, by = 0.1),
                    include.lowest = TRUE, right = TRUE,
                    labels = c("0.05", "0.15", "0.25", "0.35", "0.45",
                               "0.55", "0.65", "0.75", "0.85", "0.95"))

# 3. 汇总每个 Tau_bin 和 Category 组合下的基因数量
df_summary <- data %>%
  group_by(Tau_bin, Category) %>%
  summarise(NumberOfGenes = n(), .groups = "drop") %>%
  # 保证每个 Tau_bin 与 Category 的组合都存在，不存在的填 0
  complete(Tau_bin, Category, fill = list(NumberOfGenes = 0))

# 4. 设置因子顺序
# 对 x 轴 Tau_bin 按从小到大的顺序排序
df_summary$Tau_bin <- factor(df_summary$Tau_bin,
                             levels = c("0.05", "0.15", "0.25", "0.35", "0.45",
                                        "0.55", "0.65", "0.75", "0.85", "0.95"))
# 控制堆叠顺序：ggplot2 中，堆叠顺序与因子水平的“逆序”有关，
# 因此将 Category 设置为：c("DevAS", "AS", "NonAS")
df_summary$Category <- factor(df_summary$Category,
                              levels = c("DevAS", "AS", "NonAS"))

# 5. 自定义颜色
color_values <- c(
  "NonAS" = "#D9D9D9",
  "AS"    = "#80B1D3",
  "DevAS" = "#457B9D"
)

# 6. 绘制堆叠柱形图
p <- ggplot(df_summary, aes(x = Tau_bin, y = NumberOfGenes, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_values) +
  labs(x = "Tau", y = "Number of Genes") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 12),
    axis.title         = element_text(size = 14),
    legend.title       = element_blank()
  )

# 显示图形
print(p)

ggsave("zhuxingtu123.pdf", p, width = 8, height = 6, dpi = 300)

