
# 转换宽表为长表，拆解列名中的结构
long_data <- target_genes[["Tribolium_castaneum"]] %>%
  pivot_longer(
    cols = matches(".+__.+(__.+)"),  # 匹配符合 pattern 的列
    names_to = "Age_Time_Tissue",    # 存储原始列名（如 "Adult__9_days__Head"）
    values_to = "PSI"                # PSI值的列名
  ) %>%
  separate(
    Age_Time_Tissue,
    into = c("Age", "Time", "Tissue"),  # 拆分三部分
    sep = "__",                         # 分隔符为双下划线
    remove = TRUE                       # 移除原始合并列
  )

long_data_cleaned <- long_data[!is.na(long_data$PSI), ]
long_data_cleaned <- long_data_cleaned %>%
  rowwise() %>%
  mutate(PSI_mean = mean(as.numeric(unlist(strsplit(PSI, ","))), na.rm = TRUE)) %>%
  ungroup()

# 自定义顺序（注意：重复项将自动去除）
custom_order <- c(
  "Egg-3_hours","Egg-8-16_hours","Egg-10-11_hours", "Egg-10-12_hours","Egg-10-13_hours","Egg-10-14_hours","Egg-10-16_hours","Egg-16-24_hours","Egg-24-48_hours","Egg-3_days","Egg-48-52_hours","Egg-52-56_hours"
)

# 数据处理与整合
final_data <- long_data_cleaned %>%
  # 1. 重新编码 Tissue
  mutate(Tissue = if_else(Tissue %in% c("Egg", "Whole_body", "Pupa"), "Whole body", Tissue),
         # 合并 Age 与 Time 生成 StageTime，注意levels采用 unique(custom_order)以剔除重复项
         StageTime = factor(paste0(Age, "-", Time), levels = unique(custom_order))) %>%
  # 2. 筛选指定区域（根据需要调整）
  filter(exonStart_0base == 6555661 | exonStart_0base == 6555648 ) %>%
  # 3. 拆分 PSI 列
  mutate(PSI_list = strsplit(PSI, ",")) %>%
  unnest(PSI_list) %>%
  mutate(PSI_value = as.numeric(PSI_list)) %>%
  # 4. 在每个 StageTime 与 Tissue 分组内去除离群值（IQR方法）
  group_by(StageTime, Tissue) %>%
  filter(PSI_value >= quantile(PSI_value, 0.25) - 1.5 * IQR(PSI_value),
         PSI_value <= quantile(PSI_value, 0.75) + 1.5 * IQR(PSI_value)) %>%
  ungroup() %>%
  # 5. 为每个分组内数据生成索引
  group_by(StageTime, Tissue) %>%
  mutate(PointIndex = row_number()) %>%
  ungroup() %>%
  mutate(StageTimeIndex = paste0(StageTime, "-", PointIndex)) %>%
  # 6. 按 StageTime 与 PointIndex 排序，并为每个 Tissue 内生成连续索引 NumericIndex
  arrange(StageTime, PointIndex) %>%
  group_by(Tissue) %>%
  mutate(NumericIndex = row_number()) %>%
  ungroup()

# ===== 新增颜色映射功能 =====

# 创建一个基础的自定义颜色集
color_pool <- c(
  "#E41A1C","#377EB8","#984EA3","#8DA0CB", "#4DAF4A",  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62",  "#A6D854", "#FFD92F", "#E5C494",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)

# 函数：更新颜色映射
update_color_mapping <- function(existing_mapping, new_tissues, color_pool) {
  # 获取已使用的颜色
  used_colors <- unname(existing_mapping)
  # 找出未使用的颜色
  available_colors <- setdiff(color_pool, used_colors)
  # 为新组织分配颜色
  new_mapping <- setNames(
    available_colors[1:length(new_tissues)],
    new_tissues
  )
  # 合并现有映射和新映射
  updated_mapping <- c(existing_mapping, new_mapping)
  return(updated_mapping)
}

# 尝试加载现有的颜色映射
if (file.exists("tissue_color_mapping.rds")) {
  tissue_color_mapping <- readRDS("tissue_color_mapping.rds")
} else {
  tissue_color_mapping <- setNames(character(), character())
}

# 获取当前数据中的所有唯一组织
current_tissues <- unique(final_data$Tissue)

# 找出新的组织
new_tissues <- setdiff(current_tissues, names(tissue_color_mapping))

# 更新颜色映射
tissue_color_mapping <- update_color_mapping(tissue_color_mapping, new_tissues, color_pool)

# 保存更新后的颜色映射
saveRDS(tissue_color_mapping, "tissue_color_mapping.rds")

# 使用一致的颜色映射绘图
p <- ggplot(final_data, aes(x = NumericIndex, y = PSI_value, color = Tissue, group = Tissue)) +
  geom_point(alpha = 0.7, size = 5) +  # 增加点的大小
  geom_smooth(se = FALSE, method = "loess", linewidth = 3) +  # 增加线段的粗细
  theme_bw() +
  labs(x = "Development", y = "PSI") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = tissue_color_mapping[current_tissues]) +  # 使用一致的颜色映射
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"  # 移除图例
  )

print(p)

ggsave(
  filename = "Tribolium_castaneum_Marf1.pdf",  # 文件路径和名称
  plot = p,                                   # 需要保存的 ggplot 对象
  width = 6,                                  # 图像宽度（英寸）
  height = 6,                                 # 图像高度（英寸）
  dpi = 300,                                  # 图像分辨率
  device = "pdf",                             # 文件格式
  limitsize = TRUE                            # 是否限制图像大小
)
