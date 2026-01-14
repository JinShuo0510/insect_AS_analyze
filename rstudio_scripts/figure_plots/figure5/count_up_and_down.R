# 加载必要的包
library(dplyr)

# 定义辅助函数：计算外显子长度是否小于28 nt
check_exon_length <- function(exonStart, exonEnd) {
  exon_length <- abs(exonEnd - exonStart) + 1  # 计算外显子长度
  return(exon_length < 28)
}

# 假设 Bombyx_mori_event_fits_all_exon 为列表，各元素为不同组织或条件下的tibble数据
# 对每个tibble进行过滤，仅保留外显子长度小于28 nt的记录
micro_exon_list <- lapply(Bombyx_mori_event_fits_all_exon, function(df) {
  df %>% filter(check_exon_length(exonStart_0base, exonEnd))
})

# 如果需要将所有筛选后的结果整合到一个数据框中，可以使用 bind_rows，并添加一个标识列（如组织名称）
micro_exon_all <- bind_rows(micro_exon_list, .id = "Tissue")

# 查看结果（例如：前六行）
print(head(micro_exon_all))
