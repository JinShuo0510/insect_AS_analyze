###############################################################################
# 0. 依赖包 ── 若缺请先安装：
# install.packages(c("dplyr","tidyr","purrr","ggplot2","binom","scales","patchwork"))
###############################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(binom)
library(scales)
library(patchwork)

###############################################################################
# 1. 全局设置
###############################################################################
stages   <- c("Cell","Egg","Larva","Pupa","Adult")   # 阶段顺序
datasets <- list(
  "New exon"    = result_table_detail_result,
  "Alternative" = result_table_alternative_detail_result
)
palette  <- c("New exon"="#E41A1C",
              "Alternative"="#377EB8")

###############################################################################
# 2. 基础工具
###############################################################################
is_yes <- function(x){
  tolower(trimws(as.character(x))) %in% c("yes","true","1")
}

###############################################################################
# 3. 统计函数
###############################################################################
## 3.1 devAS 比例
devAS_stat <- function(df, label){
  map_dfr(stages, function(s){
    ex_col <- paste0(s, "_exist")
    is_col <- paste0(s, "_is_devas")
    mask   <- is_yes(df[[ex_col]])
    n      <- sum(mask)
    k      <- sum(is_yes(df[[is_col]]) & mask, na.rm = TRUE)
    if (n == 0) return(tibble())
    ci <- binom.confint(k, n, methods = "wilson")
    tibble(dataset=label, stage=s,
           p=k/n, lower=ci$lower, upper=ci$upper, n=n)
  })
}

## 3.2 平均 PSI（每个外显子在同一阶段只计一次）
psi_stat <- function(df, label){
  psi_cols <- grep("_PSI$", names(df), value = TRUE)
  df_yes   <- df %>% mutate(across(ends_with("_exist"), is_yes))
  
  map_dfr(stages, function(s){
    ex_col <- paste0(s, "_exist")
    cols_s <- grep(paste0("^", s, "-"), psi_cols, value = TRUE)
    
    df_stage <- df_yes %>% filter(.data[[ex_col]])
    if (nrow(df_stage) == 0) return(tibble())
    
    row_mean <- rowMeans(df_stage[, cols_s], na.rm = TRUE)
    
    n   <- length(row_mean)
    m   <- mean(row_mean, na.rm = TRUE)
    se  <- sd(row_mean, na.rm = TRUE) / sqrt(n)
    tibble(dataset=label, stage=s, n=n,
           mean=m, lower=m-1.96*se, upper=m+1.96*se)
  })
}

## 3.3 通用布尔列比例
bool_prop_stat <- function(df, label, bool_col){
  map_dfr(stages, function(s){
    ex_col <- paste0(s, "_exist")
    mask   <- is_yes(df[[ex_col]])
    n      <- sum(mask)
    k      <- sum(is_yes(df[[bool_col]]) & mask, na.rm = TRUE)
    if (n == 0) return(tibble())
    ci <- binom.confint(k, n, methods = "wilson")
    tibble(dataset=label, stage=s,
           p=k/n, lower=ci$lower, upper=ci$upper, n=n)
  })
}

###############################################################################
# 4. 数据汇总
###############################################################################
plot_devAS <- imap_dfr(datasets, devAS_stat)                                    %>%
  mutate(stage=factor(stage, levels=stages))
plot_PSI   <- imap_dfr(datasets, psi_stat)                                      %>%
  mutate(stage=factor(stage, levels=stages))
plot_div3  <- imap_dfr(datasets, bool_prop_stat, "is_divide_by3")               %>%
  mutate(stage=factor(stage, levels=stages))
plot_cdsOL <- imap_dfr(datasets, bool_prop_stat, "is_cds_overleap")             %>%
  mutate(stage=factor(stage, levels=stages))

###############################################################################
# 5. 绘图
###############################################################################
make_prop_plot <- function(df, ylab){
  ggplot(df, aes(stage, p, group=dataset, color=dataset, fill=dataset)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.25, colour=NA) +
    geom_line(size=1) +
    geom_point(size=2) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    scale_y_continuous(labels=percent_format(1)) +
    labs(x=NULL, y=ylab, color=NULL, fill=NULL) +
    theme_classic(base_size=13) +
    theme(legend.position="none")
}

p1 <- make_prop_plot(plot_devAS, "Percentage of DevAS")
p3 <- make_prop_plot(plot_div3,  "Percentage of 3N")
p4 <- make_prop_plot(plot_cdsOL, "Percentage of CDS overleap")

p2 <- ggplot(plot_PSI, aes(stage, mean, group=dataset, color=dataset, fill=dataset)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.25, colour=NA) +
  geom_line(size=1) + geom_point(size=2) +
  scale_color_manual(values=palette) +
  scale_fill_manual(values=palette) +
  labs(x=NULL, y="Mean PSI", color=NULL, fill=NULL) +
  theme_classic(base_size=13) +
  theme(legend.position="none")

###############################################################################
# 6. 多面板输出（2×2）
###############################################################################
fig6 <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("d-f.pdf",fig6, height=6,width=6)
###############################################################################
# 7. 可选：分母一致性自检
###############################################################################
check_n <- function(df, label){
  map_dfr(stages, function(s){
    ex_col <- paste0(s, "_exist")
    tibble(dataset=label, stage=s, n_exist_yes=sum(is_yes(df[[ex_col]])))
  })
}
true_denoms <- bind_rows(check_n(result_table_detail_result,"New exon"),
                         check_n(result_table_alternative_detail_result,"Alternative"))
reported_denoms <- bind_rows(plot_devAS, plot_div3, plot_cdsOL) %>% select(dataset, stage, n) %>% distinct()
inner_join(true_denoms, reported_denoms, by=c("dataset","stage")) %>%
  mutate(match = n_exist_yes == n) %>% print()
