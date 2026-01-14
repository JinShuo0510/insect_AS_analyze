# Load required libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Define folders and parameters
folders <- list(Head = "Testis", Ovary = "Ovary")
title_map <- c(Head = "C. Head", Ovary = "D. Ovary")
x_map     <- c(Head = "log2(fold change)", Ovary = "F‑statistic")
use_x     <- c(Head = "logFC",      Ovary = "F_stat")
thr_map   <- list(
  Head  = list(h = -log10(0.05), v = c(-1,1)),
  Ovary = list(h = -log10(0.05), v = NULL)
)

plots <- list()
for (nm in names(folders)) {
  dir <- folders[[nm]]
  
  # 1. Read data
  dev <- read.delim(file.path(dir, "devDTx_nonDEG_with_exons.tsv"),
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  std <- read.delim(file.path(dir, "DTx_nonDEG_with_exons.tsv"),
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # 2. Ensure fold-change named logFC
  if (nm=="Head" && !"logFC" %in% colnames(dev)) {
    dev <- rename(dev, logFC = foldChange)
    std <- rename(std, logFC = foldChange)
  }
  
  # 3. Merge + annotate significance
  dev <- dev %>% mutate(DevSig = TRUE)
  std <- std %>% mutate(DevSig = FALSE)
  df  <- bind_rows(dev, std) %>% 
    mutate(
      negLogFDR = -log10(FDR),
      xval      = .data[[ use_x[nm] ]]
    )
  
  # 4. Top 10 transcripts by negLogFDR
  top10 <- df %>% 
    filter(DevSig) %>% 
    arrange(desc(negLogFDR)) %>% 
    distinct(id, .keep_all=TRUE) %>% 
    slice_head(n=10)
  
  # 5. Build volcano
  p <- ggplot(df, aes(x=xval, y=negLogFDR)) +
    # background points
    geom_point(aes(fill=DevSig), shape=21, color="white", size=2.5, alpha=0.8) +
    scale_fill_manual(values=c("FALSE"="grey80","TRUE"="#E64B35")) +
    # thresholds
    geom_hline(yintercept=thr_map[[nm]]$h, linetype="dashed", size=0.4) +
    { if (!is.null(thr_map[[nm]]$v))
      geom_vline(xintercept=thr_map[[nm]]$v, linetype="dashed", size=0.4)
    } +
    # labels
    geom_text_repel(data=top10, aes(label=id),
                    size=3, max.overlaps=15, segment.size=0.3) +
    # titles and axes
    labs(
      title = title_map[nm],
      x     = x_map[nm],
      y     = "-log10(FDR)"
    ) +
    # style to match bubble charts on left
    theme_minimal(base_size=12) +
    theme(
      plot.title        = element_text(hjust=0, vjust=1, face="bold"),
      axis.title        = element_text(face="bold"),
      panel.grid.major  = element_line(color="grey90", size=0.3),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill="white", color=NA),
      legend.position   = "none"
    )
  
  plots[[nm]] <- p
}

# Combine vertically
combined <- plots$Head / plots$Ovary + plot_layout(heights = c(1,1))

# Display
print(combined)

# Save to file
ggsave("volcano_plots_styled_testis.pdf", plots$Head,
       width=8, height=6, dpi=300)
ggsave("volcano_plots_styled_ovary.pdf", plots$Ovary,
       width=8, height=6, dpi=300)
