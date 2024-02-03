rm(list = ls())

library(stringi)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

# KEGG
KEGG_data <- read.delim("data/gene_data/KEGG.txt")

colnames(KEGG_data)
KEGG_data <- KEGG_data[, c(2, 3, 5, 10)]

N <- 40
KEGG_data_topN <- KEGG_data[order(KEGG_data$PValue), ][1:N, ]

ggplot(KEGG_data_topN, aes(x = Count, y = reorder(Term, Count))) + 
  geom_point(aes(color = -log10(PValue), size = Count)) +
  scale_color_gradient(low = 'lightblue', high = 'red3') +  # 调整颜色梯度以匹配生存分析图
  scale_size(range = c(3, 5)) +  # 调整点的大小范围
  labs(title = "KEGG Pathway Enrichment Analysis",
       x = "Count",
       y = "Pathway",
       color = "-log10(PValue)",
       size = "Count") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),  # 调整图例标题的大小
    legend.text = element_text(size = 8)  # 调整图例文本的大小
  )

# install.packages("ggsci")
library(ggsci)
library(colorspace)
#
# red_color <- rgb(255, 165, 165, maxColorValue = 255) 
# white_color <- rgb(255, 255, 255, maxColorValue = 255) 
# blended_color <- blendColors(c(red_color, white_color), c(0.7, 0.3))

ggplot(KEGG_data_topN, aes(x = reorder(Term, -Count), y = Count, fill = -log10(PValue))) + 
  geom_col() +
  scale_fill_gradient(low = 'lightblue', high = 'red3')  +
  labs(title = "KEGG Pathway Enrichment Analysis", 
       x = "Pathway", 
       y = "Count",
       fill = "-log10(PValue)") +
  coord_flip() +  
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

KEGG_data_topN$Term
