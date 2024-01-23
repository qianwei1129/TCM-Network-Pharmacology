# 聚类分析

load("gene_clinical_data")


library(pheatmap)
library(cluster)
library(ggplot2)



data[is.na(data)] <- 0
gene_data <- as.data.frame(data[, (6:ncol(data))])
clinical_data <- as.data.frame(data[, (1:5)])


