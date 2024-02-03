rm(list = ls())

if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("eulerr", quietly = TRUE)) install.packages("eulerr")

library(eulerr)
library(UpSetR)
library(ggvenn)
library(limma)
library(VennDiagram)
library(venn)

load("files/surv_final_data")

# 绘制韦恩图

aging_data_1 <- readxl::read_xlsx("data/washed_data/HAGR_data.xlsx")
aging_data_2 <- readxl::read_xlsx("data/washed_data/msigdb_data.xlsx")
aging_gene <- as.data.frame(c(aging_data_1$Gene, aging_data_2$gene))
colnames(aging_gene) <- c("gene")

tcm_data <- readxl::read_xlsx("data/washed_data/TCM_Type.xlsx")
tcm_gene <- as.data.frame(tcm_data[(85:nrow(tcm_data)), 1])
colnames(tcm_gene) <- c("gene")

HCC_data <- readxl::read_xlsx("data/washed_data/exp_clinical")
HCC_gene <- HCC_data$...1
HCC_gene <- HCC_gene[1:(length(HCC_gene)-5)]
HCC_gene <- as.data.frame(HCC_gene)
colnames(HCC_gene) <- c("gene")

tcm_gene <- data.frame(tcm_gene)
HCC_gene <- data.frame(HCC_gene)
aging_gene <- data.frame(aging_gene)

typeof(tcm_gene[1, 1])
typeof(aging_gene[1, 1])
typeof(HCC_gene[1, 1])

HCC_gene <- as.list(HCC_gene)
tcm_gene <- as.list(tcm_gene)
aging_gene <- as.list(aging_gene)

HCC_gene <- readxl::read_xlsx("data/gene_data/HCC_gene.xlsx")
aging_gene <- readxl::read_xlsx("data/gene_data/aging_gene.xlsx")
tcm_gene <- readxl::read_xlsx("data/gene_data/tcm_gene.xlsx")

temp <- unique(intersect(HCC_gene$gene, tcm_gene$gene))
intersect_gene <- as.data.frame(unique(intersect(temp, aging_gene$gene)))

writexl::write_xlsx(intersect_gene, "data/gene_data/intersect_gene.xlsx")

input <- list("HCC" = HCC_gene[, 1], 
              "TCM" = tcm_gene[, 1],
              "Aging" = aging_gene[, 1])

venn.plot <- venn.diagram(
  x = input,
  category.names = c("HCC", "TCM", "Aging"),
  filename = "figures/VennDiagram.png",
  imagetype = "png",
  height = 3500,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 5,
  col = "black",
  fill = c("cornflowerblue", "green", "red"),
  alpha = 0.3,
  label.col = "black",
  cex = 3,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "red"),
  cat.cex = 3,
  category.col = c("darkblue", "darkgreen", "red")
)

writexl::write_xlsx(aging_gene, "data/gene_data/aging_gene.xlsx")
writexl::write_xlsx(tcm_gene, "data/gene_data/tcm_gene.xlsx")
writexl::write_xlsx(HCC_gene, "data/gene_data/HCC_gene.xlsx")



# 绘制t-SNE图
library(Rtsne)
library(ggplot2)
library(limma)

set.seed(1234)

gene_risk_data <- cbind(gene_data, risk_group)

metadata <- as.data.frame(rownames(gene_risk_data))

tsne_out <- Rtsne(
  gene_data,
  dims = 2,
  pca = TRUE,
  max_iter= 1000,
  theta = 0.4,
  perplexity = 10,
  check_duplicates = FALSE,
  verbose = F
)

str(tsne_out)

tsne_result <- as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("tSNE1", "tSNE2")

ggplot(tsne_result,
       aes(tSNE1, tSNE2, color = risk_group)) +
  geom_point()

