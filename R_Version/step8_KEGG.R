rm(list = ls())

library(stringi)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
# install.packages("ggsci")
library(ggsci)
library(colorspace)

# KEGG
KEGG_data <- read.delim("data/gene_data/KEGG.txt")

gene_list <- strsplit(KEGG_data$Genes, ",")
num_genes <- sapply(gene_list, length)

GeneRatio <- num_genes / 20

KEGG_data$GeneRatio <- GeneRatio

colnames(KEGG_data)
KEGG_data <- KEGG_data[, c(2, 3, 5, 6, 14)]

colnames(KEGG_data) <- c("Description", "Count", "pvalue", "geneID", "GeneRatio")
KEGG_data <- KEGG_data[, c("Description"	, "GeneRatio",	"pvalue",	"geneID",	"Count")]

KEGG_data$Description <- gsub(".*:", "", KEGG_data$Description)
KEGG_data$geneID <- gsub(",", "/", KEGG_data$geneID)

writexl::write_xlsx(KEGG_data, "data/gene_data/KEGG_Data.xlsx")


N <- 36
KEGG_data_topN <- KEGG_data[order(KEGG_data$PValue), ][1:N, ]
KEGG_data_topN$Term <- gsub(".*:", "", KEGG_data_topN$Term)

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

-log10(KEGG_data_topN$PValue)



# GO
GO_BP <- read.delim("data/gene_data/BP.txt")
GO_CC <- read.delim("data/gene_data/CC.txt")
GO_MF <- read.delim("data/gene_data/MF.txt")

GO_Data <- rbind(GO_BP, GO_CC, GO_MF)

# writexl::write_xlsx(GO_Data, "data/gene_data/GO_Data.xlsx")

GO_Data <- GO_Data[GO_Data$PValue < 0.05, ]

GO_Data <- GO_Data[, c(1, 2, 10)]
colnames(GO_Data) <- c("subgroup", "GOterm", "Enrichment Score")

GO_Data$GOterm <- gsub(".*~", "", GO_Data$GOterm)
colnames(GO_Data)
GO_Data$subgroup <- gsub(".*_(.*?)_.*", "\\1", GO_Data$subgroup)

# 假设数据框名为df，列名为column_name
GO_Data$subgroup <- gsub("BP", "Biological process", GO_Data$subgroup)
GO_Data$subgroup <- gsub("CC", "Cellular component", GO_Data$subgroup)
GO_Data$subgroup <- gsub("MF", "Molecular function", GO_Data$subgroup)

colnames(GO_Data)

GO_Data <- GO_Data[, c("GOterm", "subgroup", "Enrichment Score")]

writexl::write_xlsx(GO_Data, "data/gene_data/GO_Data.xlsx")

GO_Data <- readxl::read_xlsx("data/gene_data/GO_Data.xlsx")
GO_Data$`Enrichment Score` <- log10(GO_Data$`Enrichment Score`)


barplot(GO_Data, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
        showCategory =10, #只显示前10
        split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图

a <- ggplot(ego, aes(x = reorder(Term, -Fold.Enrichment), y = Fold.Enrichment)) +
  geom_col() +
  coord_flip() + # 翻转坐标轴，使术语在y轴上显示
  labs(x = "GO Term", y = "Fold Enrichment", title = "GO Enrichment Analysis Results") +
  theme_minimal()

print(a)

library(ggplot2)
library(org.Hs.eg.db)
library(enrichplot)


load("files/final_data")


gene_name <- colnames(final_gene_data)

ego <- enrichGO(gene = gene_name,
                OrgDb = org.Hs.eg.db, # 选择对应你研究物种的数据库
                keyType = "SYMBOL", # 这里的类型应该与你的基因列表匹配
                pvalueCutoff = 0.05, # p值截断
                qvalueCutoff = 0.05,
                readable = T) 

temp_ego <- as.data.frame(ego)

barplot(ego)
ego <- (as.data.frame(ego))
summary(ego)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

p <- barplot(ego, showCategory=10, split="ONTOLOGY", measure="GeneRatio", colorBy="p.adjust") + 
  facet_wrap(~ONTOLOGY, scales="free") + 
  theme_minimal()
