# 衰老相关基因与HCC基因、中药基因求交集（降维）
rm(list = ls())

load("files/exp_clinical")

library(dplyr)


# aging HCC
aging_data_1 <- readxl::read_xlsx("data/washed_data/HAGR_data.xlsx")
aging_data_2 <- readxl::read_xlsx("data/washed_data/msigdb_data.xlsx")

any(is.na(aging_data_1))
any(is.na(aging_data_2))

aging_gene <- as.data.frame(c(aging_data_1$Gene, aging_data_2$gene))
colnames(aging_gene) <- c("gene")


load("files/gene_clinical_data")

HCC_gene <- as.data.frame(rownames(exp_clinical))
colnames(HCC_gene) <- c("gene")

aging_HCC_gene <- unique(intersect(aging_gene, HCC_gene))


# aging HCC tcm
tcm_data <- readxl::read_xlsx("data/washed_data/TCM_Type.xlsx")
tcm_gene <- as.data.frame(tcm_data[(85:nrow(tcm_data)), 1])

colnames(tcm_gene) <- c("gene")
aging_HCC_tcm_gene <- unique(intersect(aging_HCC_gene, tcm_gene))

save(aging_HCC_tcm_gene, file = "files/agging_HCC_tcm_gene")
