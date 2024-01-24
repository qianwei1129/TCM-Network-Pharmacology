load("files/gene_clinical_data")
load("files/agging_HCC_tcm_gene")


library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

gene_data <- data[, unlist(aging_HCC_tcm_gene)]
data <- cbind(data[, c("age", "time", "flag")], gene_data)

# 预处理
data <- data.frame(lapply(data, as.numeric))

# 计算风险评分
gene_model <- coxph(Surv(age, time, flag) ~ ., data = data)
coefficients <- as.data.frame(gene_model$coefficients)

coefficients[is.na(coefficients)] <- 0
gene_var <- as.list(colnames(gene_data))

gene_var[2
