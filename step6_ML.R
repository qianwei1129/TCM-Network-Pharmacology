load("files/gene_clinical_data")
load("files/agging_HCC_tcm_gene")


library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(glmnet)
library(mice)


# 预处理
data <- data.frame(lapply(data, as.numeric))

gene_data <- data[, unlist(aging_HCC_tcm_gene)]
clinical_data <- as.matrix(data[, c("time", "flag")])

data <- cbind(data[, c("time", "flag")], gene_data)


# 风险回归
surv_model <- coxph(Surv(time, flag) ~ ., data = data)
coefficients <- as.data.frame(surv_model$coefficients)
coefficients[is.na(coefficients)] <- 0

gene_coef_data <- rbind(gene_data, t(coefficients))
gene_coef_data <- gene_coef_data[-197, ]
rownames(gene_coef_data)[197] <- c("coef")

surv_results <- sweep(gene_data, 2, t(coefficients), `*`)
surv_row_sums <- rowSums(surv_results)

surv_row_sums

surv_median_value <- median(surv_row_sums)
surv_median_value


# lasso回归
gene_data <- mice(gene_data, m = 5, method = 'pmm', maxit = 10)
gene_data_completed <- complete(gene_data, 1)

class(gene_data_completed)
class(clinical_data)

dim(gene_data_completed)
dim(clinical_data)

lasso_model <-
  glmnet(gene_data_completed, 
         clinical_data, 
         alpha = 1)
