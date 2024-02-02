rm(list = ls())

load("files/exp_clinical")
load("files/agging_HCC_tcm_gene")

library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(glmnet)
library(mice)
library(dplyr)
library(car)
library(randomForestSRC)
library(caret)


## 预处理
rownames(exp_clinical) <- data_rownames
exp_clinical <- replace(exp_clinical, is.na(exp_clinical), 0)
temp <- as.data.frame(t(exp_clinical))

temp <- data.frame(lapply(temp, as.numeric))
temp <- replace(temp, is.na(temp), 0)

colnames(temp) <- data_rownames
colnames(temp)

gene_data <- as.data.frame(temp[, unlist(aging_HCC_tcm_gene)])
temp_gene_data <- as.data.frame(lapply(gene_data, as.numeric))

str(temp_gene_data)

temp_gene_data <- scale(temp_gene_data)
gene_data <- as.data.frame(temp_gene_data)
# temp <- scale(gene_data)

rownames(gene_data) <- t(colnames(exp_clinical))

exp_clinical <- as.data.frame(t(exp_clinical))
clinical_data <- as.data.frame(cbind(exp_clinical$time, exp_clinical$flag))

colnames(clinical_data) <- c("time", "flag")

clinical_data <- as.data.frame(lapply(clinical_data, as.numeric))
clinical_data[is.na(clinical_data)] <- 0

sum(clinical_data$flag)

rownames(clinical_data) <- rownames(gene_data)

rows_to_remove <- rownames(clinical_data[clinical_data$time == 0, ])

gene_data <- gene_data[!rownames(gene_data) %in% rows_to_remove, ]
clinical_data <- clinical_data[!rownames(clinical_data) %in% rows_to_remove, ]
# clinical_data <- clinical_data[!rownames(clinical_data) %in% rows_to_remove, ]

# rownames(gene_data) <- t(colnames(exp_clinical))

# clinical_data <- as.matrix(temp[, c("time", "flag")])
# clinical_data <- as.data.frame(clinical_data)
# 
# clinical_data <- as.data.frame(lapply(clinical_data, as.numeric))
# sum(clinical_data$flag)
# 
data <- cbind(clinical_data, gene_data)
# # temp_data <- cbind(clinical_data, temp)
# 
# sum(is.na(clinical_data))
# sum(is.na(gene_data))

# gene_data <- preProcess(gene_data, method = "range")

save(data, gene_data, clinical_data, file = "files/surv_pre_data")


## 风险回归
rm(list = ls())
load("files/surv_pre_data")

str(clinical_data)
str(gene_data)

custom_control <- coxph.control(eps = 1e-09, 
                                iter.max = 10000, 
                                toler.chol = 1e-10)

surv_model <- coxph(Surv(time, flag) ~ ., 
                    data = data, 
                    control = custom_control)

# temp_surv_model <- coxph(
#  Surv(time, flag) ~ ., 
#  data = temp_data, 
#  control = coxph.control(iter.max = 10000)
#)

surv_model$coefficients

# surv_model <- temp_surv_model

coefficients <- as.data.frame(surv_model$coefficients)
coefficients[is.na(coefficients)] <- 0
coefficients

gene_data <- data.frame(gene_data)
gene_coef_data <- rbind(gene_data, t(coefficients))
rownames(gene_coef_data)[nrow(gene_coef_data)] <- c("coef")

surv_results <- sweep(gene_data, 2, t(coefficients), `*`)
surv_row_sums <- rowSums(surv_results)

surv_row_sums[is.na(surv_row_sums)] <- 0
surv_row_sums

surv_median_value <- median(surv_row_sums)
surv_median_value


# 分类并添加分组
risk_group <- ifelse(surv_row_sums > surv_median_value,
                     "High Risk",
                     "Low Risk")
surv_results$risk_group <- risk_group
final_surv_results <- surv_results[, !apply(surv_results, 
                                      2,
                                      function(x) all(x == 0))]

save(data, 
     risk_group,
     gene_data, 
     clinical_data, 
     coefficients, 
     final_surv_results, 
     file = "files/surv_final_data")

# 绘制生存曲线
rm(list = ls())

load("files/surv_final_data")

colnames(final_surv_results)


temp_data <- data
temp_data$time <- ceiling(temp_data$time/30)
# temp_data <- lapply(temp_data, as.numeric)

# 剔除大于5年的
rows_to_remove <- rownames(temp_data[temp_data$time > 60, ])
temp_data <- temp_data[!rownames(temp_data) %in% rows_to_remove, ]
clinical_data <- clinical_data[!rownames(clinical_data) %in% rows_to_remove, ]
risk_group <- as.data.frame(risk_group)
risk_group <- risk_group[!rownames(risk_group) %in% rows_to_remove, ]


temp_surv_obj <- Surv(temp_data$time, temp_data$flag)
temp_data <- cbind(temp_data, risk_group)
temp_surv_fit <- survfit(temp_surv_obj ~ risk_group, 
                    data = temp_data)
temp_max <- max(temp_data$time)
ggsurvplot(temp_surv_fit, 
           data = temp_data,
           risk.table = TRUE, # 显示风险表
           pval = TRUE, # 显示P值
           conf.int = TRUE, # 显示置信区间
           xlim = c(0, temp_max), # 可选：设置X轴的范围
           xlab = "Times(months)", # 设置X轴标签
           ylab = "Survival probability", # 设置Y轴标签
)

coefficients$gene_name <- rownames(coefficients)
colnames(coefficients)[1] <- c("coef")
writexl::write_xlsx(coefficients, "data/washed_data/coef.xlsx")

surv_obj <- Surv(data$time, data$flag)
surv_fit <- survfit(surv_obj ~ risk_group, 
                    data = data)
max <- max(data$time)
ggsurvplot(surv_fit, 
           data = final_surv_results,
           risk.table = TRUE, # 显示风险表
           pval = TRUE, # 显示P值
           conf.int = TRUE, # 显示置信区间
           xlim = c(0, max), # 可选：设置X轴的范围
           xlab = "Times(days)", # 设置X轴标签
           ylab = "Survival probability", # 设置Y轴标签
           title = "Kaplan-Meier 生存曲线" # 设置图标题
)


