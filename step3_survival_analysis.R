load("exp_clinical")

install.packages("survival") 
install.packages("survminer")
install.packages("glmnet")
install.packages("caret")

library(survival)
library(survminer)
library(glmnet)
library(dplyr)
library(caret)

# 绘制Kaplan-Meier曲线
exp_clinical <- as.data.frame(t(exp_clinical))
exp_clinical$time <- as.numeric(exp_clinical$time)
exp_clinical <- exp_clinical[!is.na(exp_clinical$time), ]

# 观察小范围内数据框
View(exp_clinical[(1:50), (1:50)])
View(gene_data[(1:50), (1:50)])
View(normalized_data[(1:50), (1:50)])

gene_data <- exp_clinical[(1:(nrow(exp_clinical)-5)), ]
clinical_data <- exp_clinical[((nrow(exp_clinical)-4) : nrow(exp_clinical)), ]


# 选择至少在一个样本中表达量高于该基因90%分位数的基因
sapply(gene_data, class)

gene_data[is.na(gene_data)] <- 0
gene_data[] <- lapply(gene_data, function(x) as.numeric(as.character(x)))

gene_data_scale <- as.data.frame(scale(gene_data))

normalize_minmax <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

gene_data <- as.data.frame(lapply(gene_data_scale, normalize_minmax))

names(gene_data) <- names(clinical_data)
data <- rbind(clinical_data, gene_data)
data <- as.data.frame(t(data))

rownames(data)[6:nrow(data)] <- rownames(exp_clinical)[1:(nrow(exp_clinical)-5)]
save(data, file = "gene_clinical_data")


# 做相关性分析
load("gene_clinical_data")

View(data[(1:50), (1:50)])
