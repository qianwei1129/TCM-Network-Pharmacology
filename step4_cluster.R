# 聚类分析
load("files/gene_clinical_data")


library(pheatmap)
library(cluster)
library(ggplot2)
library(factoextra)


# 数据预处理
data[is.na(data)] <- 0

gene_data <- as.data.frame(t(data[, (6:ncol(data))]))
clinical_data <- as.data.frame(t(data[, (1:5)]))

gene_data[] <- lapply(gene_data, function(x) {
  as.numeric(unlist(x))
})

clinical_data[] <- lapply(clinical_data, function(x) {
  as.numeric(unlist(x))
})

gene_data <- as.matrix(gene_data)
clinical_data <- as.matrix(clinical_data)
typeof(gene_data)


# 计算Hopkins统计量
res <- get_clust_tendency(gene_data, 40, graph = TRUE)
res$hopkins_stat


# PCA降维处理
pca_result <- prcomp(gene_data, scale. = FALSE)
summary(pca_result)

plot(cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2), 
     xlab = "Number of Principal Components", 
     ylab = "Cumulative Proportion of Variance Explained", type = "b")

threshold <- 0.8
cum_var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
num_components <- which(cum_var_explained >= threshold)[1]

num_components

gene_data_pca <- pca_result$x[, 1:num_components]
pca_result$rotation

data <- cbind(gene_data_pca, clinical_data)

# 将pca载荷数据应用于临床数据
clinical_data_scale <- scale(clinical_data)

typeof(clinical_data_scale)
typeof(pca_result$rotation)

clinical_data_scale[is.na(clinical_data_scale)] <- 0
matrix2[is.na(matrix2)] <- 0

temp <- pca_result$rotation
clinical_data_prcomp <- (clinical_data_scale) %*% (temp)

clinical_data_pca <- clinical_data_prcomp[, 1:num_components]
data <- rbind(as.data.frame(clinical_data_pca), gene_data_pca)
data <- t(data)

colnames(data)

save(data, file = "files/data_pca")

# 估计聚合簇数
load("data_pca")

set.seed(1234)

data_data_frame <- as.data.frame(data)
View(data_data_frame[(1:3), (1:50)])

clusters <- kmeans(data, centers = 2)
km.res <- kmeans(data, 2, nstart = 50)
fviz_cluster(km.res, data)
