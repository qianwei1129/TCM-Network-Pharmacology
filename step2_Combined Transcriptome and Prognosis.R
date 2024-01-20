library(limma)
library(dplyr)




# 读取文件夹
data_path <- "data/washed_data/TCGA/gdc_download_20240120_110416.041068/"
data_files <- list.files(data_path, full.names = TRUE, recursive = TRUE)


# 剔除txt文件格式
txt_files <- list.files(data_path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
file.remove(txt_files)
data_files <- list.files(data_path, full.names = TRUE, recursive = TRUE)


# 读取临床信息
data_clinical <- read.table("data/washed_data/TCGA/clinical.cart.2024-01-20/clinical.tsv", 
                            header = TRUE, sep = "\t", fill = TRUE)


# 循环读取文件夹内容
# 转换为矩阵格式并清洗
# 用第一个文件创建标准
data <- read.table(data_files[1], sep = "\t", header = TRUE, check.names = FALSE)
data <- data[(5:nrow(data)), ]
data <- distinct(data, gene_name, .keep_all = TRUE)
rownames(data) <- data$gene_name

data <- as.matrix(data)
data_exp <- data[, (4:ncol(data))] 

dimnames <- list(rownames(data_exp), colnames(data_exp))
data_exp <- matrix(as.numeric(as.matrix(data_exp)), nrow=nrow(data_exp), dimnames=dimnames)

data_exp <- avereps(data_exp)
data_exp <- data_exp[rowMeans(data_exp)>0,]

data_exp <- subset(data_exp, select = "stranded_first")
colnames(data_exp) <- data_clinical$case_submitter_id[1]

exp <- as.data.frame(t(data_exp))


# 用循环针对每个文件
row_flag = 1

for (data_file in data_files) {
  data <- read.table(data_file, sep = "\t", header = TRUE, check.names = FALSE)
  data <- data[(5:nrow(data)), ]
  data <- distinct(data, gene_name, .keep_all = TRUE)
  rownames(data) <- data$gene_name
  
  data <- as.matrix(data)
  data_exp <- data[, (4:ncol(data))] 
  
  dimnames <- list(rownames(data_exp), colnames(data_exp))
  data_exp <- matrix(as.numeric(as.matrix(data_exp)), nrow=nrow(data_exp), dimnames=dimnames)
  
  data_exp <- avereps(data_exp)
  data_exp <- data_exp[rowMeans(data_exp)>0,]
  
  data_exp <- subset(data_exp, select = "stranded_first")
  colnames(data_exp) <- data_clinical$case_submitter_id[row_flag]
  
  row_flag = row_flag + 1
  
  data_exp <- as.data.frame(t(data_exp))
  
  print(data_file)
  exp <- bind_rows(exp, data_exp)
}

exp <- as.data.frame(t(exp))

writexl::write_xlsx(exp, "data/washed_data/TCGA/exp.xlsx")
View(exp)
