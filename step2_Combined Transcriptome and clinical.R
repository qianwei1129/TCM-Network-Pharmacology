library(limma)
library(dplyr)
library(openxlsx)
library(sets)
library(xlsx)



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
exp$gene_name <- row.names(exp)

write.xlsx(exp, "data/washed_data/TCGA/exp.xlsx", rowNames = TRUE)

data <- readxl::read_xlsx("data/washed_data/TCGA/exp.xlsx")
data <- as.data.frame(data)
rownames(data) <- data$gene_name
write.xlsx(data, "data/washed_data/TCGA/exp.xlsx", rowNames = TRUE)


## 添加临床信息
colnames(data_clinical)
clinical <- data_clinical[, c("case_submitter_id", # 编号
                              "age_at_index", # 年龄（是得病日期还是死亡年龄）
                              "premature_at_birth", # 出生与死亡日期
                              "days_to_death", # 出生与死亡日期第二方案
                              "occupation_duration_years", # 出生与死亡日期第三方案、是否死亡
                              "country_of_residence_at_enrollment",
                              "percent_tumor_invasion")]

for (i in 1:nrow(clinical)) {
  if (clinical$premature_at_birth[i] == "--	--") {
    if (clinical$days_to_death[i] != "--	--") {
      clinical$premature_at_birth[i] <- clinical$days_to_death[i]
    } else {
      clinical$premature_at_birth[i] <- clinical$occupation_duration_years[i]
    }
  }
}

clinical$occupation_duration_years <- ifelse(clinical$occupation_duration_years == "Dead", 1, 0)

clinical$country_of_residence_at_enrollment <- as.numeric(clinical$country_of_residence_at_enrollment)

extract_numbers <- function(s) {
  # 提取所有数字并分配到两列
  nums <- regmatches(s, gregexpr("[0-9]+", s))[[1]]
  if (length(nums) == 0) {
    return(c(NA, NA))
  } else if (length(nums) == 1) {
    return(c(nums, nums))
  } else {
    return(nums[1:2])
  }
}

clinical[c("first", "second")] <- t(apply(clinical, 1, function(x) extract_numbers(x["premature_at_birth"])))

clinical <- clinical[, c("case_submitter_id", # 编号
                         "age_at_index", # 年龄
                         "first", 
                         "second", 
                         "country_of_residence_at_enrollment",
                         "percent_tumor_invasion",
                         "occupation_duration_years")] # 出生与死亡日期第三方案、是否死亡

clinical <- clinical[complete.cases(clinical["first"]), ]
clinical <- clinical[clinical$first != 0, ]

colnames(clinical) <- c("id", "age", "first", "second", "to_death", "to_follow_up", "flag")

clinical$to_follow_up <- as.numeric(clinical$to_follow_up)
clinical$time <- ifelse(is.na(clinical$to_follow_up), clinical$to_death, clinical$to_follow_up)

clinical <- as.data.frame(t(clinical[, c("id", # 编号
                                         "age", # 年龄
                                         "first", 
                                         "second", 
                                         "time",
                                         "flag")])) # 出生与死亡日期第三方案、是否死亡

exp <- readxl::read_xlsx("data/washed_data/TCGA/exp.xlsx")
gene <- readxl::read_xlsx("data/washed_data/TCM_Type.xlsx")

colnames(clinical) <- clinical[1, ]
clinical <- clinical[-1, ]

colnames(exp)
data_rownames <- c(as.list(exp$...1), "age", "first", "second", "time", "flag")

common_columns <- intersect(colnames(exp), colnames(clinical))
exp <- exp[common_columns]
clinical <- clinical[common_columns]

exp_clinical <- rbind(exp, clinical)
exp_clinical$id <- data_rownames

rownames(exp_clinical) <- exp_clinical$id
exp_clinical <- as.data.frame(exp_clinical)

gene_name <- exp_clinical$id
save(exp_clinical, gene_name, file = "files/exp_clinical")



