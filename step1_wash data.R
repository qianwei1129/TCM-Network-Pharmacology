# step 1:
# 在利用各个网站找出数据后，对数据进行清洗以获得type.csv


## data input
tcm_data <- c("黄芪", "甘草")
tcm_flag <- c("tcm", "tcm")
tcm_data <- data.frame(tcm = tcm_data, type = tcm_flag)
names(tcm_data) <- c("name", "type")

chem_data <- readxl::read_xlsx("data/pre/step_1/chemicals_gc.xlsx")[, c(1, 2)]
chem_data <- rbind(chem_data, readxl::read_xlsx("data/pre/step_1/chemicals_hq.xlsx")[, c(1, 2)])

Network_tcm_chem <- chem_data
names(Network_tcm_chem) <- c("source", "target")

chem_data <- unique(chem_data[, 2])
chem_data <- cbind(chem_data, rep("chem", nrow(chem_data)))
names(chem_data) <- c("name", "type")

pro_data <- readxl::read_xlsx("data/pre/step_1/proteins_gc.xlsx")
pro_data <- unique(rbind(pro_data, readxl::read_xlsx("data/pre/step_1/proteins_hq.xlsx")))

Network_chem_pro <- pro_data
names(Network_chem_pro) <- c("source", "target")

pro_data <- pro_data[, 2]
pro_data <- cbind(pro_data, rep("gene", nrow(pro_data)))
names(pro_data) <- c("name", "type") 


Network <- rbind(Network_chem_pro, Network_tcm_chem)


## data output
data_out <- rbind(tcm_data, chem_data, pro_data)
writexl::write_xlsx(data_out, "data/pre/step_1/washed_data/TCM_Type.xlsx")
writexl::write_xlsx(Network, "data/pre/step_1/washed_data/TCM_Network.xlsx")


# step 2:
# 清洗衰老相关基因
msigdb_data <- clusterProfiler::read.gmt(gmtfile = "data/pre/step_2/msigdb_data.gmt") 
View(msigdb_data)

msigdb_data_up <- msigdb_data[grepl("UP$", msigdb_data$term), ]
msigdb_data_up$term <- "up"
rownames(msigdb_data_up) <- NULL

msigdb_data_down <- msigdb_data[grepl("DN", msigdb_data$term), ]
msigdb_data_down$term <- "down"
rownames(msigdb_data_down) <- NULL

View(msigdb_data_up)



HAGR_data_over <- readxl::read_xlsx("data/pre/step_2/HAGR_data.xlsx", sheet = 2)
HAGR_data_under <- readxl::read_xlsx("data/pre/step_2/HAGR_data.xlsx", sheet = 3)
View(HAGR_data_over)
View(HAGR_data_under)

colnames(HAGR_data_over)[1] <- "temp"
colnames(HAGR_data_under)[1] <- "temp"

HAGR_data <- rbind(HAGR_data_over, HAGR_data_under)
colnames(HAGR_data) <- HAGR_data[1, ]
HAGR_data <- HAGR_data[(2:nrow(HAGR_data)), ]

HAGR_data$Entrez[(1:449)] <- "UP"
HAGR_data$Entrez[(450:nrow(HAGR_data))] <- "down"
colnames(HAGR_data)[2] <- "flag"

View(HAGR_data)

writexl::write_xlsx(HAGR_data, "data/washed_data/HAGR_data.xlsx")
writexl::write_xlsx(msigdb_data, "data/washed_data/msigdb_data.xlsx")