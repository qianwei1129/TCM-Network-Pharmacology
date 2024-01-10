# step 1:
# 在利用各个网站找出数据后，对数据进行清洗以获得type.csv


## data input
tcm_data <- c("黄芪", "甘草")
tcm_flag <- c("tcm", "tcm")
tcm_data <- data.frame(tcm = tcm_data, type = tcm_flag)
names(tcm_data) <- c("name", "type")

chem_data <- readxl::read_xlsx("data/pre/chemicals_gc.xlsx")[, c(1, 2)]
chem_data <- rbind(chem_data, readxl::read_xlsx("data/pre/chemicals_hq.xlsx")[, c(1, 2)])

Network_tcm_chem <- chem_data
names(Network_tcm_chem) <- c("source", "target")

chem_data <- unique(chem_data[, 2])
chem_data <- cbind(chem_data, rep("chem", nrow(chem_data)))
names(chem_data) <- c("name", "type")

pro_data <- readxl::read_xlsx("data/pre/proteins_gc.xlsx")
pro_data <- unique(rbind(pro_data, readxl::read_xlsx("data/pre/proteins_hq.xlsx")))

Network_chem_pro <- pro_data
names(Network_chem_pro) <- c("source", "target")

pro_data <- pro_data[, 2]
pro_data <- cbind(pro_data, rep("gene", nrow(pro_data)))
names(pro_data) <- c("name", "type") 


Network <- rbind(Network_chem_pro, Network_tcm_chem)


## data output
data_out <- rbind(tcm_data, chem_data, pro_data)
writexl::write_xlsx(data_out, "data/washed_data/Type.xlsx")
writexl::write_xlsx(Network, "data/washed_data/Network.xlsx")
