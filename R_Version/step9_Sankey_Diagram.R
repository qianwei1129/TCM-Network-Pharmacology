# install.packages("ggplot2")
# install.packages("ggalluvial")
# install.packages("networkD3")


rm(list = ls())

load("files/surv_final_data")


library(networkD3)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(jsonlite)


gene_name <- colnames(gene_data)


data <- readxl::read_xlsx("data/washed_data/TCM_Network.xlsx")


data <- data[data[[2]] %in% gene_name, ]
colnames(data) <- c("chem", "protein")


data$freq <- 1
data <- as.data.frame(data)


temp_data <- data %>%
  group_by(protein) %>%
  mutate(group_protein = cur_group_id()) %>%
  ungroup() 


temp_data <- temp_data %>%
  group_by(chem) %>%
  mutate(group_chem = cur_group_id()) %>%
  ungroup() 


temp_data <- temp_data %>%
  group_by(group_chem) %>%
  mutate(freq = n()) %>%
  ungroup()  


# df <- df %>%
#   mutate(group_chem = (group_chem - 1) %% 12 + 1)
# 
# 
# df <- df %>%
#   mutate(cycle_group = (group - 1) %% 6 + 1)


num_to_letter <- function(num) {
  letters <- c(LETTERS, paste0(rep(LETTERS, each=26), LETTERS))
  if (num > length(letters)) {
    stop("Number is too large to convert")
  }
  return(letters[num])
}


temp_data <- temp_data %>%
  mutate(group_chem = sapply(group_chem, num_to_letter))

temp_data <- temp_data %>%
  mutate(group_protein = sapply(group_protein, num_to_letter))


ggplot(temp_data, aes(y = freq, axis1 = chem, axis2 = protein)) +
  geom_alluvium(aes(fill = group_protein), width = 1/2) +
  geom_stratum(aes(fill = group_chem), width = 1/2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, vjust = -0.5) +
  scale_x_discrete(limits = c("chem", "protein"), expand = c(0, 0)) +
  scale_fill_manual(values = rep_len(rainbow(n = length(unique(temp_data$group_chem))), length(unique(temp_data$group_chem)))) +
  ggtitle("Your Title Here") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )


ggsave("plot.png", width = 10, height = 50, units = "cm")
