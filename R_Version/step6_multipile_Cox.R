rm(list = ls())
load("files/final_data")


library(survival)
library(survminer)
#install.packages("forestplot")
library(forestplot)
library(pathview)


final_gene_data <- as.data.frame(t(final_gene_data))

p <- pathview(gene.data = gse16873.d[, 1], pathway.id = "05200", species = "hsa",
              out.suffix = "figures/pathway_1", kegg.native = T)

print(p)

pathview(gene.data = final_gene_data, limit = list(gene = 2),# limit调整颜色bar的上下值
         pathway.id = "05200", 
         species = "hsa", 
         kegg.native = F,sign.pos= "bottomleft",#sign.pos更改签名的位置
         out.suffix = "4")


custom_control <- coxph.control(eps = 1e-09, 
                                iter.max = 10000, 
                                toler.chol = 1e-10)

surv_model <- coxph(Surv(time, flag) ~ ., 
                    data = data, 
                    control = custom_control)

cox_model <- coxph(Surv(time, flag) ~ ., 
                   data = data)

risk_scores <- predict(cox_model, newdata = data, type = "risk")


library(pROC)

roc_obj <- roc(data$flag, risk_scores)


roc_data <- data.frame(
  specificity = 1 - roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)


auc_value <- auc(roc_obj)


ggroc <- ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line(size = 1, color = "red") +  # 将线条颜色改为淡蓝色
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC =", round(auc_value, 2)), color = "black", size = 6) +
  labs(title = "ROC Curve", x = "1-Specificity", y = "Sensitivity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.background = element_rect(fill = "white", colour = "black")   # 白色背景，黑色边框
  ) +
  scale_x_continuous(breaks = seq(0, 2, by = 0.2), limits = c(0, 1)) +  # X轴刻度和限制
  scale_y_continuous(breaks = seq(0, 2, by = 0.2), limits = c(0, 1)) +  # Y轴刻度和限制
  theme(axis.line = element_line(color = "black")) + # 黑色的X轴和Y轴线条
  geom_abline(linetype = "dotted", color = "black", size = 1)

# 打印绘制的ROC曲线
print(ggroc)








# # 提取 Cox 回归结果中的系数、标准误差和p值
# coefficients <- coef(cox_model)
# se <- sqrt(diag(vcov(cox_model)))
# p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
# 
# # 创建一个数据框，包含系数、标准误差和p值
# forest_data <- data.frame(
#   HR = exp(coefficients),
#   lower = exp(coefficients - 1.96 * se),
#   upper = exp(coefficients + 1.96 * se),
#   p_value = p_values
# )
# 
# # 使用 forestplot() 函数生成森林图
# forestplot(forest_data, 
#            is.summary=c(TRUE, FALSE, FALSE, TRUE),
#            col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))








surv_model$coefficients

coefficients <- as.data.frame(surv_model$coefficients)
coefficients[is.na(coefficients)] <- 0
coefficients

gene_data <- data.frame(final_gene_data)
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

clinical_data <- data[, c(1, 2)]

save(data, 
     risk_group,
     gene_data, 
     coefficients, 
     clinical_data,
     final_surv_results, 
     file = "files/surv_final_data")

# 绘制生存曲线
rm(list = ls())
load("files/surv_final_data")

colnames(final_surv_results)

sum(data$flag)

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


for (i in (3:22)) {
  
  surv_median_value <- median(temp_data[[i]])
  surv_median_value
  
  # 分类并添加分组
  risk_group <- ifelse(temp_data[[i]] > surv_median_value,
                       "High Risk",
                       "Low Risk")
  
  # surv_results$risk_group 
  
  temp_data_1 <- as.data.frame(cbind(temp_data[[i]], risk_group))
  temp_surv_fit <- survfit(temp_surv_obj ~ risk_group, 
                           data = temp_data_1)
  temp_max <- max(temp_data$time)
  figure <- ggsurvplot(temp_surv_fit,
                       data = temp_data_1,
                       risk.table = TRUE, # 显示风险表
                       pval = TRUE, # 显示P值
                       conf.int = TRUE, # 显示置信区间
                       xlim = c(0, temp_max), # 可选：设置X轴的范围
                       xlab = "Times(months)", # 设置X轴标签
                       ylab = colnames(temp_data)[i], # 设置Y轴标签
  )
  print(figure)
}


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
