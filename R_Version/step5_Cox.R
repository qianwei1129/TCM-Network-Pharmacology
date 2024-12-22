rm(list = ls())
load("files/surv_pre_data")

library(survival)
library(glmnet)
library(randomForest)

significant_genes_pearson <- data.frame(Gene=character(), 
                                        Correlation=numeric(), 
                                        PValue=numeric(), 
                                        stringsAsFactors=FALSE)

significant_genes_spearman <- data.frame(Gene=character(), 
                                         Correlation=numeric(), 
                                         PValue=numeric(), 
                                         stringsAsFactors=FALSE)

significant_genes <- data.frame(Gene=character(), 
                                HR=numeric(), 
                                LowerCI=numeric(), 
                                UpperCI=numeric(), 
                                PValue=numeric(), 
                                stringsAsFactors=FALSE)


significant_genes_spearman_partial <- data.frame(Gene=character(), 
                                                 Partial_Correlation=numeric(), 
                                                 PValue=numeric(), 
                                                 stringsAsFactors=FALSE)


for(gene in colnames(gene_data)) {
  
  # Pearson correlation analysis
  pearson_corr <- cor.test(gene_data[[gene]], data$time, method = "pearson")
  if(!is.na(pearson_corr$estimate)) {
    significant_genes_pearson <- rbind(significant_genes_pearson, 
                                       data.frame(Gene=gene, 
                                                  Correlation=pearson_corr$estimate, 
                                                  PValue=pearson_corr$p.value))
  }
  
  # Spearman correlation analysis
  spearman_corr <- cor.test(gene_data[[gene]], data$time, method = "spearman")
  if(!is.na(spearman_corr$estimate)) {
    significant_genes_spearman <- rbind(significant_genes_spearman, 
                                        data.frame(Gene=gene, 
                                                   Correlation=spearman_corr$estimate, 
                                                   PValue=spearman_corr$p.value))
  }
}


for(gene in colnames(gene_data)) {
  
  formula <- as.formula(paste("Surv(time, flag) ~", gene))
  
  fit <- coxph(formula, data=data)
  
  summary_fit <- summary(fit)
  
  if(summary_fit$coefficients[5] < 0.1) {
    significant_genes <- rbind(significant_genes, 
                               data.frame(Gene=gene, 
                                          HR=summary_fit$coefficients[2], 
                                          LowerCI=summary_fit$conf.int[,"lower .95"], 
                                          UpperCI=summary_fit$conf.int[,"upper .95"],
                                          PValue=summary_fit$coefficients[5]))
  }
}

for(gene in colnames(gene_data)) {
  # Spearman partial correlation analysis
  spearman_partial_corr <- cor.test(gene_data[[gene]], data$time, 
                                    method = "spearman", partial = TRUE)
  if(!is.na(spearman_partial_corr$estimate)) {
    significant_genes_spearman_partial <- rbind(significant_genes_spearman_partial, 
                                                data.frame(Gene=gene, 
                                                           Partial_Correlation=spearman_partial_corr$estimate, 
                                                           PValue=spearman_partial_corr$p.value))
  }
}

print("Cox单变量相关")
print(significant_genes)

print("Pearson相关分析")
temp_pearson <- significant_genes_pearson[significant_genes_pearson$PValue < 0.1, ]
print(temp_pearson)

print("Spearman相关分析")
temp_spearman <- significant_genes_spearman[significant_genes_spearman$PValue < 0.1, ]
print(temp_spearman)

print("偏相关分析")
temp_cor <- significant_genes_spearman_partial[significant_genes_spearman_partial$PValue < 0.1, ]
print(temp_cor)

final_data <- c(significant_genes$Gene, 
                temp_spearman$Gene,
                temp_spearman$Gene,
                temp_cor$Gene)

final_data <- unique(final_data)
final_data

final_gene_data <- data[colnames(gene_data) %in% final_data]
final_clinical_data <- clinical_data

temp_final <- c(final_data, "time", "flag")
final_data <- data[colnames(data) %in% temp_final]

data <- cbind(final_clinical_data,
              final_gene_data)

save(final_gene_data, final_clinical_data, data, file = "files/final_data")
