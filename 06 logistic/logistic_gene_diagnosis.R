# 加载必要的包  
library(readr)  



# 读取CSV文件  

data = read.csv("D:\\xqm2\\xqm\\xqm\\2024\\9月\\0905\\2 对筛出来的26个gene进行logistic分析，看是不是与疾病显著相关\\26 gene_exp_symbol_diagnosis_no_yes.csv")

data$DiseaseStatus <- as.factor(data$DiseaseStatus)

model <- glm(DiseaseStatus ~ BTG2, data = data, family = binomial,control=list(maxit=100))

# 查看模型摘要
summary(model)

# 可以进一步提取模型的系数、p值等
coefficients(model)

# 如果想要得到每个基因与疾病状态的Logistic回归结果，可以使用循环遍历所有基因  
single_gene_results <- list()  
for (gene in names(data)[-ncol(data)]) {  
  # 拟合模型  
  model <- glm(DiseaseStatus ~ ., data = data[, c(gene, "DiseaseStatus")], family = binomial,control=list(maxit=100))  
  
  # 提取系数和对应的统计量  
  coef_summary <- summary(model)$coefficients[2, ] # 排除截距项  
  single_gene_results[[gene]] <- coef_summary  
}  

# 将结果转换为数据框  
single_gene_results_df <- do.call(rbind, lapply(single_gene_results, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))  

# 为数据框添加基因名称作为行名  
rownames(single_gene_results_df) <- names(single_gene_results)  

# 查看结果  
print(single_gene_results_df)  

# 如果需要，可以将结果保存为CSV文件  
write.csv(single_gene_results_df, file = "logistic_results2.csv", row.names = TRUE)
