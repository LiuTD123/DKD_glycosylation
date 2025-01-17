###按照contral组和STZ组做差异表达  limma包默认是排序靠后的 vs 排序靠前的
setwd("D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\3 糖尿病肾病正向分析-DEGs\\3 control和case差异表达分析")
df <- read.csv("DKD_exp_control_case.csv")
group_list <- c(rep('contral',46),rep('DKD',50))
group_list
library(limma)
design <- model.matrix(~ group_list)
design
fit <- lmFit(df, design)
fit <- eBayes(fit) 
allDiff <- topTable(fit, number = Inf) 
head(allDiff,12)
write.csv(allDiff,"DKD_exp_control_case_limma_diff.csv")
