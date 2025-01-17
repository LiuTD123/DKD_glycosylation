library(sva)

setwd("D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-合并去批次化\\输入数据 样本按照case-control排")
GSE30122 <- read.table("GSE30122 gene id_exp result_mean.csv",header=T,row.names=1,sep=",")
#GSE30528 <- read.table("GSE30528 gene id_exp result_mean.csv",header=T,row.names=1,sep=",")
GSE96804 <- read.table("GSE96804 gene id_exp result_mean.csv",header=T,row.names=1,sep=",")
#GSE25906 <- read.table("GSE25906_result.csv",header=T,row.names=1,sep=",")

##intersect函数取交集
#mrna_names <- intersect(rownames(GSE30122),rownames(GSE30528),rownames(GSE96804))
mrna_names <- intersect(rownames(GSE30122),rownames(GSE96804))
#mrna_names1 <- intersect(rownames(GSE30528),rownames(GSE96804))
#mrna_intersect <- intersect(mrna_names,rownames(GSE96804))
mrna_intersect <- mrna_names
#mrna_intersect <- intersect(mrna_names,GSE96804)
##cbind进行合并
expr <- cbind(GSE30122[mrna_intersect,],GSE96804[mrna_intersect,])
write.table(expr,"gebing_tissue_shunxu.txt")
##43,157,94,60分别是提取出的每个数据集的样本数量，数据集1命名为batch1，数据集2命名为batch2，数据集3命名为batch3，数据集4命名为batch4
batch <- paste0("batch",rep(c(1,2),c(69,61)))
##17,26代表数据集1中有17个实验组，26个对照组，以此类推
tissue <- rep(c("case","control","case","control"),c(19,50,41,20))
table(batch,tissue)
mod <- model.matrix(~tissue)

expr_batch<-ComBat(dat=expr,batch=batch,mod=mod)
write.table(expr_batch,"mrna_expr_batch.txt")
write.table(expr_batch[,1:69],"GSE30122_after_batch.txt")
write.table(expr_batch[,70:130],"GSE96804_after_batch.txt")
#write.table(expr_batch[,92:152],"GSE96804_after_batch.txt")
#write.table(expr_batch[,295:354],"GSE25906_after_batch.txt")


library(gplots)
heatmap.2(as.matrix(expr),cexrow=0.8,cexcol=1.0)
heatmap.2(as.matrix(expr_batch),cexrow=0.8,cexcol=1.0)

library(ggbiplot)
library(devtools)
install_github("vqv/ggbiplot")
 
library(scatterplot3d)
pca.plot(expr)
pca.plot(expr_batch) 

ggbiplot(expr)
pca_res <- prcomp(expr, scale = TRUE)  

# 现在pca_res是一个prcomp对象，我们可以使用ggbiplot来绘制它  
g <- ggbiplot(pca_res, obs.scale = 1, var.scale = 1, groups = NULL)
filename <- "pca_plot.png"  
ggsave(filename, g, width = 8, height = 6, dpi = 300)

pca_res <- prcomp(expr, scale = TRUE)  


pca_res1 <- prcomp(expr_batch, scale = TRUE)  
# 现在pca_res是一个prcomp对象，我们可以使用ggbiplot来绘制它  
g <- ggbiplot(pca_res1, obs.scale = 1, var.scale = 1, groups = NULL)
filename <- "pca_plot_after.png"  
ggsave(filename, g, width = 8, height = 6, dpi = 300)

 
