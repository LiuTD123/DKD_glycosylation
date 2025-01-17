#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
rm(list=ls())
rm(datExpr0)
# Load the WGCNA package

library(dynamicTreeCut)
library(fastcluster)
library(xfun)
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the female liver data set

femData = read.csv("D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\1 数据清洗\\输入数据\\mrna_expr_batch.csv");
dim(femData);#行为基因，列为不同样本的基因表达量或其他信息
names(femData);


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
#提取出表达量的数据 ，删去不需要的数据重新生成矩阵;datExpr0 就是一个以每行为样本，一列为一个基因的数据框
datExpr0 = as.data.frame(t(femData[, -c(1:1)]));#提取加转置
names(datExpr0) = femData$symbol; #基因名字
rownames(datExpr0) = names(femData)[-c(1:1)];#样品名字


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
#检查缺失值和异常值
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#如果gsg$allOK的结果为TRUE，证明没有缺失值，可以直接下一步。如果为FALSE，则需要用以下函数进行删除缺失值。

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
#聚类所有样本，观察是否有离群值或异常值
sampleTree = hclust(dist(datExpr0), method = "average");

# 绘制示例树:打开一个尺寸为12 × 9英寸的图形输出窗口
# 如果窗口太大或太小，用户应该改变尺寸
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#dev.off()

#=====================================================================================
#
#  Code chunk 6
#如果有离群值，则要删去离群的样本，如果没有则跳过下一步。
#=====================================================================================
#删除离群样本
# 划定需要剪切的枝长
abline(h = 100, col = "red");

#这时候会从高度为15这里横切，把离群样本分开
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)

#保留非离群(clust==1)的样本
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ] #去除离群值后的数据
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
#载入表型数据
traitData = read.csv("D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\1 数据清洗\\输入数据\\ClinicalTraits.csv");
dim(traitData)#每行是一个样本，每列是一种信息
names(traitData)
traitData[1:4,1:4]

#删除我们不需要的数据
#allTraits = traitData[, -c(31, 16)];
allTraits = traitData[, c(2:4) ];
dim(allTraits)
names(allTraits)

#因为表型数据和基因表达量数据中并不是全部信息都能匹配上，例如表型数据中的样本可能在基因表达量数据中不存在，这时候需要将两者进行匹配，找共有的数据才能分析。
#将临床表征数据和表达数据进行匹配（用样本名字进行匹配）
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$symbol);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================
# 可视化表型数据与基因表达量数据的联系，重构样本聚类树
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
#pdf(file = "Plots/sampleClustering_addTraits.pdf", width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
#dev.off()


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================
#图片结果解释了临床数据和基因表达量的关联程度，保存数据。
#颜色越深，代表这个表型数据与这个样本的基因表达量关系越密切。
save(datExpr, datTraits, file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\1 数据清洗\\结果\\FemaleLiver-01-dataInput.RData")


