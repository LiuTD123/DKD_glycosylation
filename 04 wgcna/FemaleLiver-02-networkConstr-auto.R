#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

rm(list=ls())

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the data saved in the first part
lnames = load(file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\2 网络构建\\输入数据\\FemaleLiver-01-dataInput.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# 选择软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# 无标度拓扑拟合指数
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##系统推荐的软阈值
sft$powerEstimate
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
# 报错：https://blog.csdn.net/liyunfan00/article/details/91686840
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)
# 改回去
cor<-stats::cor

class(net)
names(net)
table(net$colors)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\2 网络构建\\结果\\FemaleLiver-02-networkConstruction-auto.RData")


