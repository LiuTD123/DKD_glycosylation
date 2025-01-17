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

# Load the expression and trait data saved in the first part
lnames = load(file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\3 模块和性状相联系\\输入数据\\FemaleLiver-01-dataInput.RData");
lnames

# Load network data saved in the second part.
lnames = load(file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\3 模块和性状相联系\\输入数据\\FemaleLiver-02-networkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs0[1:6,1:6]

MEs = orderMEs(MEs0)
MEs[1:6,1:6]

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor[1:2,1:2]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue[1:2,1:2]


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
textMatrix[1:2,1:2]

par(mar = c(9, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Define variable weight containing the weight column of datTrait
DKD = as.data.frame(datTraits$DKD);
names(DKD) = "DKD"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, DKD, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(DKD), sep="");
names(GSPvalue) = paste("p.GS.", names(DKD), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

module = "salmon"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DKD",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

names(datExpr)[moduleColors=="salmon"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

annot = read.csv(file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\3 模块和性状相联系\\输入数据\\GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$symbol)

# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

# Create the starting data frame
#geneInfo0 = data.frame(substanceBXH = probes,
                      #geneSymbol = annot$gene_symbol[probes2annot],
#LocusLinkID = annot$LocusLinkID[probes2annot],
#moduleColor = moduleColors,
#geneTraitSignificance,
#GSPvalue)
geneInfo0 = data.frame(
                       geneSymbol = annot$symbol[probes2annot],
                      
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, DKD, use = "p")));

# Add module membership information in the chosen order
for(mod in 1:ncol(geneModuleMembership)){
   oldNames = names(geneInfo0)
   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]]);
   names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),paste("p.MM.", modNames[modOrder[mod]], sep=""))
   print(mod)  
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.DKD));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

write.csv(geneInfo, file = "D:\\xqm2\\xqm\\xqm\\2024\\6月\\0628\\1 糖尿病肾病正向分析-WCGNA\\3 模块和性状相联系\\\\结果\\geneInfo.csv")


