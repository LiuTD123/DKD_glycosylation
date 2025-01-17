# 加载必要的库  
library(easyTCGA)  
library(dplyr)  
library(clusterProfiler)  
library(msigdbr)  
library(enrichplot)  
library(ggplot2)  
library(edgeR)

setwd = ("D:\\xqm2\\xqm\\xqm\\2024\\9月\\0903\\1 Z008-L002 GSEA结果整理\\前后5个通路可视化\\输入数据")
# 读取表达数据  
expr_data <- read.csv("DKD_exp_control_case.csv", row.names = 1)  # 假设第一列是gene id  

# 读取样本分组数据  
sample_groups <- read.csv("group.csv", stringsAsFactors = FALSE)  
names(sample_groups) <- c("sample", "status")  # 确保列名正确  

# 将样本分组转换为因子  
sample_groups$status <- factor(sample_groups$status, levels = c("control", "DKD"))  

# 确保表达数据中的样本名与分组数据中的一致  
expr_data <- expr_data[, names(expr_data) %in% sample_groups$sample]  

# 检查数据维度  
dim(expr_data)  
head(sample_groups)

# 如果easyTCGA不适用，我们可以使用limma包  
library(limma)  

if (any(expr_data < 0)) {  
  warning("表达数据中包含负值，将负值转换为0")  
  expr_data[expr_data < 0] <- 0  
}

# 创建设计矩阵  
design <- model.matrix(~ 0 + factor(sample_groups$status))  
colnames(design) <- levels(factor(sample_groups$status))  

# 创建DGEList对象  
dgelist <- DGEList(counts = as.matrix(expr_data))  

# 使用limma进行差异表达分析  
v <- voom(dgelist, design)  
fit <- lmFit(v, design)  
cont.matrix <- makeContrasts(DKDvsControl = DKD - control, levels = design)  
fit2 <- contrasts.fit(fit, cont.matrix)  
fit2 <- eBayes(fit2)  

# 提取结果  
deg_res <- topTable(fit2, adjust = "fdr", number = Inf, sort.by = "p")  

# 查看前几行结果  
head(deg_res)  

# 保存结果  
#write.csv(deg_res, "DKD_vs_Control_deg_res.csv", row.names = FALSE)


# 提取结果  
deg_res <- topTable(fit2, adjust = "fdr", number = Inf, sort.by = "p")  

# 添加gene_symbol列，内容为行名  
deg_res$gene_symbol <- rownames(deg_res)  

# 查看前几行结果以确认更改  
head(deg_res)  

# 保存结果，此时gene_symbol将是最后一列（因为我们在最后添加了它）  
write.csv(deg_res, "DKD_vs_Control_deg_res1.csv", row.names = FALSE)

suppressMessages(library(clusterProfiler))

gene_entrezid <- bitr(geneID = deg_res$gene_symbol
                      , fromType = "SYMBOL" # 从symbol
                      , toTyp = "ENTREZID" # 转成ENTREZID
                      , OrgDb = "org.Hs.eg.db"
)

head(gene_entrezid)

gene_entrezid <- merge(gene_entrezid,deg_res,by.x = "SYMBOL", by.y = "gene_symbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist <- sort(genelist,decreasing = T)

head(genelist)

library(msigdbr)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

gsea_res <- GSEA(genelist,
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH"
)

gsea_res_symbol <- setReadable(gsea_res,"org.Hs.eg.db",keyType = "ENTREZID")
write.csv(gsea_res_symbol, "gsea_res_symbol.csv")

library(GseaVis)
gseaNb(gsea_res, geneSetID = "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT")

####多个通路拼在一起
terms <- gsea_res@result$ID[1:4]##第一条到第四条通路

gseaplot_list <- lapply(terms, function(x){
  gseaNb(object = gsea_res,
         geneSetID = x,
         termWidth = 30,
         addPval = T,
         pvalX = 0.75,
         pvalY = 0.6
  )
})

# 可以直接拼
cowplot::plot_grid(plotlist=gseaplot_list, ncol = 2)

####多个通路拼在一起
# 定义你想要展示的通路的索引  
terms_indices <- c(248, 279, 227, 229)  

# 使用索引从gsea_res@result$ID中提取这些通路的ID  
terms <- gsea_res@result$ID[terms_indices]  

# 使用lapply为每个选定的通路生成GSEA图  
gseaplot_list <- lapply(terms, function(x){  
  gseaNb(object = gsea_res,  
         geneSetID = x,  
         termWidth = 30,  
         addPval = TRUE,  # 注意这里使用TRUE而不是T  
         pvalX = 0.75,  
         pvalY = 0.6  
  )  
})  

# 使用cowplot的plot_grid函数来拼接这些图，这里假设你想把它们放在两列  
# 根据通路的数量，你可能需要调整ncol的值  
# 如果有4个通路，ncol=2是合适的  
# 如果有更多或更少的通路，你可能需要调整ncol的值  
cowplot::plot_grid(plotlist = gseaplot_list, ncol = 2)

####可视化多个通路
terms <- gsea_res@result$ID[1:3]

gseaNb(object = gsea_res,
       geneSetID = terms,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05,pvalY = 0.05)