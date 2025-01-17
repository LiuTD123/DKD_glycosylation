options(timeout = Inf)

rm(list = ls())
setwd("D:\\xqm2\\xqm\\xqm\\2024\\9月\\0910\\2 高低分组免疫浸润\\成功")

library(psych)
library(ggcorrplot)
library(corrplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(reshape2)
library(magrittr)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(immunedeconv)
library(tidyr)
library(IOBR)
library(reshape2)
library(remotes)
library(EPIC)
install.packages("immunedeconv")
remotes::install_github("omnideconv/immunedeconv")
install_github('icbi-lab/immunedeconv')
install.packages("purrr")
install.packages("C:\\Users\\Administrator\\AppData\\Local\\Temp\\RtmpkV5KEV\\file57b051af221c\\immunedeconv_2.1.0.tar.gz", repos = NULL, type = "source")
update.packages("purrr")
remove.packages("purrr")
install.packages("purrr", version = "1.0.2")



library(testit)
library(data.tree)
library(limSolve)
library(ComICS)
library(pracma)
library(sva)
library(singscore)
library(quantiseqr)
library(purrr)
library(EstimateLibraryComplexity)


# Estimate----
gene_expr <- readRDS("D:\\xqm2\\xqm\\xqm\\2024\\9月\\0910\\2 高低分组免疫浸润\\实验\\all gene_exp DKD group.rds")
riskScore <- read.table("D:\\xqm2\\xqm\\xqm\\2024\\9月\\0910\\2 高低分组免疫浸润\\实验\\risk group.txt",head = T)

# data_group <- riskScore[,c(258,259,260,261)]
# dat_count <- dat_count[,rownames(data_group)]
#
# risk <- readRDS("../00_rawdat/00_temp/risk_train.rds")
risk_sample<- rownames(riskScore)

head(gene_expr)
colnames(gene_expr)
rownames(gene_expr)
eset_stad <- as.data.frame(gene_expr[,risk_sample])
# 免疫浸润分析
# IOBR包继承了8种免疫浸润分析的方法，使用方法如下
t1 <- proc.time()

estimate<-deconvo_tme(eset = eset_stad, method = "estimate")

t2 <- proc.time()
t <- t2-t1
print(paste0('执行时间：',t[3][[1]],'秒'))
estimate_dat <- estimate[,2:4]
estimate_puri <- estimate[,5]

rownames(estimate_dat) <- estimate$ID
rownames(estimate_puri) <- estimate$ID

boxDat <- cbind(estimate_dat,riskScore[,7,drop = F])
boxpuri <- cbind(estimate_puri,riskScore[,7,drop = F])

boxDat_melt <- melt(boxDat,by = c("risk"))
boxpuri_melt <- melt(boxpuri,by = c("risk"))

boxDat_melt$risk %<>% factor(.,levels = c("low","high"))
boxpuri_melt$risk %<>% factor(.,levels = c("low","high"))

boxDat_melt<-transform(boxDat_melt,pos=as.numeric(variable),
                       pos_adj=ifelse(risk == "low", -0.66,0.66))
boxpuri_melt<-transform(boxpuri_melt,pos=as.numeric(variable),
                        pos_adj=ifelse(risk == "low", -0.66,0.66))

stat_res <- boxDat_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ risk) %>%
  adjust_pvalue(method = "BH") %>%    # method BH == fdr
  add_significance("p")

statpuri_res <- boxpuri_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ risk) %>%
  adjust_pvalue(method = "BH") %>%    # method BH == fdr
  add_significance("p")

statpuri_res

openxlsx::write.xlsx(stat_res, file = "02.estimate_sig_0.66.xlsx",rowNames = T)
openxlsx::write.xlsx(estimate_dat, file = "01.estimate_dat_0.66.xlsx",rowNames = T)

boxDat_melt$risk <- factor(boxDat_melt$risk,levels = c("low","high"))
boxpuri_melt$risk <- factor(boxpuri_melt$risk,levels = c("low","high"))

p <-  ggplot(boxDat_melt,
             aes(x=variable,y=value,
                 fill=risk,
                 outlier.shape = NA,
                 bxp.errorbar = T)) +
  geom_violin()+
  geom_boxplot(width=0.3,
               alpha=0.8,
               position = position_dodge(0.9))+
  annotate(geom = "text", x = stat_res$variable, y = 5000, size = 5, family = "Times",
           label =as.character(stat_res$p.signif)) +
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.title.x =element_text(size=16,family = "Times", face = "bold"),
        axis.text.x =element_text(angle=45,size=16,hjust = 1,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"))+
  theme(
    legend.title=element_text(size=15) , legend.text=element_text(size=14))+
  labs(x="",y="Score")+
  ggsci::scale_fill_jama()+
  coord_flip()

ggsave(filename = '04.box_plot_0.66.pdf',p,w=8,h=6)
ggsave(filename = '04.box_plot_0.66.png',p,w=8,h=6)

p_puri <-  ggplot(boxpuri_melt,
                  aes(x=variable,y=value,
                      fill=risk,
                      outlier.shape = NA,
                      bxp.errorbar = T)) +
  geom_violin()+
  geom_boxplot(width=0.3,
               alpha=0.8,
               position = position_dodge(0.9))+
  annotate(geom = "text", x = statpuri_res$variable, y = 1.2, size = 5, family = "Times",
           label =as.character(statpuri_res$p.signif)) +
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.title.x =element_text(size=16,family = "Times", face = "bold"),
        axis.text.x =element_text(angle=45,size=16,hjust = 1,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"))+
  theme(
    legend.title=element_text(size=15) , legend.text=element_text(size=14))+
  labs(x="",y="Score")+
  ggsci::scale_fill_jama()+
  coord_flip()

ggsave(filename = '04.puri_box_plot_0.66.pdf',p_puri,w=8,h=6)
ggsave(filename = '04.puri_box_plot_0.66.png',p_puri,w=8,h=6)

# ssGSEA----

gene_expr <- readRDS("all gene_exp DKD group.rds")
riskScore <- read.table("risk group.txt",head = T)

risk <- riskScore
# risk <- readRDS("../00_rawdat/00_temp/risk_train.rds") 
risk_sample<- rownames(risk)

TrainExp <- readRDS("all gene_exp DKD group.rds")
TrainGroup <-group <- risk[,7,drop = F] %>% rownames_to_column(.,var = "sample")
multi <- read.table("Lasso_Coefficients1.xls", header = T) %>% as.data.frame()#cox_result_step2.rds;mul_cox_result.rds

modelgene <- read.csv("hub_gene.csv")
multi <- multi[multi$gene %in% modelgene$x,]
intersect_gene <- multi$gene

TrainExp <- TrainExp[,TrainGroup$sample]
gene_set <- read.table("mmc3.txt", header = T, sep ="\t")
gene_list <- split(as.matrix(gene_set)[,1], gene_set[,2])
expr <- TrainExp
dat.final <- as.matrix(expr)
# 
# param	
# A parameter object of one of the following classes:
#   
#   A gsvaParam object built using the constructor function gsvaParam. This object will trigger gsva() to use the GSVA algorithm by Hänzelmann et al. (2013).
# 
# A plageParam object built using the constructor function plageParam. This object will trigger gsva() to use the PLAGE algorithm by Tomfohr et al. (2005).
# 
# A zscoreParam object built using the constructor function zscoreParam. This object will trigger gsva() to use the combined z-score algorithm by Lee et al. (2008).
# 
# A ssgseaParam object built using the constructor function ssgseaParam. This object will trigger gsva() to use the ssGSEA algorithm by Barbie et al. (2009).
gsvapar <- ssgseaParam(dat.final, gene_list)

ssgsea_score = gsva(gsvapar,
                    verbose = TRUE)

write.csv(ssgsea_score,
          file = "01.ssgsea_result_cell 1.csv",
          quote = F)

group <- risk[,7,drop = F] %>% rownames_to_column(.,var = "sample")
colnames(group) <- c('sample', 'group')
group_case<-group[group$group=='high',]
group_case<-as.character(group_case$sample)
group_control<-group[group$group=='low',]
group_control<-as.character(group_control$sample)

tiics_result <- read.csv('01.ssgsea_result_cell 1.csv',check.names = F, row.names = 1)

tiics_result <- tiics_result[,group$sample]%>% as.matrix()

pvalue = padj <- matrix(0, nrow(tiics_result), 1)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_case],
                                       tiics_result[i, group_control])$p.value
  # log2FoldChange[i, 1] = mean(tiics_result[i, group_case]) - 
  #   mean(tiics_result[i, group_control])
}

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_case],
                                       tiics_result[i, group_control])$p.value
  
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(#log2FoldChange, 
  pvalue, 
  padj,
  row.names = rownames(tiics_result))

rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$pvalue<0.05),]
dim(diff_Table)
# 21 4

write.csv(rTable,
          file = "02.cell_all_wilcox_test 1.csv",
          quote = F,
          row.names = F)
write.csv(diff_Table,
          file = "02.cell_diff_wilcox_test 1.csv",
          quote = F,
          row.names = F)
# 箱线图----

# library(Ipaper)

cell.data <- data.frame(Immune_Cell=rownames(tiics_result), 
                        tiics_result, 
                        pvalue=rTable$pvalue)
plot.cell <- cell.data[which(cell.data$pvalue<0.05),]
#plot.cell<-cell.data
diff_tiics <- rownames(plot.cell)
violin_dat <- gather(plot.cell, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_control,
                           "low", "high") 
violin_dat$Group <- factor(violin_dat$Group, levels = c("low", "high"))
violin_dat <- violin_dat[,-2]
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  # geom_violin(trim=T,color=alpha('black',alpha = 0.5)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  # stat_boxplot(geom="errorbar",
  #              width=0.1,
  #              position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = 21,
               outlier.fill = "black",
               outlier.size = 0.5,
               color=ggplot2::alpha('black',alpha = 0.5))+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  # geom_point(aes(fill = Group),
  #            size = 0.05,
  #            position = position_dodge(0.9))+
  scale_fill_manual(values= c("#45a9b8","#f76a56"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     method = "wilcox.test", #没有直接用差异分析的结果，但检验方法是一样的
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        text=element_text(family = 'Times'),
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs
ggsave('02.ssgsea_Box 1.pdf',boxplot_diff_TIICs,w=14,h=6)
ggsave('02.ssgsea_Box 1.png',boxplot_diff_TIICs,w=14,h=6)

# 
# # ----------------------------cibersort评估细胞丰度-免疫细胞比例--------------------------
# library(CIBERSORT)
# library(ggplot2)
# 
# a = gene_expr
# library(stringr)
# 
# k = !duplicated(rownames(a));table(k)
# exp = a[k,]
# 
# rownames(exp) = unique(rownames(a))
# colnames(exp) = str_remove(colnames(exp),"TPM")
# 
# exp[1:4,1:4]
# exp2 = as.data.frame(exp)
# exp2 = rownames_to_column(exp2)
# write.table(exp2,file = "exp.txt",row.names = F,quote = F,sep = "\t")
# 
# # 内置数据LM22.txt，记录了22种免疫细胞的基因表达特征数据。
# lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
# 
# exp2 <- na.omit(exp2)
# TME.results = cibersort(lm22f, 
#                         "exp.txt",
#                         QN = F
# )
# 
# save(TME.results,file = "TME.result.RData")
# TME.results[1:4,1:4]
# #             B cells naive B cells memory Plasma cells T cells CD8
# # GSM404005    0.00000000    0.118715823    0.5064640 0.002298742
# # GSM404006    0.02075989    0.000000000    0.7821873 0.007487980
# # GSM404007    0.00000000    0.224421334    0.2433118 0.001972307
# # GSM404008    0.03454889    0.002089211    0.5739203 0.000000000
# 
# re <- TME.results[,-(23:25)]
# 
# # 堆积柱状图
# library(RColorBrewer)
# mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
# Group = str_sub(colnames(exp),1,str_length(colnames(exp)))
# 
# group <- readRDS("../00data/tcga_group.rds")
# dat <- re %>% 
#   as.data.frame() %>%
#   rownames_to_column("Sample") %>% 
#   mutate(group = group$group) %>% 
#   gather(key = Cell_type,value = Proportion,-Sample,-group) %>% 
#   arrange(group)
# 
# dat$Sample = factor(dat$Sample,ordered = T,levels = unique(dat$Sample)) #定横坐标顺序
# 
# write.csv(dat,file = "cell_prop.csv")
# 
# 
# # 先把group排序，然后将sample设为了因子，确定排序后的顺序为水平，所以两图的顺序是对应的。
# dat2 = data.frame(a = 1:ncol(exp),
#                   b = 1,
#                   group = sort(group$group)) 
# 
# p1 = ggplot(dat2,aes(x = a, y = b)) + 
#   geom_tile(aes(fill = group)) + 
#   scale_fill_manual(values = mypalette(22)[1:length(unique(Group))]) +
#   theme(panel.grid = element_blank(), 
#         panel.background = element_blank(), 
#         axis.line = element_blank(), 
#         axis.ticks = element_blank(), 
#         axis.text = element_blank(), 
#         axis.title = element_blank()) + 
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(fill = "Group")
# 
# p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) + 
#   geom_bar(stat = "identity") +
#   labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
#   theme_bw() +
#   theme(#axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   ) + 
#   scale_y_continuous(expand = c(0.01,0)) +
#   scale_fill_manual(values = mypalette(22))
# 
# library(patchwork)
# 
# p <- p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
#   theme(legend.position = "bottom")
# 
# ggsave(filename = '04.histon_plot.pdf',p,w=15,h=8)
# ggsave(filename = '04.histon_plot.png',p,w=15,h=8)
# 
# # #------------------------- 箱线图--------------------------
# # # 全是0的行去掉
# # k = colSums(re)>0;table(k)
# # 
# # ## k
# # ## FALSE  TRUE 
# # ##     1    21
# # 
# # re = re[,k]
# # library(tinyarray)
# # 
# # group <- readRDS("../00data/tcga_group.rds")
# # # group <- t(group)
# # # factor(group[2,])
# # 
# # p <- draw_boxplot(t(re)%>%as.data.frame(),factor(group$group),
# #                   drop = T,
# #                   method = "wilcox.test",
# #                   color = mypalette(length(unique(Group))))+
# #   scale_fill_manual(values = c("#5F9EA0","#FF7F24"))+
# #   labs(x = "Cell Type", y = "Estimated Proportion") 
# # 
# # ggsave(filename = '05.box_plot2.pdf',p,w=10,h=6)
# # ggsave(filename = '05.box_plot2.png',p,w=15,h=8)
# 

# ---------------------------堆积柱状图-------------------------------------
# --------------------细胞相关性------------------------
tiics_result <- read.csv('01.ssgsea_result_cell 1.csv',check.names = F, row.names = 1)
# %>% lc.tableToNum

# tiics_result <- t(tiics_result) %>% as.data.frame()

# diff <- read.csv('./02.cell_diff_wilcox_test.csv')
# tiics_result <- tiics_result[,diff$immune_cell]
diff_Table

diffcell <- rownames(diff_Table)

res3 <- tiics_result[diffcell,]
res3 <- t(res3)%>%as.data.frame()

# #过滤掉表达为0的
# res3 <- res3[,which(colSums(res3) > 0)]
# res3 <- res3[,order(colnames(res3),decreasing = F)]

cor_data <- cor(res3,method="spearman")
corp <- cor_pmat(res3)

write.csv(cor_data,'03.cell_cor_r.csv',quote=F)
write.csv(corp,'03.cell_cor_p.csv',quote=F)
env.cor <- round(cor((res3),method="spearman"), 3)
# env.p <-round(cor_pmat((gene_dat),method = "spearman"),3) 
cor_p <- WGCNA::corPvalueStudent(env.cor,nrow(res3))

pdf("cor.pdf", width = 7, height = 7)
cor.plot<-corrplot(corr =env.cor,p.mat = cor_p,type="upper",
                   col = colorRampPalette(c("blue", "white", "red"))(50),
                   tl.pos="lt",tl.col="black", 
                   insig = "label_sig", sig.level = c(.001,.01, .05),
                   pch.cex=1,pch.col = "black",order = "AOE")
cor.plot<-corrplot(corr = env.cor,type="lower",add=TRUE,method="number",
                   col = colorRampPalette(c("blue", "white", "red"))(50),
                   tl.pos="n",tl.col="black",tl.cex=1.2,
                   diag=FALSE, cl.pos="n",pch.col = "black",
                   number.cex = 0.7,order = "AOE")
dev.off()

#  生物标志物与关键免疫细胞相关性----
exprlog <- gene_expr
group <- group
multi <- read.csv("hub_gene.csv")

gene <- multi$x

expr <- t(exprlog) %>% as.data.frame()
expr <- expr[,gene]
group <- risk[,7,drop = F] %>% rownames_to_column(.,var = "sample")
colnames(group) <- c('sample', 'group')
# library(lance)
tiics_result <- read.csv('01.ssgsea_result_cell 1.csv',check.names = F, row.names = 1)
# tiics_result <- t(tiics_result) %>% as.data.frame()
colnames(tiics_result)
tiics_result <- tiics_result[,group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[rownames(expr),]

diff <- diffcell
tiics_result <- tiics_result[,diff]


tem <- intersect(rownames(expr), rownames(tiics_result))
expr <- expr[tem,]
tiics_result <- tiics_result[tem,]
identical(rownames(expr), rownames(tiics_result))
cor_r <- cor(expr,tiics_result,method = "spearman") 
#cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(dat_exp_diff)))

d <- corr.test(expr,tiics_result,use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"04.correlation_cor.csv")

# 相关性热图带显著性----
data <- read.csv('04.correlation_cor.csv',row.names = 1,check.names = F)
data <- data %>%
  mutate(text = case_when( #设置label，并加入判断，当P值符合特定条件就显示"\n"外加特定数量的*号
    Pvalue <= 0.001 ~ "\n***", #P<0.001就显示回车加三个星号
    between(Pvalue, 0.001, 0.01) ~ "\n**", #P为0.001-0.01 显示回车加两个*号
    between(Pvalue, 0.01, 0.05) ~ "\n*",  #P为0.01-0.05 显示回车加一个星号
    T ~ ""))
# data <- data %>%
#   mutate(text = case_when(  # 一定要 get 到 case_when() 函数奥秘
#     Pvalue > 0 ~ paste(round(Pvalue, 2), "+"), # round() 只保留两位小数
#     Pvalue < 0 ~ paste(round(Pvalue, 2), "-")))
p <- 
  ggplot(data, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + # 这里可以用 windowns 小工具 takecolor 取色，看中哪个文章就吸哪个文章
  # 比如这篇 https://www.nature.com/articles/nmeth.1902 
  geom_text(aes(label = text),col ="black",size = 5) +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(hjust = 1, size = 14, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold")) + #调整y轴文字
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) +   # 修改 legend 内容
  scale_x_discrete(position = "top") #
p
ggsave(file=paste0('correlation_biomarker_','.png'), height = 8, width = 12, p)
ggsave(file=paste0('correlation_biomarker_','.pdf'), height = 8, width = 12, p)

