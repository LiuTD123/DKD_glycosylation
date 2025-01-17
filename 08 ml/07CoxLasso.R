options(timeout = Inf)

rm(list = ls())
setwd("C:\\Users\\Administrator\\Desktop\\workdir\\01_Z007-L001\\07CoxLasso")

library(ggthemes)
library(survivalROC)
library(survminer)
library(grid)
library(ggvenn)
library(pROC)
library(xgboost)
library(randomForest)
library(caret)
library(DALEX)
library(stats)
library(e1071)
library(glmnet)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(scales)
library(survminer)
library(glmnet)
library(readxl)
library(readr)
library(tidyverse)
library(plyr)
library(glmnet)
library(survminer)
library(gridExtra)
library(patchwork)
library(forestplot)
library(survival)
library(dplyr)
library(plyr)
library(tibble)
library(tidyr)

fpkm <- readRDS("../00data/tcga_fpkm.rds")
group <- readRDS("../00data/tcga_group.rds")
# keygene <- read.csv("/data/nas1/liuhouyan/01-BJTC-577-2/04_candidate_gene/DEGs-PRGs.csv")
# load("/data/nas1/liuhouyan/01-BJTC-577-2/15_wgcna/cluster_and_expression.RData")
keygene <- read.csv("../03venn/DEGs_glygenes.csv")

# ppi_top20 <- read.csv("/data/nas1/liuhouyan/01-BJTC-577-2/04_candidate_gene/ppi/10.Venn_Con.csv")

survival_dat <- readRDS("../00data/survival_dat.rds")
colnames(survival_dat) <- c("id","OS.time","OS")
# survival_dat$id <- substr(survival_dat$id, 0, 15)

fpkm_1 <- fpkm[keygene$x,] %>%t %>% as.data.frame() %>% rownames_to_column(var = "id")
# fpkm_1 <- fpkm[ppi_top20$x,] %>%t %>% as.data.frame() %>% rownames_to_column(var = "id")
fpkm_1$id <- gsub('\\.','\\-',fpkm_1$id)

fpkm_tem <- merge(fpkm_1,survival_dat,by = "id")

rownames(fpkm_tem) <- fpkm_tem$id
fpkm_tem <- fpkm_tem[,-1]
# fpkm_tem <- fpkm_tem[,-11]
df_merge <- fpkm_tem

# unicox ------------------------------------------------------------------

diff_expr_clinical <- df_merge
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum

covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))
univ_models <- lapply(univ_formulas,
                      function(x) coxph(x, data = diff_expr_clinical))

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
write.table(res, file = "03.uni_cox_OS_res.txt", sep = "\t", quote = F, row.names = T)

# 选择单因素P值 0.05
res_results<- na.omit(res[which(as.numeric(res$p.value) < 0.05),])
write.table(res_results, file = "04.uni_cox_OS_sig0.05.txt", sep = "\t", quote = F, row.names = T)
dim(res_results)
# [1] 14 2

res_results_all <- res_results
res_results1 <- res_results_all[! rownames(res_results_all)%in%c("CHST8","HS6ST3","EXTL1"),]
res_results2 <- res_results_all[c("CHST8","HS6ST3","EXTL1"),]
# 单因素cox森林图----
res_results<- res_results_all[order(res_results_all$p.value),]
# res_results <- res_results[-5,]
# res_results <- res_results[-3,]
# res_results <- res_results[-9,]
res_results_plot <- res_results
unix_res <- tidyr::separate(res_results_plot, "HR (95% CI for HR)", into = c("HR", "HR.95L", "HR.95H"), sep = " ") %>% 
  tidyr::separate("HR.95L", into = c("HR.95L", "HR.95H"), sep = "\\-") %>% 
  tibble::rownames_to_column(var="GeneName")
unix_res$HR.95L <- gsub("\\(", "", unix_res$HR.95L)
unix_res$HR.95H <- gsub("\\)", "", unix_res$HR.95H)
unix_res[, 2:ncol(unix_res)] <- as.data.frame(apply(unix_res[, 2:ncol(unix_res)], 2, as.numeric))
unix_res <- unix_res[order(unix_res$p.value),]

hz <- paste(round(unix_res$HR,4),
            "(",round(unix_res$HR.95L,3),
            "-",round(unix_res$HR.95H,3),")",sep = "")
tabletext <- cbind(c(NA,"GeneName", unix_res$GeneName),
                   c(NA,"P value", ifelse(unix_res$p.value<0.001,
                                          "< 0.001",
                                          round(unix_res$p.value,4))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
# 17

pdf(file = "13.univariate_cox_forest.pdf", height = 7, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE, rep(FALSE, 54)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,unix_res$HR),
           lower=c(NA,NA,unix_res$HR.95L), #95%置信区间下限
           upper=c(NA,NA,unix_res$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=1.2, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface = "bold", fontfamily="Times"),
                          summary = gpar(cex=1.2, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.2, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily="Times")),
           xlab="Hazard Ratio",
           cex.lab=2,
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()

png(filename = "13.univariate_cox_forest.png", height = 7, width = 10, units = "in", res = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE, rep(FALSE, 54)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,unix_res$HR),
           lower=c(NA,NA,unix_res$HR.95L), #95%置信区间下限
           upper=c(NA,NA,unix_res$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=1.2, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface = "bold", fontfamily="Times"),
                          summary = gpar(cex=1.2, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.2, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily="Times")),
           xlab="Hazard Ratio",
           cex.lab=2,
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()

save.image("unicox.Rdata")
# ---------------------------------------------------------------------
# # PH检验

for (gene in rownames(res_results)){
  # gene <- cox.zph(univ_models$gene)
  # print(gene)
  print(paste0(gene,"<-","cox.zph(univ_models$",gene,")"))
  print(paste0(gene, "<- ggcoxzph(",gene, ",caption = '",gene,"')"))
  print(paste0("ggsave(","'./ph/ggcoxzph_",gene,".pdf'",",arrangeGrob(grobs = ",gene,"))"))
  print(paste0("ggsave(","'./ph/ggcoxzph_",gene,".png'",",arrangeGrob(grobs = ",gene,"))"))
}

HS6ST3<-cox.zph(univ_models$HS6ST3)
HS6ST3<- ggcoxzph(HS6ST3,caption = 'HS6ST3')
ggsave('./ph/ggcoxzph_HS6ST3.pdf',arrangeGrob(grobs = HS6ST3))
ggsave('./ph/ggcoxzph_HS6ST3.png',arrangeGrob(grobs = HS6ST3))
MFNG<-cox.zph(univ_models$MFNG)
MFNG<- ggcoxzph(MFNG,caption = 'MFNG')
ggsave('./ph/ggcoxzph_MFNG.pdf',arrangeGrob(grobs = MFNG))
ggsave('./ph/ggcoxzph_MFNG.png',arrangeGrob(grobs = MFNG))
SLC35D1<-cox.zph(univ_models$SLC35D1)
SLC35D1<- ggcoxzph(SLC35D1,caption = 'SLC35D1')
ggsave('./ph/ggcoxzph_SLC35D1.pdf',arrangeGrob(grobs = SLC35D1))
ggsave('./ph/ggcoxzph_SLC35D1.png',arrangeGrob(grobs = SLC35D1))
GALNT7<-cox.zph(univ_models$GALNT7)
GALNT7<- ggcoxzph(GALNT7,caption = 'GALNT7')
ggsave('./ph/ggcoxzph_GALNT7.pdf',arrangeGrob(grobs = GALNT7))
ggsave('./ph/ggcoxzph_GALNT7.png',arrangeGrob(grobs = GALNT7))
HS2ST1<-cox.zph(univ_models$HS2ST1)
HS2ST1<- ggcoxzph(HS2ST1,caption = 'HS2ST1')
ggsave('./ph/ggcoxzph_HS2ST1.pdf',arrangeGrob(grobs = HS2ST1))
ggsave('./ph/ggcoxzph_HS2ST1.png',arrangeGrob(grobs = HS2ST1))
B4GALT6<-cox.zph(univ_models$B4GALT6)
B4GALT6<- ggcoxzph(B4GALT6,caption = 'B4GALT6')
ggsave('./ph/ggcoxzph_B4GALT6.pdf',arrangeGrob(grobs = B4GALT6))
ggsave('./ph/ggcoxzph_B4GALT6.png',arrangeGrob(grobs = B4GALT6))
UST<-cox.zph(univ_models$UST)
UST<- ggcoxzph(UST,caption = 'UST')
ggsave('./ph/ggcoxzph_UST.pdf',arrangeGrob(grobs = UST))
ggsave('./ph/ggcoxzph_UST.png',arrangeGrob(grobs = UST))
EXTL1<-cox.zph(univ_models$EXTL1)
EXTL1<- ggcoxzph(EXTL1,caption = 'EXTL1')
ggsave('./ph/ggcoxzph_EXTL1.pdf',arrangeGrob(grobs = EXTL1))
ggsave('./ph/ggcoxzph_EXTL1.png',arrangeGrob(grobs = EXTL1))
CHST1<-cox.zph(univ_models$CHST1)
CHST1<- ggcoxzph(CHST1,caption = 'CHST1')
ggsave('./ph/ggcoxzph_CHST1.pdf',arrangeGrob(grobs = CHST1))
ggsave('./ph/ggcoxzph_CHST1.png',arrangeGrob(grobs = CHST1))
CHPF<-cox.zph(univ_models$CHPF)
CHPF<- ggcoxzph(CHPF,caption = 'CHPF')
ggsave('./ph/ggcoxzph_CHPF.pdf',arrangeGrob(grobs = CHPF))
ggsave('./ph/ggcoxzph_CHPF.png',arrangeGrob(grobs = CHPF))
DPM2<-cox.zph(univ_models$DPM2)
DPM2<- ggcoxzph(DPM2,caption = 'DPM2')
ggsave('./ph/ggcoxzph_DPM2.pdf',arrangeGrob(grobs = DPM2))
ggsave('./ph/ggcoxzph_DPM2.png',arrangeGrob(grobs = DPM2))
ALG3<-cox.zph(univ_models$ALG3)
ALG3<- ggcoxzph(ALG3,caption = 'ALG3')
ggsave('./ph/ggcoxzph_ALG3.pdf',arrangeGrob(grobs = ALG3))
ggsave('./ph/ggcoxzph_ALG3.png',arrangeGrob(grobs = ALG3))
B3GNT6<-cox.zph(univ_models$B3GNT6)
B3GNT6<- ggcoxzph(B3GNT6,caption = 'B3GNT6')
ggsave('./ph/ggcoxzph_B3GNT6.pdf',arrangeGrob(grobs = B3GNT6))
ggsave('./ph/ggcoxzph_B3GNT6.png',arrangeGrob(grobs = B3GNT6))
CHST8<-cox.zph(univ_models$CHST8)
CHST8<- ggcoxzph(CHST8,caption = 'CHST8')
ggsave('./ph/ggcoxzph_CHST8.pdf',arrangeGrob(grobs = CHST8))
ggsave('./ph/ggcoxzph_CHST8.png',arrangeGrob(grobs = CHST8))

a <- paste(rownames(res_results),collapse = ",")

ggsave("ggcoxzph_all.pdf",w = 17,h = 15,
       arrangeGrob(grobs = c(HS6ST3,MFNG,SLC35D1,GALNT7,HS2ST1,B4GALT6,UST,EXTL1,CHST1,CHPF,DPM2,ALG3,B3GNT6,CHST8)))

# FOS<-cox.zph(univ_models$FOS)
# FOS<- ggcoxzph(FOS,title = "FOS")
# ggsave("ggcoxzph_FOS.pdf",arrangeGrob(grobs = FOS))
# ggsave("ggcoxzph_FOS.png",arrangeGrob(grobs = FOS))
# 
# STK40<-cox.zph(univ_models$STK40)
# STK40<- ggcoxzph(STK40,title = "STK40")
# ggsave("ggcoxzph_STK40.pdf",arrangeGrob(grobs = STK40))
# ggsave("ggcoxzph_STK40.png",arrangeGrob(grobs = STK40))

# ----------------------------
# PH检验
x <- df_merge[,c(112,111)]
uniSigExp = df_merge
uniSigExp = uniSigExp[,rownames(res_results)]
uniSigExp <- merge(uniSigExp,x, by = "row.names")
rownames(uniSigExp) <- uniSigExp$Row.names
uniSigExp <- uniSigExp[,-1]
dat <- uniSigExp
outPH=data.frame()
for(i in colnames(dat[,1:16])){
  cox <- coxph(Surv(OS.time, OS) ~ dat[,i], data = dat)
  test.ph <- cox.zph(cox)
  #coxP=test.ph$coefficients[,"Pr(>|z|)"]
  outPH=rbind(outPH,
              cbind(id=i,
                    p=test.ph$table[1,"p"])
  )
}

sigPH=outPH[as.numeric(as.vector(outPH$p))>0.05,]  # 46
write.table(sigPH,file="02.PH.Sig.txt",sep="\t",row.names=F,quote=F)
# -----------------------------------------

# Lasso -------------------------------------------------------------------

# diff_expr_clinical <- read.table("./02.TCGA-BLCA_logtpm_cli.xls", header = T)
# res_results_0.05 <- read.table("./04.uni_cox_OS_sig0.05.txt", header = T, sep = "\t", check.names = F)
x_all <- subset(diff_expr_clinical, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results)]
y_all <- subset(diff_expr_clinical, select = c(OS, OS.time))

fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS),
              family = "cox")
plot(fit, xvar = "lambda",label = TRUE, las=1)

png(filename = "05.lasso_model.png", height = 450, width = 600)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "05.lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()

set.seed(1)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),
                  nfold=10,
                  family = "cox")
plot(cvfit, las =1)

png(filename = "06.lasso_verify.png", height = 450, width = 600)
plot(cvfit, las =1)
dev.off()
pdf(file = "06.lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()

# 1.Lasso模型   ----------------
# library(lance）

fpkm <- readRDS("..\\00data\\tcga_fpkm.rds")
group <- readRDS("../00data/tcga_group.rds")
## 01.获取数据集 -----------------------------------------------------------  
dat <- fpkm
# dat <- log2(dat+1)
keygene <- read.csv("../03venn/DEGs_glygenes.csv")

survival <- readRDS("../00data/survival_dat.rds")

gene <- res_results


## 02.合并生存数据
survival_dat <- t(dat[rownames(gene),colnames(dat) %in% survival$sample])
train_dat <- survival_dat %>% data.frame()
train_dat$sample <- rownames(train_dat)
train_dat <- merge(survival,train_dat,by='sample')
rownames(train_dat) <- train_dat$sample
train_dat<-train_dat[,c(-1)]
colnames(train_dat)


### 03.LASSO

train_data <- train_dat
x_all <- subset(train_data, select = -c(fustat, futime))
y_all <- subset(train_data, select = c(fustat, futime))


##  04.拟合模型 ----------------------------------------------------------------
fit <- glmnet(as.matrix(x_all), Surv(y_all$futime,y_all$fustat), 
              family = "cox") 


## 05.交叉验证 -----------------------------------------------------------------
set.seed(1)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$futime,y_all$fustat),nfold=10,
                  family = "cox") 
##画图

coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1]  0.02011


# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
length(coef.min)
# 14
coef.min
# HS6ST3   1.26063722
# MFNG     0.24022013
# SLC35D1  .         
# GALNT7   .         
# HS2ST1  -0.12715515
# B4GALT6 -0.11150171
# UST      0.15696477
# EXTL1    .         
# CHST1    .         
# CHPF     0.02560628
# DPM2     0.13624936
# ALG3     0.13105151
# B3GNT6  -0.10084413
# CHST8    0.52714065

# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids

saveRDS(lasso_geneids, file= "lasso_geneids.rds")

saveRDS(coef.min, file="coef.min.rds")


lasso_geneids <- data.frame(symbol=lasso_geneids)
write.csv(lasso_geneids, "01.lasso_genes.csv",quote = F,row.names = F)

x <- coef(fit) 
tmp <- as.data.frame(as.matrix(x)) 
tmp$coef <- row.names(tmp) 
tmp <- reshape::melt(tmp, id = "coef") 
tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
tmp$coef <- gsub('_','-',tmp$coef) 
tmp$lambda <- fit$lambda[tmp$variable+1] 
# extract the lambda values 
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] 
# compute L1 norm

#图片美化
head(tmp)

pdf("02.lasso_model_name.pdf",height = 5, width = 7,family='Times')
ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Log(Lambda)") + 
  ylab('Coefficients')+ 
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values= c('#ff9898','#dedb8e','#99e0ab','#D94F04','#007172',"#E9967A","#FA8072",
                               '#025259','#c49d93','#aec6e8','#F2C6C2','#86A69D',"black",  "#FFA07A"))+

  scale_x_continuous(expand = c(0.01,0.01))+ 
  scale_y_continuous(expand = c(0.01,0.01))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=15,color='black'), 
        axis.text = element_text(size=12,color='black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size=12,color='black'), 
        legend.position = 'right')+ 
  annotate('text',x = -4.3,y=0.18,
           label=paste0('Optimal Lambda =', cvfit$lambda.min),
           color='black')+ 
  guides(col=guide_legend(ncol = 1))
dev.off()

#--------------------------- 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = cvfit$lambda.min)  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1] 0.02011
df.coef = cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame
write.table(df.coef, "Lasso_Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)  #6个

save.image("lasso_over.Rdata")
load("lasso_over.Rdata")



# --------------------机器学习-------------------------
options(timeout = Inf)

dat <- readRDS("..\\00data\\tcga_count.rds")
group <- readRDS("..\\00data\\tcga_group.rds")
rownames(group) <- group$sample
group <- group[,-1] %>% data.frame() 
colnames(group) <- c("group")
# group <- select(group,-1)
candidate_gene <- lasso_geneids$symbol


train_mat <- as.data.frame(t(dat[candidate_gene,]))
train_mat$group <- group$group
train_mat$group<-ifelse(train_mat$group=="Tumor",1,0) 

# train.dat[,1]
# rownames(train.dat)<-train.dat[,1]
# 
# train.dat2 <- train.dat[,-1]
# dat<-train.dat2[DE_migenes$x,group$sample]%>%t%>%as.data.frame()
# 
# dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
# 
# table(group$sample)
# dat$group<-factor(dat$group,levels = c('PD','Control'))

input<-train_mat

# ------------------------------------------------


# ----------------------------------------------------------------

ctrl <- trainControl(method="cv", number=10)
# classif_knn <- train(group~., data = input,
# method = "knn", trControl=ctrl, tuneLength=10)


classif_pls <- train(group~., data = input,
                     method = "pls",validation = 'CV')

# 支持向量机模型
classif_svm <- train(group~., data = input, 
                     method = "svmRadial")

# 随机森林模型
classif_rf <- train(group~., data = input, 
                    method = "rf",
                    ntree = 20)
# 广义线性模型
classif_glm <- train(group~., data = input,
                     method = 'glm')


# # 极限梯度提升模型
# x = model.matrix(group~.,input)[,-310]
# model_martix_train<-model.matrix(group~., input)[,-310]
# data_train <- xgb.DMatrix(x , label =as.numeric(input$group))
# 
# params <- list(
#   objective = "reg:squarederror"
# )
# 
# classif_xgboost <- xgb.train(params, data_train, nrounds = 100)
# 
# save(classif_xgboost,classif_glm,classif_rf,classif_svm,file="fourModel_classif.RData")

# explainer_knn<-explain(classif_knn,label = "knn",
#                        data = input,
#                        y = input$group)

explainer_pls<-explain(classif_pls,label = "PLS",
                       data = input,
                       y = input$group)

explainer_svm<-explain(classif_svm,label = "SVM",
                       data = input,
                       y = input$group)
explainer_rf<-explain(classif_rf,label = "RF",
                      data = input,
                      y = input$group)
explainer_glm<-explain(classif_glm,label = "GLM",
                       data = input,
                       y = input$group)

# -------------------------------xgboost
# predict_logit <- function(model,x){
#   raw_x <-predict(model,x)
#   exp(raw_x)/(1+exp(raw_x))
# }
# 
# logit <- function(x){
#   exp(x)/(1+exp(x))
# }
# 
# explainer_xgboost<-explain(classif_xgboost,
#                            label = "xgboost",
#                        data = x,
#                        y = as.numeric(input$group),
#                        predict_function = predict_logit,
#                        link = logit
#                        )

# save(explainer_xgboost,explainer_glm,explainer_rf,explainer_svm,file="fourModel_explainer.RData")

# model performance
# per_knn<-model_performance(explainer_knn)

per_pls<-model_performance(explainer_pls, measure = "auc")

per_svm<-model_performance(explainer_svm, measure = "auc")
per_rf<-model_performance(explainer_rf, measure = "auc")
per_glm<-model_performance(explainer_glm, measure = "auc")
# per_xgboost<-model_performance(explainer_xgboost)

auc_pls <- auc(per_pls$residuals$observed,per_pls$residuals$predicted)
auc_glm <- auc(per_glm$residuals$observed,per_glm$residuals$predicted)
auc_rf <- auc(per_rf$residuals$observed,per_rf$residuals$predicted)
auc_svm<- auc(per_svm$residuals$observed,per_svm$residuals$predicted)

# ROC曲线
pdf("01.svm_ROC.pdf",height = 4,width = 6)
plot(per_svm, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_svm, 5)))
dev.off()

# pdf("01.knn_ROC.pdf",height = 4,width = 6)
# plot(per_knn, geom = "roc")
# dev.off()

pdf("01.rf_ROC.pdf",height = 4,width = 6)
plot(per_rf, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_rf, 5)))
dev.off()

pdf("01.glm_ROC.pdf",height = 4,width = 6)
plot(per_glm, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_glm, 5)))
dev.off()

pdf("01.pls_ROC.pdf",height = 4,width = 6)
plot(per_pls, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_pls, 5)))
dev.off()
# pdf("01.xgboost_ROC.pdf",height = 4,width = 6)
# plot(per_xgboost, geom = "roc")
# dev.off()
# ------------------------------------------------------

# plot model performance
pdf("04.machine_learning_residuals_line.pdf",height = 4,width = 6)
plot(per_svm,per_rf,per_glm,per_pls)
dev.off()

pdf("04.machine_learning_residuals_box.pdf",height = 4,width = 6)
plot(per_svm,per_rf,per_glm,per_pls,geom = "boxplot")
dev.off()

# importance
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_pls<-variable_importance(
  explainer_pls,
  loss_function = loss_root_mean_square
)
pdf("04.machine_learning_importance_plot2.pdf",height = 20,width = 8)
plot(importance_svm,importance_rf,importance_glm,importance_pls)
dev.off()

# # 筛选四个模型中RSME值前50%的基因的交集，确定为特征基因----
# quantile(importance_glm$dropout_loss, 0.5)
glm <- importance_glm[importance_glm$dropout_loss < quantile(importance_glm$dropout_loss, 0.6),]
pls <- importance_pls[importance_pls$dropout_loss < quantile(importance_pls$dropout_loss, 0.6),]
svm <- importance_svm[importance_svm$dropout_loss < quantile(importance_svm$dropout_loss, 0.6),]
rf <- importance_rf[importance_rf$dropout_loss < quantile(importance_rf$dropout_loss, 0.6),]

# xgboost <- importance_xgboost[importance_xgboost$dropout_loss < 0.345,]

# knn_gene <- unique(knn$variable) %>% as.data.frame()
# knn_gene <- subset(knn_gene,knn_gene$. != "_full_model_" & knn_gene$. != "group")

glm_gene <- unique(glm$variable) %>% as.data.frame()
glm_gene <- subset(glm_gene,glm_gene$. != "_full_model_" & glm_gene$. != "group")
# 
pls_gene <- unique(pls$variable) %>% as.data.frame()
pls_gene <- subset(pls_gene,pls_gene$. != "_full_model_" & pls_gene$. != "group")
# 
svm_gene <- unique(svm$variable) %>% as.data.frame()
svm_gene <-subset(svm_gene,svm_gene$. != "_full_model_" & svm_gene$. != "group")

rf_gene <- unique(rf$variable) %>% as.data.frame()
rf_gene <- subset(rf_gene,rf_gene$. != "_full_model_" & rf_gene$. != "group")

# rf_gene <- subset(rf_gene,rf_gene$. != "_full_model_" & rf_gene$. != "group")

# xgboost_gene <- unique(xgboost$variable) %>% as.data.frame()
# xgboost_gene <- subset(xgboost_gene,xgboost_gene$. != "_full_model_" & xgboost_gene$. != "group")

# 交集基因
hub_gene <- Reduce(intersect, list(glm_gene$.,rf_gene$.,pls_gene$.,svm_gene$.)) # 5个
write.csv(hub_gene,'hub_gene.csv',row.names = T)

# 绘图

mydata<-list('GLM'= glm_gene$.,"PLS"= pls_gene$.,"SVM" = svm_gene$., "RF" = rf_gene$.)
pdf('02.venn.pdf',w=7,h=7)
ggvenn(mydata,
       c('GLM',"PLS","SVM","RF"),
       fill_color = c("#ffb2b2","#b2e7cb","#7570B3","orange"),
       show_percentage = T,
       fill_alpha = 0.5,
       stroke_alpha = 1,
       stroke_size = 0.4,
       text_size = 5,
       stroke_color=NA,
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#7570B3","orange"),
       set_name_size = 8,
       text_color = 'black'
)
dev.off()

save.image("ml.Rdata")


# -----------------------训练集验证------------
load("ml.Rdata")

df_merge <- df_merge[,c(hub_gene,"OS.time","OS")]
diff_expr_clinical <- df_merge
# riskScore <- predict(step_Cox,type="lp",newdata=diff_expr_clinical)
riskScore <- predict(cvfit, newx = as.matrix(x_all), s=cvfit$lambda.min)
riskScore <- data.frame(riskScore)

identical(rownames(riskScore), rownames(diff_expr_clinical))
diff_expr_clinical$riskScore <- riskScore$X1

# 中位数
# diff_expr_clinical$risk <- ifelse(diff_expr_clinical$riskScore>median(diff_expr_clinical$riskScore), "High", "Low")
# 最佳截断值
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical, #数据集
                                   minprop = 0.25,
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "riskScore"  #需要计算的数据列名
)
OS.cut
#           cutpoint statistic
# riskScore 1.240599  4.960954

res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
diff_expr_clinical$risk <- ifelse(diff_expr_clinical[,"riskScore"] > res["riskScore",1],'High','Low' )

write.table(diff_expr_clinical, file = "10.all_train_riskScore.txt", sep = "\t", row.names = T, col.names = T, quote = F)

diff_expr_clinical <- read.table("10.all_train_riskScore.txt")
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = diff_expr_clinical)

customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}
train_km <- ggsurvplot(kmfit,
                       pval = TRUE, 
                       pval.method = T,
                       conf.int = F,
                       legend.labs=c("High risk","Low risk" ),
                       legend.title="Risk Score",
                       xlab = "Overall Survival",
                       title="Train KM",
                       font.main = c(15,"bold"),
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       linetype = "strata", 
                       # surv.median.line = "hv",
                       font.family = "Times",
                       risk.table.y.text.col = T,
                       risk.table.y.text = T,
                       risk.table.height = 0.35,
                       ggtheme = theme_bw(), 
                       palette = c("#A73030FF", "#0073C2FF"))
train_km
train_km$plot <- train_km$plot + labs(
  title    = "Survival Curves",
  subtitle = "Based on Kaplan-Meier Estimates"
)
train_km$table <- train_km$table + labs(
  caption  = "Created with Train Data"
)
train_km <- customize_labels(
  train_km,
  font.title    = c(16, "bold"),
  font.subtitle = c(15, "bold.italic"),
  font.caption  = c(14, "plain", "orange"),
  font.x        = c(14, "bold.italic"),
  font.y        = c(14, "bold.italic"),
  font.xtickslab = c(12, "plain")
)

train_km
pvalue <- stringr::str_extract(train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]],
                               "\\d.*")
train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(pvalue)))
train_km
pdf(file = "09.train_km.pdf", family = "Times", height = 6, width = 6, onefile = F)
print(train_km)
dev.off()
png(filename = "09.train_km.png", family = "Times", height = 6, width = 6, units = "in", res = 600)
print(train_km)
dev.off()


# ROC----

rt = subset(diff_expr_clinical, select = c(OS, OS.time, riskScore))
# rt$OS.time <- rt$OS.time / 365
colnames(rt) <- c("fustat","futime","riskscore")
ROC <- rt
cutoff_1 <- 1
cutoff_2 <- 3
cutoff_3 <- 5
year_1= survivalROC(Stime=ROC$futime,##生存时间
                    status=ROC$fustat,## 终止事件
                    marker = ROC$riskscore, ## marker value
                    predict.time = cutoff_1,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_2= survivalROC(Stime=ROC$futime,##生存时间
                    status=ROC$fustat,## 终止事件
                    marker = ROC$riskscore, ## marker value
                    predict.time = cutoff_2,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_3= survivalROC(Stime=ROC$futime,##生存时间
                    status=ROC$fustat,## 终止事件
                    marker = ROC$riskscore, ## marker value
                    predict.time = cutoff_3,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
if(T){
  pdf(file = paste0("01.ModelROC_train.pdf"),width = 8,height = 8)
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("01.ModelROC_train.png"),width= 8, height= 8, units="in", res=300)
  dev.control("enable")
  par(mar = c(5,5,5,2))
  plot(year_1$FP, year_1$TP,
       type="l",col="red",xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Fraction",
       ylab="True Positive Fraction",
       main="COAD-train, Method=KM\n Year = 3,5,7",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5
  )
  abline(0,1,col="gray",lty=2)
  lines(year_2$FP, year_2$TP, type="l",col="#EB4B17",xlim=c(0,1), ylim=c(0,1))
  lines(year_1$FP, year_1$TP, type="l",col="#2775AB",xlim=c(0,1), ylim=c(0,1))
  lines(year_3$FP, year_3$TP, type="l",col="#4C8045",xlim=c(0,1), ylim=c(0,1))
  legend(0.6,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,3)),
                   paste("AUC of 3 year =",round(year_2$AUC,3)),
                   paste("AUC of 5 year =",round(year_3$AUC,3))),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,col=c("#2775AB","#EB4B17",'#4C8045'),
         bty = "n",# bty框的类型
         seg.len=1,cex=1.2)#
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
}

# riskScore dis----

riskScore <- riskScore[,1]
# median(riskScore)
# [1] -0.5101641
risk_dis <- ggplot(diff_expr_clinical, aes(x=reorder(rownames(diff_expr_clinical), riskScore), y=riskScore, color = risk)) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = rownames(diff_expr_clinical)[order(diff_expr_clinical$riskScore)][c(1,50,100,150,200,250,300,350,400)],
                   labels = c(1,50,100,150,200,250,300,350,400),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(diff_expr_clinical[which(diff_expr_clinical$risk=="Low"),]) + 0.5, lty = 2) +
  geom_hline(yintercept = res$cutpoint, lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train Dataset Risk Score Distribution") + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        # legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15),
        text = element_text(family = "Times", face = "bold"))
risk_dis
ggsave(filename = "14.riskScore_dis.pdf", height = 3, width = 5, risk_dis)
ggsave(filename = "14.riskScore_dis.png", height = 3, width = 5, risk_dis)

surv_stat <- ggplot(diff_expr_clinical, aes(x=reorder(rownames(diff_expr_clinical), riskScore),
                                            y=OS.time,
                                            color = factor(OS,
                                                           levels = c(0,1),
                                                           labels = c("Alive", 
                                                                      "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = rownames(diff_expr_clinical)[order(diff_expr_clinical$riskScore)][c(1,50,100,150,200,250,300,350,400)],
                   labels = c(1,50,100,150,200,250,300,350,400),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(diff_expr_clinical[which(diff_expr_clinical$risk=="Low"),]) + 0.5, lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Overall Survival (days)",
       title = "Train Dataset Overall Survival (days) Distribution") + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        # legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15),
        text = element_text(family = "Times", face = "bold"))
surv_stat
ggsave(filename = "15.OS_dis.pdf", height = 5, width = 7, surv_stat)
ggsave(filename = "15.OS_dis.png", height = 5, width = 7, surv_stat)

# 
# # # 预后基因的表达热图----
# train_count <- df_merge
# 
# dim(train_count)# 166 229
# multi <- read.table("Lasso_Coefficients.xls", header = T) %>% as.data.frame()#cox_result_step2.rds;mul_cox_result.rds
# 
# gene_list <-  multi$gene
# model_gene <- unique(gene_list) #目的基因
# length(model_gene) # 3
# exprSet <- train_count
# # exprSet=log(edgeR::cpm(count)+1)
# # exprSet <- edgeR::cpm(count)
# mat <- exprSet
# mat <- t(scale(t(mat)))#对原始count进行标准化
# mat[mat > 2] <- 2
# mat[mat < -2] <- -2
# 
# 
# diff_expr_clinical <- read.table("10.all_train_riskScore.txt")
# # risk <- diff_expr_clinical[,c(1888,1889,1890,1891)]
# risk <- diff_expr_clinical[,c(9,8,7,6)]
# 
# group <- rownames_to_column(risk,var = "sample")
# group <- group[order(group$risk),]
# group$risk <- factor(group$risk,levels = c("Low","High"))
# annotation_col <- group[,2,drop = F]
# 
# group_color <- c('Low' = "#0073C2FF",  'High' = "#A73030B2")
# names(group_color) <- levels(group$risk)
# ann_colors <- list(
#   group = group_color)
# rownames(annotation_col) <- group$sample
# mat<- t(mat)
# mat <- mat[model_gene,group$sample]
# newnames <- lapply(rownames(mat),function(x) bquote(italic(.(x))))
# 
# library(pheatmap)
# pheatplot <- pheatmap(mat=mat,
#                       annotation_col = annotation_col,
#                       annotation_colors = ann_colors,
#                       annotation_legend = T,#
#                       color = colorRampPalette(c("#0073C2FF", "white", "firebrick3"))(50),
#                       fontsize = 8,
#                       fontsize_row = 10,
#                       # cellheight=20, cellwidth=1.1,
#                       fontfamily = "Times",
#                       labels_row = as.expression(newnames),
#                       show_colnames = FALSE,
#                       cluster_cols = F,
#                       cluster_rows = T,
#                       annotation_names_row = F,
#                       treeheight_row = 15
# )
# 
# # dev.off()
# # ggsave(width=8,height=3,'05.heatmap.pdf',pheatplot)
# # ggsave(width=8,height=3,'05.heatmap.png',pheatplot)
# pheatplot
# png(filename = "train_heatmap.png", height = 3, width = 8, units = 'in',res = 600,family='Times')
# pheatplot
# dev.off()
# pdf(file = "train_heatmap.pdf", height = 3,width = 8,family='Times')
# pheatplot
# dev.off()
# 
# 
