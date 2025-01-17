options(timeout = Inf)

rm(list = ls())
setwd("D:\\xqm2\\xqm\\xqm\\2024\\9月\\0905\\4 对lasso分析筛选出来的gene进行四个模型构建 残差选择\\实验")

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
install.packages("DALEX")


# --------------------机器学习-------------------------
options(timeout = Inf)

dat <- readRDS("26 gene_exp.rds")
group <- readRDS("group.rds")
rownames(group) <- group$sample
group <- group[,-1] %>% data.frame() 
colnames(group) <- c("group")
# group <- select(group,-1)
lasso_geneids <- read.csv("lasso_genes.csv")
candidate_gene <- lasso_geneids$symbol


train_mat <- as.data.frame(t(dat[candidate_gene,]))
train_mat$group <- group$group
train_mat$group<-ifelse(train_mat$group=="DKD",1,0) 

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



