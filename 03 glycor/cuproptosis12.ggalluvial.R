######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("ggplot2")
#install.packages("ggalluvial")


#???ð?
library(dplyr)
library(ggalluvial)
library(ggplot2)

inputFile="corResult.txt"         #???????Ľ????ļ?
setwd("C:\\Users\\86150\\Desktop\\Fer+m6A+LIHC\\12.ggalluvial")     #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#??????ɫ
mycol=rainbow(length(unique(rt[,"Cuproptosis"])), s=0.8, v=0.8)
#mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
#????ɣ??ͼ
p1<-ggplot(data=rt, aes(axis1 =lncRNA , axis2 =Cuproptosis, y = 1))+
  geom_alluvium(aes(fill =Cuproptosis), width = 0.1, knot.pos = 0.1, reverse = F)+ 
  geom_stratum(fill=NA, color=NA, alpha= 0.5, width = 0.1)+
  geom_text(stat = 'stratum', size =1.5, color='black', label.strata = T)+
  scale_fill_manual(values = mycol) +
  scale_x_discrete(limits = c('lncRNA','mFRG'), expand=c(0, 0))+
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(),axis.text.x = element_blank()) + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  coord_flip()+ggtitle("")   

#????ͼ??
pdf(file="ggalluvial.pdf", width=18, height=6)
print(p1)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

