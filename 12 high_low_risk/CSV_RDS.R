library(readr)
# 假设你的CSV文件名为"data.csv"，位于工作目录中
data <- read.csv("all gene_exp DKD group.csv")
# 将数据保存为RDS文件，文件名为"data.rds"
saveRDS(data, "all gene_exp DKD group.rds")
