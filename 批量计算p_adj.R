library(plyr)
library(dplyr)
library(ggplot2)


df <- read.csv('1.csv')
result_list <- list()


for(i in colnames(df)){

    if( i == 'sample'){
        next()
    }
    data = TukeyHSD(aov(as.formula(paste(i, " ~ sample")), data = df))$sample ## as.formula(paste(i, " ~ sample") 可以调用外部变量

    filename = paste(i,'_padj.csv',sep = "")

    result_list[[filename]] <- data

    write.csv(data,file = filename)
}