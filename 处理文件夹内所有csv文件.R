library(plyr)
library(dplyr)
library(ggplot2)
library(ggsignif)



a <-  system('ls',inter = TRUE)
result_list <- list()

for (i in a) {
  if(endsWith(i,'.csv')){
  data <- read.csv(i)
  names(data) = c('sample','type','mean','se')

  ## 将sample，type 转为因子，并按照给定的顺序排列  使ggplot2绘图时按其排序
  data$sample <- factor(data$sample, level = c("开花前", "初花期", "盛花期", "幼果期"))
  data$type <- factor(data$type, level = c('叶', '花', '果'))
  
  ## 使用 paste0 来拼接文件名 paste0 不需要设定sep
  filename <- gsub('.csv','',i)
  
  ## 将结果保存到列表中
  result_list[[filename]] <- data

  fig_color = c('#8CCDBF','#CDE0A5','#F9DB95','#EF8476')


  figname = gsub('.csv','.png',i)

  title_name = gsub('.csv','',i)
  ## ylab的名称
  yname = '相对表达量'

## 绘图
  p = ggplot(data, aes(x = type, y = mean, fill = sample)) +
      geom_bar(stat = "identity", position = 'dodge', width = 0.5) + # 簇状图，position = 'dodge'
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.5), width = 0.25) + # position = position_dodge(0.5) 使误差线与bar对应
      theme_classic() +
      labs(title = title_name,x = NULL, y = yname) +  # yname 是上面的变量
      scale_fill_manual(values = fig_color) +
      theme(legend.title = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 2,size = 24,family = 'Times New Roman'))

  ggsave(p,file = figname,dpi = 1200)

}else{
    next()
  }

}