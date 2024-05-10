library(plyr)
library(dplyr)
library(ggplot2)
library(ggsignif)


df <- read.csv('1.csv')
result_list <- list()

for (i in colnames(df)) {
  data <- ddply(df, 'sample', summarise, 
                mean = mean(!!sym(i)),  # 使用 !!sym(i) 插入循环变量作为列名
                sd = sd(!!sym(i)), 
                n1 = sum(!is.na(!!sym(i))), 
                se = sd/sqrt(n1))
  

  ## 添加type
  ## 使用gsub保留最后一个字
  result <- gsub(".*(.{1})", "\\1", data$sample)
  data$type <- result
  
  ## 更改sample的名称，为簇状图分类做准备
  sample <- c('开花前', '初花期', '盛花期', '幼果期', '开花前', '初花期', '盛花期', '幼果期', '开花前', '初花期', '盛花期', '幼果期')
  data$sample <- sample
  
  ## 将sample，type 转为因子，并按照给定的顺序排列  使ggplot2绘图时按其排序
  data$sample <- factor(data$sample, level = c("开花前", "初花期", "盛花期", "幼果期"))
  data$type <- factor(data$type, level = c('叶', '花', '果'))
  
  ## 使用 paste0 来拼接文件名
  filename <- paste0(i, '.csv')
  
  ## 将结果保存到列表中
  result_list[[filename]] <- data
  
  ## 写入CSV文件
  write.csv(data, file = filename, row.names = FALSE, quote = FALSE)

  fig_color = c('#8CCDBF','#CDE0A5','#F9DB95','#EF8476')


## ylab的名称
yname = bquote(.(i) ~ '含量/ng' ~ '·' ~ g^'-1')

## 绘图
p = ggplot(data, aes(x = type, y = mean, fill = sample)) +
    geom_bar(stat = "identity", position = 'dodge', width = 0.5) + # 簇状图，position = 'dodge'
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.5), width = 0.25) + # position = position_dodge(0.5) 使误差线与bar对应
    theme_classic() +
    labs(x = NULL, y = yname) +  # yname 是上面的变量
    scale_fill_manual(values = fig_color) +
    theme(legend.title = element_blank())+
    geom_signif(
      xmin = c(1, 0.8, 1, 1.2, 1.8, 2, 2.2), ## 分组设置，一一对应
      xmax = c(2, 1.0, 1.2, 0.8, 2.0, 2.2, 1.8), ## 同上
      annotations = c("***", "***", "NS", "***", "***", "***", "***"),## 同上
      y_position = c(46, 28, 32, 38, 25, 35, 40), ## 标识（横向）的位置
      textsize = 3, 
      vjust = 0, ## 调整annotions的位置
      tip_length = 0 ## 有没有向下的线
  ) 


  figname =  paste0(i, '.png')
  ggsave(p,file = figname,dpi = 1200)
}

## result_list 中保存了每一列的处理结果，你可以根据需要进一步操作这些结果。
