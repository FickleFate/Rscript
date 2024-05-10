library(plyr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(optparse)

#// 描述参数的解析方式
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = FALSE,
              help = "选择文件"),
  make_option(c("-r", "--reference"), type = "character", default = FALSE,
              help = "内参基因名称"),
  make_option(c("-p", "--path"), type = "character", default = FALSE,
              help = "文件保存路径")
)


#// 解析参数
args = parse_args(OptionParser(option_list = option_list, usage = "定量数据分析并作图"))


#// 开始处理数据
data = read.xlsx(args$file)

levels = unique(data$Sample)

data$Sample = factor(data$Sample,level = levels)

gene_name = sort(unique(data$Target))

orign_data_list = list()
AFD_data_list = list()

ref = data[grepl(args$reference,data$Target),]

for (i in gene_name){

    ## 跳过内参基因
    if (i == args$reference){
        next()
    }

    data_name = i

    ## 获得目标基因的信息
    gene_data = data[grepl(i,data$Target),]
    gene_data_deal = gene_data
    
    ## 计算基因与内参Cq值的差值，并添加到数据中
    GSRm = gene_data$Cq - ref$Cq.Mean
    
    gene_data$GSRm = GSRm


    ## 计算对照与内参的差值平均值  并添加到数据框
    GSRm_mean = mean(gene_data[1:3,]$GSRm)

    gene_data$GSRm_mean = GSRm_mean

    ## 计算基因与内参的差值 减去 对照与内参的差值平均值 并添加到数据框
    FSG = gene_data$GSRm-gene_data$GSRm_mean

    gene_data$FSG = FSG

    ## 计算power值
    powers = 2^-gene_data$FSG
    
    gene_data$powers = powers


    orign_data_list[[data_name]] = gene_data

    ## 调用ddply函数，计算power 的均值，标准差，标准误
    gene_data_deal= ddply(gene_data, 'Sample', summarise, 
                mean = mean(powers),  # 使用 !!sym(i) 插入循环变量作为列名
                sd = sd(powers), 
                n = sum(!is.na(Sample)), 
                se = sd/sqrt(n))

    AFD_data_list[[data_name]] = gene_data_deal





    p = ggplot(gene_data_deal, aes(x = Sample, y = mean)) +
    geom_bar(stat = "identity", fill = '#8CCDBF',width = 0.5) + # fill 更改柱子颜色
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.5), width = 0.25) + # position = position_dodge(0.5) 使误差线与bar对应
    theme_classic() +
    labs(title = data_name,x = NULL, y = '相对表达量')+
    theme(plot.title = element_text(hjust = 0.5, vjust = 2,size = 24,family = 'Times New Roman'),
          axis.text.x = element_text(size = 14),  # 更改横坐标标签字体大小
          axis.text.y = element_text(size = 14),  # 更改纵坐标标签字体大小
          axis.title.x = element_text(size = 18), # 更改横坐标标题字体大小
          axis.title.y = element_text(size = 18))   # yname 是上面的变量
    paths = args$path
    filename = paste(paths,data_name,'.png',sep = '')

    ggsave(p,file = filename,dpi = 1200)
 }