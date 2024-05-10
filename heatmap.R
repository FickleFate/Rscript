# 加载包
library(pheatmap)
library(RColorBrewer)


# 读取数据
df <- read.csv('/home/fate/OneDrive/postgrad/试验相关/开花通路基因表达量/甲基化相关/闽楠甲基化相关/PbMET1_exp.csv',row.name = 1)

# 读取注释文件
# 第一列为对应上面表格的gene名  第二列为注释信息 
# 第一列不要表头，第二列需要s


p1 = pheatmap(df, 
    # border_color = "black",
    scale = 'row',
    # display_numbers = TRUE,
    cluster_cols= FALSE,
    # cluster_rows = FALSE,
    cellwidth =30,
    cellhight =30,
    width = 15,
    hight = 50,
    fontsize = 8,
    # 添加注释
    # annotation_row = anno_for_MADS,
    color = colorRampPalette(colors = c("blue","white","red"))(100),
    angle_col = 45,
    # 保存文件
    filename = '/home/fate/OneDrive/postgrad/试验相关/开花通路基因表达量/甲基化相关/闽楠甲基化相关/MET1.pdf')
