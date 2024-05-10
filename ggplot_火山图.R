library(ggplot2)


df = read.table('开题/1.DEGs_files/control_treat.DESeq2.select.txt')

p <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
        geom_point(size = 1) +                                                                                      #绘制散点图
        scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +                  #自定义点的颜色
        labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '') +         #坐标轴标题
        theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),                      #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
        geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +                                               #添加阈值线
        geom_hline(yintercept = 2, lty = 3, color = 'black') +
        xlim(-12, 12) + ylim(0, 35)                                                                                 #定义刻度边界