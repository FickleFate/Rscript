library(DESeq2)



dat = read.csv('/home/fate/OneDrive/postgrad/转录组分析/DEG_way/old/gene_counts.csv',row.name = 1)

coldata = read.csv('/home/fate/OneDrive/postgrad/转录组分析/DEG_way/old/sample_info.csv',row.name = 1)

dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~ type)

dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

#注意，需将 treat 在前，control 在后，意为 treat 相较于 control 中哪些基因上调/下调
res <- results(dds1, contrast = c('type', 'treat', 'CK'))

#输出表格至本地
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'M12h_treat_CK_all.txt', col.names = NA, sep = '\t', quote = FALSE)


##筛选差异表达基因
#首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
#log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'

#输出选择的差异基因总表
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'M12h_treat_CK_selected.table', sep = '\t', col.names = NA, quote = FALSE)

#根据 up 和 down 分开输出
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'M12h_treat_CK.up.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'M12h_treat_CK.down.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1, file = 'M12h_treat_CK_NA.txt', sep = '\t', col.names = NA, quote = FALSE)

##ggplot2 差异火山图
library(ggplot2)

#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
geom_point(size = 1) +  #绘制散点图


































































scale_color_manual(values = c('steelblue','gray','brown'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '') +  #坐标轴标题
theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
    panel.background = element_rect(color = 'black', fill = 'transparent'), 
    legend.key = element_rect(fill = 'transparent')) +
geom_vline(xintercept = c(-1, 1), lty = 12, color = 'black') +  #添加阈值线
geom_hline(yintercept = 2, lty = 12, color = 'black') +
xlim(-12, 12) + ylim(0, 125)  #定义刻度边界

# ggsave('DEGS.pdf',p)

