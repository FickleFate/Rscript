#! /usr/bin/Rscript

library(ggtree)
library(ggplot2)
library(gggenes)
library(treeio)
library(patchwork)


# 进化树源文件
tree <- read.tree("/home/fate/Dropbox/wsl/SPL/gene_structure/image_files/data/1.treefile")
group <- read.csv("/home/fate/Dropbox/wsl/SPL/gene_structure/image_files/data/Pb_SPL_class.csv")
high <- read.csv("/home/fate/Dropbox/wsl/SPL/gene_structure/image_files/data/high1.csv")

# 合并进化树相关信息
spl <- full_join(tree,group,by = "label")

# 绘制进化树
tree <- ggtree(spl, branch.length='none')+
    geom_tiplab()+
    geom_highlight(high,aes(node = id,fill = type),align = "both")+
    theme(legend.position = "none")+
    xlim(NA,12)

# 基因结构信息
feature <- read.csv('/home/fate/Dropbox/wsl/SPL/gene_structure/image_files/data/Rename.structure.csv')

# Y轴按进化树排序
count = c("PbSPL7","PbSPL11","PbSPL17",
    "PbSPL3","PbSPL18","PbSPL6",
    "PbSPL4","PbSPL12","PbSPL8",
    "PbSPL1","PbSPL5","PbSPL15","PbSPL20",
    "PbSPL2","PbSPL19","PbSPL13",
    "PbSPL10","PbSPL16","PbSPL14","PbSPL9")
feature$label <- factor(feature$label, levels = count)

# 绘制基因结构图
structure <- ggplot(feature, aes(xmin = start, xmax = end, y = label, fill = feature)) +
     geom_gene_arrow() +
     facet_wrap(~ label, scales = "free", ncol = 1) +
     scale_fill_brewer(palette = "Set3")+
     theme_genes()+
     labs(x = NULL, y = "")

# # 保存图像
# ggsave("gene_structure.jpg",plot = structure,width = 15,height = 12,units = 'in',dpi = 1000)
# ggsave("tree.jpg",plot = tree,width = 15,height = 12,units = 'in', dpi = 1000)

print(structures)
