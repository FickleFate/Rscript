library(ggplot2)
library(ggthemes)

# miRNA作用位点文件
target_info <- read.csv('/home/fate/Dropbox/wsl/miRNA_data/miR156/miR156_target_prediction/miR_target_info.csv',
    header = TRUE)

# mRNA特征文件 如CDS UTR 等
mRNA_featrure <- read.csv('/home/fate/Dropbox/wsl/miRNA_data/miR156/target_mRNA_feature.csv',
    header = TRUE)

sort = c("Pb_SPL_mRNA_2",
"Pb_SPL_mRNA_3",
"Pb_SPL_mRNA_4",
"Pb_SPL_mRNA_7",
"Pb_SPL_mRNA_8",
"Pb_SPL_mRNA_10",
"Pb_SPL_mRNA_11",
"Pb_SPL_mRNA_12",
"Pb_SPL_mRNA_13",
"Pb_SPL_mRNA_17",
"Pb_SPL_mRNA_20"
)

sort_1 = rev(sort)

mRNA_featrure$gene <- factor(mRNA_featrure$gene,levels = sort_1)
target_info$Target <- factor(target_info$Target,levels = sort_1)

# 作图

fig1 <- ggplot()+
    # 画出mRNA长度
    geom_linerange(data = mRNA_featrure, mapping = aes(xmin = st, xmax = ed, y = gene,color = feature),linewidth = 2)+
    # 划出miRNA靶位点
    geom_linerange(data = target_info,mapping = aes(xmin = Target_start, xmax = Target_end, y = Target ,color = 'miR_binding_site'),linewidth = 11)+
    labs(x = NULL,y=NULL)+
    scale_color_manual(values = c("#FFBE7A","#FA7F6F","#8ECFC9"))+
    theme_classic()



print(fig1)

# ggsave('/home/fate/Dropbox/wsl/miRNA_data/miR156binding_site.jpg',plot = fig1,dpi = 900)

