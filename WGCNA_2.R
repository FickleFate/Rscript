### 导入需要的包
library(tidyverse)
library(magrittr)
library(WGCNA)
library(pheatmap)
library(reshape2)
options(stringsAsFactors = FALSE)
library(ggplot2)
library(RColorBrewer)
enableWGCNAThreads()


### 读取数据
expr_normalized <- read.csv('/home/fate/Music/OneDrive/postgrad/转录组数据/WGCNA_result_finally/data_used/diff_gene.csv',row.names = 1)
sample_info <- read.csv('/home/fate/Music/OneDrive/postgrad/转录组数据/WGCNA_result_finally/data_used/sample_info.csv',row.names  = 1 )


### WGCNA
# 转置矩阵，符合WGCNA的要求
datExpr = t(expr_normalized)

# 挑选软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 画出结果
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")

abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 选择第一个到达红线的数字作为软阈值
# 矫正

# sft$powerEstimate即为最佳软阈值
picked_power = sft$powerEstimate





# 一步构建TOM表格
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
net <- blockwiseModules(datExpr,                # <= input here

                            # == Adjacency Function ==
                            power = picked_power,                # <= power here
                            networkType = "signed",

                            # == Tree and Block Options ==
                            deepSplit = 2,
                            pamRespectsDendro = F,
                            # detectCutHeight = 0.75,
                            minModuleSize = 30,
                            maxBlockSize = 4000,

                            # == Module Adjustments ==
                            reassignThreshold = 0,
                            mergeCutHeight = 0.1,

                            # == TOM == Archive the run results in TOM file (saves time)
                            saveTOMs = F,
                            saveTOMFileBase = "ER",

                            # == Output Options
                            numericLabels = T,
                            verbose = 3)

cor <- temp_cor     # Return cor function to original namespace



# 保存图片
pdf("cluster.pdf")

# 绘制聚类图
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
    net$dendrograms[[1]],
    mergedColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05 )

dev.off()

# # 保存module信息
module_df <- data.frame(
    gene_id = names(net$colors),
    colors = labels2colors(net$colors)
)

module_df[1:5,]

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")




#// 模块与性状的关系

nSamples = nrow(datExpr)

design=model.matrix(~0+ sample_info$type)

nnames = gsub("sample_info\\$type","",colnames(design))

colnames(design) = nnames

moduleColors <- labels2colors(net$colors)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

MEs = orderMEs(MEs0)

# 计算相关性矩阵和p值矩阵
moduleTraitCor <- cor(MEs, design, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


# 计算R2
# moduleTraitR2 <- round(moduleTraitCor^2, digits = 2)



# 将相关性矩阵和p值矩阵转换为长格式以便绘图
cor_df <- melt(moduleTraitCor)
pvalue_df <- melt(moduleTraitPvalue)

# 合并相关性和p值
cor_df$pvalue <- pvalue_df$value



pvalue_formatted <- format(cor_df$pvalue, digits = 2)
colors = brewer.pal(11, "RdYlBu")

# 绘制ggplot2热图
p1 =ggplot(cor_df, aes(x = Var2, y = Var1, fill = value)) +
  	geom_tile() +
  	scale_fill_gradient2(
    	low = colors[11], mid = colors[6], high = colors[1],
    	midpoint = 0,
    	limits = c(-1, 1),
    	name = "Correlation") +
  	geom_text(aes(label = paste(round(value, 2), "\n (",pvalue_formatted , ")", sep = "")),
        size = 3, color = "black") +
  	theme_minimal() +
  	theme(
   		axis.text.x = element_text(size = 20),  # x轴标签字体大小
    	axis.text.y = element_text(size = 20))+
  	labs(
    	x = "Trait",
    	y = "Module",
    	title = "Module-trait relationships") +
  	theme(
    	plot.title = element_text(size = 20))
    	

ggsave(p1,file = "/home/fate/Music/OneDrive/postgrad/转录组数据/WGCNA_result_finally/result_need_save/heatmap.png",dpi = 1200)