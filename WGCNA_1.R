
### 导入需要的包


library(tidyverse)
library(magrittr)
library(WGCNA)
library(DESeq2)  
library(genefilter)
library(pheatmap)
library(reshape2)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()


## 读取数据


expr_normalized <- read.csv('/home/fate/Documents/trns_analyze/WGCNA_by_FPKM/data_used/diff_gene.csv',row.names = 1)
sample_info <- read.csv('/home/fate/Documents/trns_analyze/WGCNA_by_FPKM/data_used/sample_info.csv',row.names  = 1 )


## 标准化


RS_data_0 = as.matrix(RS_data[,-1])
row.names(RS_data_0) =RS_data$gene_id


# sample_info  为 样品信息矩阵
# 新建变量，将data中名称一行转化为列，并命名为 sample
# %>% 即  向右操作符，forward-pipe operator）是最常用的一种操作符，就是把左侧准备的数据或表达式，传递给右侧的函数调用或表达式进行运行，可以连续操作就像一个链条一样
# mutate 即 添加新的变量
# gsub 即 替换字符 
sample_info <- read.csv('/home/fate/Documents/trns_analyze/WGCNA_by_FPKM/sample_info.csv',row.names  = 1 )


# 构建dds矩阵
# design 为 差异比较矩阵
dds <- DESeqDataSetFromMatrix(countData = RS_data_0, colData = sample_info, design= ~ type)

#对原始dds进行normalize
dds <- DESeq(dds)

# 获取标准化之后的表格
expr_normalized <- getVarianceStabilizedData(dds)


## WGCNA


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



# 性状与表型关联图

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(datExpr, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = sample_info$type

# tidy & plot data
mME = MEs0 %>%
    pivot_longer(-treatment) %>%
    mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
    )



p1 <- mME %>% ggplot(., aes(x=treatment, y=name, fill=value,label = value)) +
    geom_text() +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(
    low = '#001219',
    high = '#9B2226',
    mid = 'white',
    midpoint = 0,
    limit = c(-1,1)) +
    theme(axis.text.x = element_text(angle=90)) +
    labs(title = "Module-trait Relationships", y = "Modules", fill="corr")









## 基因与模块的关联程度
## 这个sample——info 中要用数字表示性状


moduleLabels = net$colors

moduleColors = labels2colors(net$colors)

MEs = net$MEs

geneTree = net$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

express =as.data.frame(datTraits$type);
names(express) = "express"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, express, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(express), sep="");
names(GSPvalue) = paste("p.GS.", names(express), sep="");



module = "greenyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body express",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



# 挑选一个或多个感兴趣的模块
# 这个变量导出时也要用到
modules_of_interest = "greenyellow"
# pick out a few modules of interest here
# modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
    subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
    mutate(
    gene_id = row.names(.)
    ) %>%
    pivot_longer(-gene_id) %>%
    mutate(
    module = module_df[gene_id,]$colors
    )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
    geom_line(aes(color = module),
            alpha = 0.2) +
    theme_bw() +
    theme(
    axis.text.x = element_text(angle = 90)
    ) +
    facet_grid(rows = vars(module)) +
    labs(x = "treatment",
        y = "normalized expression")






# 生成并导出网络
# 选择感兴趣的模块
modules_of_interest = "greenyellow"
genes_of_interest = module_df %>%
    subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
    mutate(
    gene1 = row.names(.)
    ) %>%
    pivot_longer(-gene1) %>%
    dplyr::rename(gene2 = name, correlation = value) %>%
    unique() %>%
    subset(!(gene1==gene2)) %>%
    mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
    )

head(edge_list)

write_delim(edge_list,
            file = "greenyellow.tsv",
            delim = "\t")

