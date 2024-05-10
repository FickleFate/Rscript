#! /usr/bin/R

setwd("/home/fate/Dropbox/wsl/SPL/others/gene_location")

library(RIdeogram)
library(rsvg)



# 获得基因密度
gene_density <- GFFex(input = "/home/fate/Dropbox/wsl/SPL/others/gene_location/Pbournei.v2.genome.gff",
    karyotype = "/home/fate/Dropbox/wsl/SPL/others/gene_location/ChrLen.txt",
    feature = "gene", window = 1000000)

# 读取Chr文件和feature文件

chr_info <- read.table('/home/fate/Dropbox/wsl/SPL/others/gene_location/ChrLen.txt',header = TRUE)

chr_marker <- read.csv('/home/fate/Dropbox/wsl/SPL/others/gene_location/location_info.table',header = TRUE)

ideogram(karyotype = chr_info, 
    overlaid = gene_density,label = chr_marker, 
    label_type = "marker",
    colorset1 = c("#FFF0F3","#FF4D6D","#590D22"))



# rsvg_png("/home/fate/Dropbox/wsl/SPL/others/gene_location/chromosome.svg", file = "seqlog.png", width = 1200, height = 700)
# print(info)
# print(gene_density)

