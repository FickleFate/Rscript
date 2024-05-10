library(ggplot2)

protein_len <- read.table("/home/fate/Dropbox/wsl/SPL/others/SBP/protein_len.txt",
    header = TRUE, sep = ",")
SBP_info <- read.csv("/home/fate/Dropbox/wsl/SPL/others/SBP/SBP_info.csv",
    header = TRUE)

or <- c("PbSPL1",
"PbSPL2",
"PbSPL3",
"PbSPL4",
"PbSPL5",
"PbSPL6",
"PbSPL7",
"PbSPL8",
"PbSPL9",
"PbSPL10",
"PbSPL11",
"PbSPL12",
"PbSPL13",
"PbSPL14",
"PbSPL15",
"PbSPL16",
"PbSPL17",
"PbSPL18",
"PbSPL19",
"PbSPL20")

a <- rev(or)

protein_len$Sequence <- factor(protein_len$Sequence,levels = a)
# SBP_info$Sequence <- factor(SBP_info$Sequence,levels = a)



p1 <- ggplot() +
    geom_linerange(data = protein_len,
        mapping = aes(xmin = 0, xmax = Len, y = Sequence),
        size = 0.5) +
    geom_linerange(data = SBP_info,
        mapping = aes(xmin = from, xmax = to, y = Sequence),
        size = 6,
        color = "#a4133c") +
        labs(x = NULL,y=NULL) +
    theme_classic()


ggsave("/home/fate/Dropbox/wsl/SPL/others/SBP/SBP_domain.png",
    p1,dpi = 900)

# print(p1)