library(ggplot2)

df = read.csv('/home/fate/Dropbox/wsl/meme/infomation.csv')

p1 = ggplot()+
    geom_linerange(data = df, mapping = aes(xmin =0, xmax = 2000, y = geneid), color = "#000000",linewidth = 0.5)+
    geom_linerange(data = df,mapping = aes(xmin = st, xmax = st+20, y = geneid ,color = type),linewidth = 6)+
    labs(x = NULL,y=NULL)+
    # scale_color_manual(values = c("#FFBE7A","#FA7F6F","#8ECFC9"))+
    theme_classic()

print(p1)