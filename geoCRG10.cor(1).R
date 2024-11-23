library(corrplot)
library(circlize)

inputFile="diffGeneExp.txt"  
setwd("C:\\biowolf\\geoCRG\\10.cor")     

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)

cor1=cor(rt)

col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)

colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()


pdf(file="corrplot.pdf", width=7, height=7)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
         )
dev.off()