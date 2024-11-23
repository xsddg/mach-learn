library(reshape2)
library(ggpubr)
inputFile="CIBERSORT-Results.txt"     
setwd("C:\\Users\\Administrator\\Desktop\\PCD在脓毒症中的作用\\备选\\14.核心基因免疫浸润(有问题)")   

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))

pdf(file="barplot.pdf", width=26, height=12)
col=rainbow(nrow(data), s=0.6, v=1.0)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.8)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="#0073C2")
text(a1[conNum]/2,-0.035,"Control",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="#E64B35")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")

group=levels(factor(data$Type))
bioCol=c("#0073C2","#E64B35","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Type",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Type",
				  add="point",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()
