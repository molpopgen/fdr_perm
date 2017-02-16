library(ggplot2)
library(reshape)
x=read.table("mean_fdr.txt",header=T)
xr=melt(x,id=c("upreg","tpr","noise"))
print(head(xr))
p = ggplot(xr,aes(tpr,value)) + geom_line(aes(color=variable,linetype=variable)) + facet_wrap(upreg~noise) + theme_bw() + xlab("True positive rate") + ylab("False discovery rate") + theme(axis.text.x=element_text(angle=75,hjust=1))

ggsave("fdr.pdf",p)
