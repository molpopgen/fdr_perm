library(ggplot2)
library(reshape)
x=read.table("mean_fdr.txt",header=T)
xr=melt(x,id=c("upreg","tpr","noise"))
print(head(xr))
p = ggplot(xr,aes(tpr,value)) + geom_line(aes(linetype=variable)) + facet_wrap(upreg~noise) 

ggsave("fdr.pdf",p)
