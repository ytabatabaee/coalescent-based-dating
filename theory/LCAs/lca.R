lg = read.csv('lcas-pairs.txt.xz',sep=" ",h=F)
nrow(lg)
head(lg)

ls = read.csv('lcas-st-pairs.txt', sep=" ", h=F)
nrow(ls)
head(ls)

l = merge(lg,ls,by=1:2,all=T)
nrow(l)
head(l)

l$GS = l$V3.x - l$V3.y 
l = l[l$GS!=0,]
head(l)
require(ggplot2)
ggplot(aes(x=GS),data=l)+stat_density()+scale_x_continuous(trans = "log10",name="GS (generations)")+
  theme_bw()+geom_vline(color="red",xintercept = 400000)+
  coord_cartesian(xlim=c(100,10^7))
ggsave("GSdist.pdf",width=3,height = 1.2)
