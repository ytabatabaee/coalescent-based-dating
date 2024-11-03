require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern);library(stringr);library(dplyr)
library(tidyr)

df=read.csv('genera_mdcat_castlespro_concat.csv')
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=genera,color=method),data=df)+
  geom_point()+
  scale_color_brewer(palette = "Dark2",name=element_blank())+
  geom_segment(data = df2, aes(y=as.numeric(factor(genera)),yend=as.numeric(factor(genera)),x=`MD-CAT+CASTLES-Pro`,xend=`MD-CAT+Concat(RAxML)`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_classic()+
  theme(legend.position = 'bottom')+ 
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         legend.title=element_blank())
ggsave("mdcat_median_genera.pdf",width=3,height = 5)

corr=read.csv('castles_caml_corr.csv')


ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.6)+
  scale_x_continuous(trans="log10",name="CASTLES-Pro length")+
  scale_y_continuous(trans="log10",name="CAML length")+
  facet_wrap(~Method,ncol=4)+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = c(.3,.25)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("median_corr_r.pdf",width=9.5,height = 3.5)



corr=read.csv('mdcat_median_castles_caml_corr.csv')


ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.7)+
  scale_x_continuous(trans="log10",name="MD-CAT+CASTLES-Pro length")+
  scale_y_continuous(trans="log10",name="MD-CAT+CAML length")+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = c(.2,.8)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("mdcat_median_corr_r.pdf",width=4,height = 3.5)


corr=read.csv('mdcat_median_castles_caml_corr.csv')
a=read.csv('annotations.txt',s="\t")
a[a$clade_mag7=="Palaeognathae","color_mag7_hex"]="#9D6D24"
a$taxa <- gsub('_',' ',a$taxa)
corr_per_family=(merge(corr,a[,c(1,3,4,5,8)],by.x="Taxon",by.y="taxa"))

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=clade_mag7, y=1/as.numeric(l1.l2), fill=clade_mag7)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  #scale_color_brewer(palette = "Spectral",name="")+
  scale_color_manual(name = "", values = unique(corr_per_family$clade_mag7), labels = unique(corr_per_family$clade_mag7)) +
  scale_y_continuous(name="Branch length ratio" )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(0,2))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("corr_violin_higher_clade.pdf",width=6,height=4)

newd <-  corr_per_family %>% group_by(order) %>% filter(n()>2)
ggplot(newd[newd$Branch.Type=='terminal',], aes(x=order, y=1/as.numeric(l1.l2), fill=order)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  #scale_color_manual(name = "", values = unique(corr_per_family$order), labels = unique(corr_per_family$order)) +
  scale_y_continuous(name="Branch length ratio" )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(0,5))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("corr_violin_order.pdf",width=8,height=5.5)

newd <-  corr_per_family %>% group_by(family) %>% filter(n()>2)
ggplot(newd[newd$Branch.Type=='terminal',], aes(x=family, y=1/as.numeric(l1.l2), fill=family)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  #scale_color_manual(name = "", values = unique(corr_per_family$order), labels = unique(corr_per_family$order)) +
  scale_y_continuous(name="Branch length ratio" )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(0,8))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("corr_violin_family.pdf",width=9,height=7)
