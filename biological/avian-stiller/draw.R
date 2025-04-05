require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern);library(stringr);library(dplyr)
library(tidyr)

df=read.csv('genera_castlespro_concat.csv')
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C", "#FDBF6F","#FF7F00", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("median_genera.pdf",width=3,height = 5.5)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  #geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
  #             colour = "#FB9A90", 
  #             arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C", "#FDBF6F","#FF7F00", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'bottom')+ 
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         legend.title=element_blank())
ggsave("legend_avian.pdf",width=8,height = 5.5)

df=read.csv('families_castlespro_concat.csv')
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2)+
  #geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
  #             colour = "#7C8385", 
  #             arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C", "#FDBF6F","#FF7F00", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("median_families.pdf",width=3.1,height = 5.5)

df=read.csv('orders_castlespro_concat.csv')
#df = df[df$method!="MD-CAT+CASTLES-Pro" & df$method!="MD-CAT+ConBL",]
#df = df[df$method!="TreePL+CASTLES-Pro" & df$method!="TreePL+ConBL",]
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=order,color=method),data=df)+
  geom_point(size=2)+
  #geom_segment(data = df2, aes(y=as.numeric(factor(order)),yend=as.numeric(factor(order)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+Concat(RAxML)`), 
  #             colour = "#7C8385", 
  #             arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C", "#FDBF6F","#FF7F00", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("median_orders.pdf",width=3,height = 5.5)


df=read.csv('genera_castlespro_concat.csv')
df = df[df$method!="TreePL+CASTLES-Pro" & df$method!="MD-CAT+CASTLES-Pro",]
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A","#E31A1C", "#FF7F00"), name="")+
  theme(legend.position = 'bottom')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("median_genera.pdf",width=3,height = 5.5)

df=read.csv('families_castlespro_concat.csv')
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C", "#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("median_families.pdf",width=3.1,height = 5.5)


df=read.csv('genera_castlespro_concat.csv')
df = df[df$method!="MD-CAT+CASTLES-Pro" & df$method!="MD-CAT+ConBL",]
#df = df[df$method!="TreePL+CASTLES-Pro" & df$method!="TreePL+ConBL",]
df = df[df$method!="wLogDate+CASTLES-Pro" & df$method!="wLogDate+ConBL",]
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("mdcat_median_genera.pdf",width=3,height = 5.5)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`TreePL+CASTLES-Pro`,x=`TreePL+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("treepl_median_genera.pdf",width=3,height = 5.5)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`wLogDate+CASTLES-Pro`,x=`wLogDate+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("wlogdate_median_genera.pdf",width=3,height = 5.5)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2)+
  #geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+ConBL`), 
  #             colour = "#FB9A90", 
  #             arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'bottom')+ 
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         legend.title=element_blank())
ggsave("legend_wlogdate_mcmctree.pdf",width=8,height = 5.5)

ggplot(aes(x=age,y=genera,color=method),data=df[df$method!='MD-CAT+Concat(RAxML)',])+
  geom_point(size=2)+
  scale_color_brewer(palette = "Dark2",name=element_blank())+
  geom_segment(data = df2, aes(y=as.numeric(factor(genera)),yend=as.numeric(factor(genera)),xend=`MD-CAT+CASTLES-Pro`,x=`MCMCtree`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_classic()+
  theme(legend.position = 'bottom')+ 
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         legend.title=element_blank())
ggsave("mcmctree_genera.pdf",width=3,height = 5.5)


df=read.csv('families_castlespro_concat.csv')
df = df[df$method!="MD-CAT+CASTLES-Pro" & df$method!="MD-CAT+ConBL",]
#df = df[df$method!="TreePL+CASTLES-Pro" & df$method!="TreePL+ConBL",]
#df = df[df$method!="wLogDate+CASTLES-Pro" & df$method!="wLogDate+ConBL",]
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+Concat(RAxML)`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("mdcat_median_families.pdf",width=3.1,height = 5.5)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`TreePL+CASTLES-Pro`,x=`TreePL+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("treepl_median_families.pdf",width=3.1,height = 5.5)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`wLogDate+CASTLES-Pro`,x=`wLogDate+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("wlogdate_median_families.pdf",width=3.1,height = 5.5)

df=read.csv('orders_castlespro_concat.csv')
df = df[df$method!="MD-CAT+CASTLES-Pro" & df$method!="MD-CAT+ConBL",]
#df = df[df$method!="TreePL+CASTLES-Pro" & df$method!="TreePL+ConBL",]
#df = df[df$method!="wLogDate+CASTLES-Pro" & df$method!="wLogDate+ConBL",]
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=order,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(order)),yend=as.numeric(factor(order)),xend=`MD-CAT+CASTLES-Pro`,x=`MD-CAT+Concat(RAxML)`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FB9A99","#E31A1C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("mdcat_median_orders.pdf",width=3,height = 5.5)


ggplot(aes(x=age,y=order,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(order)),yend=as.numeric(factor(order)),xend=`TreePL+CASTLES-Pro`,x=`TreePL+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", "#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("treepl_median_orders.pdf",width=3,height = 5.5)

ggplot(aes(x=age,y=order,color=method),data=df)+
  geom_point(size=2)+
  geom_segment(data = df2, aes(y=as.numeric(factor(order)),yend=as.numeric(factor(order)),xend=`wLogDate+CASTLES-Pro`,x=`wLogDate+ConBL`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#6A3D9A", '#B2DF8A',"#33A02C"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("wlogdate_median_orders.pdf",width=3,height = 5.5)



corr=read.csv('castles_caml_corr.csv')

ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.6)+
  scale_x_continuous(trans="log10",name="Branch length (CASTLES-Pro)")+
  scale_y_continuous(trans="log10",name="Branch length (ConBL)")+
  facet_wrap(~Method,ncol=4)+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = c(.3,.25)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("median_corr_r.pdf",width=9.5,height = 3.5)

corr=read.csv('node_age_corr.csv')

ggplot(aes(x=l1,y=l2,color=Node.Type),data=corr)+
  geom_point(alpha=0.6)+
  scale_x_continuous(trans="log10",name="Node age (CASTLES-Pro)")+
  scale_y_continuous(trans="log10",name="Node age (ConBL)")+
  facet_wrap(~Method,ncol=4)+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = 'none') + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("node_age_corr.pdf",width=9.5,height = 3.5)



corr=read.csv('mdcat_median_castles_caml_corr.csv')

ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.7)+
  scale_x_continuous(trans="log10",name="MD-CAT+CASTLES-Pro length")+
  scale_y_continuous(trans="log10",name="MD-CAT+ConBL length")+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = c(.2,.8)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("mdcat_median_corr_r.pdf",width=4,height = 3.5)

corr=read.csv('wlogdate_median_castles_caml_corr.csv')

ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.7)+
  scale_x_continuous(trans="log10",name="wLogDate+CASTLES-Pro length")+
  scale_y_continuous(trans="log10",name="wLogDate+ConBL length")+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = c(.2,.8)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("wlogdate_median_corr_r.pdf",width=4,height = 3.5)

corr=read.csv('treepl_median_castles_caml_corr.csv')

ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.7)+
  scale_x_continuous(trans="log10",name="TreePL+CASTLES-Pro length")+
  scale_y_continuous(trans="log10",name="TreePL+ConBL length")+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2,100),ylim=c(10^-2,100))+
  theme_classic()+
  theme(legend.position = c(.2,.8)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("treepl_median_corr_r.pdf",width=4,height = 3.5)

corr=read.csv('mdcat_median_castles_caml_corr.csv')
a=read.csv('annotations.txt',s="\t")
a[a$clade_mag7=="Palaeognathae","color_mag7_hex"]="#9D6D24"
a$taxa <- gsub('_',' ',a$taxa)
corr_per_family=(merge(corr,a[,c(1,3,4,5,8)],by.x="Taxon",by.y="taxa"))

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=clade_mag7, y=1/as.numeric(l1.l2), fill=clade_mag7)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  #scale_color_brewer(palette = "Spectral",name="")+
  scale_color_manual(name = "", values = unique(corr_per_family$clade_mag7), labels = unique(corr_per_family$clade_mag7)) +
  scale_y_continuous(trans="log2",labels=c(expression(1/3),1,expression(3)),limit = c(-1, 4),breaks=c(1/3,1,3),name="Branch length ratio (log scale)")+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  #coord_cartesian(ylim=c(0,3.8))+
  coord_cartesian(ylim=c(-1,4.5))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  #annotate(geom="text",label="a)", x = 1, y = 3.1, size = 6, width=2)+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("mdcat_corr_violin_higher_clade_log.pdf",width=6,height=4)

newd <-  corr_per_family %>% group_by(order) %>% filter(n()>2)
ggplot(newd[newd$Branch.Type=='terminal',], aes(x=order, y=1/as.numeric(l1.l2), fill=order)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  #scale_color_manual(name = "", values = unique(corr_per_family$order), labels = unique(corr_per_family$order)) +
  scale_y_continuous(trans="log2",labels=c(expression(1/6),1,expression(6)),limit = c(-1, 6),breaks=c(1/6,1,6),name="Branch length ratio (log scale)")+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  #coord_cartesian(ylim=c(0,5))+
  coord_cartesian(ylim=c(-1,8))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("corr_violin_order_log.pdf",width=8,height=5.5)

newd <-  corr_per_family %>% group_by(family) %>% filter(n()>2)
ggplot(newd[newd$Branch.Type=='terminal',], aes(x=family, y=1/as.numeric(l1.l2), fill=family)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  #scale_color_manual(name = "", values = unique(corr_per_family$order), labels = unique(corr_per_family$order)) +
  scale_y_continuous(trans="log2",labels=c(expression(1/8),1,expression(8)),limit = c(-1, 8),breaks=c(1/8,1,8),name="Branch length ratio (log scale)")+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  #coord_cartesian(ylim=c(0,8))+
  theme_classic()+
  coord_cartesian(ylim=c(-1,10))+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("corr_violin_family_log.pdf",width=9,height=7)

corr=read.csv('treepl_median_castles_caml_corr.csv')
a=read.csv('annotations.txt',s="\t")
a[a$clade_mag7=="Palaeognathae","color_mag7_hex"]="#9D6D24"
a$taxa <- gsub('_',' ',a$taxa)
corr_per_family=(merge(corr,a[,c(1,3,4,5,8)],by.x="Taxon",by.y="taxa"))

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=clade_mag7, y=1/as.numeric(l1.l2), fill=clade_mag7)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  #scale_color_brewer(palette = "Spectral",name="")+
  scale_color_manual(name = "", values = unique(corr_per_family$clade_mag7), labels = unique(corr_per_family$clade_mag7)) +
  scale_y_continuous(trans="log2",labels=c(expression(1/2),1,expression(2)),breaks=c(1/2,1,2),name="Branch length ratio (log scale)")+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(-1,2.5))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  #annotate(geom="text",label="a)", x = 1, y = 3.1, size = 6, width=2)+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("treepl_corr_violin_higher_clade_log.pdf",width=6,height=4)

newd <-  corr_per_family %>% group_by(order) %>% filter(n()>2)
ggplot(newd[newd$Branch.Type=='terminal',], aes(x=order, y=1/as.numeric(l1.l2), fill=order)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
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
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  #scale_color_manual(name = "", values = unique(corr_per_family$order), labels = unique(corr_per_family$order)) +
  scale_y_continuous(name="Branch length ratio" )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  #coord_cartesian(ylim=c(0,8))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("corr_violin_family.pdf",width=9,height=7)

corr=read.csv('wlogdate_median_castles_caml_corr.csv')
a=read.csv('annotations.txt',s="\t")
a[a$clade_mag7=="Palaeognathae","color_mag7_hex"]="#9D6D24"
a$taxa <- gsub('_',' ',a$taxa)
corr_per_family=(merge(corr,a[,c(1,3,4,5,8)],by.x="Taxon",by.y="taxa"))

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=clade_mag7, y=1/as.numeric(l1.l2), fill=clade_mag7)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  #scale_color_brewer(palette = "Spectral",name="")+
  scale_color_manual(name = "", values = unique(corr_per_family$clade_mag7), labels = unique(corr_per_family$clade_mag7)) +
  scale_y_continuous(trans="log2",labels=c(expression(1/2),1,expression(2)),breaks=c(1/2,1,2),name="Branch length ratio (log scale)")+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(-1,2.5))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  #annotate(geom="text",label="a)", x = 1, y = 3.1, size = 6, width=2)+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("wlogdate_corr_violin_higher_clade_log.pdf",width=6,height=4)

corr=read.csv('median_castles_caml_corr.csv')
a=read.csv('annotations.txt',s="\t")
a[a$clade_mag7=="Palaeognathae","color_mag7_hex"]="#9D6D24"
a$taxa <- gsub('_',' ',a$taxa)
corr_per_family=(merge(corr,a[,c(1,3,4,5,8)],by.x="Taxon",by.y="taxa"))

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=clade_mag7, y=1/as.numeric(l1.l2), fill=clade_mag7)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  facet_wrap(~Method)+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  #scale_color_brewer(palette = "Spectral",name="")+
  scale_color_manual(name = "", values = unique(corr_per_family$clade_mag7), labels = unique(corr_per_family$clade_mag7)) +
  scale_y_continuous(trans="log2",labels=c(expression(1/3),1,expression(3)),breaks=c(1/3,1,3),name="Branch length ratio (log scale)\n (ConBL / CoalBL)")+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(0.25,4))+
  theme_classic()+
  theme(legend.position = c(0.5,0.18), 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  #annotate(geom="text",label="a)", x = 1, y = 3.1, size = 6, width=2)+
  guides(fill=guide_legend(nrow=3),
         linetype=guide_legend(nrow=3))
ggsave("corr_violin_higher_clade_log.pdf",width=8,height=3)


