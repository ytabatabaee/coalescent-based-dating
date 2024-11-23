require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern);library(stringr);library(dplyr)
library(tidyr)

corr=read.csv('treepl_castles_caml_corr.csv')

ggplot(aes(x=l1,y=l2,color=Branch.Type),data=corr)+
  geom_point(alpha=0.09)+
  scale_x_continuous(trans="log10",name="Branch length (CASTLES-Pro)")+
  scale_y_continuous(trans="log10",name="Branch length (CAML)")+
  facet_wrap(~Topology,ncol=4)+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-3,100),ylim=c(10^-3,100))+
  theme_classic()+
  theme(legend.position = c(.92,.25)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("suboscines_treepl_corr_main.pdf",width=6,height = 3)

corr=read.csv('node_age_corr.csv')

ggplot(aes(x=l1,y=l2,color=Node.Type),data=corr)+
  geom_point(alpha=0.1)+
  scale_x_continuous(trans="log10",name="Node age (CASTLES-Pro)")+
  scale_y_continuous(trans="log10",name="Node age (CAML)")+
  facet_wrap(~Topology,ncol=4)+
  stat_smooth(se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="black",linetype=2)+
  coord_cartesian(xlim=c(10^-2.3,100),ylim=c(10^-2.3,100))+
  theme_classic()+
  theme(legend.position = 'none') + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("suboscines_treepl_node_age_corr_main.pdf",width=6,height = 3)

corr=read.csv('treepl_castles_caml_corr.csv')
a=read.csv('Species_name_map_uids.csv')
a$tipnamecodes <- gsub('_',' ',a$tipnamecodes)
corr_per_family=(merge(corr,a[,c(3,9)],by.x="Taxon",by.y="tipnamecodes"))
corr_per_family=corr_per_family[!is.na(corr_per_family$howardmoore.family),]

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=howardmoore.family, y=1/as.numeric(l1.l2), fill=howardmoore.family)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  scale_y_continuous(name="Branch length ratio" )+
  facet_wrap(~Topology,ncol=1)+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(0,25))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("subsocines_ratio_family.pdf",width=8.5,height=8.5)

df=read.csv('treeness.csv')
df$comb = paste(df$Topology,df$Method)

ggplot(data=df, aes(x=Method, y=Treeness, fill=comb)) +
  facet_wrap(~Topology,ncol=3)+
  geom_bar(stat="identity",colour="black")+
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","red")) +
  scale_x_discrete(name="Branch length")+
  coord_cartesian(ylim=c(0,0.48))+
  geom_text(aes(label=Treeness), vjust=-0.3, size=3.5)+
  theme_classic()+
  theme(legend.position = "none")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("subsocines_treeness.pdf",width=5,height=2.5)


df=read.csv('genera_treepl_castlespro_concat.csv')
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=genus,color=method),data=df)+
  geom_point(size=2.5)+
  geom_segment(data = df2, aes(y=as.numeric(factor(genus)),yend=as.numeric(factor(genus)),xend=`TreePL+CASTLES-Pro`,x=`TreePL+Concat(RAxML)`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("suboscines_genera.pdf",width=3.7,height = 6.5)

df=read.csv('families_treepl_castlespro_concat.csv')
df2 <- df %>% spread(method,age)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2.5)+
  geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`TreePL+CASTLES-Pro`,x=`TreePL+Concat(RAxML)`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("suboscines_family.pdf",width=3.7,height = 6.5)

ggplot(aes(x=age,y=family,color=method),data=df)+
  geom_point(size=2.5)+
  geom_segment(data = df2, aes(y=as.numeric(factor(family)),yend=as.numeric(factor(family)),xend=`TreePL+CASTLES-Pro`,x=`TreePL+Concat(RAxML)`), 
               colour = "#7C8385", 
               arrow = arrow(length = unit(6,'pt')))+
  theme_bw()+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'bottom')+ 
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         legend.title=element_blank())
ggsave("legend_treepl.pdf",width=5,height = 6.5)


df=read.csv('caml_ltt.csv')

ggplot(df[df$time<62.042784,], aes(x =time,y=lineages,color=method)) +
  geom_step(size=0.8) +
  theme_classic()+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  geom_vline(color="grey50",linetype=2,xintercept = 60)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         legend.title=element_blank())
ggsave("caml_ltt_log.pdf",width=4,height = 3)

ggplot(df[df$time<62.042784 & df$time>60,], aes(x =time,y=lineages,color=method)) +
  geom_step(size=0.8) +
  theme_classic()+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         legend.title=element_blank())
ggsave("caml_ltt_log_end.pdf",width=4,height = 3)

df=read.csv('suboscines_ltt.csv')

ggplot(df[df$time<64.6951439999999,], aes(x =time,y=lineages,color=method)) +
  geom_step(size=0.8) +
  theme_classic()+
  facet_wrap(~factor(topology,c('CAML', 'ASTRAL', 'wASTRAL')))+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  geom_vline(color="grey50",linetype=2,xintercept = 62)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         legend.title=element_blank())
ggsave("stiller_ltt_log.pdf",width=9,height = 3)

ggplot(df[df$time<64.6951439999999 & df$time>62,], aes(x =time,y=lineages,color=method)) +
  geom_step(size=0.8) +
  theme_classic()+
  facet_wrap(~factor(topology,c('CAML', 'ASTRAL', 'wASTRAL')))+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_manual(values=c("#FDBF6F","#FF7F00"), name="")+
  theme(legend.position = 'none')+ 
  geom_vline(color="grey50",linetype=2,xintercept = 62)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         legend.title=element_blank())
ggsave("stiller_ltt_log_end.pdf",width=9,height = 3)

ggplot(df[df$time<64.6951439999999,], aes(x=time,y=lineages,color=interaction(method,topology,sep='+')))+
  geom_step(size=0.8)+
  theme_classic()+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_brewer(palette = "Paired",name="")+
  theme(legend.position = 'none')+ 
  geom_vline(color="grey50",linetype=2,xintercept = 62)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("stiller_ltt_log_all.pdf",width=4,height = 3)

ggplot(df[df$time<64.6951439999999 & df$time>62,], aes(x=time,y=lineages,color=interaction(method,topology,sep='+')))+
  geom_step(size=0.8)+
  theme_classic()+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_brewer(palette = "Paired",name="")+
  theme(legend.position = 'none')+ 
  geom_vline(color="grey50",linetype=2,xintercept = 62)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         legend.title=element_blank())
ggsave("stiller_ltt_log_all_end.pdf",width=4,height = 3)


ggplot(df[df$time<64.6951439999999,], aes(x=time,y=lineages,color=interaction(method,topology,sep='+')))+
  geom_step(size=0.8)+
  theme_classic()+
  scale_x_continuous(name="Time (million years)")+
  scale_y_continuous(name="Number of lineages",trans = 'log10')+
  scale_color_brewer(palette = "Paired",name="")+
  theme(legend.position = 'bottom')+ 
  geom_vline(color="grey50",linetype=2,xintercept = 62)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         legend.title=element_blank())
ggsave("legend_stiller_5.pdf",width=8,height = 3)

# lineage-through-time plots
library(phytools)

tree1 <- read.tree(file="./treepl_castles_T400F.examl.rooted.tre")
tree1<-force.ultrametric(tree1)

tree2 <- read.tree(file="./treepl_concat_T400F.examl.rooted.tre")
tree2<-force.ultrametric(tree2)

is.ultrametric(tree1)
is.ultrametric(tree2)

trees <- c(tree1, tree2)
trees<-lapply(trees,force.ultrametric)
class(trees)<-"multiPhylo"
ltt(trees,log=TRUE)#,xlim=c(60,65))

#ltt(tree2,log=TRUE,plot=TRUE)#,xlim=c(60,65))


