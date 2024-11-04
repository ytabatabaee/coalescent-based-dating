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

ggplot(corr_per_family[corr_per_family$Branch.Type=='terminal',], aes(x=howardmoore.family, y=1/as.numeric(l1.l2), fill=howardmoore.family)) + 
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #scale_fill_brewer(palette = clad,name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  scale_y_continuous(name="Branch length ratio" )+
  facet_wrap(~Topology,ncol=1)+
  geom_hline(yintercept=1, linetype="dashed", color = "grey40")+
  coord_cartesian(ylim=c(0,15))+
  theme_classic()+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_x_discrete(name="", labels = "")+
  guides(color=guide_legend(ncol=3, byrow=TRUE),
         linetype=guide_legend(ncol=3, byrow=TRUE))
ggsave("subsocines_ratio_family.pdf",width=8.5,height=8.5)

