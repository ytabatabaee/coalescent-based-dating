require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern);require(tidyverse)

## S100 dataset

### node age
s0=read.csv('s100_dating_unit_node_age.csv')
s3=read.csv('s100_dating_normalized_n3_node_age.csv')
s10=read.csv('s100_dating_normalized_n10_node_age.csv')

s0$Calibrations <- 0
s3$Calibrations <- 3
s10$Calibrations <- 10
s <- rbind(s0,s3,s10)


s = s[s$Taxon.Type!="terminal",]
s$Method = factor(s$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$isconcat = factor(grepl("Concat", s$Method))
s$age.est = ifelse(s$age.est <=0, 1e-6, s$age.est)
s$log10err = log10(s$age.est / s$age.true )
s$abserr = abs(s$age.true - s$age.est)
s$se = (s$age.est - s$age.true)^2 
s$bias = s$age.est - s$age.true
s <- s |> mutate(conditionAD = case_when(
  AD >= 0.3 & AD < 0.4 ~ '[0.3,0.4)',
  AD >= 0.4 & AD <= 0.5 ~ '[0.4,0.5)',
  AD >= 0.5 & AD <= 0.6 ~ '[0.5,0.6)',
))

dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+Taxon.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE), shape='none')
ggsave("S100-node-age_abserr-dating_calib_pro-arrow.pdf",width=5,height = 2.5)


dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+Taxon.Type+isconcat~'se' ,value.var = "se",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~Condition,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE), shape='none')
ggsave("S100-node-age_rmse-dating_calib_pro-arrow.pdf",width=5,height = 2.5)

ggplot(aes(x=as.factor(Calibrations),y=abserr,color=Method),
       data=dcast(data=s, Condition+Calibrations+Method+replicate+Taxon.Type~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_pro-node-age.pdf",width=5,height =2.5)

ggplot(aes(x=as.factor(Calibrations),y=sqrt(se),color=Method),
       data=dcast(data=s, Condition+Calibrations+Method+replicate+Taxon.Type~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-rmse-dating_calib_pro-node-age.pdf",width=5,height =2.5)

ggplot(aes(x=as.factor(Calibrations),y=age.est-age.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Estimated" - "true length (bias)"))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_wrap(~Condition,ncol=4)+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib_broken-line_node_age.pdf",width=5,height = 2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+Taxon.Type+isconcat~'bias' ,value.var = "bias",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression(atop("Estimated" - "true length (bias)", paste("(node age)"))))+
  facet_wrap(~Condition,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE), shape='none')
ggsave("S100-node-age-bias_dating_calib_pro-arrow_main.pdf",width=5,height = 2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+Taxon.Type+isconcat~'bias' ,value.var = "bias",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_wrap(~Condition,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE), shape='none')
ggsave("S100-node-age-bias_dating_calib_pro-arrow.pdf",width=5,height = 2.5)

ggplot(data=dcast(data=s,Condition+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=as.factor(Calibrations),y=log10err,color=Method))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "none", legend.direction = "horizontal")
ggsave("S100-dating_error-logabs-calib_node_age.pdf",width=5,height = 2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+isconcat~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x)))
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(log10err = mean(log10err)) %>% pivot_wider(names_from = isconcat,values_from = log10err) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~Condition,ncol=4)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-2)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-dating_error-logerror_node_age_arrow.pdf",width=5,height = 2.5)


## TMRCA for root-unfixed, S100
s=read.csv('s100_dating_tmrca.csv')
s$Method = factor(s$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$isconcat = factor(grepl("Concat", s$Method))
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 

ggplot(aes(x=Condition, y=l.true/l.est,color=Method,shape=isconcat),data=s[s$Outgroup=='True'&s$Calibrations==3,])+
  scale_y_continuous(trans="identity",name=expression("True / estimated tMRCA"))+
  scale_x_discrete(name="Sequence length")+
  coord_cartesian(ylim=c(0,3))+
  #facet_grid(~Calibrations)+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-tmrca-bias_dating_calib_root-unfixed_main.pdf",width=3.5,height = 3)

ggplot(aes(x=Condition, y=l.true/l.est,color=Method,shape=isconcat),data=s[s$Outgroup=='True',])+
  scale_y_continuous(trans="identity",name=expression("True / estimated tMRCA"))+
  scale_x_discrete(name="Sequence length")+
  coord_cartesian(ylim=c(0,3))+
  facet_grid(~Calibrations)+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-tmrca-bias_dating_calib_root-unfixed.pdf",width=7,height = 3)


### TREENESS for root-unfixed
## TMRCA for root-unfixed, S100
s=read.csv('s100_dating_treeness.csv')
s$Method = factor(s$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$isconcat = factor(grepl("Concat", s$Method))
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 

ggplot(aes(x=Condition, y=l.true/l.est,color=Method,shape=isconcat),data=s[s$Outgroup=='True',])+
  scale_y_continuous(trans="identity",name=expression("True / estimated treeness"))+
  scale_x_discrete(name="Sequence length")+
  facet_grid(~Calibrations)+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-treeness_dating_root-unfixed.pdf",width=7,height = 3)


## Varying number of genes, S100
s=read.csv('s100_dating_normalized_n3_genes.csv')
s$Calibrations <- 3

s$Method = factor(s$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$isconcat = factor(grepl("Concat", s$Method))
s$l.est = ifelse(s$l.est <=0, 1e-6, s$l.est)
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 
s$bias = s$l.est - s$l.true

ggplot(data=dcast(data=s,Method+replicate+isconcat+genes+Condition~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=as.factor(genes),y=log10err,color=Method,shape=isconcat))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~Condition,ncol=4)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_x_discrete(name="Number of genes")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "none", legend.direction = "horizontal")
ggsave("S100-dating_error-logerror_n3_genes.pdf",width=6,height = 2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+genes+Branch.Type+isconcat~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x)))
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,genes,isconcat,datingMethod) %>%
  summarise(log10err = mean(log10err)) %>% pivot_wider(names_from = isconcat,values_from = log10err) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~Condition,ncol=4)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of genes",breaks = c(1, 2, 3, 4),label = c("50", "200", "500", "1000"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-dating_error-logerror_n3_genes_arrow.pdf",width=6,height = 2.5)


ggplot(aes(x=as.factor(genes),y=abserr,color=Method,shape=isconcat),
       data=dcast(data=s, Condition+genes+Method+replicate+Branch.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  scale_x_discrete(name="Number of genes")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-abserr_dating_genes_n3_genes.pdf",width=6,height =2.5)

ggplot(aes(x=as.factor(genes),y=bias,color=Method,shape=isconcat),
       data=dcast(data=s, Condition+genes+Method+replicate+Branch.Type+isconcat~'bias' ,value.var = "bias",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_wrap(~Condition,ncol=4)+
  scale_x_discrete(name="Number of genes")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-bias_dating_genes_n3_genes.pdf",width=6,height =2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+genes+Branch.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,genes,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of genes",breaks = c(1, 2, 3, 4),label = c("50", "200", "500", "1000"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-abserr_dating_genes_n3-genes_arrow_main.pdf",width=7,height = 2.5)


dtemp=dcast(data=s,
            Condition+Method+replicate+genes+Branch.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,genes,isconcat,Branch.Type,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(Branch.Type~Condition)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of genes",breaks = c(1, 2, 3, 4),label = c("50", "200", "500", "1000"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-abserr_dating_genes_n3-genes_arrow_broken.pdf",width=8,height = 3.5)


dtemp=dcast(data=s,
            Condition+Method+replicate+genes+Branch.Type+isconcat~'bias' ,value.var = "bias",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,genes,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_wrap(~Condition,ncol=4)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of genes",breaks = c(1, 2, 3, 4),label = c("50", "200", "500", "1000"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-bias_dating_genes_n3-genes_arrow.pdf",width=6,height = 2.5)

dtemp=dcast(data=s[s$Condition=='800bp',],
            Condition+Method+replicate+genes+Branch.Type+isconcat~'bias' ,value.var = "bias",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,genes,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(trans="identity",name=expression(atop("Estimated" - "true length (bias)", paste("(branch length)"))))+
  facet_wrap(~Condition,ncol=4)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of genes",breaks = c(1, 2, 3, 4),label = c("50", "200", "500", "1000"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-bias_dating_genes_n3-genes_arrow_800bp.pdf",width=2.8,height = 2.5)

ggplot(aes(x=as.factor(genes),y=sqrt(se),color=Method,shape=isconcat),
       data=dcast(data=s, Condition+genes+Method+replicate+Branch.Type+isconcat~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~Condition,ncol=4)+
  scale_x_discrete(name="Number of genes")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-rmse_dating_genes_n3_genes.pdf",width=6,height =2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+genes+Branch.Type+isconcat~'se' ,value.var = "se",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,genes,isconcat,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~Condition,ncol=4)+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(genes)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')))+
  scale_x_continuous(name="Number of genes",breaks = c(1, 2, 3, 4),label = c("50", "200", "500", "1000"))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),shape='none')
ggsave("S100-rmse_dating_genes_n3-genes_arrow.pdf",width=6,height = 2.5)


## S100 branch length

s0=read.csv('s100_dating_unit.csv')
s3=read.csv('s100_dating_normalized_n3.csv')
s10=read.csv('s100_dating_normalized_n10.csv')

s0$Calibrations <- 0
s3$Calibrations <- 3
s10$Calibrations <- 10
s <- rbind(s0,s3,s10)

s$Method = factor(s$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$isconcat = factor(grepl("Concat", s$Method))
s$l.est = ifelse(s$l.est <=0, 1e-6, s$l.est)
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 
s$bias = s$l.est - s$l.true
s <- s |> mutate(conditionAD = case_when(
  AD >= 0.3 & AD < 0.4 ~ '[0.3,0.4)',
  AD >= 0.4 & AD <= 0.5 ~ '[0.4,0.5)',
  AD >= 0.5 & AD <= 0.6 ~ '[0.5,0.6)',
))


ggplot(data=dcast(data=s[s$Calibrations==3,],Condition+Method+replicate+isconcat~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=Condition,y=log10err,color=Method,shape=isconcat))+
  scale_y_continuous(trans="identity",name="Mean log10 error\n(branch length)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary()+
  #stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_error-logerror_pro_n3.pdf",width=3.5,height = 3)

ggplot(data=dcast(data=s[s$Calibrations==3,],Condition+Method+replicate+isconcat~'se' ,value.var = "se",fun.aggregate = mean), aes(x=Condition,y=sqrt(se),color=Method,shape=isconcat))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)\n(branch length)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary()+
  #stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_rmse-overall_pro.pdf",width=3.5,height = 3)

ggplot(data=dcast(data=s,Condition+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=as.factor(Calibrations),y=log10err,color=Method))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_error-logabs-calib.pdf",width=6,height = 2.5)


ggplot(data=dcast(data=s[s$Condition=='800bp',],Condition+Method+replicate+Calibrations+conditionAD~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=as.factor(Calibrations),y=log10err,color=Method))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~conditionAD,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_error-logabs-calib_ils.pdf",width=4.5,height = 2.5)

ggplot(aes(x=l.true,y=l.est,color=Branch.Type,linetype),
       data=s[!s$Condition=='true gene trees'  & s$Calibrations==3,])+
  facet_grid(Method~Condition)+
  scale_x_continuous(trans="log10",name="True length")+
  scale_y_continuous(trans="log10",name="Estimated length")+
  scale_color_brewer(palette = "Dark2")+
  geom_abline(color="grey30",linetype=2)+
  geom_point(alpha=0.1,size=0.5)+
  stat_smooth(se=F,alpha=1,size=0.7,method="glm",formula=y ~ poly(x, 2))+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave("S100-dating-correlation.png",width=8,height =8)



ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est-true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib_overall.pdf",width=3.5,height =  4)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est-true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
    guides(color=guide_legend(nrow=2, byrow=TRUE),
           fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_overall.pdf",width=7,height =  5)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method,shape=isconcat), data=s[s$Calibrations==3,])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_n3-pro.pdf",width=4,height =  3)


ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_wrap(~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib-line.pdf",width=9,height =  4)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  facet_grid(Branch.Type~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib_broken-2.pdf",width=9,height = 6)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  facet_wrap(~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib-2-pro.pdf",width=9,height = 6)

ggplot(aes(x=as.factor(Calibrations),y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_grid(Branch.Type~Condition)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib_broken-line.pdf",width=9,height = 4)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  facet_wrap(~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib.pdf",width=9,height = 4)

ggplot(aes(x=as.factor(genes),y=abserr,color=Method),
       data=dcast(data=s, Condition+genes+Method+replicate+Branch.Type~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_pro.pdf",width=7,height =2.5)


ggplot(aes(x=as.factor(Calibrations),y=abserr,color=Method),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_no_OG.pdf",width=7,height =2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+Branch.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error\n(branch length)")+
  facet_wrap(~Condition,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("S100-error-perrep_dating_calib-arrow-main.pdf",width=5,height = 2.5)

dtemp=dcast(data=s,
            Condition+Method+replicate+Calibrations+Branch.Type+isconcat~'bias' ,value.var = "bias",fun.aggregate = mean)
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Condition,Calibrations,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(trans="identity",name=expression(atop("Estimated" - "true length (bias)", paste("(branch length)"))))+
  facet_wrap(~Condition,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "10"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE), shape='none')
ggsave("S100-bias_dating_calib-arrow_main.pdf",width=5,height = 2.5)

ggplot(aes(x=as.factor(Calibrations),y=abserr,color=Method),
       data=dcast(data=s[s$Condition=='800bp',], conditionAD+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~conditionAD,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  #scale_x_discrete(name="True gene tree discordance (ILS)")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = 'none', legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_ils.pdf",width=4.5,height =2.5)

ggplot(aes(x=as.factor(Calibrations),y=abserr,color=Method),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = 'none', legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib.pdf",width=6,height =2.5)


ggplot(aes(x=as.factor(Calibrations),y=sqrt(se),color=Method),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = 'none', legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-rmse-perrep_dating_calib.pdf",width=6,height =2.5)

ggplot(aes(x=as.factor(Calibrations),y=sqrt(se),color=Method),
       data=dcast(data=s[s$Condition=='800bp',], Condition+Method+replicate+Branch.Type+Calibrations+conditionAD~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~conditionAD,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = 'none', legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-rmse-perrep_dating_calib_ils.pdf",width=4.5,height =2.5)


ggplot(aes(x=as.factor(Calibrations),y=abserr,color=Method),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(Branch.Type~Condition)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = c(0.5,0.85), legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_broken.pdf",width=10,height =4)

ggplot(aes(x=as.factor(Calibrations),y=abserr,color=Method),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Condition,ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = 'bottom', legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib.pdf",width=9,height =3)

ggplot(aes(x=Condition,y=abserr,color=Method),
       data=dcast(data=subset(s, Calibrations %in% c('3') & l.true <= 1), Condition+Method+replicate+Branch.Type+Calibrations+l.true~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~Condition,ncol=4)+
  facet_wrap(~cut(l.true,c(0,0.01,0.1,1),right=F),ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_bl.pdf",width=7,height =3.5)

ggplot(aes(x=Condition,y=log10err,color=Method),
       data=dcast(data=subset(s, Calibrations %in% c('3') & l.true <= 1), Condition+Method+replicate+Branch.Type+Calibrations+l.true~'log10err' ,value.var = "log10err", fun.aggregate = function(x) mean(abs(x))))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  #facet_wrap(~Condition,ncol=4)+
  facet_wrap(~cut(l.true,c(0,0.01,0.1,1),right=F),ncol=4)+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  #stat_summary()+
  #stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-logerror-perrep_dating_calib_bl.pdf",width=7,height =3.5)



ggplot(aes(x=Condition,y=l.est-l.true,color=Method),
       data=dcast(data=subset(s, Calibrations %in% c('3') & l.true <= 1)))+
  scale_y_continuous(trans="identity",name="Bias")+
  #facet_wrap(~Condition,ncol=4)+
  facet_wrap(~cut(l.true,c(0,0.01,0.1,1),right=F),ncol=4)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '10'))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_bl.pdf",width=7,height =3.5)



ggplot(aes(x=Condition,y=sqrt(se),color=Method),
       data=dcast(data=s,Condition+Method+replicate+Branch.Type~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #stat_summary()+
  #stat_summary(aes(group=Method),geom="line")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  guides(fill=guide_legend(title="Method"))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-rmse-perrep_dating.pdf",width=7,height = 5)


## MVROOT dataset

### TMRCA

m=read.csv('mvroot_estgt_dating_tmrca.csv')
m$outgroup = factor(grepl("outgroup.1", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)', 'test', 'MCMCtree'))
m = m[m$Method!="MCMCtree",]
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$isconcat = factor(grepl("Concat", m$Method))
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$ratio = m$l.true/m$l.est
m$se = (m$l.est - m$l.true)^2 
outgroup.labs <- c("With outgroup","No outgroup")
names(outgroup.labs) <- c(TRUE, FALSE)
m <- m |> mutate(conditionAD = case_when(
  AD >= 0 & AD < 0.25 ~ '[0,0.25)',
  AD >= 0.25 & AD <= 0.5 ~ '[0.25,0.5)',
  AD >= 0.5 & AD <= 0.75 ~ '[0.5,0.75)',
  AD >= 0.75 & AD <= 1 ~ '[0.75,1)',
))

ggplot(aes(x=ratevar, y=l.est/l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==FALSE & m$Calibrations==3,])+
  scale_y_continuous(trans="identity",name=expression("True / estimated tMRCA"))+
  scale_x_discrete(labels=c("High","Medium","Low"),name="Clock deviation")+
  coord_cartesian(ylim=c(0,5))+
  #facet_wrap(~Calibrations)+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-tmrca-bias_dating_calib_root-unfixed_ratevar.pdf",width=4,height = 3)


ggplot(aes(x=cut(AD,4), y=(l.true/l.est),color=Method,shape=isconcat),
       data=m[m$outgroup ==FALSE & m$Calibrations==3,])+
  scale_y_continuous(trans="identity",name=expression("True / estimated tMRCA"))+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  coord_cartesian(ylim=c(0,5))+
  #facet_wrap(~Calibrations)+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-tmrca_dating.pdf",width=5,height = 3)

dtemp=merge(
  dcast(data=m[m$outgroup ==FALSE,],
        outgroup+Method+replicate+Calibrations+isconcat+conditionAD~'ratio' ,value.var = "ratio",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==FALSE,], replicate+outgroup~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod,outgroup) %>%
  summarise(ratio = mean(ratio)) %>% pivot_wider(names_from = isconcat,values_from = ratio) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="True / estimated tMRCA")+
  facet_grid(~conditionAD)+#,labeller = labeller(outgroup = outgroup.labs))+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2),label = c("3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C", "#6A3D9A"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-tmrca-arrow.pdf",width=5,height = 2.5)

ggplot(aes(x=cut(AD,4), y=(l.true)/l.est,color=Method,shape=isconcat), data=m)+
  scale_y_continuous(trans="identity",name=expression("True / estimated tMRCA"))+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  coord_cartesian(ylim=c(0,4.5))+
  facet_grid(outgroup~Calibrations, labeller = labeller(outgroup = outgroup.labs))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape="none")
ggsave("MV-tmrca_dating_calib.pdf",width=8,height = 6.2)

ggplot(aes(x=ratevar, y=(l.true)/l.est,color=Method,shape=isconcat), data=m)+
  scale_y_continuous(trans="identity",name=expression("True / estimated tMRCA"))+
  scale_x_discrete(labels=c("High","Medium","Low"),name="Clock deviation")+
  coord_cartesian(ylim=c(0,4.5))+
  facet_grid(outgroup~Calibrations, labeller = labeller(outgroup = outgroup.labs))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape="none")
ggsave("MV-tmrca_dating_calib_ratevar.pdf",width=8,height = 4.5)

### TREENESS

m=read.csv('mvroot_estgt_dating_treeness.csv')
m$outgroup = factor(grepl("outgroup.1", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)', 'test', 'MCMCtree'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$isconcat = factor(grepl("Concat", m$Method))
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$bias = m$l.est - m$l.true
m = m[m$Method!="MCMCtree",]
m$ratio = m$l.true/m$l.est
m$se = (m$l.est - m$l.true)^2 
outgroup.labs <- c("With outgroup","No outgroup")
names(outgroup.labs) <- c(TRUE, FALSE)
m <- m |> mutate(conditionAD = case_when(
  AD >= 0 & AD < 0.25 ~ '[0,0.25)',
  AD >= 0.25 & AD <= 0.5 ~ '[0.25,0.5)',
  AD >= 0.5 & AD <= 0.75 ~ '[0.5,0.75)',
  AD >= 0.75 & AD <= 1 ~ '[0.75,1)',
))

dtemp=merge(
  dcast(data=m[m$outgroup ==FALSE,],
        outgroup+Method+replicate+Calibrations+isconcat+conditionAD~'ratio' ,value.var = "ratio",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==FALSE,], replicate+outgroup~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod,outgroup) %>%
  summarise(ratio = mean(ratio)) %>% pivot_wider(names_from = isconcat,values_from = ratio) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="True / estimated treeness")+
  facet_grid(~conditionAD)+#,labeller = labeller(outgroup = outgroup.labs))+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2),label = c("3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')))+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C", "#6A3D9A"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-treeness-arrow.pdf",width=5,height = 2.5)

ggplot(aes(x=ratevar, y=(l.est)/l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==FALSE & m$Calibrations==3,])+
  scale_y_continuous(trans="identity",name=expression("Estimated / true treeness"))+
  scale_x_discrete(labels=c("High","Medium","Low"),name="Clock deviation")+
  #facet_wrap(~Calibrations)+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-treeness-dating_ratevar.pdf",width=4,height = 3)

ggplot(aes(x=conditionAD, y=l.true/l.est,color=Method,shape=isconcat),
       data=m[m$outgroup ==FALSE & m$Calibrations==3,])+
  scale_y_continuous(trans="identity",name=expression("True / estimated treeness"))+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  #facet_wrap(~Calibrations)+
  #coord_cartesian(ylim=c(0,3))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-treeness_dating_main.pdf",width=5,height = 3)


ggplot(aes(x=conditionAD, y=l.est-l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==FALSE & m$Calibrations==3,])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true treeness (bias)"))+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  #facet_wrap(~Calibrations)+
  #coord_cartesian(ylim=c(0,3))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-treeness_dating.pdf",width=5,height = 3)

ggplot(aes(x=conditionAD, y=l.est-l.true,color=Method,shape=isconcat), data=m)+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true treeness (bias)"))+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_grid(outgroup~Calibrations, labeller = labeller(outgroup = outgroup.labs))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape="none")
ggsave("MV-treeness_dating_calib.pdf",width=8,height = 6.2)

ggplot(aes(x=ratevar, y=l.est-l.true,color=Method,shape=isconcat), data=m)+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true treeness (bias)"))+
  scale_x_discrete(labels=c("High","Medium","Low"),name="Clock deviation")+
  facet_grid(outgroup~Calibrations, labeller = labeller(outgroup = outgroup.labs))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape="none")
ggsave("MV-treeness_dating_calib_ratevar.pdf",width=8,height = 4.5)

### BRANCH LENGTH

m0 = read.csv('mvroot_estgt_dating_unit.csv')
m3 = read.csv('mvroot_estgt_dating_n3_normalized.csv')
m5 = read.csv('mvroot_estgt_dating_n5_normalized.csv')
m0$Calibrations <- 0
m3$Calibrations <- 3
m5$Calibrations <- 5
m <- rbind(m0,m3,m5)
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m = m[m$Method!="MCMCtree",]
#m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))#, 'MCMCtree'))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)', 'MCMCtree'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$isconcat = factor(grepl("Concat", m$Method))

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-3, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$bias = m$l.est - m$l.true
m$se = (m$l.est - m$l.true)^2 
ratevar.labs <- c("Low","Medium", "High")
names(ratevar.labs) <- c("5","1.5", "0.15")
m <- m |> mutate(conditionAD = case_when(
  AD >= 0 & AD < 0.25 ~ '[0,0.25)',
  AD >= 0.25 & AD <= 0.5 ~ '[0.25,0.5)',
  AD >= 0.5 & AD <= 0.75 ~ '[0.5,0.75)',
  AD >= 0.75 & AD <= 1 ~ '[0.75,1)',
))


dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'bias' ,value.var = "bias",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,Branch.Type,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression(atop("Estimated" - "true length (bias)", paste("(branch length)"))))+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')), size=0.7)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  c(0.5,0.8), legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        legend.text=element_text(size=11))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         shape='none')
ggsave("MV-bias_dating_bycalib-pro_broken.pdf",width=10,height = 4)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'bias' ,value.var = "bias",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,Branch.Type,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        legend.text=element_text(size=11))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         shape='none')
ggsave("MV-bias_dating_bycalib-pro_broken-arrow.pdf",width=6,height = 3.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+ratevar~'bias' ,value.var = "bias",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,Branch.Type,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        legend.text=element_text(size=11))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         shape='none')
ggsave("MV-ratevar-bias_dating_bycalib-pro_broken-arrow.pdf",width=4.5,height = 3.5)


ggplot(aes(x=as.factor(Calibrations), y=l.est-l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="Number of calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_grid(Branch.Type~conditionAD)+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  #scale_shape_manual(guide="none")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         shape="none")
ggsave("MV-bias_dating_bycalib-pro.pdf",width=6,height = 3.5)


ggplot(aes(x=as.factor(Calibrations), y=l.est-l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(name="Number of calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),shape="none")
ggsave("MV-bias-ratevar_dating_bycalib-pro.pdf",width=4.5,height = 3.5)

ggplot(aes(x=as.factor(Calibrations), y=l.est-l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(name="Number of calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  theme_classic()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  #scale_shape_manual(guide="none")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),shape="none")
ggsave("legend.pdf",width=8,height = 3.5)


### ABS

ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-abserr_dating_bycalib_line-pro.pdf",width=6,height = 2.5)

ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         shape="none")
ggsave("MV-abserr_dating_bycalib_line-pro_broken.pdf",width=6,height = 3.5)

ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+ratevar~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         shape="none")
ggsave("MV-abserr_ratevar_dating_bycalib_line-pro_broken.pdf",width=4.5,height = 3.5)


ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-abserr_ratevar-dating_bycalib_line-pro.pdf",width=4.5,height = 2.5)


dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-abserr_ratevar-dating_bycalib_line-pro_arrow.pdf",width=4.5,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,Branch.Type,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-abserr_ratevar-dating_bycalib_line-pro_arrow_broken.pdf",width=4.5,height = 3.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-abserr-dating_bycalib_line-pro_arrow.pdf",width=6,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,Branch.Type,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-abserr-dating_bycalib_line-pro_arrow_broken.pdf",width=6,height = 3.5)

### RMSE

ggplot(aes(color=Method, y=sqrt(se),x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'se' ,value.var = "se",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-rmse_dating_bycalib_line-pro.pdf",width=6,height = 2.5)

ggplot(aes(color=Method, y=sqrt(se),x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'se' ,value.var = "se",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-rmse_ratevar-dating_bycalib_line-pro.pdf",width=4.5,height = 2.5)

ggplot(aes(color=Method, y=sqrt(se),x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'se' ,value.var = "se",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         shape="none")
ggsave("MV-rmse_dating_bycalib_line-pro_broken.pdf",width=6,height = 3.5)

ggplot(aes(color=Method, y=sqrt(se),x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+ratevar~'se' ,value.var = "se",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE), shape="none")
ggsave("MV-rmse_ratevar_dating_bycalib_line-pro_broken.pdf",width=4.5,height = 3.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'se' ,value.var = "se",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-rmse_ratevar-dating_bycalib_line-pro_arrow.pdf",width=4.5,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'se' ,value.var = "se",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,Branch.Type,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-rmse_ratevar-dating_bycalib_line-pro_arrow_broken.pdf",width=4.5,height = 3.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'se' ,value.var = "se",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-rmse-dating_bycalib_line-pro_arrow.pdf",width=6,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'se' ,value.var = "se",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,Branch.Type,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-rmse-dating_bycalib_line-pro_arrow_broken.pdf",width=6,height = 3.5)

### LOGERROR

ggplot(aes(color=Method, y=log10err,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-log10err_ratevar-dating_bycalib_line-pro.pdf",width=4.5,height = 2.5)

ggplot(aes(color=Method, y=log10err,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-log10err_dating_bycalib_line-pro.pdf",width=6,height = 2.5)

ggplot(aes(color=Method, y=log10err,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_grid(Branch.Type~conditionAD)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-log10err_dating_bycalib_line-pro_broken.pdf",width=6,height = 3.5)

ggplot(aes(color=Method, y=log10err,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+ratevar~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-ratevar-log10err_dating_bycalib_line-pro_broken.pdf",width=4.5,height = 3.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,datingMethod) %>%
  summarise(log10err = mean(log10err)) %>% pivot_wider(names_from = isconcat,values_from = log10err) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-logerr_ratevar-dating_bycalib_line-pro_arrow.pdf",width=4.5,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+ratevar+isconcat~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,Branch.Type,datingMethod) %>%
  summarise(log10err = mean(log10err)) %>% pivot_wider(names_from = isconcat,values_from = log10err) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_grid(Branch.Type~ratevar,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-logerr_ratevar-dating_bycalib_line-pro_arrow_broken.pdf",width=4.5,height = 3.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod) %>%
  summarise(log10err = mean(log10err)) %>% pivot_wider(names_from = isconcat,values_from = log10err) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-logerr-dating_bycalib_line-pro_arrow.pdf",width=6,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Branch.Type+isconcat+conditionAD~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,Branch.Type,isconcat,datingMethod) %>%
  summarise(log10err = mean(log10err)) %>% pivot_wider(names_from = isconcat,values_from = log10err) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_grid(Branch.Type~conditionAD,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-logerr-dating_bycalib_line-pro_arrow_broken.pdf",width=6,height = 3.5)

### MCMCtree figures

m0 = read.csv('mvroot_estgt_dating_unit.csv')
m3 = read.csv('mvroot_estgt_dating_n3_normalized.csv')
m5 = read.csv('mvroot_estgt_dating_n5_normalized.csv')
m0$Calibrations <- 0
m3$Calibrations <- 3
m5$Calibrations <- 5
m <- rbind(m0,m3,m5)
m$outgroup = factor(grepl("outgroup.0", m$Condition))
#m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))#, 'MCMCtree'))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)', 'MCMCtree'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$isconcat = factor(grepl("Concat", m$Method))

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-3, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$bias = m$l.est - m$l.true
m$se = (m$l.est - m$l.true)^2 
ratevar.labs <- c("Low","Medium", "High")
names(ratevar.labs) <- c("5","1.5", "0.15")
m <- m |> mutate(conditionAD = case_when(
  AD >= 0 & AD < 0.25 ~ '[0,0.25)',
  AD >= 0.25 & AD <= 0.5 ~ '[0.25,0.5)',
  AD >= 0.5 & AD <= 0.75 ~ '[0.5,0.75)',
  AD >= 0.75 & AD <= 1 ~ '[0.75,1)',
))

ggplot(aes(color=Method, y=abserr,x=cut(AD,c(0,0.25,0.5,0.75,1))),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Calibrations==5 & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+Branch.Type+ratevar~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Calibrations==5 & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error\n(branch length)")+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary()+
  facet_grid(~Branch.Type)+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  theme(legend.position ='none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        legend.text = element_text(angle=15,size=8),
        axis.text.x = element_text(angle=15,size=8))+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-abserr-ILS-line-no-outgroup_broken.pdf",width=4.5,height = 2.5)

ggplot(aes(color=Method, y=abserr,x=cut(AD,c(0,0.25,0.5,0.75,1)),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Calibrations==5 & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+Branch.Type+ratevar+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Calibrations==5 & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error\n(branch length)")+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary()+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  theme(legend.position ="none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        legend.text = element_text(angle=15,size=8),
        axis.text.x = element_text(angle=15,size=8))+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-abserr-ratevar-ILS-line-no-outgroup.pdf",width=6.5,height = 2.5)

ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+Branch.Type+ratevar+Calibrations+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error\n(branch length)")+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '5'))+
  stat_summary()+
  facet_grid(ratevar~conditionAD,,labeller = labeller(ratevar = ratevar.labs))+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  theme(legend.position ='none', legend.direction = "horizontal",
        legend.box.margin = margin(0))+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-abserr-ILS-mcmctree.pdf",width=5,height = 4.5)

ggplot(aes(color=Method, y=bias,x=as.factor(Calibrations)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+Branch.Type+ratevar+Calibrations+conditionAD~'bias' ,value.var = "bias",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '5'))+
  stat_summary()+
  facet_grid(Branch.Type~conditionAD,,labeller = labeller(ratevar = ratevar.labs))+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position ='none', legend.direction = "horizontal",
        legend.box.margin = margin(0))+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-bias-ILS-mcmctree.pdf",width=6,height = 3.5)

ggplot(aes(color=Method, y=bias,x=as.factor(Calibrations)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+Branch.Type+ratevar+Calibrations~'bias' ,value.var = "bias",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '5'))+
  stat_summary()+
  facet_grid(Branch.Type~ratevar,,labeller = labeller(ratevar = ratevar.labs))+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position ='none', legend.direction = "horizontal",
        legend.box.margin = margin(0))+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-bias-ratevar-mcmctree.pdf",width=4.5,height = 3.5)

ggplot(aes(color=Method, y=bias,x=as.factor(Calibrations)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+Branch.Type+ratevar+Calibrations~'bias' ,value.var = "bias",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '5'))+
  stat_summary()+
  facet_grid(Branch.Type~ratevar,,labeller = labeller(ratevar = ratevar.labs))+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position ='bottom', legend.direction = "horizontal",
        legend.box.margin = margin(0))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("legend-mcmctree.pdf",width=7,height = 3.5)



### MVroot node age

m0=read.csv('mvroot_estgt_dating_unit_node_age.csv')
m3=read.csv('mvroot_estgt_dating_n3_node_age_normalized.csv')
m5=read.csv('mvroot_estgt_dating_n5_node_age_normalized.csv')

m0$Calibrations <- 0
m3$Calibrations <- 3
m5$Calibrations <- 5
m <- rbind(m0,m3,m5)

m$outgroup = factor(grepl("outgroup.0", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)', 'MCMCtree'))
m$ratevar = as.factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$ratevar = factor(m$ratevar)
levels(m$ratevar) = list("0.15" = "High", "1.5" = "Med", "5" = "Low")
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m = m[m$Method!="MCMCtree" & m$Taxon.Type=="internal",]
m$isconcat = factor(grepl("Concat", m$Method))

### Comment out to include negative branch lengths.
m$log10err = log10(m$age.est / m$age.true )
m$abserr = abs(m$age.true - m$age.est)
m$se = (m$age.est - m$age.true)^2 
m$bias = m$age.est - m$age.true
ratevar.labs <- c("Low","Medium", "High")
names(ratevar.labs) <- c("5","1.5", "0.15")
m <- m |> mutate(conditionAD = case_when(
  AD >= 0 & AD < 0.25 ~ '[0,0.25)',
  AD >= 0.25 & AD <= 0.5 ~ '[0.25,0.5)',
  AD >= 0.5 & AD <= 0.75 ~ '[0.5,0.75)',
  AD >= 0.75 & AD <= 1 ~ '[0.75,1)',
))

ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),],
               outgroup+Method+replicate+ratevar+Calibrations+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE & m$Method %in% c('MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'MCMCtree'),], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error\n(node age)")+
  scale_x_discrete(name="Number of Calibrations",labels=c('0', '3', '5'))+
  stat_summary()+
  facet_grid(ratevar~conditionAD,,labeller = labeller(ratevar = ratevar.labs))+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_manual(values=c("#FB9A99","#E31A1C", "#6A3D9A"), name="")+
  theme_classic()+
  theme(legend.position ='none', legend.direction = "horizontal",
        legend.box.margin = margin(0))+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-abserr-ILS-mcmctree_node_age.pdf",width=5,height = 4.5)


dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Taxon.Type+isconcat+conditionAD~'bias' ,value.var = "bias",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,Taxon.Type,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0","3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        legend.text=element_text(size=11))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         shape='none')
ggsave("MV-bias_dating_bycalib-pro_broken-arrow.pdf",width=6,height = 3.5)


ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+isconcat+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-node-age-abserr_dating_bycalib_line-pro.pdf",width=6,height = 2.5)


ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+ratevar+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-node-age-abserr_ratevar-dating_bycalib_line-pro.pdf",width=4.5,height = 2.5)

ggplot(aes(color=Method, y=sqrt(se),x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+isconcat+conditionAD~'se' ,value.var = "se",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-node_age_rmse_dating_bycalib_line-pro.pdf",width=6,height = 2.5)

ggplot(aes(color=Method, y=sqrt(se),x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+ratevar+isconcat~'se' ,value.var = "se",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_discrete(name="Number of Calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_rmse_ratevar-dating_bycalib_line-pro.pdf",width=4.5,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+Taxon.Type+isconcat+conditionAD~'abserr' ,value.var = "abserr",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_abserr-dating_bycalib_line-pro_arrow.pdf",width=6,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+ratevar+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,datingMethod) %>%
  summarise(abserr = mean(abserr)) %>% pivot_wider(names_from = isconcat,values_from = abserr) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Mean absolute error\n(node age)")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_abserr_ratevar-dating_bycalib_line-pro_arrow_main.pdf",width=4.5,height = 2.5)


dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+ratevar+isconcat~'se' ,value.var = "se",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~ratevar,ncol=4,labeller = as_labeller(c('0.15'='High','1.5'='Medium','5'='Low')))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_rmse_ratevar-dating_bycalib_line-pro_arrow.pdf",width=4.5,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+isconcat+conditionAD~'se' ,value.var = "se",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod) %>%
  summarise(rmse = mean(sqrt(se))) %>% pivot_wider(names_from = isconcat,values_from = rmse) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name="Root mean square error (RMSE)")+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(4,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        theme(legend.text=element_text(size=15)))+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_rmse-dating_bycalib_line-pro_arrow.pdf",width=6,height = 2.5)

ggplot(aes(x=as.factor(Calibrations), y=age.est-age.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="Number of calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_wrap(~conditionAD,ncol=4)+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  #scale_shape_manual(guide="none")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         shape="none")
ggsave("MV-node-age-bias_dating_bycalib-pro.pdf",width=6,height = 2.5)


ggplot(aes(x=as.factor(Calibrations), y=age.est-age.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(name="Number of calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_wrap(~ratevar,ncol=4,labeller = labeller(ratevar = ratevar.labs))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),shape="none")
ggsave("MV-node-age-bias-ratevar_dating_bycalib-pro.pdf",width=4.5,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+isconcat+conditionAD~'bias' ,value.var = "bias",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,conditionAD,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_wrap(~conditionAD,ncol=4)+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        legend.text=element_text(size=11))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_bias_dating_bycalib-pro_broken-arrow.pdf",width=6,height = 2.5)

dtemp=merge(
  dcast(data=m[m$outgroup ==TRUE,],
        outgroup+Method+replicate+Calibrations+isconcat+ratevar~'bias' ,value.var = "bias",fun.aggregate = mean),
  dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)) 
dtemp$datingMethod = sub("\\+.*","",dtemp$Method)
dtemp %>% group_by(Calibrations,ratevar,isconcat,datingMethod) %>%
  summarise(bias = mean(bias)) %>% pivot_wider(names_from = isconcat,values_from = bias) %>%
  ggplot(aes(color=datingMethod))+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  facet_wrap(~ratevar,ncol=4,labeller = labeller(ratevar = ratevar.labs))+
  scale_x_continuous(name="Number of Calibrations",breaks = c(1, 2, 3),label = c("0", "3", "5"))+
  geom_segment(aes(yend=`FALSE`,                   y=`TRUE`,
                   x=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8,
                   xend=(as.numeric(factor(Calibrations)))+(as.numeric(factor(datingMethod))-4)/8),
               arrow = arrow(length = unit(5,'pt')), size=0.6)+
  scale_color_manual(values=c("#1F78B4","#E31A1C", "#FF7F00", "#33A02C"), name="")+
  theme_classic()+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  theme(legend.position =  'none', legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0),
        legend.text=element_text(size=11))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE),
         shape='none')
ggsave("MV-node_age_ratevar-bias_dating_bycalib-pro_broken-arrow.pdf",width=4.5,height = 2.5)
