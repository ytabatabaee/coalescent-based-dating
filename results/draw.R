require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern)


s=read.csv('s100_dating_all_branches.csv')
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$l.est = ifelse(s$l.est <=0, 1e-6, s$l.est)
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 

s3=read.csv('s100_dating_n3_all_branches.csv')
s3$Calibrations <- 3
s10=read.csv('s100_dating_n10_all_branches.csv')
s10$Calibrations <- 10
s <- rbind(s3,s10)
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$l.est = ifelse(s$l.est <=0, 1e-6, s$l.est)
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 



ggplot(data=dcast(data=s[!s$Condition=='true gene trees',],Condition+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=Condition,y=log10err,fill=interaction(Method,Calibrations)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  coord_cartesian(ylim=c(0,1))+
  guides(fill=guide_legend(title="Method, #calibrations"))
ggsave("S100-dating_calib_error-logabs-overall.pdf",width=7,height = 5)


ggplot(aes(x=l.true,y=l.est,color=Branch.Type,linetype),
       data=s[!s$Condition=='true gene trees',])+
  facet_grid(Method~Condition)+
  scale_x_continuous(trans="log10",name="True length")+
  scale_y_continuous(trans="log10",name="Estimated length")+
  scale_color_brewer(palette = "Dark2")+
  coord_cartesian(xlim=c(10^-4,0.9),ylim=c(10^-4,0.9))+
  geom_abline(color="grey30",linetype=2)+
  geom_point(alpha=0.1,size=0.5)+
  stat_smooth(se=F,alpha=1,size=0.7,method="glm",formula=y ~ poly(x, 2))+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave("S100-dating-correlation.png",width=8,height =8)

ggplot(aes(x= Condition,y=l.est-l.true,color=interaction(Method,Calibrations), pattern=Calibrations), data=s[!s$Condition=='true gene trees',])+
  facet_wrap(~reorder(Branch.Type,-l.true),ncol=2)+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est-true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  coord_cartesian(ylim=c(-2,2))
ggsave("S100-bias_dating_calib.pdf",width=9,height =  4)

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
ggsave("S100-bias_dating_n3_overall.pdf",width=3.5,height =  4)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est-true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()
  #theme(legend.position = "bottom", legend.direction = "horizontal",
  #      axis.title.x = element_blank())
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_overall.pdf",width=7,height =  4)

ggplot(aes(x=Condition,y=abserr,fill=interaction(Method,Calibrations),pattern=Calibrations),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method, #calibrations"))
  #coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib.pdf",width=7,height = 5)


ggplot(aes(x=Condition,y=sqrt(se),fill=interaction(Method,Calibrations),pattern=Calibrations),
       data=dcast(data=s,Condition+Method+replicate+Branch.Type+Calibrations~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method, #calibrations"))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-rmse-perrep_dating_calib.pdf",width=7,height = 5)


m = read.csv('mvroot_all_branches_estgt_dating.csv')
head(m)
nrow(m)
unique(m$Method)
m$se = (m$l.est - m$l.true)^2 
mvariants = m$Method %in% c("Naive")
names(m) =  c(names(s)[1:4],"AD", "GTEE", names(s)[5:7])
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat+RAxML', 'wLogDate+CASTLES', 'wLogDate+Concat+RAxML'))
#m$ratevar =  unique(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))

### Comment out to include negative branch lengths.
summary(with(m[m$Method =="CASTLES" ,],l.est < 0))
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 

ggplot(aes(color=Method, y=log10err,x=cut(AD,4)),
       data=merge(
         dcast(data=m[!mvariants & m$outgroup ==TRUE,],
               outgroup+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[!mvariants  & m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  coord_cartesian(ylim=c(0.6,2.1))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-logerr-perrep-ILS-bymethod_dating.pdf",width=6.2*0.95,height = 4.3*0.95)

ggplot(aes(color=Method, y=log10err,x=cut(AD,c(0,25,35,50,60,70,85)/100)),
       data=merge(
         dcast(data=m[!mvariants & m$outgroup ==TRUE & !grepl("Pat",m$Method),],
               outgroup+Method+replicate~'log10err' ,
               value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[!mvariants  & m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_color_brewer(palette = "Paired",name="")+
  #geom_boxplot(outlier.size = 0)+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  theme_bw()+
  theme(legend.position =  c(.2,.8), 
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  coord_cartesian(ylim=c(0.8,2.1))
ggsave("MV-logerr-ILS-dating.pdf",width=6,height = 4)

ggplot(aes(x=Method, y=l.true-l.est,color=Method),
       data=m[!mvariants,])+
  scale_y_continuous(trans="identity",name=expression("True" - "Estimated length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("MV-bias_bymethod_pc3_dating.pdf",width=6.4,height = 5)

ggplot(aes(x=ratevar, y=l.true-l.est,color=Method,shape=outgroup),
       data=m[!mvariants,])+
  scale_y_continuous(trans="identity",name=expression("True" - "Estimated length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  #scale_fill_manual(values=c("white","grey70"),name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_x_discrete(labels=c("High","Med","Low"),name="Clock deviation")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-bias_pc3.pdf",width=6.4,height = 5)
