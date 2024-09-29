require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern)



s0=read.csv('s100_dating_unit_castles_pro.csv')
s0$Calibrations <- 0
s3=read.csv('s100_dating_normalized_n3_castles_pro.csv')
s3$Calibrations <- 3
s10=read.csv('s100_dating_normalized_n10_castles_pro.csv')
s10$Calibrations <- 10
s <- rbind(s0,s3,s10)
s$Method = factor(s$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))
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


ggplot(data=dcast(data=s[!s$Condition=='true gene trees',],Condition+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=Condition,y=log10err,color=Method))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
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
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_error-logabs-overall_pro.pdf",width=7,height = 5)

ggplot(data=dcast(data=s[!s$Condition=='true gene trees',],Condition+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=Condition,y=log10err,color=Method))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~Calibrations,ncol=2)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_error-logabs-calib.pdf",width=8,height = 4.5)


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


ggplot(aes(x= Condition,y=l.est-l.true,color=Method, pattern=Calibrations), data=s[!s$Condition=='true gene trees',])+
  facet_wrap(~reorder(Branch.Type,-l.true),ncol=2)+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est-true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-bias_dating_calib-2.pdf",width=9,height =  4)


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

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
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
ggsave("S100-bias_dating_calib-2.pdf",width=8,height =  4)


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

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
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
ggsave("S100-bias_dating_calib_broken-line.pdf",width=9,height = 3)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
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
ggsave("S100-bias_dating_calib.pdf",width=9,height = 3.5)


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
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib.pdf",width=9,height =4)

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
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib_broken.pdf",width=9,height =5)

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


m3=read.csv('mvroot_estgt_dating_n3_root_unfixed_tmrca.csv')
m5=read.csv('mvroot_estgt_dating_n5_root_unfixed_tmrca.csv')
m3$Calibrations <- 3
m5$Calibrations <- 5
m <- rbind(m3,m5)
m$outgroup = factor(grepl("outgroup.0", m$Condition))
#m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))#, 'MCMCtree'))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))#, 'MCMCtree'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$isconcat = factor(grepl("Concat", m$Method))
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 


ggplot(aes(color=Method, shape=isconcat,y=log10err,x=Calibrations),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_wrap(~cut(AD,4))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-logerr-tmrca_dating_calib_normalized_bycalib.pdf",width=12,height = 4.5)



ggplot(aes(x=cut(AD,4), y=(l.est)/l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==FALSE & m$Calibrations==5,])+
  scale_y_continuous(trans="identity",name=expression("Estimated / true tMRCA"))+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  coord_cartesian(ylim=c(0,3))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-tmrca-bias_dating_calib_root-unfixed.pdf",width=6,height = 3)

ggplot(aes(x=ratevar, y=(l.est)/l.true,color=Method),
       data=m[m$outgroup ==FALSE & m$Calibrations==5,])+
  scale_y_continuous(trans="identity",name=expression("Estimated / true tMRCA"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  coord_cartesian(ylim=c(0,3))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 1)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-tmrca-ratevar-bias_dating_calib_root-unfixed.pdf",width=6,height = 4)



m0 = read.csv('mvroot_estgt_dating_castlespro.csv')
m3 = read.csv('mvroot_estgt_dating_n3_normalized_castlespro.csv')
m5 = read.csv('mvroot_estgt_dating_n5_normalized_castlespro.csv')
m0$Calibrations <- 0
m3$Calibrations <- 3
m5$Calibrations <- 5
m <- rbind(m0,m3,m5)
m$outgroup = factor(grepl("outgroup.0", m$Condition))
#m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))#, 'MCMCtree'))
m$Method = factor(m$Method, levels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'))#, 'MCMCtree'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$isconcat = factor(grepl("Concat", m$Method))

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 




ggplot(aes(color=Method, y=log10err,x=Calibrations),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_wrap(~cut(AD,4))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  #stat_summary()+
  #stat_summary(aes(group=Method),geom="line")+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-logerr-perrep-ILS-bymethod_dating_calib_normalized_bycalib.pdf",width=12,height = 4.5)

ggplot(aes(x=ratevar, y=l.est/l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  facet_wrap(~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  #scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  #scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-bias_dating_calib_root-unfixed.pdf",width=8,height = 5)

ggplot(aes(x=ratevar, y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  #scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  #scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=4, byrow=TRUE),
         fill=guide_legend(nrow=4, byrow=TRUE))
ggsave("MV-bias_dating.pdf",width=5.5,height = 5)

ggplot(aes(x=cut(AD,5), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  #facet_wrap(~Branch.Type,ncol=2)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  #scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  #scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        # axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-bias_dating_ILS_root_unfixed.pdf",width=7,height = 5)

ggplot(aes(x=cut(AD,4), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_grid(Branch.Type~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  #scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  #scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        # axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-bias_dating_ILS_calib_broken_normalized.pdf",width=12,height = 8)

ggplot(aes(x=as.factor(Calibrations), y=l.est-l.true,color=Method,shape=isconcat),
       data=m[m$outgroup ==TRUE & m$Method %in% c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),])+
  scale_y_continuous(trans="identity",name=expression("Estimated" - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="Number of calibrations")+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  facet_grid(Branch.Type~cut(AD,c(0,0.25,0.5,0.75,1)))+
  theme_classic()+
  theme(legend.position =  c(0.5,0.85), legend.direction = "horizontal",
        # axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  #scale_color_brewer(palette = "Paired",name="")+
  scale_colour_manual(name="", labels=c('LSD+CASTLES-Pro', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES-Pro', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES-Pro', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)'),values=brewer.pal(n=8,"Paired"))+
  #scale_shape_manual(guide="none")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         shape="none")
ggsave("MV-bias_dating_bycalib-pro_broken.pdf",width=10,height = 4)

ggplot(aes(x=cut(AD,4), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==TRUE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_wrap(~Calibrations,ncol=2)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
       # axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-bias_dating_ILS_broken_normalized-pro.pdf",width=10,height = 5)


ggplot(aes(color=Method, y=abserr,x=as.factor(Calibrations),shape=isconcat),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~cut(AD,4),ncol=4)+
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
         dcast(data=m,
               outgroup+Method+replicate+Calibrations+Branch.Type+isconcat~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m, outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_grid(outgroup~cut(AD,4))+
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
ggsave("MV-abserr_dating_bycalib_line-pro-outgroup.pdf",width=6,height = 2.5)


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

ggplot(aes(color=Method, y=log10err,x=as.factor(Calibrations)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+Calibrations+Branch.Type~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  facet_wrap(~cut(AD,4),ncol=4)+
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

ggplot(aes(color=Method, y=abserr,x=Calibrations),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  facet_wrap(~cut(AD,4))+
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
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-abserr-perrep-ILS-bymethod_dating_calib_normalized_bycalib.pdf",width=12,height = 4.5)


ggplot(aes(color=Method, y=abserr,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate+Calibrations~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  facet_wrap(~Calibrations)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=8))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-abserr-ILS-line.pdf",width=10,height = 4)


m = read.csv('mvroot_estgt_dating_n5_normalized.csv')
m$Calibrations <- 5
m$se = (m$l.est - m$l.true)^2 
m$outgroup = factor(grepl("outgroup.0", m$Condition))
#m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)', 'MCMCtree'))
m$Method = factor(m$Method, levels=c('MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'MCMCtree'))
m$ratevar = as.factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$ratevar = factor(m$ratevar)
levels(m$ratevar) = list("0.15" = "High", "1.5" = "Med", "5" = "Low")

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 

ggplot(aes(color=Method, y=abserr,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==TRUE,],
               outgroup+Method+replicate+ratevar~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==TRUE,], outgroup+replicate+ratevar~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  #facet_wrap(~ratevar)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=8))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-abserr-ILS-line-no-outgroup.pdf",width=4,height = 4)

ggplot(aes(x=cut(AD,5), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==TRUE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_wrap(ratevar~Branch.Type)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  #scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  #scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        # axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-bias_dating_ILS_no_outgroup.pdf",width=10,height = 5)


