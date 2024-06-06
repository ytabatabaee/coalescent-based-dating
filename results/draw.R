require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr);require(ggpattern)


s=read.csv('s100_dating_unit.csv')
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")
s$Method = factor(s$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))
s$l.est = ifelse(s$l.est <=0, 1e-6, s$l.est)
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 

s3=read.csv('s100_dating_n3.csv')
s3$Calibrations <- 3
s10=read.csv('s100_dating_n10.csv')
s10$Calibrations <- 10
s <- rbind(s3,s10)
s$Method = factor(s$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))
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



ggplot(data=dcast(data=s[!s$Condition=='true gene trees',],Condition+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=Condition,y=log10err,color=interaction(Method,Calibrations)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  stat_summary()+
  stat_summary(aes(group=interaction(Method,Calibrations)),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method, #calibrations"))
ggsave("S100-dating_error-logabs-calib.pdf",width=7,height = 5)

ggplot(data=dcast(data=s[!s$Condition=='true gene trees',],Condition+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))), aes(x=Condition,y=log10err,color=Method))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary()+
  #stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())
ggsave("S100-dating_error-logabs-overall.pdf",width=7,height = 5)

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
       data=s[!s$Condition=='true gene trees',])+
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
  facet_wrap(~Calibrations,ncol=2)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib.pdf",width=8,height =  4)

ggplot(aes(x= Condition,y=l.est-l.true,color=Method), data=s[!s$Condition=='true gene trees',])+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  scale_y_continuous(name=expression("Bias (est. - true length)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  facet_grid(Branch.Type~Calibrations)+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
#coord_cartesian(xlim=c(1,5),clip="off",ylim=c(-0.06,0.125))
ggsave("S100-bias_dating_calib_broken.pdf",width=8,height =  8)


ggplot(aes(x=Condition,y=abserr,color=interaction(Method,Calibrations),pattern=Calibrations),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  stat_summary()+
  stat_summary(aes(group=interaction(Method,Calibrations)),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method, #calibrations"))
  #coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib.pdf",width=7,height = 5)


ggplot(aes(x=Condition,y=abserr,color=Method),
       data=dcast(data=s, Condition+Method+replicate+Branch.Type+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  facet_wrap(~Calibrations,ncol=2)+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
  guides(fill=guide_legend(title="Method"))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-error-perrep_dating_calib-line.pdf",width=8,height = 4.5)


ggplot(aes(x=Condition,y=sqrt(se),color=interaction(Method,Calibrations),pattern=Calibrations),
       data=dcast(data=s,Condition+Method+replicate+Branch.Type+Calibrations~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  #stat_summary(position = position_dodge(width=0.86))+
  stat_summary()+
  stat_summary(aes(group=interaction(Method,Calibrations)),geom="line")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired",name="")+
  theme_bw()+
  guides(fill=guide_legend(title="Method, #calibrations"))
#coord_cartesian(ylim=c(0,0.2),xlim=c(1,5),clip="off")
ggsave("S100-rmse-perrep_dating_calib.pdf",width=7,height = 5)


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


m = read.csv('mvroot_estgt_dating.csv')
head(m)
nrow(m)
unique(m$Method)
m$se = (m$l.est - m$l.true)^2 
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 


m3 = read.csv('mvroot_estgt_dating_n3.csv')
m5 = read.csv('mvroot_estgt_dating_n5.csv')
m3$Calibrations <- 3
m5$Calibrations <- 5
m <- rbind(m3,m5)
m$se = (m$l.est - m$l.true)^2 
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 


ggplot(aes(color=Method, y=log10err,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate+Calibrations~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  facet_wrap(~Calibrations,ncol=2)+
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
ggsave("MV-logerr-perrep-ILS-bymethod_dating_calib.pdf",width=9,height = 4.5)

ggplot(aes(x=ratevar, y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
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
ggsave("MV-bias_dating_calib.pdf",width=8,height = 5)

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
ggsave("MV-bias_dating_ILS.pdf",width=7,height = 5)

ggplot(aes(x=cut(AD,4), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_wrap(~Calibrations)+
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
ggsave("MV-bias_dating_ILS_calib_broken.pdf",width=8,height = 4.5)

ggplot(aes(x=cut(AD,5), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  scale_x_discrete(name="True gene tree discordance (ILS)")+
  facet_wrap(~Branch.Type,ncol=2)+
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
ggsave("MV-bias_dating_ILS_broken.pdf",width=10,height = 5)


ggplot(aes(color=Method, y=abserr,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate+Calibrations~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  facet_wrap(~Calibrations,ncol=2)+
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
ggsave("MV-abserr-perrep-ILS-bymethod_dating_calib.pdf",width=9,height = 4.5)

ggplot(aes(color=Method, y=abserr,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  #geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary()+
  stat_summary(aes(group=Method),geom="line")+
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
ggsave("MV-abserr-ILS-line.pdf",width=6.2*0.95,height = 4.3*0.95)


m = read.csv('mvroot_estgt_dating_n3.csv')
head(m)
nrow(m)
unique(m$Method)
m$se = (m$l.est - m$l.true)^2 
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m$Method = factor(m$Method, levels=c('LSD+CASTLES', 'LSD+Concat(RAxML)', 'wLogDate+CASTLES', 'wLogDate+Concat(RAxML)', 'MD-Cat+CASTLES', 'MD-Cat+Concat(RAxML)', 'TreePL+CASTLES', 'TreePL+Concat(RAxML)'))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))

### Comment out to include negative branch lengths.
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 

ggplot(aes(color=Method, y=abserr,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate+Branch.Type~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~Branch.Type,ncol=2)+
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
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-abserr-n3-perrep-ILS-bymethod_dating.pdf",width=6.2*0.95,height = 4.3*0.95)


ggplot(aes(color=Method, y=log10err,x=cut(AD,4)),
       data=merge(
         dcast(data=m[m$outgroup ==FALSE,],
               outgroup+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[m$outgroup ==FALSE,], outgroup+replicate~'AD' ,value.var = "AD",fun.aggregate = mean)))+
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
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-logerr-perrep-ILS-bymethod_dating_n3.pdf",width=6.2*0.95,height = 4.3*0.95)

ggplot(aes(x=cut(AD,5), y=l.est-l.true,color=Method),
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
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-bias_n3_dating.pdf",width=5,height = 5)


ggplot(aes(x=cut(AD,5), y=l.est-l.true,color=Method),
       data=m[m$outgroup ==FALSE,])+
  scale_y_continuous(trans="identity",name=expression("Est." - "true length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  facet_wrap(~Branch.Type,ncol=2)+
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
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("MV-bias_n3_dating_broken.pdf",width=8,height = 5)
