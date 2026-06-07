rm(list=ls())
gc()
pacman::p_load(char = c("dplyr","visreg","lmerTest","effectsize","cowplot","ggplot2","scico"))

figure.path <- 'figures/'
load('data_new_3_GIFT.Rdata')

mydata$newOBJIDsic=as.character(mydata$newOBJIDsic)


#Fig.a###############################
library(effectsize)
data.lme=mydata
lm.mpd = lmer(MPD ~ Climate.PC1*HM+ Climate.PC2*HM + (1 |TDWG),data=data.lme)

lm.deltampd.full.nonna.wt = lmer(delta.MPD.full.nonna ~ Climate.PC1*HM+ Climate.PC2*HM + 
                                   (1 |TDWG),weights = 1/(data.lme$sd.MPD.full.nonna)^2,data=data.lme)

lm.deltampd.full.wt = lmer(delta.MPD.full ~ Climate.PC1*HM+ Climate.PC2*HM + 
                             (1 |TDWG),weights = 1/(data.lme$sd.MPD.full)^2,data=data.lme)

lm.deltampd.eco.wt = lmer(delta.MPD.eco ~ Climate.PC1*HM+ Climate.PC2*HM + 
                            (1 |TDWG),weights = 1/(data.lme$sd.MPD.eco)^2,data=data.lme)

lm.deltampd.sdm.wt = lmer(delta.MPD.sdm ~ Climate.PC1*HM+ Climate.PC2*HM + 
                            (1 |TDWG),weights = 1/(data.lme$sd.MPD.sdm)^2,data=data.lme)
data.temp=data.lme[!is.na(data.lme$delta.MPD.eco.sdm),]
lm.deltampd.eco.sdm.wt = lmer(delta.MPD.eco.sdm ~ Climate.PC1*HM+ Climate.PC2*HM + 
                            (1 |TDWG),weights = 1/(data.temp$sd.MPD.eco.sdm)^2,data=data.temp)
lm.deltampd.useful.wt = lmer(delta.MPD.useful ~ Climate.PC1*HM+ Climate.PC2*HM + 
                            (1 |TDWG),weights = 1/(data.lme$sd.MPD.useful)^2,data=data.lme)
##plot mixed effect###
mix.res=rbind(effectsize(lm.mpd),effectsize(lm.deltampd.full.wt),effectsize(lm.deltampd.full.nonna.wt),
              effectsize(lm.deltampd.eco.wt),effectsize(lm.deltampd.sdm.wt),
              effectsize(lm.deltampd.eco.sdm.wt),effectsize(lm.deltampd.useful.wt))
mix.res$sig=19
mix.res$sig[mix.res$CI_low*mix.res$CI_high < 0]=1
mix.res$var=c(rep("MPD",6),rep("dGlobal nat.",6),rep("fGlobal nonnative flora",6),rep("Continent nat.",6),
              rep("bClimate nat.",6),rep("aClimate continent nat.",6),rep("Econ. use flora",6))
mix.res=mix.res[-which(mix.res$Parameter == "(Intercept)"),]
mix.res$Parameter=gsub(pattern = "HM:Climate.PC2",replacement = "Climate.PC2:HM",mix.res$Parameter)


pd <- position_dodge(width = 0.75)

p.a = ggplot(mix.res, aes(x=Parameter, y=Std_Coefficient,colour = var, group = var))+
  geom_hline(yintercept=0, linetype="dashed",color="gray",size=0.8) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high,colour = var),
                width = 0.5,size=0.7, position = pd) +
  geom_point(aes(colour = var),lwd=2.4, shape = mix.res$sig,position = pd) +
  
  scale_x_discrete(limits=c('Climate.PC2:HM','Climate.PC1:HM','HM','Climate.PC2','Climate.PC1'),labels=c('PCPrec:HM','PCTemp:HM','HM','PCPrec','PCTemp'))+
  scale_color_manual(name = "",breaks = c("MPD","fGlobal nonnative flora","Econ. use flora","dGlobal nat.", "Continent nat.", "bClimate nat.","aClimate continent nat."),
                     values = c("black","#E69F00","#CC6600","gray65","#336633","dodgerblue4","slateblue4"))+
  theme_classic() +
  labs(x=NULL, y = "Effect size (95% confidence interval)") +
  theme(panel.grid =element_blank(),
        plot.title = element_text(size = 9,face = "bold",hjust = 0),
        axis.ticks.x=element_line(size=1),
        axis.title.x =element_text(size=9), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 7,colour ="black"),
        axis.text.y = element_text(size = 9,colour ="black"),
        legend.title=element_blank(),
        legend.text =element_text(size=8),
        legend.key.width = unit(1,'line'),
        legend.key.height = unit(0.8,'line'),
        legend.position = c(0.8, 0.15),
        legend.background=element_blank()) +
  coord_flip()
p.a

#Fig.b###############################

results=mydata
interaction.MPD = lmer(MPD ~ Climate.PC1*HM+ Climate.PC2*HM + (1 |TDWG),data=results)

n=3
cols=scale_color_manual( values = scico(n = n, palette = "bamako"),name="HM",labels=c("0.016 (10% quantile)","0.221 (50% quantile)","0.529 (90% quantile)"))

p.1.mpd = visreg(interaction.MPD,"Climate.PC1", by=c("HM"), breaks=n, 
                 overlay = TRUE, partial = FALSE, rug = FALSE,
                 plot=T)

p.b=ggplot(p.1.mpd$fit, aes(Climate.PC1, visregFit, color=factor(HM))) +
  cols+
  geom_line(size=1.4) +
  labs(x=expression("PC"[Temp]), y = "MPD (Mya)") +
  theme_classic()+
  theme(axis.ticks.x=element_line(size=1),
        axis.title=element_text(size=9), 
        axis.text= element_text(size = 7,colour ="black"),
        legend.text =element_text(size=7),
        legend.title=element_text(size=7),
        legend.title.align = 0.5,
        legend.key.width = unit(0.8,'line'),
        legend.key.height = unit(0.8,'line'),
        legend.justification=c(0,0), legend.position=c(0.0,0.6),
        legend.background=element_blank())

p.2.mpd = visreg(interaction.MPD,"Climate.PC2", by=c("HM"), breaks=n,
                 overlay = TRUE, partial = FALSE, rug = FALSE,
                 plot=FALSE)

p.c=ggplot(p.2.mpd$fit, aes(Climate.PC2, visregFit, color=factor(HM))) +
  cols+
  geom_line(size=1.4) +
  labs(x=expression("PC"[Prec]), y = "MPD (Mya)") +
  theme_classic()+
  theme(axis.ticks.x=element_line(size=1),
        axis.title=element_text(size=9), 
        axis.text= element_text(size = 7,colour ="black"),
        legend.position = "none")

#delatMPD
interaction.delta.MPD.full = lmer(delta.MPD.full ~ Climate.PC1*HM+ Climate.PC2*HM + (1 |TDWG),weights = 1/(results$sd.MPD.full)^2,data=results)
p.1 = visreg(interaction.delta.MPD.full,"Climate.PC1", by=c("HM"), breaks=n, 
             overlay = TRUE, partial = FALSE, rug = FALSE,
             plot=F)

p.d=ggplot(p.1$fit, aes(Climate.PC1, visregFit, color=factor(HM))) +
  cols+
  geom_line(size=1.4) +
  labs(x=expression("PC"[Temp]), y = expression(Delta*"MPD (Mya)")) +
  theme_classic()+
  theme(axis.ticks.x=element_line(size=1),
        axis.title=element_text(size=9), 
        axis.text= element_text(size = 7,colour ="black"),
        legend.position = "none")

p.2 = visreg(interaction.delta.MPD.full,"Climate.PC2", by=c("HM"), breaks=n,
             overlay = TRUE, partial = FALSE, rug = FALSE,
             plot=FALSE)

p.e=ggplot(p.2$fit, aes(Climate.PC2, visregFit, color=factor(HM))) +
  cols+
  geom_line(size=1.4) +
  labs(x=expression("PC"[Prec]), y = expression(Delta*"MPD (Mya)")) +
  theme_classic()+
  theme(axis.ticks.x=element_line(size=1),
        axis.title=element_text(size=9), 
        axis.text= element_text(size = 7,colour ="black"),
        legend.position = "none")

#merge 5 plots####
pdf(paste(figure.path,"/Fig.4.pdf",sep=""), height=5,width=10)

ggdraw() +
  draw_plot(p.a, x = 0, y = 0, width = 0.48, height = 0.98) +
  draw_plot(p.b, x = 0.5, y = 0.5, width = 0.25, height = 0.48) +
  draw_plot(p.c, x = 0.5, y = 0.0, width = 0.25, height = 0.48) +
  draw_plot(p.d, x = 0.75, y = 0.5, width = 0.25, height = 0.48) +
  draw_plot(p.e, x = 0.75, y = 0.0, width = 0.25, height = 0.48) +
  
  draw_plot_label(label = c("a","b","c","d","e"), size = 12,
                  x = c(0,0.5, 0.5, 0.75, 0.75), y = c(1,1,0.5,1,0.5))
dev.off()
