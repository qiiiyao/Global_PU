rm(list=ls())
gc()
pacman::p_load(char = c("tidyverse", "magrittr", "ggpubr",
                        "raster", "rgdal", "viridis",
                        "scico", "ggpmisc", "sf",
                        "rnaturalearth","gridExtra",
                        "dplyr","visreg","lmerTest","cowplot"))

figure.path <- 'figures/'
load('Sample_code/data_new_3_GIFT.Rdata')
mydata$newOBJIDsic=as.character(mydata$newOBJIDsic)

#Fig.2a,b###############################

# read the dataset that include the data for Fig 1
load("Sample_code/data.for.shp.plots.4.Rdata")
sf.polygons.analysed=merge(sf.polygons.analysed, mydata, by.x='OBJIDsic',by.y='newOBJIDsic',all.y=T)


## bin MPD and delta.MPD.full to groups
breaks_MPD <- c(min(sf.polygons.analysed$MPD), seq(235, 250, 2.5), max(sf.polygons.analysed$MPD))
breaks_delta.MPD.full <- c(min(sf.polygons.analysed$delta.MPD.full), seq(-8, 4, 2), max(sf.polygons.analysed$delta.MPD.full))
labels_MPD <- round(breaks_MPD[-1], 2)
labels_MPD_full <- round(breaks_MPD, 2)
labels_delta.MPD.full <- round(breaks_delta.MPD.full[-1], 2)
labels_delta.MPD.full_full <- round(breaks_delta.MPD.full, 2)
sf.polygons.analysed %<>%
  mutate(MPD.group = cut(MPD, breaks = breaks_MPD, include.lowest = TRUE, labels = labels_MPD)) %>%
  mutate(delta.MPD.full.group = cut(delta.MPD.full, breaks = breaks_delta.MPD.full, include.lowest = TRUE, labels = labels_delta.MPD.full))

# define color gradients
colors1 <- scico::scico(n=8, palette = "lajolla") 
colors2 <- scico::scico(n=10, palette = "vik")[1:8] 
#plot(1:8, col = colors2, pch = 19, cex = 5)

# define a function to make legend of color gradients
legend.func <- function(mycolors, mylabels) {
  group <- rep("cc", 8)
  condition <- letters[1:8]
  value <- rep(1, 8)
  df.legend <- data.frame(group, condition, value)
  mycolors.corrected <- rev(mycolors)
  ggplot(df.legend, aes(fill = condition, y = value, x = group)) +
    geom_bar(position = "stack", stat = "identity", color = "white") +
    scale_fill_manual(values = mycolors.corrected) +
    theme_classic() +
    theme(
      legend.position = "none", aspect.ratio = 0.03,
      axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(size = 7, color = "black")
    ) +
    scale_y_continuous(breaks = 0:8, labels = mylabels) +
    coord_flip() +
    xlab("")
}

# make the map for MPD
scale_manual_MPD <- list(
  scale_fill_manual(values = colors1, drop = FALSE),
  scale_color_manual(values = colors1, drop = FALSE)
)

p.MPD <-
  ggplot() +
  geom_sf(data = countries, color = NA, fill = "gray") +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = sf.polygons.analysed, aes(fill = MPD.group), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12)
  ) +
  scale_manual_MPD
p.MPD

legend.MPD <- legend.func(mycolors = colors1, mylabels = labels_MPD_full) +
  ggtitle("MPD (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
p.MPD <- ggplotGrob(p.MPD)
legend.MPD <- ggplotGrob(legend.MPD)
p.a <- arrangeGrob(p.MPD, legend.MPD, ncol = 1, layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(p.a)
# make the map for delta.MPD.full
scale_manual_delta.MPD.full <- list(
  scale_fill_manual(values = colors2, drop = FALSE),
  scale_color_manual(values = colors2, drop = FALSE)
)

p.delta.MPD.full <-
  ggplot() +
  geom_sf(data = countries, color = NA, fill = "gray") +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = sf.polygons.analysed, aes(fill = delta.MPD.full.group), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12)
  ) +
  scale_manual_delta.MPD.full
p.delta.MPD.full
legend.delta.MPD.full <- legend.func(mycolors = colors2, mylabels = labels_delta.MPD.full_full) +
  ggtitle(expression(Delta*"MPD (Mya)")) +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
p.delta.MPD.full <- ggplotGrob(p.delta.MPD.full)
legend.delta.MPD.full <- ggplotGrob(legend.delta.MPD.full)
p.b <- arrangeGrob(p.delta.MPD.full, legend.delta.MPD.full, ncol = 1, layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(p.b)


#Fig.2c,d###############################
#MPD~LAT#
p.c=ggplot(mydata, aes(x = LAT,y = MPD)) +
  geom_point(shape=19,alpha = 1,size=0.7,color="gray90") +
  #scale_color_manual(values = c("#601200","#013E7D","gray60"), drop = FALSE,
  #                  breaks = c("Positive","Negative","NA"), labels = c("Significantly positive","Significantly negative","Non-significant"))+
  geom_smooth(method = "lm", formula = y ~ x, se = T,size=1.5,color="gray23")+
  stat_regline_equation(aes(x = LAT,y = MPD,label = ..adj.rr.label..),
           label.x = 5,label.y = 228,size=2.3)+
  theme_classic() +
  labs(x = "Absolute latitude (°)",y = "MPD (Mya)")+
  theme(axis.title =element_text(size=9),
        axis.text = element_text(size = 7,colour = 'black'),
        legend.position = c(0.7,0.9),
        legend.background = element_blank(),
        legend.text = element_text(size=7),
        legend.key.width = unit(1,'line'),
        legend.key.height = unit(1,'line'),
        legend.title = element_blank()
  )
p.c

#delta MPD~LAT#
cbPalette <- c("gray65","#E69F00","#CC6600","#336633","dodgerblue4","slateblue4")#
p.d=
  ggplot(mydata) +
  geom_point(aes(x = LAT, y = delta.MPD.full),shape=19,alpha = 1,size=0.7,color="gray90") +
  stat_smooth(aes(x = LAT, y = delta.MPD.useful, colour ="Useful",fill="Useful"),
              method = "lm", formula =y~x, se = T,lwd=1.3,alpha=0.35) + 
  stat_smooth(aes(x = LAT, y = delta.MPD.full.nonna, colour ="Full.nonna",fill = "Full.nonna"),
              method = "lm", formula = y~x,se = T,lwd=1.5,alpha=0.35) + 
  stat_smooth(aes(x = LAT, y = delta.MPD.sdm, colour ="MAXENT",fill="MAXENT"),
              method = "lm", formula =y~x,se = T,lwd=1.5,alpha=0.45)+
  stat_smooth(aes(x = LAT, y = delta.MPD.eco, colour ="TDWG",fill="TDWG"),
              method = "lm", formula = y~x, se = T,lwd=1.5,alpha=0.35) + 
  stat_smooth(aes(x = LAT, y = delta.MPD.eco.sdm, colour ="Eco.Sdm",fill="Eco.Sdm"),
              method = "lm", formula =y~x, se = T,lwd=1.3,alpha=0.45) + 
  stat_smooth(aes(x = LAT, y = delta.MPD.full, colour ="Full",fill = "Full"),
              method = "lm", formula = y~x,se = T,lwd=1.5,alpha=0.4) + 
  scale_color_manual(name = "",breaks = c("Full.nonna","Useful","Full", "TDWG", "MAXENT","Eco.Sdm"),
                     values = c("Full"=cbPalette[1], "Full.nonna"=cbPalette[2], "TDWG"=cbPalette[4],"MAXENT"=cbPalette[5],"Eco.Sdm"=cbPalette[6],"Useful"=cbPalette[3]))+
  scale_fill_manual(name = "",breaks = c("Full.nonna","Useful","Full", "TDWG", "MAXENT","Eco.Sdm"),
                    values = c("Full"=cbPalette[1], "Full.nonna"=cbPalette[2], "TDWG"=cbPalette[4],"MAXENT"=cbPalette[5],"Eco.Sdm"=cbPalette[6],"Useful"=cbPalette[3])) +
  stat_regline_equation(aes(x = LAT,y = delta.MPD.eco.sdm,label = paste("Climate","continent","nat.",..adj.rr.label.., sep = "*\"  \"*")),
           label.x = 8,label.y = -16.275,size=2.3)+
  stat_regline_equation(aes(x = LAT,y = delta.MPD.sdm,label = paste("Climate","nat.",..adj.rr.label.., sep = "*\"  \"*")),
           label.x = 8,label.y = -14.64,size=2.3)+
  stat_regline_equation(aes(x = LAT,y = delta.MPD.eco,label = paste("Continent","nat.",..adj.rr.label.., sep = "*\"  \"*")),
           label.x = 8,label.y = -13.005,size=2.3)+
  stat_regline_equation(aes(x = LAT,y = delta.MPD.full,label = paste("Global","nat.",..adj.rr.label.., sep = "*\"  \"*")),
           label.x = 8,label.y = -11.37,size=2.3)+
  stat_regline_equation(aes(x = LAT,y = delta.MPD.useful,label = paste("Econ.","use", "flora",..adj.rr.label.., sep = "*\"  \"*")),
           label.x = 8,label.y = -9.735,size=2.3)+
  stat_regline_equation(aes(x = LAT,y = delta.MPD.full.nonna,label = paste("Global","flora",..adj.rr.label.., sep = "*\"  \"*")),
           label.x = 8,label.y = -8.11,size=2.3)+
  theme_classic() +
  labs(x = "Absolute latitude (°)",y = expression(Delta*"MPD (Mya)"))+
  theme(axis.title =element_text(size=9),
        axis.text = element_text(size = 7,colour = 'black'),
        legend.title=element_blank(),
        legend.text =element_blank(),
        legend.key.width = unit(0.8,'line'),
        legend.key.height = unit(0.56,'line'),
        legend.position = c(0.09,0.21),
        legend.background=element_blank())

p.d

#Fig. 2e,f####

n=3
cols=scale_color_manual( values = scico(n = n, palette = "bamako"),name="HM",labels=c("0.016 (10% quantile)","0.221 (50% quantile)","0.529 (90% quantile)"))
interaction.hm.lat.mpd = lmer(MPD ~ LAT*HM + (1 |TDWG),data=mydata)
v.hm.lat.mpd = visreg(interaction.hm.lat.mpd,"LAT", by=c("HM"), breaks=n, 
                      overlay = TRUE, partial = FALSE, rug = FALSE,plot=T)
p.e=ggplot(v.hm.lat.mpd$fit, aes(LAT, visregFit, color=factor(HM))) +
  cols+
  geom_line(size=1.5) +
  labs(x="Absolute latitude (°)", y = "MPD (Mya)") +
  theme_classic()+
  theme(axis.ticks.x=element_line(size=1),
        axis.title=element_text(size=9), 
        axis.text= element_text(size = 7,colour ="black"),
        legend.text =element_text(size=7),
        legend.title=element_text(size=7),
        legend.title.align = 0.3,
        legend.key.width = unit(0.8,'line'),
        legend.key.height = unit(0.8,'line'),
        legend.justification=c(0,0), legend.position=c(0.0,0.0),
        legend.background=element_blank())
p.e
interaction.hm.lat.deltampd.full = lmer(delta.MPD.full ~ LAT*HM + (1 |TDWG),data=mydata)
v.hm.lat.deltampd.full = visreg(interaction.hm.lat.deltampd.full,"LAT", by=c("HM"), breaks=n, 
                                overlay = TRUE, partial = FALSE, rug = FALSE,plot=T)
p.f=ggplot(v.hm.lat.deltampd.full$fit, aes(LAT, visregFit, color=factor(HM))) +
  cols+
  geom_line(size=1.5) +
  labs(x="Absolute latitude (°)", y = expression(Delta*"MPD (Mya)")) +
  theme_classic()+
  theme(axis.ticks.x=element_line(size=1),
        axis.title=element_text(size=9), 
        axis.text= element_text(size = 7,colour ="black"),
        legend.position = "none")
p.f
#merge all plots####
pdf(paste(figure.path,"/Fig.2.pdf",sep=""), height=6,width=11)
ggdraw() +
  draw_plot(p.a, x = -0.03, y = .5, width = 0.58, height = 0.49) +
  draw_plot(p.b, x = -0.03, y = 0, width = 0.58, height = 0.49) +
  draw_plot(p.c, x = 0.5, y = 0.52, width = 0.25, height = 0.45) +
  draw_plot(p.d, x = 0.5, y = 0.02, width = 0.25, height = 0.45) +
  draw_plot(p.e, x = 0.75, y = 0.52, width = 0.25, height = 0.45) +
  draw_plot(p.f, x = 0.75, y = 0.02, width = 0.25, height = 0.45) +
  draw_plot_label(label = c("a","b","c","d","e","f"), size = 12,
                  x = c(0.05,0.05,0.5, 0.5, 0.75, 0.75), y = c(1,0.5,1,0.5,1,0.5))
dev.off()


