rm(list=ls())
gc()
library("dplyr")
library("tidyverse")
library("lme4")
library("lmerTest")
library("scico")

figure.path <- 'figures/'
load('data_new_3_GIFT.Rdata')
df.MPD.each.Sp=df.MPD.each.Sp.487

#Fig. 3 ########################
#GLMM#
# get data

df <- as_tibble(df.MPD.each.Sp)

df <- df %>% 
  mutate(LAT_abs_scaled = scale(abs(LAT)),
         MPD_scaled = scale(MPD))
# the region column is numeric, change it to character
df <- df %>% 
  mutate(region = as.character(newOBJIDsic))
set.seed(1) 
# the model was performed using a 70GB node, and taking about 16 hours to finish
mod.glmm <- glmer(ys ~ LAT_abs_scaled * MPD_scaled + (LAT_abs_scaled * MPD_scaled|species) + (1|region), 
                  data = df,
                  family = binomial(link = "cloglog"),
                  nAGQ=0, 
                  control=glmerControl(optimizer = "nloptwrap"))

# make data for predictions using the coefficient table
range(df$LAT_abs_scaled) 
range(df$MPD_scaled) 

LAT_abs_scaled <- seq(-1.90, 2.6, 0.5)
MPD_scaled <- seq(-7.7, 2.9, 0.5)

new_dt <- crossing(LAT_abs_scaled, MPD_scaled)
df_eff <- round(summary(mod.glmm)$coefficients, digits = 3) %>% as_tibble()

# manually calculate the predictions without considering the random effect from species
inverse.cloglog.f <- function(x) 1-exp(-exp(x))
new_dt <- new_dt %>% 
  # the linear part of the GLMM
  mutate(mod.glmm_X = LAT_abs_scaled * df_eff$Estimate[2] + MPD_scaled * df_eff$Estimate[3] + LAT_abs_scaled * MPD_scaled * df_eff$Estimate[4]) %>% 
  # transform it to the logistic equation
  mutate(predicted = inverse.cloglog.f(mod.glmm_X))

# organize data
new_dt$LAT_abs=(new_dt$LAT_abs_scaled * sd(df$LAT))+mean(df$LAT)
new_dt$MPD=(new_dt$MPD_scaled * sd(df$MPD))+mean(df$MPD)

p.a = ggplot(data = new_dt, aes(x =MPD , y = predicted, color = LAT_abs, group = as.factor(LAT_abs))) +
  geom_line() +
  theme_bw() +
  scale_color_scico(palette = 'imola')+
  labs(x='MPD (Mya)', y = "Predicted naturalization success") +
  ylim(0.18,1)+
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_line(size=1),
        axis.title.x = element_text(size=11), 
        axis.title.y = element_text(size=11), 
        axis.text.x = element_text(size = 9,colour ="black"),
        axis.text.y = element_text(size = 9,colour ="black"),
        legend.position = 'top', legend.box = "horizontal") 

pdf(paste(figure.path,"Fig.3.pdf",sep=""), height=3.6,width=3.5)
p.a
dev.off()
