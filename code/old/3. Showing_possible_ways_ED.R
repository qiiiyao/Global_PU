#Qi Yao
#Email:qiyao.eco@gmail.com

#0 loading----------------------------
rm(list = ls())
library(reshape2)
library(plyr)
library(betareg)
library(MASS)
#library(gdm)
library(vegan)
library(rgeos)
library(RColorBrewer)
library(scico)
library(sjPlot)
library(dplyr)
library(purrr)
library(nlme)
library(spdep)
library(stringr)
library(ggplot2)

colors_4d = c('#91989f', '#efbb24', '#1b813e', '#398fb7')
colors_10d = c('#91989f',
               colorRampPalette(c('#FFF9B0', '#FFE066', '#FFC300', '#FF9E00', '#FF7A00', '#E85C00'))(6),
               colorRampPalette(c('#9BD4E4', '#398FB7', '#004A78'))(3))

setwd("D:/R projects/Global_ED")
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_clear = df
rm(df)




f.hex = function(dat) {
  dat %>%
    ggplot(aes(x = native,
               y = native.and.naturalized)) +
    geom_hex(#bins = 30,
             aes(fill = stat(count))) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "black",
      size = 1,
      linetype = "dashed"
    ) +
    coord_fixed(ratio = 1,
                xlim = range(c(dat$native, dat$native.and.naturalized)),
                ylim = range(c(dat$native, dat$native.and.naturalized)))
}


#1. Species ED----------------------------
#load species level ED
load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')
load('results/primary_results/distances_beta/ED_mammal_all.rdata')
load('results/primary_results/distances_beta/ED_mammal_native.rdata')

load('results/primary_results/distances_beta/ED_bird_extant.rdata')
load('results/primary_results/distances_beta/ED_bird_native.rdata')

load('results/primary_results/distances_beta/ED_plant_extant.rdata')
load('results/primary_results/distances_beta/ED_plant_native.rdata')
load("results/primary_results/predictors_GIFT/pd_plant_native.RDATA")

load('results/primary_results/distances_beta/ED_fish_extant.rdata')
load('results/primary_results/distances_beta/ED_fish_native.rdata')
load('results/primary_results/distances_beta/LCBD_fish_extant.rdata')
load("results/primary_results/predictors_basins/pd_fish_native.RDATA")

### Mammal
ED_mammal_extant_1 = cbind(RegionID = colnames(phy_turn_mammal_native),
                           ED_mammal_all$df)
colnames(ED_mammal_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_mammal_native_1 = cbind(RegionID = colnames(phy_turn_mammal_native),
                           ED_mammal_native$df)
colnames(ED_mammal_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_mammal = ED_mammal_native_1 %>% left_join(ED_mammal_extant_1,
                                                   by = 'RegionID')


ggthemr::ggthemr(palette = "fresh", layout = "clean")
sp_ED_mammal = Delta_ED_mammal %>%
   dplyr::select(mean_native_ED,
          mean_all_ED) %>%
  dplyr::rename(
    native = mean_native_ED,
     native.and.naturalized = mean_all_ED
     ) %>%
   f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(sp_ED[native])) +
  ylab(expression(sp_ED[native+naturalized]))+ 
  annotate("text", x = 6, y = 8.5,
           label = "Pathway 1")+
  annotate("text", x = 8.5, y = 6,
           label = "Pathway 2")+
  #ggtitle('Mammals') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5)) 

### Bird
ED_bird_extant_1 = cbind(RegionID = df_clear$RegionID,
                         ED_bird_extant$df)
colnames(ED_bird_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_bird_native_1 = cbind(RegionID = df_clear$RegionID,
                         ED_bird_native$df)
colnames(ED_bird_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_bird = ED_bird_native_1 %>% left_join(ED_bird_extant_1,
                                               by = 'RegionID')


range(c(Delta_ED_bird$mean_native_ED, Delta_ED_bird$mean_all_ED))
sp_ED_bird = Delta_ED_bird %>%
  dplyr::select(mean_native_ED,
                mean_all_ED) %>%
  dplyr::rename(
    native = mean_native_ED,
    native.and.naturalized = mean_all_ED
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(sp_ED[native])) +
  ylab(expression(sp_ED[native+naturalized]))+ 
  annotate("text", x = 10, y = 15,
           label = "Pathway 1")+
  annotate("text", x = 15, y = 10,
           label = "Pathway 2")+ 
  #ggtitle('Birds') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5)) 



### Plants
ED_plant_extant_1 = cbind(Region_id = sort(unique(pd_plant_native$Region_id)),
                          ED_plant_extant$df)
colnames(ED_plant_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_plant_native_1 = cbind(Region_id = sort(unique(pd_plant_native$Region_id)),
                          ED_plant_native$df)
colnames(ED_plant_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_plant = ED_plant_native_1 %>% left_join(ED_plant_extant_1,
                                                 by = 'Region_id')


sp_ED_plant = Delta_ED_plant %>%
  dplyr::select(mean_native_ED,
                mean_all_ED) %>%
  dplyr::rename(
    native = mean_native_ED,
    native.and.naturalized = mean_all_ED
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(sp_ED[native])) +
  ylab(expression(sp_ED[native+naturalized]))+ 
  annotate("text", x = 5, y = 15,
           label = "Pathway 1")+
  annotate("text", x = 15, y = 5,
           label = "Pathway 2")+ 
  #ggtitle('Plants')+
  theme(plot.title = element_text(face = "bold", size = 12, 
                            hjust = 0.5)) 


### Fishes
turnover_fish_extant_mat = LCBD_fish_extant$beta_mat

ED_fish_extant_1 = cbind(Basin.name = sort(colnames(turnover_fish_extant_mat)),
                         ED_fish_extant$df)
colnames(ED_fish_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_fish_native_1 = cbind(Basin.name = sort(unique(pd_fish_native$X1.Basin.Name)),
                         ED_fish_native$df)
colnames(ED_fish_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_fish = ED_fish_native_1 %>% left_join(ED_fish_extant_1,
                                               by = 'Basin.name')


sp_ED_fish = Delta_ED_fish %>%
  dplyr::select(mean_native_ED,
                mean_all_ED) %>%
  dplyr::rename(
    native = mean_native_ED,
    native.and.naturalized = mean_all_ED
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(sp_ED[native])) +
  ylab(expression(sp_ED[native+naturalized]))+ 
  annotate("text", x = 20, y = 75,
           label = "Pathway 1")+
  annotate("text", x = 75, y = 20,
           label = "Pathway 2")+ 
  #ggtitle('Fishes') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5)) 








#2. Regional ED----------------------------
#load regional ED
load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')
load('results/primary_results/distances_beta/phy_turn_mammal_all.rdata')
load("results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path3_5_mat.rdata")
load("results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path4_7_mat.rdata")
load("results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path6_mat.rdata")

load('results/primary_results/distances_beta/LCBD_bird_extant.rdata')
load('results/primary_results/distances_beta/LCBD_bird_native.rdata')


##2.1 Mammal -----
turnover_mammal_extant_mat = phy_turn_mammal_all
turnover_mammal_native_mat = phy_turn_mammal_native

idx = which(upper.tri(turnover_mammal_native_mat),
           arr.ind = TRUE)

df_turnover_mammal_native = data.frame(
  Region_pair = paste(rownames(turnover_mammal_native_mat)[idx[,1]],
                      colnames(turnover_mammal_native_mat)[idx[,2]],
                      sep = "."),
  native_turnover = turnover_mammal_native_mat[idx])

df_turnover_mammal_extant = data.frame(
  Region_pair = paste(rownames(turnover_mammal_extant_mat)[idx[,1]],
                      colnames(turnover_mammal_extant_mat)[idx[,2]],
                      sep = "."),
  extant_turnover = turnover_mammal_extant_mat[idx])
Delta_turnover_mammal_extant = df_turnover_mammal_native %>% 
  left_join(df_turnover_mammal_extant,
            by = 'Region_pair')

df_turnover_mammal_path6 = data.frame(
  Region_pair = paste(rownames(turnover_mammal_path6_mat)[idx[,1]],
    colnames(turnover_mammal_path6_mat)[idx[,2]],
    sep = "."),
  extant_turnover = turnover_mammal_path6_mat[idx])
Delta_turnover_mammal_path6 = df_turnover_mammal_native %>% 
  left_join(df_turnover_mammal_path6,
            by = 'Region_pair')

df_turnover_mammal_path3_5 = data.frame(
  Region_pair = paste(rownames(turnover_mammal_path3_5_mat)[idx[,1]],
                      colnames(turnover_mammal_path3_5_mat)[idx[,2]],
                      sep = "."),
  extant_turnover = turnover_mammal_path3_5_mat[idx])
Delta_turnover_mammal_path3_5 = df_turnover_mammal_native %>% 
  left_join(df_turnover_mammal_path3_5,
            by = 'Region_pair')

df_turnover_mammal_path4_7 = data.frame(
  Region_pair = paste(rownames(turnover_mammal_path4_7_mat)[idx[,1]],
                      colnames(turnover_mammal_path4_7_mat)[idx[,2]],
                      sep = "."),
  extant_turnover = turnover_mammal_path4_7_mat[idx])
Delta_turnover_mammal_path4_7 = df_turnover_mammal_native %>% 
  left_join(df_turnover_mammal_path4_7,
            by = 'Region_pair')


ggthemr::ggthemr(palette = "fresh", layout = "clean")
regional_ED_mammal_extant = Delta_turnover_mammal_extant %>%
  dplyr::select(native_turnover,
                extant_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = extant_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  scale_x_continuous(breaks = c(0,1))+
  scale_y_continuous(breaks = c(0,1))+
  #ggtitle('Mammals') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_mammal_path6 = Delta_turnover_mammal_path6 %>%
  dplyr::select(native_turnover,
                path6_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path6_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 6")+
  scale_x_continuous(breaks = c(0,1))+
  scale_y_continuous(breaks = c(0,1))+
  #ggtitle('Mammals') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


regional_ED_mammal_path3_5 = Delta_turnover_mammal_path3_5 %>%
  dplyr::select(native_turnover,
                path3_5_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path3_5_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  scale_x_continuous(breaks = c(0,1))+
  scale_y_continuous(breaks = c(0,1))+
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 3")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 5")+
  #ggtitle('Mammals') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_mammal_path4_7 = Delta_turnover_mammal_path4_7 %>%
  dplyr::select(native_turnover,
                path4_7_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path4_7_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  scale_x_continuous(breaks = c(0,1))+
  scale_y_continuous(breaks = c(0,1))+
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 4")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 7")+
  #ggtitle('Mammals') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


##2.2 Bird ----
turnover_bird_extant_mat = phy_turn_bird_all
turnover_bird_native_mat = phy_turn_bird_native

idx = which(upper.tri(turnover_bird_native_mat),
            arr.ind = TRUE)

df_turnover_bird_native = data.frame(
  Region_pair = paste(rownames(turnover_bird_native_mat)[idx[,1]],
                      colnames(turnover_bird_native_mat)[idx[,2]],
                      sep = "."),
  native_turnover = turnover_bird_native_mat[idx])

df_turnover_bird_extant = data.frame(
  Region_pair = paste(rownames(turnover_bird_extant_mat)[idx[,1]],
                      colnames(turnover_bird_extant_mat)[idx[,2]],
                      sep = "."),
  extant_turnover = turnover_bird_extant_mat[idx])
Delta_turnover_bird_extant = df_turnover_bird_native %>% 
  left_join(df_turnover_bird_extant,
            by = 'Region_pair')

df_turnover_bird_path6 = data.frame(
  Region_pair = paste(rownames(turnover_bird_path6_mat)[idx[,1]],
                      colnames(turnover_bird_path6_mat)[idx[,2]],
                      sep = "."),
  path6_turnover = turnover_bird_path6_mat[idx])
Delta_turnover_bird_path6 = df_turnover_bird_native %>% 
  left_join(df_turnover_bird_path6,
            by = 'Region_pair')

df_turnover_bird_path3_5 = data.frame(
  Region_pair = paste(rownames(turnover_bird_path3_5_mat)[idx[,1]],
                      colnames(turnover_bird_path3_5_mat)[idx[,2]],
                      sep = "."),
  path3_5_turnover = turnover_bird_path3_5_mat[idx])
Delta_turnover_bird_path3_5 = df_turnover_bird_native %>% 
  left_join(df_turnover_bird_path3_5,
            by = 'Region_pair')

df_turnover_bird_path4_7 = data.frame(
  Region_pair = paste(rownames(turnover_bird_path4_7_mat)[idx[,1]],
                      colnames(turnover_bird_path4_7_mat)[idx[,2]],
                      sep = "."),
  path4_7_turnover = turnover_bird_path4_7_mat[idx])
Delta_turnover_bird_path4_7 = df_turnover_bird_native %>% 
  left_join(df_turnover_bird_path4_7,
            by = 'Region_pair')


ggthemr::ggthemr(palette = "fresh", layout = "clean")
regional_ED_bird_extant = Delta_turnover_bird_extant %>%
  dplyr::select(native_turnover,
                extant_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = extant_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_bird_path6 = Delta_turnover_bird_path6 %>%
  dplyr::select(native_turnover,
                path6_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path6_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 6")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


regional_ED_bird_path3_5 = Delta_turnover_bird_path3_5 %>%
  dplyr::select(native_turnover,
                path3_5_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path3_5_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 3")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 5")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_bird_path4_7 = Delta_turnover_bird_path4_7 %>%
  dplyr::select(native_turnover,
                path4_7_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path4_7_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 4")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 7")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


##2.3 Plant -----
turnover_plant_extant_mat = phy_turn_plant_all
turnover_plant_native_mat = phy_turn_plant_native

idx = which(upper.tri(turnover_plant_native_mat),
            arr.ind = TRUE)

df_turnover_plant_native = data.frame(
  Region_pair = paste(rownames(turnover_plant_native_mat)[idx[,1]],
                      colnames(turnover_plant_native_mat)[idx[,2]],
                      sep = "."),
  native_turnover = turnover_plant_native_mat[idx])

df_turnover_plant_extant = data.frame(
  Region_pair = paste(rownames(turnover_plant_extant_mat)[idx[,1]],
                      colnames(turnover_plant_extant_mat)[idx[,2]],
                      sep = "."),
  extant_turnover = turnover_plant_extant_mat[idx])
Delta_turnover_plant_extant = df_turnover_plant_native %>% 
  left_join(df_turnover_plant_extant,
            by = 'Region_pair')

df_turnover_plant_path6 = data.frame(
  Region_pair = paste(rownames(turnover_plant_path6_mat)[idx[,1]],
                      colnames(turnover_plant_path6_mat)[idx[,2]],
                      sep = "."),
  path6_turnover = turnover_plant_path6_mat[idx])
Delta_turnover_plant_path6 = df_turnover_plant_native %>% 
  left_join(df_turnover_plant_path6,
            by = 'Region_pair')

df_turnover_plant_path3_5 = data.frame(
  Region_pair = paste(rownames(turnover_plant_path3_5_mat)[idx[,1]],
                      colnames(turnover_plant_path3_5_mat)[idx[,2]],
                      sep = "."),
  path3_5_turnover = turnover_plant_path3_5_mat[idx])
Delta_turnover_plant_path3_5 = df_turnover_plant_native %>% 
  left_join(df_turnover_plant_path3_5,
            by = 'Region_pair')

df_turnover_plant_path4_7 = data.frame(
  Region_pair = paste(rownames(turnover_plant_path4_7_mat)[idx[,1]],
                      colnames(turnover_plant_path4_7_mat)[idx[,2]],
                      sep = "."),
  path4_7_turnover = turnover_plant_path4_7_mat[idx])
Delta_turnover_plant_path4_7 = df_turnover_plant_native %>% 
  left_join(df_turnover_plant_path4_7,
            by = 'Region_pair')


ggthemr::ggthemr(palette = "fresh", layout = "clean")
regional_ED_plant_extant = Delta_turnover_plant_extant %>%
  dplyr::select(native_turnover,
                extant_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = extant_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_plant_path6 = Delta_turnover_plant_path6 %>%
  dplyr::select(native_turnover,
                path6_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path6_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 6")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


regional_ED_plant_path3_5 = Delta_turnover_plant_path3_5 %>%
  dplyr::select(native_turnover,
                path3_5_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path3_5_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 3")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 5")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_plant_path4_7 = Delta_turnover_plant_path4_7 %>%
  dplyr::select(native_turnover,
                path4_7_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path4_7_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 4")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 7")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


## 2.4 Fish ----
turnover_fish_extant_mat = phy_turn_fish_all
turnover_fish_native_mat = phy_turn_fish_native

idx = which(upper.tri(turnover_fish_native_mat),
            arr.ind = TRUE)

df_turnover_fish_native = data.frame(
  Region_pair = paste(rownames(turnover_fish_native_mat)[idx[,1]],
                      colnames(turnover_fish_native_mat)[idx[,2]],
                      sep = "."),
  native_turnover = turnover_fish_native_mat[idx])

df_turnover_fish_extant = data.frame(
  Region_pair = paste(rownames(turnover_fish_extant_mat)[idx[,1]],
                      colnames(turnover_fish_extant_mat)[idx[,2]],
                      sep = "."),
  extant_turnover = turnover_fish_extant_mat[idx])
Delta_turnover_fish_extant = df_turnover_fish_native %>% 
  left_join(df_turnover_fish_extant,
            by = 'Region_pair')

df_turnover_fish_path6 = data.frame(
  Region_pair = paste(rownames(turnover_fish_path6_mat)[idx[,1]],
                      colnames(turnover_fish_path6_mat)[idx[,2]],
                      sep = "."),
  path6_turnover = turnover_fish_path6_mat[idx])
Delta_turnover_fish_path6 = df_turnover_fish_native %>% 
  left_join(df_turnover_fish_path6,
            by = 'Region_pair')

df_turnover_fish_path3_5 = data.frame(
  Region_pair = paste(rownames(turnover_fish_path3_5_mat)[idx[,1]],
                      colnames(turnover_fish_path3_5_mat)[idx[,2]],
                      sep = "."),
  path3_5_turnover = turnover_fish_path3_5_mat[idx])
Delta_turnover_fish_path3_5 = df_turnover_fish_native %>% 
  left_join(df_turnover_fish_path3_5,
            by = 'Region_pair')

df_turnover_fish_path4_7 = data.frame(
  Region_pair = paste(rownames(turnover_fish_path4_7_mat)[idx[,1]],
                      colnames(turnover_fish_path4_7_mat)[idx[,2]],
                      sep = "."),
  path4_7_turnover = turnover_fish_path4_7_mat[idx])
Delta_turnover_fish_path4_7 = df_turnover_fish_native %>% 
  left_join(df_turnover_fish_path4_7,
            by = 'Region_pair')


ggthemr::ggthemr(palette = "fresh", layout = "clean")
regional_ED_fish_extant = Delta_turnover_fish_extant %>%
  dplyr::select(native_turnover,
                extant_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = extant_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_fish_path6 = Delta_turnover_fish_path6 %>%
  dplyr::select(native_turnover,
                path6_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path6_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 6")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))


regional_ED_fish_path3_5 = Delta_turnover_fish_path3_5 %>%
  dplyr::select(native_turnover,
                path3_5_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path3_5_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 3")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 5")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))

regional_ED_fish_path4_7 = Delta_turnover_fish_path4_7 %>%
  dplyr::select(native_turnover,
                path4_7_turnover) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = path4_7_turnover
  ) %>%
  f.hex() +
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  annotate("text", x = 0.2, y = 0.8,
           label = "Pathway 4")+
  annotate("text", x = 0.8, y = 0.2,
           label = "Pathway 7")+
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))



#3. Exporting all taxon column figures----------------------------
library(ggpubr)
plot_spED_pathways = ggarrange(sp_ED_plant, sp_ED_bird,
                               sp_ED_mammal, sp_ED_fish,
                               ncol = 2, nrow = 2,
                               labels = c('a', 'b', 'c', 'd'))

plot_regionED_mammal_pathways = ggarrange(regional_ED_mammal_extant, regional_ED_mammal_path3_5,
                                   regional_ED_mammal_path4_7, regional_ED_mammal_path6,
                                   nrow = 4,
                                   common.legend = T, 
                                   legend = 'bottom')
null_plot = ggplot(NULL)
plot_ED_pathways = ggarrange(sp_ED_plant, sp_ED_bird, sp_ED_mammal, sp_ED_fish,
                             null_plot, null_plot, regional_ED_mammal_path3_5, null_plot,
                             null_plot, null_plot, regional_ED_mammal_path4_7, null_plot,
                             null_plot, null_plot, regional_ED_mammal_path6, null_plot,
                             ncol = 4, nrow = 4)



Delta_turnover_mammal_extant$status = 'Extant'
Delta_turnover_mammal_path3_5$status = 'Path: 3&5'
Delta_turnover_mammal_path4_7$status = 'Path: 4&7'
Delta_turnover_mammal_path6$status = 'Path: 6'

colnames(Delta_turnover_mammal_extant)
colnames(Delta_turnover_mammal_path3_5)
Delta_turnover_mammal_all = rbind(Delta_turnover_mammal_extant, 
                     Delta_turnover_mammal_path3_5,
                     Delta_turnover_mammal_path4_7,
                     Delta_turnover_mammal_path6) %>%
  as_tibble()

Delta_turnover_mammal_all$status = factor(Delta_turnover_mammal_all$status,
                           levels = c('Extant',
                                      'Path: 3&5',
                                      'Path: 4&7',
                                      'Path: 6'))

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_regionED_mammal_pathways = Delta_turnover_mammal_all %>%
  dplyr::select(native_turnover,
                extant_turnover,
                status) %>%
  dplyr::rename(
    native = native_turnover,
    native.and.naturalized = extant_turnover
  ) %>%
  f.hex() +
  facet_wrap(vars(status), ncol = 4)+
  scale_fill_scico(palette = "bilbao",
                   trans = 'log10',
                   direction = -1,
                   begin = 0, 
                   end = 0.7) +
  xlab(expression(Phy_turnover[native])) +
  ylab(expression(Phy_turnover[native+naturalized]))+ 
  scale_x_continuous(breaks = c(0,1))+
  scale_y_continuous(breaks = c(0,1))+
  #ggtitle('Mammals') +
  theme(plot.title = element_text(face = "bold", size = 12, 
                                  hjust = 0.5))
plot_regionED_mammal_pathways

#4. Export figures for paper----
library(devEMF)
library(cowplot)

emf('figures/plot_spED_pathways.emf',
    width = 20*0.85, height = 20*0.85, coordDPI = 600*0.8, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_spED_pathways
dev.off() #turn off device and finalize file

emf('figures/plot_regionED_mammal_pathways.emf',
    width = 20, height = 40, coordDPI = 600*0.8, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_regionED_mammal_pathways
dev.off() #turn off device and finalize file
