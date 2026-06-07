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

colors_5d = c("#6E3B12","#D3B28E","#CFE2E2","#1FA79B","#17634E")

setwd("D:/R projects/Global_ED")

# load the background data for plotting the world map plot
df_trans = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = df_trans %>% st_drop_geometry()
df_trans[is.na(df_trans$Island),]$Island = 0

shp.glonaf.new = st_read("data/Plants/shp_glonaf_new_eck4.shp")
shp.glonaf.new = shp.glonaf.new %>% st_drop_geometry()
colnames(shp.glonaf.new)[which(colnames(shp.glonaf.new) == 'Region_id')] = 'RegionID'

df_fish = st_read("data/Fishes/data/Basin042017_3119.shp")
df_fish = df_fish %>% st_drop_geometry()
df_fish[is.na(df_fish$Island),]$Island = 0
colnames(df_fish)[which(colnames(df_fish) == 'BasinName')] = 'Basin'
gc()

#2. Regional ED----------------------------
#load regional ED
load('results/primary_results/distances_beta/phy_turn_mammal_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')
load("results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path3_5_mat.rdata")
load("results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path4_7_mat.rdata")
load("results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path6_mat.rdata")

load('results/primary_results/distances_beta/phy_turn_plant_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_plant_native.rdata')
load("results/primary_results/partitioning_possible_ways/plants/turnover_plant_path3_5_mat.rdata")
load("results/primary_results/partitioning_possible_ways/plants/turnover_plant_path4_7_mat.rdata")
load("results/primary_results/partitioning_possible_ways/plants/turnover_plant_path6_mat.rdata")

load('results/primary_results/distances_beta/phy_turn_bird_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_native.rdata')
load("results/primary_results/partitioning_possible_ways/birds/turnover_bird_path3_5_mat.rdata")
load("results/primary_results/partitioning_possible_ways/birds/turnover_bird_path4_7_mat.rdata")
load("results/primary_results/partitioning_possible_ways/birds/turnover_bird_path6_mat.rdata")

load('results/primary_results/distances_beta/phy_turn_fish_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_fish_native.rdata')
load("results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path3_5_mat.rdata")
load("results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path4_7_mat.rdata")
load("results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path6_mat.rdata")


gc()

##2.1 Mammal -----
turnover_mammal_delta_mat = log((phy_turn_mammal_extant+0.001) /
                                  (phy_turn_mammal_native+0.001))

turnover_mammal_delta_mat_path4_7 = log((turnover_mammal_path4_7_mat+0.001) /
                                          (phy_turn_mammal_native+0.001))
turnover_mammal_delta_mat_path4 = turnover_mammal_delta_mat_path4_7
turnover_mammal_delta_mat_path7 = turnover_mammal_delta_mat_path4_7
turnover_mammal_delta_mat_path4[which(turnover_mammal_delta_mat_path4 < 0)] = 0
turnover_mammal_delta_mat_path7[which(turnover_mammal_delta_mat_path7 > 0)] = 0

turnover_mammal_delta_mat_path3_5 = log((turnover_mammal_path3_5_mat+0.001) /
                                          (phy_turn_mammal_native+0.001))
turnover_mammal_delta_mat_path3 = turnover_mammal_delta_mat_path3_5
turnover_mammal_delta_mat_path5 = turnover_mammal_delta_mat_path3_5
turnover_mammal_delta_mat_path3[which(turnover_mammal_delta_mat_path3 < 0)] = 0
turnover_mammal_delta_mat_path5[which(turnover_mammal_delta_mat_path5 > 0)] = 0

turnover_mammal_delta_mat_path6 = log((turnover_mammal_path6_mat+0.001) /
                                        (phy_turn_mammal_native+0.001))


region_ED_mammal_delta = data.frame(RegionID = as.integer(colnames(turnover_mammal_delta_mat)),
                                    real = colMeans(turnover_mammal_delta_mat,
                                                    na.rm = T),
                                    path3 = colMeans(turnover_mammal_delta_mat_path3,
                                                     na.rm = T),
                                    path4 = colMeans(turnover_mammal_delta_mat_path4,
                                                     na.rm = T),
                                    path5 = colMeans(turnover_mammal_delta_mat_path5,
                                                     na.rm = T),
                                    path6 = colMeans(turnover_mammal_delta_mat_path6,
                                                     na.rm = T),
                                    path7 = colMeans(turnover_mammal_delta_mat_path7,
                                                     na.rm = T)
)

region_ED_mammal_delta_2 = region_ED_mammal_delta %>% 
  mutate(path3_diff = log(1 / abs(path3 - real)),
         path4_diff = log(1 / abs(path4 - real)),
         path5_diff = log(1 / abs(path5 - real)),
         path6_diff = log(1 / abs(path6 - real)),
         path7_diff = log(1 / abs(path7 - real))) %>% 
  #mutate(path3_diff = 1/abs(path3 - real),
  #    path4_diff = 1/abs(path4 - real),
  #   path5_diff = 1/abs(path5 - real),
  # path6_diff = 1/abs(path6 - real),
  # path7_diff = 1/abs(path7 - real)) %>% 
  dplyr::select(c(RegionID, real, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-c(RegionID,real)))) %>%
  dplyr::mutate(across(-c(RegionID, real, total), ~ .x / total)) %>%
  dplyr::select(-total) %>% 
  left_join(df_trans[,c('RegionID', 'Island')], by = 'RegionID')

region_ED_mammal_delta_2 = arrange(region_ED_mammal_delta_2,
                                 region_ED_mammal_delta_2$real)

region_ED_mammal_delta_3 = region_ED_mammal_delta_2 %>% 
  mutate(f_real = 1:nrow(region_ED_mammal_delta_2))


region_ED_mammal_delta_3[region_ED_mammal_delta_3$Island == 0,]$Island = 'Mainland'
region_ED_mammal_delta_3[region_ED_mammal_delta_3$Island == 1,]$Island = 'Island'

zero_mammal_mainland = region_ED_mammal_delta_3 %>%
  filter(Island == 'Mainland') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)          

zero_mammal_island = region_ED_mammal_delta_3 %>%
  filter(Island == 'Island') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)    

region_ED_mammal_delta_long = region_ED_mammal_delta_3 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_mammal_delta_long$path = factor(region_ED_mammal_delta_long$path,
                                        levels = c("path3_diff", "path4_diff",
                                                   "path5_diff","path6_diff",
                                                   "path7_diff"))
region_ED_mammal_delta_long$Island = factor(region_ED_mammal_delta_long$Island,
                                          level = c(
                                            'Mainland', 'Island'
                                          ))

vline_data_mammal = data.frame(
  Island = c("Mainland", "Island"),
  x_val  = c(zero_mammal_mainland, zero_mammal_island)                   
)
vline_data_mammal$Island = factor(vline_data_mammal$Island,
                                  level = c(
                                    'Mainland', 'Island'
                                  ))

plot_delta_region_ED_mammal = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = f_real, y = delta_ED, fill = path),
            data = region_ED_mammal_delta_long) +
  geom_vline(data = vline_data_mammal,
             aes(xintercept = x_val), 
             linetype = 2, color = 'grey') +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  scale_y_continuous(limits = c(0,1.05), breaks = c(0, 0.5, 1))+
  labs(x = "order(log(all / native only))",
       y = "Relative path importance") +
  theme(legend.position = "right",
        legend.title = element_blank())


##2.2 Bird ----
turnover_bird_delta_mat = log((phy_turn_bird_extant+0.001) /
                                (phy_turn_bird_native+0.001))

turnover_bird_delta_mat_path4_7 = log((turnover_bird_path4_7_mat+0.001) /
                                        (phy_turn_bird_native+0.001))
turnover_bird_delta_mat_path4 = turnover_bird_delta_mat_path4_7
turnover_bird_delta_mat_path7 = turnover_bird_delta_mat_path4_7
turnover_bird_delta_mat_path4[which(turnover_bird_delta_mat_path4 < 0)] = 0
turnover_bird_delta_mat_path7[which(turnover_bird_delta_mat_path7 > 0)] = 0

turnover_bird_delta_mat_path3_5 = log((turnover_bird_path3_5_mat+0.001) /
                                        (phy_turn_bird_native+0.001))
turnover_bird_delta_mat_path3 = turnover_bird_delta_mat_path3_5
turnover_bird_delta_mat_path5 = turnover_bird_delta_mat_path3_5
turnover_bird_delta_mat_path3[which(turnover_bird_delta_mat_path3 < 0)] = 0
turnover_bird_delta_mat_path5[which(turnover_bird_delta_mat_path5 > 0)] = 0

turnover_bird_delta_mat_path6 = log((turnover_bird_path6_mat+0.001) /
                                      (phy_turn_bird_native+0.001))

region_ED_bird_delta = data.frame(RegionID = as.integer(colnames(turnover_bird_delta_mat)),
                                  real = colMeans(turnover_bird_delta_mat,
                                                  na.rm = T),
                                  path3 = colMeans(turnover_bird_delta_mat_path3,
                                                   na.rm = T),
                                  path4 = colMeans(turnover_bird_delta_mat_path4,
                                                   na.rm = T),
                                  path5 = colMeans(turnover_bird_delta_mat_path5,
                                                   na.rm = T),
                                  path6 = colMeans(turnover_bird_delta_mat_path6,
                                                   na.rm = T),
                                  path7 = colMeans(turnover_bird_delta_mat_path7,
                                                   na.rm = T)
)

region_ED_bird_delta_2 = region_ED_bird_delta %>% 
  mutate(path3_diff = log(1 / abs(path3 - real)),
         path4_diff = log(1 / abs(path4 - real)),
         path5_diff = log(1 / abs(path5 - real)),
         path6_diff = log(1 / abs(path6 - real)),
         path7_diff = log(1 / abs(path7 - real))) %>% 
  #mutate(path3_diff = 1/abs(path3 - real),
  #    path4_diff = 1/abs(path4 - real),
  #   path5_diff = 1/abs(path5 - real),
  # path6_diff = 1/abs(path6 - real),
  # path7_diff = 1/abs(path7 - real)) %>% 
  dplyr::select(c(RegionID, real, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-c(RegionID,real)))) %>%
  dplyr::mutate(across(-c(RegionID, real, total), ~ .x / total)) %>%
  dplyr::select(-total) %>% 
  left_join(df_trans[,c('RegionID', 'Island')], by = 'RegionID')

region_ED_bird_delta_2 = arrange(region_ED_bird_delta_2,
                                   region_ED_bird_delta_2$real)

region_ED_bird_delta_3 = region_ED_bird_delta_2 %>% 
  mutate(f_real = 1:nrow(region_ED_bird_delta_2))


region_ED_bird_delta_3[region_ED_bird_delta_3$Island == 0,]$Island = 'Mainland'
region_ED_bird_delta_3[region_ED_bird_delta_3$Island == 1,]$Island = 'Island'

zero_bird_mainland = region_ED_bird_delta_3 %>%
  filter(Island == 'Mainland') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)          

zero_bird_island = region_ED_bird_delta_3 %>%
  filter(Island == 'Island') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)    

region_ED_bird_delta_long = region_ED_bird_delta_3 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_bird_delta_long$path = factor(region_ED_bird_delta_long$path,
                                          levels = c("path3_diff", "path4_diff",
                                                     "path5_diff","path6_diff",
                                                     "path7_diff"))
region_ED_bird_delta_long$Island = factor(region_ED_bird_delta_long$Island,
                                            level = c(
                                              'Mainland', 'Island'
                                            ))

vline_data_bird = data.frame(
  Island = c("Mainland", "Island"),
  x_val  = c(zero_bird_mainland, zero_bird_island)                   
)
vline_data_bird$Island = factor(vline_data_bird$Island,
                                  level = c(
                                    'Mainland', 'Island'
                                  ))

plot_delta_region_ED_bird = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = f_real, y = delta_ED, fill = path),
            data = region_ED_bird_delta_long) +
  geom_vline(data = vline_data_bird,
             aes(xintercept = x_val), 
             linetype = 2, color = 'grey') +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  scale_y_continuous(limits = c(0,1.05), breaks = c(0, 0.5, 1))+
  labs(x = "order(log(all / native only))",
       y = "Relative path importance") +
  theme(legend.position = "right",
        legend.title = element_blank())



##2.3 Plant -----
turnover_plant_delta_mat = log((phy_turn_plant_extant+0.001) /
                                 (phy_turn_plant_native+0.001))

turnover_plant_delta_mat_path4_7 = log((turnover_plant_path4_7_mat+0.001) /
                                         (phy_turn_plant_native+0.001))
turnover_plant_delta_mat_path4 = turnover_plant_delta_mat_path4_7
turnover_plant_delta_mat_path7 = turnover_plant_delta_mat_path4_7
turnover_plant_delta_mat_path4[which(turnover_plant_delta_mat_path4 < 0)] = 0
turnover_plant_delta_mat_path7[which(turnover_plant_delta_mat_path7 > 0)] = 0

turnover_plant_delta_mat_path3_5 = log((turnover_plant_path3_5_mat+0.001) /
                                         (phy_turn_plant_native+0.001))
turnover_plant_delta_mat_path3 = turnover_plant_delta_mat_path3_5
turnover_plant_delta_mat_path5 = turnover_plant_delta_mat_path3_5
turnover_plant_delta_mat_path3[which(turnover_plant_delta_mat_path3 < 0)] = 0
turnover_plant_delta_mat_path5[which(turnover_plant_delta_mat_path5 > 0)] = 0

turnover_plant_delta_mat_path6 = log((turnover_plant_path6_mat+0.001) /
                                       (phy_turn_plant_native+0.001))

region_ED_plant_delta = data.frame(RegionID = as.integer(colnames(turnover_plant_delta_mat)),
                                   real = colMeans(turnover_plant_delta_mat,
                                                   na.rm = T),
                                   path3 = colMeans(turnover_plant_delta_mat_path3,
                                                    na.rm = T),
                                   path4 = colMeans(turnover_plant_delta_mat_path4,
                                                    na.rm = T),
                                   path5 = colMeans(turnover_plant_delta_mat_path5,
                                                    na.rm = T),
                                   path6 = colMeans(turnover_plant_delta_mat_path6,
                                                    na.rm = T),
                                   path7 = colMeans(turnover_plant_delta_mat_path7,
                                                    na.rm = T)
)

region_ED_plant_delta_2 = region_ED_plant_delta %>% 
  mutate(path3_diff = log(1 / abs(path3 - real)),
         path4_diff = log(1 / abs(path4 - real)),
         path5_diff = log(1 / abs(path5 - real)),
         path6_diff = log(1 / abs(path6 - real)),
         path7_diff = log(1 / abs(path7 - real))) %>% 
  #mutate(path3_diff = 1/abs(path3 - real),
  #    path4_diff = 1/abs(path4 - real),
  #   path5_diff = 1/abs(path5 - real),
  # path6_diff = 1/abs(path6 - real),
  # path7_diff = 1/abs(path7 - real)) %>% 
  dplyr::select(c(RegionID, real, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-c(RegionID,real)))) %>%
  dplyr::mutate(across(-c(RegionID, real, total), ~ .x / total)) %>%
  dplyr::select(-total) %>% 
  left_join(shp.glonaf.new[,c('RegionID', 'Island')], by = 'RegionID')

region_ED_plant_delta_2 = arrange(region_ED_plant_delta_2,
                                 region_ED_plant_delta_2$real)

region_ED_plant_delta_3 = region_ED_plant_delta_2 %>% 
  mutate(f_real = 1:nrow(region_ED_plant_delta_2))


region_ED_plant_delta_3[region_ED_plant_delta_3$Island == 0,]$Island = 'Mainland'
region_ED_plant_delta_3[region_ED_plant_delta_3$Island == 1,]$Island = 'Island'

zero_plant_mainland = region_ED_plant_delta_3 %>%
  filter(Island == 'Mainland') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)          

zero_plant_island = region_ED_plant_delta_3 %>%
  filter(Island == 'Island') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)    

region_ED_plant_delta_long = region_ED_plant_delta_3 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_plant_delta_long$path = factor(region_ED_plant_delta_long$path,
                                        levels = c("path3_diff", "path4_diff",
                                                   "path5_diff","path6_diff",
                                                   "path7_diff"))
region_ED_plant_delta_long$Island = factor(region_ED_plant_delta_long$Island,
                                          level = c(
                                            'Mainland', 'Island'
                                          ))

vline_data_plant = data.frame(
  Island = c("Mainland", "Island"),
  x_val  = c(zero_plant_mainland, zero_plant_island)                   
)
vline_data_plant$Island = factor(vline_data_plant$Island,
                                level = c(
                                  'Mainland', 'Island'
                                ))

plot_delta_region_ED_plant = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = f_real, y = delta_ED, fill = path),
            data = region_ED_plant_delta_long) +
  geom_vline(data = vline_data_plant,
             aes(xintercept = x_val), 
             linetype = 2, color = 'grey') +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  scale_y_continuous(limits = c(0,1.05), breaks = c(0, 0.5, 1))+
  labs(x = "order(log(all / native only))",
       y = "Relative path importance") +
  theme(legend.position = "right",
        legend.title = element_blank())


## 2.4 Fish ----
turnover_fish_delta_mat = log((phy_turn_fish_extant+0.001) /
                                (phy_turn_fish_native+0.001))

turnover_fish_delta_mat_path4_7 = log((turnover_fish_path4_7_mat+0.001) /
                                        (phy_turn_fish_native+0.001))
turnover_fish_delta_mat_path4 = turnover_fish_delta_mat_path4_7
turnover_fish_delta_mat_path7 = turnover_fish_delta_mat_path4_7
turnover_fish_delta_mat_path4[which(turnover_fish_delta_mat_path4 < 0)] = 0
turnover_fish_delta_mat_path7[which(turnover_fish_delta_mat_path7 > 0)] = 0

turnover_fish_delta_mat_path3_5 = log((turnover_fish_path3_5_mat+0.001) /
                                        (phy_turn_fish_native+0.001))
turnover_fish_delta_mat_path3 = turnover_fish_delta_mat_path3_5
turnover_fish_delta_mat_path5 = turnover_fish_delta_mat_path3_5
turnover_fish_delta_mat_path3[which(turnover_fish_delta_mat_path3 < 0)] = 0
turnover_fish_delta_mat_path5[which(turnover_fish_delta_mat_path5 > 0)] = 0

turnover_fish_delta_mat_path6 = log((turnover_fish_path6_mat+0.001) /
                                      (phy_turn_fish_native+0.001))

region_ED_fish_delta = data.frame(Basin = colnames(turnover_fish_delta_mat),
                                  real = colMeans(turnover_fish_delta_mat,
                                                  na.rm = T),
                                  path3 = colMeans(turnover_fish_delta_mat_path3,
                                                   na.rm = T),
                                  path4 = colMeans(turnover_fish_delta_mat_path4,
                                                   na.rm = T),
                                  path5 = colMeans(turnover_fish_delta_mat_path5,
                                                   na.rm = T),
                                  path6 = colMeans(turnover_fish_delta_mat_path6,
                                                   na.rm = T),
                                  path7 = colMeans(turnover_fish_delta_mat_path7,
                                                   na.rm = T)
)

region_ED_fish_delta_2 = region_ED_fish_delta %>% 
  mutate(path3_diff = log(1 / abs(path3 - real)),
         path4_diff = log(1 / abs(path4 - real)),
         path5_diff = log(1 / abs(path5 - real)),
         path6_diff = log(1 / abs(path6 - real)),
         path7_diff = log(1 / abs(path7 - real))) %>% 
  #mutate(path3_diff = 1/abs(path3 - real),
  #    path4_diff = 1/abs(path4 - real),
  #   path5_diff = 1/abs(path5 - real),
  # path6_diff = 1/abs(path6 - real),
  # path7_diff = 1/abs(path7 - real)) %>% 
  dplyr::select(c(Basin, real, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-c(Basin,real)))) %>%
  dplyr::mutate(across(-c(Basin, real, total), ~ .x / total)) %>%
  dplyr::select(-total) %>% 
  left_join(df_fish[,c('Basin', 'Island')], by = 'Basin')

region_ED_fish_delta_2 = arrange(region_ED_fish_delta_2,
                                  region_ED_fish_delta_2$real)

region_ED_fish_delta_3 = region_ED_fish_delta_2 %>% 
  mutate(f_real = 1:nrow(region_ED_fish_delta_2))


region_ED_fish_delta_3[region_ED_fish_delta_3$Island == 0,]$Island = 'Mainland'
region_ED_fish_delta_3[region_ED_fish_delta_3$Island == 1,]$Island = 'Island'

zero_fish_mainland = region_ED_fish_delta_3 %>%
  filter(Island == 'Mainland') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)          

zero_fish_island = region_ED_fish_delta_3 %>%
  filter(Island == 'Island') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)    

region_ED_fish_delta_long = region_ED_fish_delta_3 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_fish_delta_long$path = factor(region_ED_fish_delta_long$path,
                                         levels = c("path3_diff", "path4_diff",
                                                    "path5_diff","path6_diff",
                                                    "path7_diff"))
region_ED_fish_delta_long$Island = factor(region_ED_fish_delta_long$Island,
                                           level = c(
                                             'Mainland', 'Island'
                                           ))

vline_data_fish = data.frame(
  Island = c("Mainland", "Island"),
  x_val  = c(zero_fish_mainland, zero_fish_island)                   
)
vline_data_fish$Island = factor(vline_data_fish$Island,
                                 level = c(
                                   'Mainland', 'Island'
                                 ))

plot_delta_region_ED_fish = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = f_real, y = delta_ED, fill = path),
            data = region_ED_fish_delta_long) +
  geom_vline(data = vline_data_fish,
             aes(xintercept = x_val), 
             linetype = 2, color = 'grey') +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0, 0.5, 1))+
  labs(x = "order(log(all / native only))",
       y = "Relative path importance") +
  theme(legend.position = "right",
        legend.title = element_blank())


#3. Exporting all taxon column figures----------------------------
library(ggpubr)
plot_all_regionED_mammal_pathways = ggarrange(plot_delta_region_ED_plant, plot_delta_region_ED_bird,
                                          plot_delta_region_ED_mammal, plot_delta_region_ED_fish,
                                          ncol = 2, nrow = 2,
                                          labels = c('a', 'b', 'c', 'd'),
                                          common.legend = T)

#4. Export figures for paper----
library(devEMF)
library(cowplot)


emf('figures/plot_all_regionED_mammal_pathways.emf',
    width = 25, height = 20, coordDPI = 600*0.8, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_all_regionED_mammal_pathways
dev.off() #turn off device and finalize file
