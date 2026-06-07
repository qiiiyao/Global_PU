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

colors_5d = rev(c("#6E3B12","#D3B28E","#CFE2E2","#1FA79B","#17634E"))

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
load("results/primary_results/partitioning_possible_ways/mammals/delta_pu_mammal_dat.rdata")

load('results/primary_results/distances_beta/phy_turn_plant_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_plant_native.rdata')
load("results/primary_results/partitioning_possible_ways/plants/turnover_plant_path3_5_mat.rdata")

load('results/primary_results/distances_beta/phy_turn_bird_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_native.rdata')
load("results/primary_results/partitioning_possible_ways/birds/delta_pu_bird_dat.rdata")

load('results/primary_results/distances_beta/phy_turn_fish_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_fish_native.rdata')
load("results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path3_5_mat.rdata")



gc()

##2.1 Mammal -----
delta_pu_mammal_dat$region = as.integer(delta_pu_mammal_dat$region)
pu_mammal_extant = data.frame(RegionID = colnames(phy_turn_mammal_extant),
                                 PU = colMeans(phy_turn_mammal_extant))

pu_mammal_native = data.frame(RegionID = colnames(phy_turn_mammal_native),
                              PU = colMeans(phy_turn_mammal_native))
                                 
colnames(pu_mammal_extant)[which(colnames(pu_mammal_extant) == 'PU')] = 'all_PU'
pu_mammal_delta = pu_mammal_native %>% 
  left_join(pu_mammal_extant, by = 'RegionID')
pu_mammal_delta$RegionID = as.integer(pu_mammal_delta$RegionID)
pu_mammal_delta$delta_PU = log(pu_mammal_delta$all_PU/
                                 pu_mammal_delta$PU)

delta_pu_mammal_dat = delta_pu_mammal_dat %>% left_join(pu_mammal_delta[,c('RegionID',
                                                                           'delta_PU')],
                                                        by = join_by('region' == 'RegionID'))

delta_pu_mammal_dat2 = delta_pu_mammal_dat %>%
  dplyr::mutate(
    across(
      starts_with("delta_pu_path"),
      abs,
      .names = "{.col}_abs"
    )
  ) %>%
  mutate(
    abs_sum = delta_pu_path1_abs +
      delta_pu_path2_abs +
      delta_pu_path3_abs +
      delta_pu_path4_abs +
      delta_pu_path5_abs
  ) %>%
  dplyr::mutate(
    across(
      ends_with("_abs"),
      ~ .x / abs_sum,
      .names = "{.col}_prop")
  ) %>% 
  left_join(df_trans[,c('RegionID', 'Island')], by = join_by('region' == 'RegionID'))

colnames(delta_pu_mammal_dat2)[which(colnames(delta_pu_mammal_dat2) == 'delta_PU')] = 'real'
delta_pu_mammal_dat2 = arrange(delta_pu_mammal_dat2,
                               delta_pu_mammal_dat2$real)

delta_pu_mammal_dat3 = delta_pu_mammal_dat2 %>% 
  mutate(f_real = 1:nrow(delta_pu_mammal_dat2))

delta_pu_mammal_dat3[delta_pu_mammal_dat3$Island == 0,]$Island = 'Mainland'
delta_pu_mammal_dat3[delta_pu_mammal_dat3$Island == 1,]$Island = 'Island'

zero_mammal_mainland = delta_pu_mammal_dat3 %>%
  filter(Island == 'Mainland') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)          

zero_mammal_island = delta_pu_mammal_dat3 %>%
  filter(Island == 'Island') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)    

delta_pu_mammal_dat3_long = delta_pu_mammal_dat3 %>% 
  pivot_longer(cols = c(delta_pu_path1_abs_prop:delta_pu_path5_abs_prop),
               values_to = 'delta_PU',
               names_to = 'path')

delta_pu_mammal_dat3_long$path = factor(delta_pu_mammal_dat3_long$path,
                                          levels = c("delta_pu_path1_abs_prop", "delta_pu_path2_abs_prop",
                                                     "delta_pu_path3_abs_prop","delta_pu_path4_abs_prop",
                                                     "delta_pu_path5_abs_prop"))
delta_pu_mammal_dat3_long$Island = factor(delta_pu_mammal_dat3_long$Island,
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

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_delta_region_PU_mammal = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = f_real, y = delta_PU, fill = path),
            data = delta_pu_mammal_dat3_long) +
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

plot_delta_region_PU_mammal_2 = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = region, y = delta_PU, fill = path),
            data = delta_pu_mammal_dat3_long) +
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
delta_pu_bird_dat$region = as.integer(delta_pu_bird_dat$region)
pu_bird_extant = data.frame(RegionID = colnames(phy_turn_bird_extant),
                              PU = colMeans(phy_turn_bird_extant))

pu_bird_native = data.frame(RegionID = colnames(phy_turn_bird_native),
                              PU = colMeans(phy_turn_bird_native))

colnames(pu_bird_extant)[which(colnames(pu_bird_extant) == 'PU')] = 'all_PU'
pu_bird_delta = pu_bird_native %>% 
  left_join(pu_bird_extant, by = 'RegionID')
pu_bird_delta$RegionID =as.integer(pu_bird_delta$RegionID)

pu_bird_delta$delta_PU = log(pu_bird_delta$all_PU/
                                 pu_bird_delta$PU)

delta_pu_bird_dat = delta_pu_bird_dat %>% left_join(pu_bird_delta[,c('RegionID',
                                                                           'delta_PU')],
                                                        by = join_by('region' == 'RegionID'))

delta_pu_bird_dat2 = delta_pu_bird_dat %>%
  dplyr::mutate(
    across(
      starts_with("delta_pu_path"),
      abs,
      .names = "{.col}_abs"
    )
  ) %>%
  mutate(
    abs_sum = delta_pu_path1_abs +
      delta_pu_path2_abs +
      delta_pu_path3_abs +
      delta_pu_path4_abs +
      delta_pu_path5_abs
  ) %>%
  dplyr::mutate(
    across(
      ends_with("_abs"),
      ~ .x / abs_sum,
      .names = "{.col}_prop")
  ) %>% 
  left_join(df_trans[,c('RegionID', 'Island')], by = join_by('region' == 'RegionID'))

colnames(delta_pu_bird_dat2)[which(colnames(delta_pu_bird_dat2) == 'delta_PU')] = 'real'
delta_pu_bird_dat2 = arrange(delta_pu_bird_dat2,
                               delta_pu_bird_dat2$real)

delta_pu_bird_dat3 = delta_pu_bird_dat2 %>% 
  mutate(f_real = 1:nrow(delta_pu_bird_dat2))

delta_pu_bird_dat3[delta_pu_bird_dat3$Island == 0,]$Island = 'Mainland'
delta_pu_bird_dat3[delta_pu_bird_dat3$Island == 1,]$Island = 'Island'

zero_bird_mainland = delta_pu_bird_dat3 %>%
  filter(Island == 'Mainland') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)          

zero_bird_island = delta_pu_bird_dat3 %>%
  filter(Island == 'Island') %>%     
  arrange(real) %>%                    
  filter(real > 0) %>%                 
  slice(1) %>% pull(f_real)    

delta_pu_bird_dat3_long = delta_pu_bird_dat3 %>% 
  pivot_longer(cols = c(delta_pu_path1_abs_prop:delta_pu_path5_abs_prop),
               values_to = 'delta_PU',
               names_to = 'path')

delta_pu_bird_dat3_long$path = factor(delta_pu_bird_dat3_long$path,
                                        levels = c("delta_pu_path1_abs_prop", "delta_pu_path2_abs_prop",
                                                   "delta_pu_path3_abs_prop","delta_pu_path4_abs_prop",
                                                   "delta_pu_path5_abs_prop"))
delta_pu_bird_dat3_long$Island = factor(delta_pu_bird_dat3_long$Island,
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

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_delta_region_PU_bird = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = f_real, y = delta_PU, fill = path),
            data = delta_pu_bird_dat3_long) +
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

plot_delta_region_PU_bird_2 = ggplot() +
  facet_wrap(vars(Island), scales = "free", ncol = 2)+
  geom_area(aes(x = region, y = delta_PU, fill = path),
            data = delta_pu_bird_dat3_long) +
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
plot_all_regionPU_mammal_pathways = ggarrange(plot_delta_region_ED_plant, plot_delta_region_ED_bird,
                                              plot_delta_region_ED_mammal, plot_delta_region_ED_fish,
                                              ncol = 2, nrow = 2,
                                              labels = c('a', 'b', 'c', 'd'),
                                              common.legend = T)

plot_all_regionPU_pathways = ggarrange(plot_delta_region_PU_bird,
                                              plot_delta_region_PU_mammal,
                                              ncol = 2, nrow = 1,
                                              labels = c('a', 'b', 'c', 'd'),
                                              common.legend = T)


#4. Export figures for paper----
library(devEMF)
library(cowplot)


emf('figures/plot_all_regionPU_mammal_pathways.emf',
    width = 25, height = 20, coordDPI = 600*0.8, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_all_regionPU_mammal_pathways
dev.off() #turn off device and finalize file


emf('figures/plot_all_regionPU_pathways.emf',
    width = 25, height = 12, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_all_regionPU_pathways
dev.off() #turn off device and finalize file
