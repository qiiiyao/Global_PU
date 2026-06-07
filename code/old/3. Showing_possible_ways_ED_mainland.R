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
df_mainland = df_trans %>% filter(!(Island == 1 & 
                                      Area < 5e3))
tdwg_mainlands = as.character(sort(unique(df_mainland$RegionID)))

shp.glonaf.new = st_read("data/Plants/shp_glonaf_new_eck4.shp")
shp.glonaf.mainland = shp.glonaf.new %>% filter(!(Island == 1 & 
                                                    Area < 5e3)) 
glonaf_mainlands = as.character(sort(unique(shp.glonaf.mainland$Region_id)))


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

tdwg_mainlands_delta_mammal = intersect(colnames(turnover_mammal_delta_mat),
                                        tdwg_mainlands)

turnover_mammal_delta_mat = turnover_mammal_delta_mat[tdwg_mainlands_delta_mammal,
                                                      tdwg_mainlands_delta_mammal]
turnover_mammal_delta_mat_path3 = turnover_mammal_delta_mat_path3[tdwg_mainlands_delta_mammal,
                                                      tdwg_mainlands_delta_mammal]
turnover_mammal_delta_mat_path4 = turnover_mammal_delta_mat_path4[tdwg_mainlands_delta_mammal,
                                                      tdwg_mainlands_delta_mammal]
turnover_mammal_delta_mat_path5 = turnover_mammal_delta_mat_path5[tdwg_mainlands_delta_mammal,
                                                      tdwg_mainlands_delta_mammal]
turnover_mammal_delta_mat_path6 = turnover_mammal_delta_mat_path6[tdwg_mainlands_delta_mammal,
                                                      tdwg_mainlands_delta_mammal]
turnover_mammal_delta_mat_path7 = turnover_mammal_delta_mat_path7[tdwg_mainlands_delta_mammal,
                                                      tdwg_mainlands_delta_mammal]


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
  dplyr::select(c(RegionID, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-RegionID))) %>%
  dplyr::mutate(across(-c(RegionID, total), ~ .x / total)) %>%
  dplyr::select(-total)

region_ED_mammal_delta_long = region_ED_mammal_delta_2 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_mammal_delta_long$path = factor(region_ED_mammal_delta_long$path,
                                          levels = c("path3_diff", "path4_diff",
                                                     "path5_diff","path6_diff",
                                                     "path7_diff"))



plot_delta_region_ED_mammal = ggplot() +
  geom_area(aes(x = RegionID, y = delta_ED, fill = path),
            data = region_ED_mammal_delta_long) +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  labs(x = "Region ID",
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

tdwg_mainlands_delta_bird = intersect(colnames(turnover_bird_delta_mat),
                                        tdwg_mainlands)

turnover_bird_delta_mat = turnover_bird_delta_mat[tdwg_mainlands_delta_bird,
                                                      tdwg_mainlands_delta_bird]
turnover_bird_delta_mat_path3 = turnover_bird_delta_mat_path3[tdwg_mainlands_delta_bird,
                                                                  tdwg_mainlands_delta_bird]
turnover_bird_delta_mat_path4 = turnover_bird_delta_mat_path4[tdwg_mainlands_delta_bird,
                                                                  tdwg_mainlands_delta_bird]
turnover_bird_delta_mat_path5 = turnover_bird_delta_mat_path5[tdwg_mainlands_delta_bird,
                                                                  tdwg_mainlands_delta_bird]
turnover_bird_delta_mat_path6 = turnover_bird_delta_mat_path6[tdwg_mainlands_delta_bird,
                                                                  tdwg_mainlands_delta_bird]
turnover_bird_delta_mat_path7 = turnover_bird_delta_mat_path7[tdwg_mainlands_delta_bird,
                                                                  tdwg_mainlands_delta_bird]


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
  dplyr::select(c(RegionID, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-RegionID))) %>%
  dplyr::mutate(across(-c(RegionID, total), ~ .x / total)) %>%
  dplyr::select(-total)

region_ED_bird_delta_long = region_ED_bird_delta_2 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_bird_delta_long$path = factor(region_ED_bird_delta_long$path,
                                          levels = c("path3_diff", "path4_diff",
                                                     "path5_diff","path6_diff",
                                                     "path7_diff"))

plot_delta_region_ED_bird = ggplot() +
  geom_area(aes(x = RegionID, y = delta_ED, fill = path),
            data = region_ED_bird_delta_long) +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  labs(x = "Region ID",
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

glonaf_mainlands_delta_plant = intersect(colnames(turnover_plant_delta_mat),
                                         glonaf_mainlands)

turnover_plant_delta_mat = turnover_plant_delta_mat[glonaf_mainlands_delta_plant,
                                                  glonaf_mainlands_delta_plant]
turnover_plant_delta_mat_path3 = turnover_plant_delta_mat_path3[glonaf_mainlands_delta_plant,
                                                              glonaf_mainlands_delta_plant]
turnover_plant_delta_mat_path4 = turnover_plant_delta_mat_path4[glonaf_mainlands_delta_plant,
                                                              glonaf_mainlands_delta_plant]
turnover_plant_delta_mat_path5 = turnover_plant_delta_mat_path5[glonaf_mainlands_delta_plant,
                                                              glonaf_mainlands_delta_plant]
turnover_plant_delta_mat_path6 = turnover_plant_delta_mat_path6[glonaf_mainlands_delta_plant,
                                                              glonaf_mainlands_delta_plant]
turnover_plant_delta_mat_path7 = turnover_plant_delta_mat_path7[glonaf_mainlands_delta_plant,
                                                              glonaf_mainlands_delta_plant]


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
  dplyr::select(c(RegionID, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-RegionID))) %>%
  dplyr::mutate(across(-c(RegionID, total), ~ .x / total)) %>%
  dplyr::select(-total)

region_ED_plant_delta_long = region_ED_plant_delta_2 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_plant_delta_long$path = factor(region_ED_plant_delta_long$path,
                                        levels = c("path3_diff", "path4_diff",
                                                   "path5_diff","path6_diff",
                                                   "path7_diff"))

plot_delta_region_ED_plant = ggplot() +
  geom_area(aes(x = RegionID, y = delta_ED, fill = path),
            data = region_ED_plant_delta_long) +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  labs(x = "Region ID",
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

region_ED_fish_delta = data.frame(Basin = 1:length(colnames(turnover_fish_delta_mat)),
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
  dplyr::select(c(Basin, path3_diff:path7_diff)) %>% 
  dplyr::mutate(total = rowSums(across(-Basin))) %>%
  dplyr::mutate(across(-c(Basin, total), ~ .x / total)) #%>%
  #dplyr::select(-total)

region_ED_fish_delta_long = region_ED_fish_delta_2 %>% 
  pivot_longer(cols = c(path3_diff:path7_diff),
               values_to = 'delta_ED',
               names_to = 'path')

region_ED_fish_delta_long$path = factor(region_ED_fish_delta_long$path,
                                         levels = c("path3_diff", "path4_diff",
                                                    "path5_diff","path6_diff",
                                                    "path7_diff"))

plot_delta_region_ED_fish = ggplot() +
  geom_area(aes(x = Basin, y = delta_ED, fill = path),
            data = region_ED_fish_delta_long) +
  scale_fill_manual(values = colors_5d, 
                    labels = c('Path1', 
                               'Path2',
                               'Path3', 
                               'Path4',
                               'Path5')) +
  scale_y_continuous(limits = c(0,1), 
                     labels = c(0, 0.25, 0.5,
                                0.75, 1))+
  labs(x = "Basin ID",
       y = "Relative path importance") +
  theme(legend.position = "right",
        legend.title = element_blank())




#3. Exporting all taxon column figures----------------------------
library(ggpubr)
plot_regionED_mammal_pathways = ggarrange(plot_delta_region_ED_plant, plot_delta_region_ED_bird,
                               plot_delta_region_ED_mammal, plot_delta_region_ED_fish,
                               ncol = 2, nrow = 2,
                               labels = c('a', 'b', 'c', 'd'),
                               common.legend = T)

#4. Export figures for paper----
library(devEMF)
library(cowplot)


emf('figures/plot_regionED_mammal_pathways.emf',
    width = 20, height = 20, coordDPI = 600*0.8, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_regionED_mammal_pathways
dev.off() #turn off device and finalize file
