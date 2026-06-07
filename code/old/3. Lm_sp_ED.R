#Lirong Cai
#Email:lirong.cai18@gmail.com

#In this script, we modeled geographic patterns of phylogenetic and taxonomic turnover

#0 loading----------------------------
rm(list = ls())
library(reshape2)
library(plyr)
library(betareg)
library(MASS)
#library(gdm)
library(vegan)
library(rgeos)
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

#1. modelling for Mammals and Birds: TDWG4----------------------------
# load individual climate distances
load('results/primary_results/distances_TDWG/clim_dist.rdata')
load("results/primary_results/distances_TDWG/geodistances_grid_7.RDATA")
#load current cost distance
load("results/primary_results/distances_TDWG/graphdistances_barriers_grid_7.RDATA")
geo_dist = as.matrix(geodistances) 
geo_dist_barriers = as.matrix(graphdistances_barriers)
rm(geodistances)
rm(graphdistances_barriers)

geo_dist_df = data.frame(RegionID = as.integer(colnames(geo_dist)),
                         geo_dist = colMeans(geo_dist))
geo_dist_barriers_df = data.frame(RegionID = as.integer(colnames(geo_dist_barriers)),
                                  geo_dist_barriers = colMeans(geo_dist_barriers))

clim_dist_df = data.frame(RegionID = as.integer(colnames(clim_dist)),
                          clim_dist = colMeans(clim_dist))


# load trade flow value matrix
load('results/primary_results/distances_TDWG/tradeflow_mat_a.RDATA')
load("results/primary_results/predictors_TDWG/RegionID_direct_trade.rdata")
load("results/primary_results/predictors_TDWG/pd_bird_exotic.RDATA")
load("results/primary_results/predictors_TDWG/pd_bird_native.RDATA")
load("results/primary_results/predictors_TDWG/pd_mammal_exotic.RDATA")
load("results/primary_results/predictors_TDWG/pd_mammal_native.RDATA")

colnames(pd_bird_exotic)[which(colnames(pd_bird_exotic) == 'PD')] = 'pd_exotic'
colnames(pd_bird_native)[which(colnames(pd_bird_native) == 'PD')] = 'pd_native'
colnames(pd_mammal_exotic)[which(colnames(pd_mammal_exotic) == 'PD')] = 'pd_exotic'
colnames(pd_mammal_native)[which(colnames(pd_mammal_native) == 'PD')] = 'pd_native'

pd_mammal_native$RegionID = as.integer(pd_mammal_native$RegionID)
pd_mammal_exotic$RegionID = as.integer(pd_mammal_exotic$RegionID)
pd_bird_native$RegionID = as.integer(pd_bird_native$RegionID)
pd_bird_exotic$RegionID = as.integer(pd_bird_exotic$RegionID)

dim(tradeflow_mat_a)
## trade from i(col) to j(row)
i_mean_trade = data.frame(RegionID = as.integer(row.names(tradeflow_mat_a)),
                     i_trade_v = rowMeans(tradeflow_mat_a, na.rm = T))

#load species level ED
load("data/Mammals/Results_data/native/sp_nati.rdata")

load('results/primary_results/distances_beta/ED_mammal_all.rdata')
load('results/primary_results/distances_beta/ED_mammal_native.rdata')
load('results/primary_results/distances_beta/ED_bird_extant.rdata')
load('results/primary_results/distances_beta/ED_bird_native.rdata')

ED_mammal_extant_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_all$df)
colnames(ED_mammal_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_mammal_native_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_native$df)
colnames(ED_mammal_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_mammal = ED_mammal_native_1 %>% left_join(ED_mammal_extant_1,
                                                   by = 'RegionID')

Delta_ED_mammal$delta_mean_ED = log(Delta_ED_mammal$mean_all_ED / 
                                      Delta_ED_mammal$mean_native_ED)
Delta_ED_mammal$taxa = 'Mammal'
Delta_ED_mammal_df = Delta_ED_mammal %>%
  left_join(pd_mammal_native,  by = c('RegionID')) %>% 
  left_join(pd_mammal_exotic,  by = c('RegionID')) %>% 
  left_join(geo_dist_df,  by = c('RegionID')) %>% 
  left_join(geo_dist_barriers_df,  by = c('RegionID')) %>% 
  left_join(clim_dist_df,  by = c('RegionID')) %>%
  left_join(i_mean_trade,  by = c('RegionID'))

ED_bird_extant_1 = cbind(RegionID = df_clear$RegionID,
                           ED_bird_extant$df)
colnames(ED_bird_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_bird_native_1 = cbind(RegionID = df_clear$RegionID,
                           ED_bird_native$df)
colnames(ED_bird_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_bird = ED_bird_native_1 %>% left_join(ED_bird_extant_1,
                                                   by = 'RegionID')

Delta_ED_bird$delta_mean_ED = log(Delta_ED_bird$mean_all_ED / 
                                  Delta_ED_bird$mean_native_ED)
Delta_ED_bird$taxa = 'Bird'
Delta_ED_bird_df = Delta_ED_bird %>%
  left_join(pd_mammal_native,  by = c('RegionID')) %>% 
  left_join(pd_mammal_exotic,  by = c('RegionID')) %>% 
  left_join(geo_dist_df,  by = c('RegionID')) %>% 
  left_join(geo_dist_barriers_df,  by = c('RegionID')) %>% 
  left_join(clim_dist_df,  by = c('RegionID')) %>%
  left_join(i_mean_trade,  by = c('RegionID'))


sp_ED_tdwg_delta_df = rbind(Delta_ED_mammal_df, Delta_ED_bird_df) %>% 
  left_join(df_clear,  by = c('RegionID')) 

sp_ED_tdwg_delta_df[is.na(sp_ED_tdwg_delta_df$Island),]$Island = 0

sp_ED_tdwg_delta_df %>% filter(is.na(pd_native))
sp_ED_tdwg_delta_df %>% filter(is.na(pd_exotic))

str(sp_ED_tdwg_delta_df)


##1.1 All regions with direct and relative population * trade data----
sp_ED_tdwg_delta_df = sp_ED_tdwg_delta_df %>% filter(!is.na(pd_exotic))
sp_ED_tdwg_delta_df = sp_ED_tdwg_delta_df %>% filter(!is.na(i_trade_v))

sp_ED_tdwg_delta_df_mammal = sp_ED_tdwg_delta_df %>% filter(taxa == 'Mammal')
sp_ED_tdwg_delta_df_bird = sp_ED_tdwg_delta_df %>% filter(taxa == 'Bird')

# mammal
coords = cbind(
  sp_ED_tdwg_delta_df_mammal$Lon,
  sp_ED_tdwg_delta_df_mammal$Lat
)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(sp_ED_tdwg_delta_df_mammal$delta_mean_ED, weight_matrix)
print(moran_result)  # Result note: p-value = 0.006


sp_ED_tdwg_delta_df_mammal$modi_pd_native = scale(log(sp_ED_tdwg_delta_df_mammal$pd_native+0.001))
sp_ED_tdwg_delta_df_mammal$modi_pd_exotic = scale(log(sp_ED_tdwg_delta_df_mammal$pd_exotic+0.001))
sp_ED_tdwg_delta_df_mammal$modi_area = scale(log(sp_ED_tdwg_delta_df_mammal$Area))
sp_ED_tdwg_delta_df_mammal$modi_i_trade_v = scale(log(sp_ED_tdwg_delta_df_mammal$i_trade_v+0.001))
sp_ED_tdwg_delta_df_mammal$modi_geo_dist = scale(log(sp_ED_tdwg_delta_df_mammal$geo_dist))
sp_ED_tdwg_delta_df_mammal$modi_geo_dist_barriers = scale(log(sp_ED_tdwg_delta_df_mammal$geo_dist_barriers))
sp_ED_tdwg_delta_df_mammal$modi_clim_dist = scale(log(sp_ED_tdwg_delta_df_mammal$clim_dist))


gls_mammals = list()
gls_mammals[[1]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_mammal,
  correlation = corExp(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_mammals[[2]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_mammal,
  correlation = corGaus(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_mammals[[3]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_mammal,
  correlation = corLin(form=~Lon+Lat),
  method = "REML"
)

print(data.frame(SD=sapply(gls_mammals,
                           function(x) sd(resid(x))),AIC=sapply(gls_mammals,AIC)))

summary(gls_mammals[[1]])

lm_mammals = lm(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
                  modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
                data = sp_ED_tdwg_delta_df_mammal)
sd(resid(lm_mammals))
summary(lm_mammals)


#### bird
coords = cbind(
  sp_ED_tdwg_delta_df_bird$Lon,
  sp_ED_tdwg_delta_df_bird$Lat
)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(sp_ED_tdwg_delta_df_bird$delta_mean_ED, weight_matrix)
print(moran_result)  # Result note: p-value = 9.055e-16


sp_ED_tdwg_delta_df_bird$modi_pd_native = scale(log(sp_ED_tdwg_delta_df_bird$pd_native+0.001))
sp_ED_tdwg_delta_df_bird$modi_pd_exotic = scale(log(sp_ED_tdwg_delta_df_bird$pd_exotic+0.001))
sp_ED_tdwg_delta_df_bird$modi_area = scale(log(sp_ED_tdwg_delta_df_bird$Area))
sp_ED_tdwg_delta_df_bird$modi_i_trade_v = scale(log(sp_ED_tdwg_delta_df_bird$i_trade_v+0.001))
sp_ED_tdwg_delta_df_bird$modi_geo_dist = scale(log(sp_ED_tdwg_delta_df_bird$geo_dist))
sp_ED_tdwg_delta_df_bird$modi_geo_dist_barriers = scale(log(sp_ED_tdwg_delta_df_bird$geo_dist_barriers))
sp_ED_tdwg_delta_df_bird$modi_clim_dist = scale(log(sp_ED_tdwg_delta_df_bird$clim_dist))

plot(sp_ED_tdwg_delta_df_bird$Island, sp_ED_tdwg_delta_df_bird$delta_mean_ED)
lm_birds = lm(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
                modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
              data = sp_ED_tdwg_delta_df_bird)
summary(lm_birds)


gls_birds = list()
gls_birds[[1]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_bird,
  correlation = corExp(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_birds[[2]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_bird,
  correlation = corGaus(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_birds[[3]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_bird,
  correlation = corSpher(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_birds[[4]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_tdwg_delta_df_bird,
  correlation = corLin(form=~Lon+Lat, nugget=T),
  method = "REML"
)
sd(resid(lm_birds))
print(data.frame(SD=sapply(gls_birds,function(x) sd(resid(x))),AIC=sapply(gls_birds,AIC)))

summary(gls_birds[[1]]) # best model with corExp structure


pred_gls_birds = get_model_data(gls_birds[[1]], type = 'est') 
pred_gls_mammals = get_model_data(gls_mammals[[1]], type = 'est') 


#2. modelling for Plants: GIFT----------------------------
shp_glonaf_new = st_read("data/Plants/shp_glonaf_new_eck4.shp")

# load individual climate distances
load('results/primary_results/distances_GIFT/clim_dist.rdata')
load("results/primary_results/distances_GIFT/geodistances_grid_7.RDATA")
#load current cost distance
load("results/primary_results/distances_GIFT/graphdistances_barriers_grid_7.RDATA")
geo_dist = as.matrix(geodistances) 
geo_dist_barriers = as.matrix(graphdistances_barriers)
rm(geodistances)
rm(graphdistances_barriers)


geo_dist_df = data.frame(Region_id = as.integer(colnames(geo_dist)),
                         geo_dist = colMeans(geo_dist))
geo_dist_barriers_df = data.frame(Region_id = as.integer(colnames(geo_dist_barriers)),
                                  geo_dist_barriers = colMeans(geo_dist_barriers))

clim_dist_df = data.frame(Region_id = as.integer(colnames(clim_dist)),
                          clim_dist = colMeans(clim_dist))

# load trade flow value matrix
load('results/primary_results/distances_GIFT/tradeflow_mat_a.RDATA')
load("results/primary_results/predictors_GIFT/Region_id_direct_trade.rdata")
load("results/primary_results/predictors_GIFT/pd_plant_exotic.RDATA")
load("results/primary_results/predictors_GIFT/pd_plant_native.RDATA")

colnames(pd_plant_exotic)[which(colnames(pd_plant_exotic) == 'pd')] = 'pd_exotic'
colnames(pd_plant_native)[which(colnames(pd_plant_native) == 'pd')] = 'pd_native'
pd_plant_native$Region_id = as.integer(pd_plant_native$Region_id)
pd_plant_exotic$Region_id = as.integer(pd_plant_exotic$Region_id)

## trade from i(col) to j(row)
i_mean_trade = data.frame(Region_id = as.integer(row.names(tradeflow_mat_a)),
                          i_trade_v = rowMeans(tradeflow_mat_a, na.rm = T))

#load sp_level ED
load('results/primary_results/distances_beta/ED_plant_extant.rdata')
load('results/primary_results/distances_beta/ED_plant_native.rdata')

ED_plant_extant_1 = cbind(Region_id = sort(unique(pd_plant_native$Region_id)),
                           ED_plant_extant$df)
colnames(ED_plant_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_plant_native_1 = cbind(Region_id = sort(unique(pd_plant_native$Region_id)),
                          ED_plant_native$df)
colnames(ED_plant_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_plant = ED_plant_native_1 %>% left_join(ED_plant_extant_1,
                                                 by = 'Region_id')

Delta_ED_plant$delta_mean_ED = log(Delta_ED_plant$mean_all_ED / 
                                   Delta_ED_plant$mean_native_ED)
Delta_ED_plant$taxa = 'plant'
Delta_ED_plant_df = Delta_ED_plant %>%
  left_join(pd_plant_native,  by = c('Region_id')) %>% 
  left_join(pd_plant_exotic,  by = c('Region_id'))%>% 
  left_join(geo_dist_df,  by = c('Region_id')) %>% 
  left_join(geo_dist_barriers_df,  by = c('Region_id')) %>% 
  left_join(clim_dist_df,  by = c('Region_id')) %>%
  left_join(i_mean_trade,  by = c('Region_id'))

sp_ED_gift_delta_df = Delta_ED_plant_df %>% 
  left_join(shp_glonaf_new,  by = c('Region_id')) 


str(sp_ED_gift_delta_df)

##2.1 All regions with direct and relative population * trade data----
sp_ED_gift_delta_df = sp_ED_gift_delta_df %>% filter(!is.na(pd_exotic))
sp_ED_gift_delta_df = sp_ED_gift_delta_df %>% filter(!is.na(i_trade_v))

sp_ED_gift_delta_df$modi_pd_native = scale(log(sp_ED_gift_delta_df$pd_native+0.001))
sp_ED_gift_delta_df$modi_pd_exotic = scale(log(sp_ED_gift_delta_df$pd_exotic+0.001))
sp_ED_gift_delta_df$modi_area = scale(log(sp_ED_gift_delta_df$Area))
sp_ED_gift_delta_df$modi_i_trade_v = scale(log(sp_ED_gift_delta_df$i_trade_v+0.001))
sp_ED_gift_delta_df$modi_geo_dist = scale(log(sp_ED_gift_delta_df$geo_dist))
sp_ED_gift_delta_df$modi_geo_dist_barriers = scale(log(sp_ED_gift_delta_df$geo_dist_barriers))
sp_ED_gift_delta_df$modi_clim_dist = scale(log(sp_ED_gift_delta_df$clim_dist))

hist(sp_ED_gift_delta_df$delta_mean_ED)
#lm_plants = lm(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area 
 #              + modi_i_trade_v + Island,
  #             data = sp_ED_gift_delta_df)
#summary(lm_plants)

coords = cbind(sp_ED_gift_delta_df$Lon, sp_ED_gift_delta_df$Lat)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(sp_ED_gift_delta_df$delta_mean_ED, weight_matrix)
print(moran_result)  # Result note: p-value = < 2.2e-16

gls_plants = list()
gls_plants[[1]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_gift_delta_df,
  correlation = corExp(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_plants[[2]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_gift_delta_df,
  correlation = corGaus(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_plants[[3]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_gift_delta_df,
  correlation = corSpher(form=~Lon+Lat, nugget=T),
  method = "REML"
)

gls_plants[[4]] = gls(
  delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island + modi_i_trade_v + 
    modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
  data = sp_ED_gift_delta_df,
  correlation = corLin(form=~Lon+Lat),
  method = "REML")

print(data.frame(SD=sapply(gls_plants,function(x) sd(resid(x))),AIC=sapply(gls_plants,AIC)))

summary(gls_plants[[2]]) # best model with corGaus structure


pred_gls_plants = get_model_data(gls_plants[[2]], type = 'est') 




#3. modelling for Fishes: Basins----------------------------
df_fish = st_read("data/Fishes/data/Basin042017_3119.shp")

# load individual climate distances
load('results/primary_results/distances_basins/clim_dist.rdata')
load("results/primary_results/distances_basins/geodistances_grid_7.RDATA")
#load current cost distance
load("results/primary_results/distances_basins/graphdistances_barriers_grid_7.RDATA")
geo_dist = as.matrix(geodistances) 
geo_dist_barriers = as.matrix(graphdistances_barriers)
rm(geodistances)
rm(graphdistances_barriers)


geo_dist_df = data.frame(Basin.name = as.character(colnames(geo_dist)),
                         geo_dist = colMeans(geo_dist))
geo_dist_barriers_df = data.frame(Basin.name = as.character(colnames(geo_dist_barriers)),
                                  geo_dist_barriers = colMeans(geo_dist_barriers))

clim_dist_df = data.frame(Basin.name = as.character(colnames(clim_dist)),
                          clim_dist = colMeans(clim_dist))


# load trade flow value matrix
load('results/primary_results/distances_basins/tradeflow_mat_a.RDATA')
load("results/primary_results/predictors_basins/basin_direct_trade.rdata")
load("results/primary_results/predictors_basins/pd_fish_exotic.RDATA")
load("results/primary_results/predictors_basins/pd_fish_native.RDATA")

colnames(pd_fish_exotic)[which(colnames(pd_fish_exotic) == 'PD')] = 'pd_exotic'
colnames(pd_fish_native)[which(colnames(pd_fish_native) == 'PD')] = 'pd_native'

pd_fish_native$Basin.name = as.character(pd_fish_native$X1.Basin.Name)
pd_fish_exotic$Basin.name = as.character(pd_fish_exotic$X1.Basin.Name)

## trade from i(col) to j(row)
i_mean_trade = data.frame(Basin.name = as.character(row.names(tradeflow_mat_a)),
                          i_trade_v = rowMeans(tradeflow_mat_a, na.rm = T))

#load sp_level ED
load('results/primary_results/distances_beta/ED_fish_extant.rdata')
load('results/primary_results/distances_beta/ED_fish_native.rdata')
load('results/primary_results/distances_beta/LCBD_fish_extant.rdata')
turnover_fish_extant_mat = LCBD_fish_extant$beta_mat

ED_fish_extant_1 = cbind(Basin.name = sort(colnames(turnover_fish_extant_mat)),
                          ED_fish_extant$df)
colnames(ED_fish_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_fish_native_1 = cbind(Basin.name = sort(unique(pd_fish_native$X1.Basin.Name)),
                          ED_fish_native$df)
colnames(ED_fish_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_fish = ED_fish_native_1 %>% left_join(ED_fish_extant_1,
                                                 by = 'Basin.name')

Delta_ED_fish$delta_mean_ED = log(Delta_ED_fish$mean_all_ED / 
                                     Delta_ED_fish$mean_native_ED)
Delta_ED_fish$taxa = 'fish'
Delta_ED_fish_df = Delta_ED_fish %>%
  left_join(pd_fish_native,  by = c('Basin.name')) %>% 
  left_join(pd_fish_exotic,  by = c('Basin.name'))%>% 
  left_join(geo_dist_df,  by = c('Basin.name')) %>% 
  left_join(geo_dist_barriers_df,  by = c('Basin.name')) %>% 
  left_join(clim_dist_df,  by = c('Basin.name')) %>%
  left_join(i_mean_trade,  by = c('Basin.name'))

sp_ED_basin_delta_df = Delta_ED_fish_df %>% 
  left_join(df_fish,  by = join_by('Basin.name' == 'BasinName')) 


str(sp_ED_basin_delta_df)

##2.1 All regions with direct and relative population * trade data----
sp_ED_basin_delta_df[is.na(sp_ED_basin_delta_df$pd_exotic),
]$pd_exotic = 0
sp_ED_basin_delta_df = sp_ED_basin_delta_df %>% filter(!is.na(i_trade_v))

sp_ED_basin_delta_df$modi_pd_native = scale(log(sp_ED_basin_delta_df$pd_native+0.001))
sp_ED_basin_delta_df$modi_pd_exotic = scale(log(sp_ED_basin_delta_df$pd_exotic+0.001))
sp_ED_basin_delta_df$modi_area = scale(log(sp_ED_basin_delta_df$Surf_area))
sp_ED_basin_delta_df$modi_i_trade_v = scale(log(sp_ED_basin_delta_df$i_trade_v+0.001))
sp_ED_basin_delta_df$modi_geo_dist = scale(log(sp_ED_basin_delta_df$geo_dist))
sp_ED_basin_delta_df$modi_geo_dist_barriers = scale(log(sp_ED_basin_delta_df$geo_dist_barriers))
sp_ED_basin_delta_df$modi_clim_dist = scale(log(sp_ED_basin_delta_df$clim_dist))
sp_ED_basin_delta_df[is.na(sp_ED_basin_delta_df$Island),]$Island = 0

hist(sp_ED_basin_delta_df$delta_mean_ED)
#lm_plants = lm(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area 
#              + modi_i_trade_v + Island,
#             data = sp_ED_basin_delta_df)
#summary(lm_plants)

coords = cbind(sp_ED_basin_delta_df$Lon, sp_ED_basin_delta_df$Lat)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(sp_ED_basin_delta_df$delta_mean_ED, weight_matrix)
print(moran_result)  # Result note: p-value = < 2.2e-16

gls_fishs = list()
gls_fishs[[1]] = gls(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island +
                       modi_i_trade_v + modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
                     data = sp_ED_basin_delta_df,
                     correlation = corExp(form=~Lon+Lat, nugget=T),
                     method = "REML",
                     control = glsControl(apVar = FALSE,maxIter = 30,msMaxIter = 30,opt = "nlminb"))


gls_fishs[[2]] = gls(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island +
                       modi_i_trade_v + modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
                     data = sp_ED_basin_delta_df,
                     correlation = corGaus(form=~Lon+Lat, nugget=T),
                     method = "REML",
                     control = glsControl(apVar = FALSE,maxIter = 30,msMaxIter = 30,opt = "nlminb"))

gls_fishs[[3]] = gls(delta_mean_ED ~ modi_pd_native + modi_pd_exotic + modi_area + Island +
                       modi_i_trade_v + modi_geo_dist + modi_geo_dist_barriers + modi_clim_dist,
                     data = sp_ED_basin_delta_df,
                     correlation = corLin(form=~Lon+Lat, nugget=T),
                     method = "REML",
                     control = glsControl(apVar = FALSE,maxIter = 30,msMaxIter = 30,opt = "nlminb"))

save(gls_fishs, file = 'results/primary_results/fitted_mods/gls_sp_ED_fishes.rdata')


load('results/primary_results/fitted_mods/gls_sp_ED_fishes.rdata')
print(data.frame(SD=sapply(gls_fishs,function(x) sd(resid(x))),AIC=sapply(gls_fishs,AIC)))

summary(gls_fishs[[1]]) # best model with corGaus structure


pred_gls_fishs = get_model_data(gls_fishs[[1]], type = 'est') 





#4. Exporting all taxon column figures----------------------------
pred_gls_plants$taxa = 'Plants'
pred_gls_mammals$taxa = 'Mammals'
pred_gls_birds$taxa = 'Birds'
pred_gls_fishs$taxa = 'Fishes'

pred_gls_all = rbind(pred_gls_plants, 
                     pred_gls_mammals,
                     pred_gls_birds,
                     pred_gls_fishs) %>%
  as_tibble()

pred_gls_all$taxa = factor(pred_gls_all$taxa,
                           levels = c("Plants",
                                      "Mammals",
                                      "Birds",
                                      "Fishes"))

pred_gls_all$term = factor(pred_gls_all$term,
                           levels = c("modi_pd_exotic",
                                      "modi_pd_native",
                                      "modi_clim_dist",
                                      "modi_geo_dist_barriers",
                                      "modi_geo_dist",
                                      'modi_area',
                                      'Island',
                                      "modi_i_trade_v"))
pred_gls_all$sig = ifelse(pred_gls_all$p.value < 0.05, 1, 0)
pred_gls_all$sig = factor(pred_gls_all$sig,
                          levels = c(1,
                                     0))

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_char_of_spED = ggplot()+
  facet_wrap(vars(taxa), scales = "free", ncol = 4)+
  geom_point(mapping = aes(y = term, x = estimate, shape = sig),
             data = pred_gls_all,
             size = 2.5,
             color = 'black')+
  geom_linerange(mapping = aes(y = term, xmin = conf.low, xmax = conf.high),
                 data = pred_gls_all) + 
  geom_vline(xintercept = 0, linetype = 2, color = colors_4d[1])+
  scale_shape_manual(values = c(16, 1))+
  scale_x_continuous(limits = c(-0.052, 0.052), 
                     breaks = c(-0.05, 0, 0.05)) + 
  scale_y_discrete(labels = c('PD naturalized',
                              'PD native',
                              'Climatic distance',
                              'Geo barrier distance',
                              'Geo distance',
                              'Region size',
                              'Island',
                              'Trade import'))+
  labs(x = 'Coefficient estimate', y = '')+
  theme(#plot.title = element_blank(),
    legend.position = 'None',
    plot.title = element_text(hjust = 0.5))
plot_char_of_spED


#5. Export figures for paper----
library(devEMF)
library(cowplot)

emf('figures/plot_char_of_spED.emf',
    width = 25, height = 10*0.84, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_char_of_spED
dev.off() #turn off device and finalize file

