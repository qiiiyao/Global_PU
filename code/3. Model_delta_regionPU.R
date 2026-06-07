#Lirong Cai
#Email:lirong.cai18@gmail.com

#In this script, we modeled geographic patterns of phylogenetic and taxonomic turnover

#0 loading----------------------------
rm(list = ls())
gc()
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
library(gstat)  
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
load('results/primary_results/distances_TDWG/tradeflow_mat_s.RDATA')
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

## trade value between i(col) and j(row)
i_mean_trade = data.frame(RegionID = as.integer(row.names(tradeflow_mat)),
                          i_trade_v = rowMeans(tradeflow_mat, na.rm = T))

#load beta diversity
load('results/primary_results/distances_beta/phy_turn_mammal_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_native.rdata')

turnover_mammal_delta_mat = log((phy_turn_mammal_extant+0.001) /
                                  (phy_turn_mammal_native+0.001))
colMeans(turnover_mammal_delta_mat)

turnover_bird_delta_mat = log((phy_turn_bird_extant+0.001) /
                                (phy_turn_bird_native+0.001))

dim(turnover_mammal_delta_mat)
dim(turnover_bird_delta_mat)


turnover_mammal_delta_df = data.frame(RegionID = as.integer(colnames(turnover_mammal_delta_mat)),
                                      mean_delta_tur = colMeans(turnover_mammal_delta_mat),
                                      taxa = 'Mammal')
turnover_mammal_delta_df = turnover_mammal_delta_df %>%
  left_join(pd_mammal_native,  by = c('RegionID')) %>% 
  left_join(pd_mammal_exotic,  by = c('RegionID')) %>% 
  left_join(geo_dist_df,  by = c('RegionID')) %>% 
  left_join(geo_dist_barriers_df,  by = c('RegionID')) %>% 
  left_join(clim_dist_df,  by = c('RegionID')) %>%
  left_join(i_mean_trade,  by = c('RegionID'))
  

turnover_bird_delta_df = data.frame(RegionID = as.integer(colnames(turnover_bird_delta_mat)),
                                    mean_delta_tur = colMeans(turnover_bird_delta_mat),
                                    taxa = 'Bird')
turnover_bird_delta_df = turnover_bird_delta_df %>%
  left_join(pd_bird_native,  by = c('RegionID')) %>% 
  left_join(pd_bird_exotic,  by = c('RegionID')) %>% 
  left_join(geo_dist_df,  by = c('RegionID')) %>% 
  left_join(geo_dist_barriers_df,  by = c('RegionID')) %>% 
  left_join(clim_dist_df,  by = c('RegionID')) %>%
  left_join(i_mean_trade,  by = c('RegionID'))

turnover_tdwg_delta_df = rbind(turnover_mammal_delta_df, turnover_bird_delta_df) %>% 
  left_join(df_clear,  by = c('RegionID')) 

turnover_tdwg_delta_df[is.na(turnover_tdwg_delta_df$Island),]$Island = 0

turnover_tdwg_delta_df %>% filter(is.na(pd_native))
turnover_tdwg_delta_df %>% filter(is.na(pd_exotic))

str(turnover_tdwg_delta_df)


##1.1 All regions with direct and relative population * trade data----
turnover_tdwg_delta_df = turnover_tdwg_delta_df %>% filter(!is.na(pd_exotic))
turnover_tdwg_delta_df = turnover_tdwg_delta_df %>% filter(!is.na(i_trade_v))

turnover_tdwg_delta_df_mammal = turnover_tdwg_delta_df %>% filter(taxa == 'Mammal')
turnover_tdwg_delta_df_bird = turnover_tdwg_delta_df %>% filter(taxa == 'Bird')

# mammal
turnover_tdwg_delta_df_mammal$modi_pd_native = scale(log(turnover_tdwg_delta_df_mammal$pd_native+0.001))
turnover_tdwg_delta_df_mammal$modi_pd_exotic = scale(log(turnover_tdwg_delta_df_mammal$pd_exotic+0.001))
turnover_tdwg_delta_df_mammal$modi_area = scale(log(turnover_tdwg_delta_df_mammal$Area))
turnover_tdwg_delta_df_mammal$modi_i_trade_v = scale(log(turnover_tdwg_delta_df_mammal$i_trade_v+0.001))
turnover_tdwg_delta_df_mammal$modi_geo_dist = scale(log(turnover_tdwg_delta_df_mammal$geo_dist))
turnover_tdwg_delta_df_mammal$modi_geo_dist_barriers = scale(log(turnover_tdwg_delta_df_mammal$geo_dist_barriers))
turnover_tdwg_delta_df_mammal$modi_clim_dist = scale(log(turnover_tdwg_delta_df_mammal$clim_dist))


f1 = formula(mean_delta_tur ~ modi_pd_native + modi_pd_exotic + modi_area +
               Island+ modi_i_trade_v + 
               modi_clim_dist)

gls_mammals_1 = gls(f1, data = turnover_tdwg_delta_df_mammal)

E = residuals(gls_mammals_1)

mydata = data.frame(E, lon = turnover_tdwg_delta_df_mammal$Lon,
                    lat = turnover_tdwg_delta_df_mammal$Lat)
coordinates(mydata) = c("lon","lat")  
bubble(mydata, "E", col = c("black","grey"),
       main = "Residuals",
       xlab = "X-coordinates",  ylab = "Y-coordinates")

Vario.gls = Variogram(gls_mammals_1,
                      form =~Lon+Lat,  robust = TRUE,  resType = "pearson") 
plot(Vario.gls, smooth = TRUE)

coords = cbind(turnover_tdwg_delta_df_mammal$Lon, 
               turnover_tdwg_delta_df_mammal$Lat)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(turnover_tdwg_delta_df_mammal$mean_delta_tur, weight_matrix)
print(moran_result)  # Result note: p-value = 0.006



lme_mammals = lme(f1, random = ~1|Level2_cod,
                  correlation = corRatio(form=~Lon+Lat, nugget=T),
                  control = lmeControl(opt = 'optim'),
                  data = turnover_tdwg_delta_df_mammal)
summary(lme_mammals)

gls_mammals = list()
gls_mammals[[1]] = gls(f1, data = turnover_tdwg_delta_df_mammal)

gls_mammals[[2]] = gls(f1, data = turnover_tdwg_delta_df_mammal,
  correlation = corSpher(form=~Lon+Lat, nugget=T), 
  control = glsControl(opt = "optim"))

gls_mammals[[3]] = gls(f1, data = turnover_tdwg_delta_df_mammal,
  correlation = corLin(form=~Lon+Lat, nugget=T), 
  control = glsControl(opt = "optim"))

gls_mammals[[4]] = gls(f1, data = turnover_tdwg_delta_df_mammal,
  correlation = corRatio(form=~Lon+Lat, nugget=T))

gls_mammals[[5]] = gls(f1, data = turnover_tdwg_delta_df_mammal,
  correlation = corGaus(form=~Lon+Lat, nugget=T))

gls_mammals[[6]] = gls(f1, data = turnover_tdwg_delta_df_mammal,
  correlation = corExp(form=~Lon+Lat, nugget=T))

print(data.frame(SD=sapply(gls_mammals,function(x) sd(resid(x))),AIC=sapply(gls_mammals,AIC)))

summary(gls_mammals[[4]])


# bird
turnover_tdwg_delta_df_bird$modi_pd_native = scale(log(turnover_tdwg_delta_df_bird$pd_native+0.001))
turnover_tdwg_delta_df_bird$modi_pd_exotic = scale(log(turnover_tdwg_delta_df_bird$pd_exotic+0.001))
turnover_tdwg_delta_df_bird$modi_area = scale(log(turnover_tdwg_delta_df_bird$Area))
turnover_tdwg_delta_df_bird$modi_i_trade_v = scale(log(turnover_tdwg_delta_df_bird$i_trade_v+0.001))
turnover_tdwg_delta_df_bird$modi_geo_dist = scale(log(turnover_tdwg_delta_df_bird$geo_dist))
turnover_tdwg_delta_df_bird$modi_geo_dist_barriers = scale(log(turnover_tdwg_delta_df_bird$geo_dist_barriers))
turnover_tdwg_delta_df_bird$modi_clim_dist = scale(log(turnover_tdwg_delta_df_bird$clim_dist))


gls_birds_1 = gls(f1, data = turnover_tdwg_delta_df_bird)

E = residuals(gls_birds_1)

mydata = data.frame(E, lon = turnover_tdwg_delta_df_bird$Lon,
                    lat = turnover_tdwg_delta_df_bird$Lat)
coordinates(mydata) = c("lon","lat")  
bubble(mydata, "E", col = c("black","grey"),
       main = "Residuals",
       xlab = "X-coordinates",  ylab = "Y-coordinates")

Vario.gls = Variogram(gls_birds_1,
                      form =~Lon+Lat,  robust = TRUE,  resType = "pearson") 
plot(Vario.gls, smooth = TRUE)

coords = cbind(turnover_tdwg_delta_df_bird$Lon, 
               turnover_tdwg_delta_df_bird$Lat)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(turnover_tdwg_delta_df_bird$mean_delta_tur, weight_matrix)
print(moran_result)  # Result note: p-value = 0.006



lme_birds = lme(f1, random = ~1|Level2_cod,
                  correlation = corRatio(form=~Lon+Lat, nugget=T),
                  control = lmeControl(opt = 'optim'),
                  data = turnover_tdwg_delta_df_bird)
summary(lme_birds)

gls_birds = list()
gls_birds[[1]] = gls(f1, data = turnover_tdwg_delta_df_bird)

gls_birds[[2]] = gls(f1, data = turnover_tdwg_delta_df_bird,
                       correlation = corSpher(form=~Lon+Lat, nugget=T), 
                       control = glsControl(opt = "optim"))

gls_birds[[3]] = gls(f1, data = turnover_tdwg_delta_df_bird,
                       correlation = corLin(form=~Lon+Lat, nugget=T), 
                       control = glsControl(opt = "optim"))

gls_birds[[4]] = gls(f1, data = turnover_tdwg_delta_df_bird,
                       correlation = corRatio(form=~Lon+Lat, nugget=T))

gls_birds[[5]] = gls(f1, data = turnover_tdwg_delta_df_bird,
                       correlation = corGaus(form=~Lon+Lat, nugget=T))

gls_birds[[6]] = gls(f1, data = turnover_tdwg_delta_df_bird,
                       correlation = corExp(form=~Lon+Lat, nugget=T))

print(data.frame(SD=sapply(gls_birds,function(x) sd(resid(x))),AIC=sapply(gls_birds,AIC)))

summary(gls_birds[[4]])

pred_gls_birds = get_model_data(gls_birds[[4]], type = 'est') 
pred_gls_mammals = get_model_data(gls_mammals[[4]], type = 'est') 


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
load('results/primary_results/distances_GIFT/tradeflow_mat_s.RDATA')
load("results/primary_results/predictors_GIFT/Region_id_direct_trade.rdata")
load("results/primary_results/predictors_GIFT/pd_plant_exotic.RDATA")
load("results/primary_results/predictors_GIFT/pd_plant_native.RDATA")

colnames(pd_plant_exotic)[which(colnames(pd_plant_exotic) == 'pd')] = 'pd_exotic'
colnames(pd_plant_native)[which(colnames(pd_plant_native) == 'pd')] = 'pd_native'


pd_plant_native$Region_id = as.integer(pd_plant_native$Region_id)
pd_plant_exotic$Region_id = as.integer(pd_plant_exotic$Region_id)

## trade value between i(col) and j(row)
i_mean_trade = data.frame(Region_id = as.integer(row.names(tradeflow_mat)),
                          i_trade_v = rowMeans(tradeflow_mat, na.rm = T))

#load beta diversity
load('results/primary_results/distances_beta/phy_turn_plant_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_plant_native.rdata')

turnover_plant_delta_mat = log((phy_turn_plant_extant+0.001) /
                                 (phy_turn_plant_native+0.001))

dim(turnover_plant_delta_mat)

turnover_plant_delta_df = data.frame(Region_id = as.integer(colnames(turnover_plant_delta_mat)),
                                      mean_delta_tur = colMeans(turnover_plant_delta_mat),
                                      taxa = 'plant')
turnover_plant_delta_df = turnover_plant_delta_df %>%
  left_join(pd_plant_native,  by = c('Region_id')) %>% 
  left_join(pd_plant_exotic,  by = c('Region_id'))%>% 
  left_join(geo_dist_df,  by = c('Region_id')) %>% 
  left_join(geo_dist_barriers_df,  by = c('Region_id')) %>% 
  left_join(clim_dist_df,  by = c('Region_id')) %>%
  left_join(i_mean_trade,  by = c('Region_id'))

turnover_gift_delta_df = turnover_plant_delta_df %>% 
  left_join(shp_glonaf_new,  by = c('Region_id')) 


str(turnover_gift_delta_df)

##2.1 All regions with direct and relative population * trade data----
turnover_gift_delta_df = turnover_gift_delta_df %>% filter(!is.na(pd_exotic))
turnover_gift_delta_df = turnover_gift_delta_df %>% filter(!is.na(i_trade_v))

turnover_gift_delta_df$modi_pd_native = scale(log(turnover_gift_delta_df$pd_native+0.001))
turnover_gift_delta_df$modi_pd_exotic = scale(log(turnover_gift_delta_df$pd_exotic+0.001))
turnover_gift_delta_df$modi_area = scale(log(turnover_gift_delta_df$Area))
turnover_gift_delta_df$modi_i_trade_v = scale(log(turnover_gift_delta_df$i_trade_v+0.001))
turnover_gift_delta_df$modi_geo_dist = scale(log(turnover_gift_delta_df$geo_dist))
turnover_gift_delta_df$modi_geo_dist_barriers = scale(log(turnover_gift_delta_df$geo_dist_barriers))
turnover_gift_delta_df$modi_clim_dist = scale(log(turnover_gift_delta_df$clim_dist))

gls_plants_1 = gls(f1, data = turnover_gift_delta_df)

E = residuals(gls_plants_1)

mydata = data.frame(E, lon = turnover_gift_delta_df$Lon,
                    lat = turnover_gift_delta_df$Lat)
coordinates(mydata) = c("lon","lat")  
bubble(mydata, "E", col = c("black","grey"),
       main = "Residuals",
       xlab = "X-coordinates",  ylab = "Y-coordinates")

Vario.gls = Variogram(gls_plants_1,
                      form =~Lon+Lat,  robust = TRUE,  resType = "pearson") 
plot(Vario.gls, smooth = TRUE)

coords = cbind(turnover_gift_delta_df$Lon, 
               turnover_gift_delta_df$Lat)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(turnover_gift_delta_df$mean_delta_tur, weight_matrix)
print(moran_result)  # Result note: p-value = 0.006



lme_plants = lme(f1, random = ~1|Level2_cod,
                correlation = corRatio(form=~Lon+Lat, nugget=T),
                control = lmeControl(opt = 'optim'),
                data = turnover_gift_delta_df)
summary(lme_plants)

gls_plants = list()
gls_plants[[1]] = gls(f1, data = turnover_gift_delta_df)

gls_plants[[2]] = gls(f1, data = turnover_gift_delta_df,
                     correlation = corSpher(form=~Lon+Lat, nugget=T), 
                     control = glsControl(opt = "optim"))

gls_plants[[3]] = gls(f1, data = turnover_gift_delta_df,
                     correlation = corLin(form=~Lon+Lat, nugget=T), 
                     control = glsControl(opt = "optim"))

gls_plants[[4]] = gls(f1, data = turnover_gift_delta_df,
                     correlation = corRatio(form=~Lon+Lat, nugget=T))

gls_plants[[5]] = gls(f1, data = turnover_gift_delta_df,
                     correlation = corGaus(form=~Lon+Lat, nugget=T))

gls_plants[[6]] = gls(f1, data = turnover_gift_delta_df,
                     correlation = corExp(form=~Lon+Lat, nugget=T))

save(gls_plants, file = 'results/primary_results/gls_plants_region_ED_all.rdata')

load('results/primary_results/gls_plants_region_ED_all.rdata')

print(data.frame(SD=sapply(gls_plants,function(x) sd(resid(x))),AIC=sapply(gls_plants,AIC)))

summary(gls_plants[[3]])

pred_gls_plants = get_model_data(gls_plants[[3]], type = 'est') 



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
load('results/primary_results/distances_basins/tradeflow_mat_s.RDATA')
load("results/primary_results/predictors_basins/basin_direct_trade.rdata")
load("results/primary_results/predictors_basins/pd_fish_exotic.RDATA")
load("results/primary_results/predictors_basins/pd_fish_native.RDATA")

colnames(pd_fish_exotic)[which(colnames(pd_fish_exotic) == 'PD')] = 'pd_exotic'
colnames(pd_fish_native)[which(colnames(pd_fish_native) == 'PD')] = 'pd_native'


pd_fish_native$Basin.name = as.character(pd_fish_native$X1.Basin.Name)
pd_fish_exotic$Basin.name = as.character(pd_fish_exotic$X1.Basin.Name)

## trade value between i(col) and j(row)
i_mean_trade = data.frame(Basin.name = as.character(row.names(tradeflow_mat_s)),
                          i_trade_v = rowMeans(tradeflow_mat_s, na.rm = T))

#load beta diversity
load('results/primary_results/distances_beta/phy_turn_fish_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_fish_native.rdata')

phy_turn_fish_extant = phy_turn_fish_extant[
  colnames(phy_turn_fish_native),
  rownames(phy_turn_fish_native)]

turnover_fish_delta_mat = log((phy_turn_fish_extant+0.001) /
                                 (phy_turn_fish_native+0.001))

dim(turnover_fish_delta_mat)

turnover_fish_delta_df = data.frame(Basin.name = as.character(colnames(turnover_fish_delta_mat)),
                                     mean_delta_tur = colMeans(turnover_fish_delta_mat, na.rm = T),
                                     taxa = 'fish')

turnover_fish_delta_df = turnover_fish_delta_df %>%
  left_join(pd_fish_native,  by = c('Basin.name')) %>% 
  left_join(pd_fish_exotic,  by = c('Basin.name'))%>% 
  left_join(geo_dist_df,  by = c('Basin.name')) %>% 
  left_join(geo_dist_barriers_df,  by = c('Basin.name')) %>% 
  left_join(clim_dist_df,  by = c('Basin.name')) %>%
  left_join(i_mean_trade,  by = c('Basin.name'))

turnover_basins_delta_df = turnover_fish_delta_df %>% 
  left_join(df_fish,  by = join_by('Basin.name' == 'BasinName')) 


str(turnover_basins_delta_df)

##2.1 All regions with direct and relative population * trade data----
turnover_basins_delta_df[is.na(turnover_basins_delta_df$pd_exotic),
                         ]$pd_exotic = 0

turnover_basins_delta_df = turnover_basins_delta_df %>% 
  filter(!is.na(i_trade_v))

turnover_basins_delta_df = turnover_basins_delta_df %>% 
  filter(!is.na(mean_delta_tur))

turnover_basins_delta_df$modi_pd_native = scale(log(turnover_basins_delta_df$pd_native+0.001))
turnover_basins_delta_df$modi_pd_exotic = scale(log(turnover_basins_delta_df$pd_exotic+0.001))
turnover_basins_delta_df$modi_area = scale(log(turnover_basins_delta_df$Surf_area))
turnover_basins_delta_df$modi_i_trade_v = scale(log(turnover_basins_delta_df$i_trade_v+0.001))
turnover_basins_delta_df$modi_geo_dist = scale(log(turnover_basins_delta_df$geo_dist))
turnover_basins_delta_df$modi_geo_dist_barriers = scale(log(turnover_basins_delta_df$geo_dist_barriers))
turnover_basins_delta_df$modi_clim_dist = scale(log(turnover_basins_delta_df$clim_dist))
turnover_basins_delta_df[is.na(turnover_basins_delta_df$Island),]$Island = 0


gls_fishes_1 = gls(f1, data = turnover_basins_delta_df)

E = residuals(gls_fishes_1)

mydata = data.frame(E, lon = turnover_basins_delta_df$Lon,
                    lat = turnover_basins_delta_df$Lat)
coordinates(mydata) = c("lon","lat")  
bubble(mydata, "E", col = c("black","grey"),
       main = "Residuals",
       xlab = "X-coordinates",  ylab = "Y-coordinates")

Vario.gls = Variogram(gls_fishes_1,
                      form =~Lon+Lat,  robust = TRUE,  resType = "pearson") 
plot(Vario.gls, smooth = TRUE)

coords = cbind(turnover_basins_delta_df$Lon, turnover_basins_delta_df$Lat)

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb = knn2nb(knearneigh(coords, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix = nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result = moran.test(turnover_basins_delta_df$mean_delta_tur, weight_matrix)
print(moran_result)  # Result note: p-value = < 2.2e-16


lme_fishes = lme(mean_delta_tur ~ modi_pd_native + modi_pd_exotic + modi_area +
                  modi_i_trade_v + 
                  modi_clim_dist, random = ~1|Country,
                #correlation = corExp(form=~Lon+Lat, nugget=T),
                data = turnover_basins_delta_df)
summary(lme_fishes)

gls_fishes = list()
gls_fishes[[1]] = gls(f1, data = turnover_basins_delta_df)

### warning!!! following code are very slow to run!
gls_fishes[[2]] = gls(f1, data = turnover_basins_delta_df,
                       correlation = corSpher(form=~Lon+Lat, nugget=T), 
                      control = glsControl(opt = "optim",
                                           maxIter = 20,
                                           msMaxIter = 20))

gls_fishes[[3]] = gls(f1, data = turnover_basins_delta_df,
                       correlation = corLin(form=~Lon+Lat, nugget=T), 
                      control = glsControl(opt = "optim",
                                           maxIter = 20,
                                           msMaxIter = 20))

gls_fishes[[4]] = gls(f1, data = turnover_basins_delta_df,
                       correlation = corRatio(form=~Lon+Lat, nugget=T),
                      control = glsControl(opt = "optim",
                                           maxIter = 20,
                                           msMaxIter = 20))

gls_fishes[[5]] = gls(f1, data = turnover_basins_delta_df,
                       correlation = corGaus(form=~Lon+Lat, nugget=T),
                      control = glsControl(opt = "optim",
                                           maxIter = 20,
                                           msMaxIter = 20))

gls_fishes[[6]] = gls(f1, data = turnover_basins_delta_df,
                       correlation = corExp(form=~Lon+Lat, nugget=T),
                      control = glsControl(opt = "optim",
                                           maxIter = 20,
                                           msMaxIter = 20))

save(gls_fishes, file = 'results/primary_results/gls_fishes_region_ED_all.rdata')

load('results/primary_results/gls_fishes_region_ED_all.rdata')
print(data.frame(SD=sapply(gls_fishes,function(x) sd(resid(x))),AIC=sapply(gls_fishes,AIC)))

summary(gls_fishes[[4]]) # best model with corGaus structure

car::vif(gls_fishes[[4]])

pred_gls_fishes = get_model_data(gls_fishes[[4]], type = 'est') 



#4. Exporting all taxon column figures----------------------------
pred_gls_plants$taxa = 'Plants'
pred_gls_mammals$taxa = 'Mammals'
pred_gls_birds$taxa = 'Birds'
pred_gls_fishes$taxa = 'Fishes'

pred_gls_all = rbind(pred_gls_plants, 
                     pred_gls_mammals,
                     pred_gls_birds,
                     pred_gls_fishes) %>%
  as_tibble()

pred_gls_all$taxa = factor(pred_gls_all$taxa,
                           levels = c("Plants",
                                      "Birds",
                                      "Mammals",
                                      "Fishes"))

pred_gls_all$term = factor(pred_gls_all$term,
                           levels = c("modi_pd_exotic",
                                      "modi_pd_native",
                                      "modi_clim_dist",
                                      'modi_area',
                                      'Island',
                                      "modi_i_trade_v"))

pred_gls_all$sig = ifelse(pred_gls_all$p.value < 0.05, 1, 0)
pred_gls_all$sig = factor(pred_gls_all$sig,
                          levels = c(1,
                                    0))

library(ggh4x)
scale_list = list(
  scale_x_continuous(limits = c(-0.033, 0.033), breaks = c(-0.03, 0, 0.03)),
  scale_x_continuous(limits = c(-0.023, 0.023), breaks = c(-0.02, 0, 0.02)),
  scale_x_continuous(limits = c(-0.053, 0.053), breaks = c(-0.05, 0, 0.05)),
  scale_x_continuous(limits = c(-0.023, 0.023), breaks = c(-0.02, 0, 0.02))
)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_char_of_all_regions = ggplot()+
  facet_wrap(vars(taxa), scales = "free_x", ncol = 2)+
  geom_point(mapping = aes(y = term, x = estimate, shape = sig),
             data = pred_gls_all,
             size = 2.5,
             color = 'black')+
  geom_linerange(mapping = aes(y = term, xmin = conf.low, xmax = conf.high),
                 data = pred_gls_all) + 
  geom_vline(xintercept = 0, linetype = 2, color = colors_4d[1])+
  scale_shape_manual(values = c(16, 1))+
  scale_x_continuous(name = 'Coefficient estimate') + 
  facetted_pos_scales(
    x = scale_list
  ) + 
  scale_y_discrete(labels = c('PD naturalized',
                              'PD native',
                              'Climatic distance',
                              'Region size',
                              'Island',
                              'Trade flow'))+
  labs(x = 'Coefficient estimate', y = '')+
  theme(#plot.title = element_blank(),
    legend.position = 'None',
    plot.title = element_text(hjust = 0.5))
plot_char_of_all_regions


#5. Export figures for paper----
library(devEMF)
library(cowplot)

emf('figures/plot_char_of_all_regions.emf',
    width = 18, height = 15, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_char_of_all_regions
dev.off() #turn off device and finalize file

