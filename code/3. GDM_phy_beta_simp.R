#Lirong Cai
#Email:lirong.cai18@gmail.com

#In this script, we modeled geographic patterns of phylogenetic and taxonomic turnover

#0 loading----------------------------
rm(list = ls())
library(reshape2)
library(plyr)
library(betareg)
library(MASS)
library(gdm)
library(vegan)
library(rgeos)
library(geosphere)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

colors_4d = c('#91989f', '#efbb24', '#1b813e', '#398fb7')
colors_10d = c('#91989f',
  colorRampPalette(c('#FFF9B0', '#FFE066', '#FFC300', '#FF9E00', '#FF7A00', '#E85C00'))(6),
  colorRampPalette(c('#9BD4E4', '#398FB7', '#004A78'))(3))

setwd("D:/R projects/Global_ED")

#1. modelling for Mammals and Birds: TDWG4----------------------------
# load environmental data and distances
load("results/primary_results/predictors_TDWG/geoentities_env.rdata")
geoentities = geoentities_env[order(geoentities_env$RegionID),]
#environmental factors
env_all = geoentities %>% dplyr::select(RegionID,Lon,Lat,CECSOL,elevation,
                                        TS,PS,AI,MAT,MAP)
summary(env_all) 

#load current geographical distance
load("results/primary_results/distances_TDWG/geodistances_grid_7.RDATA")
#load current cost distance
load("results/primary_results/distances_TDWG/graphdistances_barriers_grid_7.RDATA")
#load individual climate distances
load("results/primary_results/distances_TDWG/graphdistances_clim_grid_7.RDATA")

#load beta diversity
load('results/primary_results/distances_beta/phy_turn_mammal_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_native.rdata')

dissim = list(native_mammals = phy_turn_mammal_native,
              all_mammals = phy_turn_mammal_extant,
              native_birds = phy_turn_bird_native,
              all_birds = phy_turn_bird_extant)

# Fit model
#group to compare (maximum three group)
varSet_cur = vector("list", 2)
names(varSet_cur) = c("env", "geo")
varSet_cur$env = c("CECSOL","MAP","MAT","TS", "PS","AI")
varSet_cur$geo = c("matrix_1","matrix_2","matrix_3")


gdms_phy_TDWG = list()
groupimp_phy_TDWG = list()
imp_table_phy_TDWG = list()

for (i in 1:length(dissim)){
    #i = 1
    dissimilarities = as.matrix(dissim[[i]])
    #dissimilarities = turnover_bird_native_mat
    geo_dist = as.matrix(geodistances) 
    geo_dist_barriers = as.matrix(graphdistances_barriers)
    geo_dist_clim = as.matrix(graphdistances_clim)
 
    #subset regions
    env_tmp = env_all[which(env_all$RegionID%in%as.numeric(rownames(dissimilarities))),]
    
    env_tmp = env_tmp[complete.cases(env_tmp),]
    dissimilarities = dissimilarities[which(!is.na(dissimilarities[,2]) & is.finite(dissimilarities[,2])),]
    geo_dist = geo_dist[which(!is.na(geo_dist[,2]) & is.finite(geo_dist[,2])),]
    geo_dist_barriers = geo_dist_barriers[which(!is.na(geo_dist_barriers[,2]) & is.finite(geo_dist_barriers[,2])),]
    geo_dist_clim = geo_dist_clim[which(!is.na(geo_dist_clim[,2]) & is.finite(geo_dist_clim[,2])),]

    #harmonize the various datasets:
    env_tmp = env_tmp[which(env_tmp$RegionID %in% as.numeric(row.names(dissimilarities))
                             & env_tmp$RegionID %in% as.numeric(row.names(geo_dist))
                             & env_tmp$RegionID %in% as.numeric(row.names(geo_dist_barriers))
                             & env_tmp$RegionID %in% as.numeric(row.names(geo_dist_clim))
    ),]
    
    dissimilarities = dissimilarities[which(as.numeric(rownames(dissimilarities)) %in% env_tmp$RegionID),
                                      which(as.numeric(colnames(dissimilarities)) %in% env_tmp$RegionID)]
    dim(dissimilarities)
    geo_dist = geo_dist[which(as.numeric(rownames(geo_dist)) %in% env_tmp$RegionID),
                        which(as.numeric(colnames(geo_dist)) %in% env_tmp$RegionID)]
    dim(geo_dist)
    geo_dist_barriers = geo_dist_barriers[which(as.numeric(rownames(geo_dist_barriers)) %in% env_tmp$RegionID),
                                          which(as.numeric(colnames(geo_dist_barriers)) %in% env_tmp$RegionID)]
    dim(geo_dist_barriers)
    geo_dist_clim = geo_dist_clim[which(as.numeric(rownames(geo_dist_clim)) %in% env_tmp$RegionID),
                                  which(as.numeric(colnames(geo_dist_clim)) %in% env_tmp$RegionID)]
    dim(geo_dist_clim)
    
    dissimilarities[which(dissimilarities<0)] = 0 # few cases due to rounding
    
    dissimilarities = cbind("RegionID" = as.numeric(colnames(dissimilarities)), dissimilarities)
    geo_dist = cbind("RegionID" = as.numeric(colnames(geo_dist)), geo_dist)
    geo_dist_barriers = cbind("RegionID" = as.numeric(colnames(geo_dist_barriers)), geo_dist_barriers)
    geo_dist_clim = cbind("RegionID" = as.numeric(colnames(geo_dist_clim)), geo_dist_clim)
  
    
    #select environmental data
    env_tmp_cur = env_tmp %>% dplyr::select(RegionID,Lon,Lat,
                                          MAT,MAP,TS,PS,AI)

    #Produce a GDM-formatted Site-Pair Table
    gdm_tab_tmp_cur = formatsitepair(bioData = dissimilarities, bioFormat = 3, siteColumn = "RegionID",
                                     abundance = F, XColumn = "Lon", YColumn = "Lat",
                                     predData = env_tmp_cur, 
                                     distPreds = list("Geo_Distance" = geo_dist,
                                                      "Geo_Distance_Barriers" = geo_dist_barriers,
                                                      "Geo_Distance_clim" = geo_dist_clim
                                      ))
    

    #fit gdms
    gdms_phy_TDWG[[i]] = gdm(gdm_tab_tmp_cur, geo = F)

    
    
    # relative variable importance
    spl = isplineExtract(gdms_phy_TDWG[[i]])
    
    max_height = apply(spl$y, 2, max, na.rm = TRUE)
    max_height
    
    w_rel = max_height / sum(max_height)
    w_rel 
    
    dev_pct = gdms_phy_TDWG[[i]]$explained/100         
    
    imp_scaled = w_rel * dev_pct       
    imp_scaled = c(Unexplained = 1-dev_pct, imp_scaled)
    
    imp_table_phy_TDWG[[i]] = data.frame(
      predictor        = c('Unexplained', names(max_height)),
      max_height       = c(NA, max_height),
      rel_importance   = c(NA, w_rel),     
      deviance_percent = imp_scaled,
      taxon = str_split(names(dissim)[i], '_')[[1]][2],
      status = str_split(names(dissim)[i], '_')[[1]][1]
    )
    
    
    #partition.deviance
    groupimp_phy_TDWG[[i]]=cbind(gdm.partition.deviance(sitePairTable=gdm_tab_tmp_cur,
                                             varSet_cur, partSpace=FALSE),
    taxon = str_split(names(dissim)[i], '_')[[1]][2],
    status = str_split(names(dissim)[i], '_')[[1]][1])

    
  }

names(gdms_phy_TDWG) = c('native_mammals',
                            'all_mammals',
                            'native_birds',
                            'all_birds')
names(groupimp_phy_TDWG) = c('native_mammals',
                                'all_mammals',
                                'native_birds',
                                'all_birds')
names(imp_table_phy_TDWG) = c('native_mammals',
                             'all_mammals',
                             'native_birds',
                             'all_birds')

#2. modelling for Plants: GIFT----------------------------
# load environmental data and distances
gc()
shp_glonaf_new = st_read("data/Plants/shp_glonaf_new_eck4.shp")
load("results/primary_results/predictors_GIFT/geoentities_env.rdata")
geoentities = geoentities_env[order(geoentities_env$Region_id),] %>% 
  left_join((shp_glonaf_new %>% st_drop_geometry())[,c('Region_id',
                                                       'Lon',
                                                       'Lat')],
            by = 'Region_id')
#environmental factors
env_all = geoentities %>% dplyr::select(Region_id, Lon, Lat, CECSOL, elevation,
                                        TS,PS,AI,MAT,MAP)
summary(env_all) 

#load current geographical distance
load("results/primary_results/distances_GIFT/geodistances_grid_7.RDATA")
#load current cost distance
load("results/primary_results/distances_GIFT/graphdistances_barriers_grid_7.RDATA")
#load individual climate distances
load("results/primary_results/distances_GIFT/graphdistances_clim_grid_7.RDATA")

#load beta diversity
load('results/primary_results/distances_beta/phy_turn_plant_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_plant_native.rdata')

# Fit model
dissim = list(native_plants = phy_turn_plant_native,
              all_plants = phy_turn_plant_extant)

#group to compare (maximum three group)
varSet_cur = vector("list", 2)
names(varSet_cur) = c("env", "geo")
varSet_cur$env = c("CECSOL","MAP","MAT","TS", "PS","AI")
varSet_cur$geo = c("matrix_1","matrix_2","matrix_3")

gdms_phy_GIFT = list()
groupimp_phy_GIFT = list()
imp_table_phy_GIFT = list()

for (i in 1:length(dissim)){
  #i = 1
  dissimilarities = as.matrix(dissim[[i]])
  #dissimilarities = turnover_bird_native_mat
  geo_dist = as.matrix(geodistances) 
  geo_dist_barriers = as.matrix(graphdistances_barriers)
  geo_dist_clim = as.matrix(graphdistances_clim)
  
  #subset regions
  env_tmp = env_all[which(env_all$Region_id%in%as.numeric(rownames(dissimilarities))),]
  
  env_tmp = env_tmp[complete.cases(env_tmp),]
  dissimilarities = dissimilarities[which(!is.na(dissimilarities[,2]) & is.finite(dissimilarities[,2])),]
  geo_dist = geo_dist[which(!is.na(geo_dist[,2]) & is.finite(geo_dist[,2])),]
  geo_dist_barriers = geo_dist_barriers[which(!is.na(geo_dist_barriers[,2]) & is.finite(geo_dist_barriers[,2])),]
  geo_dist_clim = geo_dist_clim[which(!is.na(geo_dist_clim[,2]) & is.finite(geo_dist_clim[,2])),]
  
  #harmonize the various datasets:
  env_tmp = env_tmp[which(env_tmp$Region_id %in% as.numeric(row.names(dissimilarities))
                          & env_tmp$Region_id %in% as.numeric(row.names(geo_dist))
                          & env_tmp$Region_id %in% as.numeric(row.names(geo_dist_barriers))
                          & env_tmp$Region_id %in% as.numeric(row.names(geo_dist_clim))
  ),]
  
  dissimilarities = dissimilarities[which(as.numeric(rownames(dissimilarities)) %in% env_tmp$Region_id),
                                    which(as.numeric(colnames(dissimilarities)) %in% env_tmp$Region_id)]
  dim(dissimilarities)
  geo_dist = geo_dist[which(as.numeric(rownames(geo_dist)) %in% env_tmp$Region_id),
                      which(as.numeric(colnames(geo_dist)) %in% env_tmp$Region_id)]
  dim(geo_dist)
  geo_dist_barriers = geo_dist_barriers[which(as.numeric(rownames(geo_dist_barriers)) %in% env_tmp$Region_id),
                                        which(as.numeric(colnames(geo_dist_barriers)) %in% env_tmp$Region_id)]
  dim(geo_dist_barriers)
  geo_dist_clim = geo_dist_clim[which(as.numeric(rownames(geo_dist_clim)) %in% env_tmp$Region_id),
                                which(as.numeric(colnames(geo_dist_clim)) %in% env_tmp$Region_id)]
  dim(geo_dist_clim)
  
  dissimilarities[which(dissimilarities<0)] = 0 # few cases due to rounding
  
  dissimilarities = cbind("Region_id" = as.numeric(colnames(dissimilarities)), dissimilarities)
  geo_dist = cbind("Region_id" = as.numeric(colnames(geo_dist)), geo_dist)
  geo_dist_barriers = cbind("Region_id" = as.numeric(colnames(geo_dist_barriers)), geo_dist_barriers)
  geo_dist_clim = cbind("Region_id" = as.numeric(colnames(geo_dist_clim)), geo_dist_clim)
  
  
  #select environmental data
  env_tmp_cur = env_tmp%>%dplyr::select(Region_id,Lon,Lat,CECSOL,
                                        MAT,MAP,TS,PS,AI)
  
  #Produce a GDM-formatted Site-Pair Table
  gdm_tab_tmp_cur = formatsitepair(bioData = dissimilarities, bioFormat = 3, siteColumn = "Region_id",
                                   abundance = F, XColumn = "Lon", YColumn = "Lat",
                                   predData = env_tmp_cur, 
                                   distPreds = list("Geo_Distance" = geo_dist,
                                                    "Geo_Distance_Barriers" = geo_dist_barriers,
                                                    "Geo_Distance_clim" = geo_dist_clim
                                   ))
  
  
  #fit gdms
  gdms_phy_GIFT[[i]] = gdm(gdm_tab_tmp_cur, geo = F)
  
  # relative variable importance
  spl = isplineExtract(gdms_phy_GIFT[[i]])

  max_height = apply(spl$y, 2, max, na.rm = TRUE)
  max_height

  w_rel = max_height / sum(max_height)
  w_rel 
    
  dev_pct = gdms_phy_GIFT[[i]]$explained/100         
  
  imp_scaled = w_rel * dev_pct       
  imp_scaled = c(Unexplained = 1-dev_pct, imp_scaled)
  
  imp_table_phy_GIFT[[i]] = data.frame(
    predictor        = c('Unexplained', names(max_height)),
    max_height       = c(NA, max_height),
    rel_importance   = c(NA, w_rel),     
    deviance_percent = imp_scaled,
    taxon = str_split(names(dissim)[i], '_')[[1]][2],
    status = str_split(names(dissim)[i], '_')[[1]][1]
  )
  
  #partition.deviance
  groupimp_phy_GIFT[[i]]=cbind(gdm.partition.deviance(sitePairTable=gdm_tab_tmp_cur,
                                                 varSet_cur, partSpace=FALSE),
                          taxon = str_split(names(dissim)[i], '_')[[1]][2],
                          status = str_split(names(dissim)[i], '_')[[1]][1])
  
}


names(gdms_phy_GIFT) = c('native_plants',
                         'all_plants')
names(groupimp_phy_GIFT) = c('native_plants',
                             'all_plants')
names(imp_table_phy_GIFT) = c('native_plants',
                             'all_plants')





#3. modelling for fishes: Basin----------------------------
# load environmental data and distances
df_fish = st_read("data/Fishes/data/Basin042017_3119.shp")
load("results/primary_results/predictors_basins/geoentities_env.rdata")

#environmental factors
colnames(geoentities_env)
env_all = geoentities_env %>% dplyr::select(BasinName, Lon, Lat,elevation,
                                        TS,PS,AI,MAT,MAP)

summary(env_all) 

#load current geographical distance
load("results/primary_results/distances_basins/geodistances_grid_7.RDATA")
#load current cost distance
load("results/primary_results/distances_basins/graphdistances_barriers_grid_7.RDATA")
#load individual climate distances
load("results/primary_results/distances_basins/graphdistances_clim_grid_7.RDATA")
geo_dist = as.matrix(geodistances) 
geo_dist_barriers = as.matrix(graphdistances_barriers)
geo_dist_clim = as.matrix(graphdistances_clim)

#load beta diversity
load('results/primary_results/distances_beta/phy_turn_fish_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_fish_native.rdata')

geo_dist = geo_dist[colnames(phy_turn_fish_extant), colnames(phy_turn_fish_extant)]
geo_dist_barriers = geo_dist_barriers[colnames(phy_turn_fish_extant), colnames(phy_turn_fish_extant)]
geo_dist_clim = geo_dist_clim[colnames(phy_turn_fish_extant), colnames(phy_turn_fish_extant)]

identical(colnames(geo_dist_barriers), colnames(geo_dist))
identical(colnames(geo_dist_clim), colnames(geo_dist))
identical(colnames(phy_turn_fish_extant), colnames(geo_dist))

env_all = env_all %>%
  slice(match(colnames(geo_dist), BasinName)) %>% 
  filter(BasinName %in% colnames(geo_dist))

identical(env_all$BasinName, colnames(geo_dist))

# add the numeric basin id for fitting models
colnames(geo_dist) = c(1:nrow(env_all))
rownames(geo_dist) = c(1:nrow(env_all))
colnames(geo_dist_barriers) = c(1:nrow(env_all))
rownames(geo_dist_barriers) = c(1:nrow(env_all))
colnames(geo_dist_clim) = c(1:nrow(env_all))
rownames(geo_dist_clim) = c(1:nrow(env_all))
env_all$BasinID = c(1:nrow(env_all))
colnames(phy_turn_fish_extant) = c(1:nrow(env_all))
rownames(phy_turn_fish_extant) = c(1:nrow(env_all))
colnames(phy_turn_fish_native) = c(1:nrow(env_all))
rownames(phy_turn_fish_native) = c(1:nrow(env_all))

# Fit model
dissim = list(native_fishes = phy_turn_fish_native,
              all_fishes = phy_turn_fish_extant)

#group to compare (maximum three group)
varSet_cur = vector("list", 2)
names(varSet_cur) = c("env", "geo")
varSet_cur$env = c("MAP","MAT","TS", "PS","AI")
varSet_cur$geo = c("matrix_1","matrix_2","matrix_3")

gdms_phy_Basin = list()
groupimp_phy_Basin = list()
imp_table_phy_Basin = list()

for (i in 1:length(dissim)){
  #i = 1
  dissimilarities = as.matrix(dissim[[i]])
  #dissimilarities = turnover_bird_native_mat

  #subset regions
  env_tmp = env_all[which(env_all$BasinID%in%as.numeric(rownames(dissimilarities))),]
  
  env_tmp = env_tmp[complete.cases(env_tmp),]
  dissimilarities = dissimilarities[which(!is.na(dissimilarities[,2]) & is.finite(dissimilarities[,2])),]
  geo_dist = geo_dist[which(!is.na(geo_dist[,2]) & is.finite(geo_dist[,2])),]
  geo_dist_barriers = geo_dist_barriers[which(!is.na(geo_dist_barriers[,2]) & is.finite(geo_dist_barriers[,2])),]
  geo_dist_clim = geo_dist_clim[which(!is.na(geo_dist_clim[,2]) & is.finite(geo_dist_clim[,2])),]
  
  #harmonize the various datasets:
  env_tmp = env_tmp[which(env_tmp$BasinID %in% as.numeric(row.names(dissimilarities))
                          & env_tmp$BasinID %in% as.numeric(row.names(geo_dist))
                          & env_tmp$BasinID %in% as.numeric(row.names(geo_dist_barriers))
                          & env_tmp$BasinID %in% as.numeric(row.names(geo_dist_clim))
  ),]
  
  dissimilarities = dissimilarities[which(as.numeric(rownames(dissimilarities)) %in% env_tmp$BasinID),
                                    which(as.numeric(colnames(dissimilarities)) %in% env_tmp$BasinID)]
  dim(dissimilarities)
  geo_dist = geo_dist[which(as.numeric(rownames(geo_dist)) %in% env_tmp$BasinID),
                      which(as.numeric(colnames(geo_dist)) %in% env_tmp$BasinID)]
  dim(geo_dist)
  geo_dist_barriers = geo_dist_barriers[which(as.numeric(rownames(geo_dist_barriers)) %in% env_tmp$BasinID),
                                        which(as.numeric(colnames(geo_dist_barriers)) %in% env_tmp$BasinID)]
  dim(geo_dist_barriers)
  geo_dist_clim = geo_dist_clim[which(as.numeric(rownames(geo_dist_clim)) %in% env_tmp$BasinID),
                                which(as.numeric(colnames(geo_dist_clim)) %in% env_tmp$BasinID)]
  dim(geo_dist_clim)
  
  dissimilarities[which(dissimilarities<0)] = 0 # few cases due to rounding
  
  dissimilarities_all = cbind("BasinID" = as.numeric(colnames(dissimilarities)), dissimilarities)
  geo_dist_all = cbind("BasinID" = as.numeric(colnames(geo_dist)), geo_dist)
  geo_dist_barriers_all = cbind("BasinID" = as.numeric(colnames(geo_dist_barriers)), geo_dist_barriers)
  geo_dist_clim_all = cbind("BasinID" = as.numeric(colnames(geo_dist_clim)), geo_dist_clim)
  
  #select environmental data
  env_tmp_cur = env_tmp %>% dplyr::select(BasinID,Lon,Lat, #elevation,
                                        MAT,MAP,TS,PS,AI)
  
  #Produce a GDM-formatted Site-Pair Table
  gdm_tab_tmp_cur = formatsitepair(bioData = dissimilarities_all, bioFormat = 3, siteColumn = "BasinID",
                                   abundance = F, XColumn = "Lon", YColumn = "Lat",
                                   predData = env_tmp_cur, 
                                   distPreds = list("Geo_Distance" = geo_dist_all,
                                                    "Geo_Distance_Barriers" = geo_dist_barriers_all,
                                                    "Geo_Distance_clim" = geo_dist_clim_all
                                   ))
  
  
  #fit gdms
  gdms_phy_Basin[[i]] = gdm(gdm_tab_tmp_cur, geo = F)
  
  # relative variable importance
  spl = isplineExtract(gdms_phy_Basin[[i]])
  
  max_height = apply(spl$y, 2, max, na.rm = TRUE)
  max_height
  
  w_rel = max_height / sum(max_height)
  w_rel 
  
  dev_pct = gdms_phy_Basin[[i]]$explained/100         
  
  imp_scaled = w_rel * dev_pct       
  imp_scaled = c(Unexplained = 1-dev_pct, imp_scaled)
  
  imp_table_phy_Basin[[i]] = data.frame(
    predictor        = c('Unexplained', names(max_height)),
    max_height       = c(NA, max_height),
    rel_importance   = c(NA, w_rel),     
    deviance_percent = imp_scaled,
    taxon = str_split(names(dissim)[i], '_')[[1]][2],
    status = str_split(names(dissim)[i], '_')[[1]][1]
  )
  
  #partition.deviance
  groupimp_phy_Basin[[i]]=cbind(gdm.partition.deviance(sitePairTable=gdm_tab_tmp_cur,
                                                      varSet_cur, partSpace=FALSE),
                               taxon = str_split(names(dissim)[i], '_')[[1]][2],
                               status = str_split(names(dissim)[i], '_')[[1]][1])
  
}


names(gdms_phy_Basin) = c('native_fishes',
                         'all_fishes')
names(groupimp_phy_Basin) = c('native_fishes',
                             'all_fishes')
names(imp_table_phy_Basin) = c('native_fishes',
                              'all_fishes')



#4. Exporting all taxon column figures----------------------------
groupimp_phy_dat = rbind(data.table::rbindlist(groupimp_phy_GIFT) %>%
  as_tibble() %>% 
  filter(!(variableSet %in% c('env',
                              'geo',
                              'ALL VARIABLES (env & geo)'))), 
  data.table::rbindlist(groupimp_phy_TDWG) %>%
    as_tibble() %>% 
    filter(!(variableSet %in% c('env',
                                'geo',
                                'ALL VARIABLES (env & geo)'))),
  data.table::rbindlist(groupimp_phy_Basin) %>%
    as_tibble() %>% 
    filter(!(variableSet %in% c('env',
                                'geo',
                                'ALL VARIABLES (env & geo)'))))
groupimp_phy_dat[groupimp_phy_dat$taxon == 'plants',]$taxon = 'Plants'
groupimp_phy_dat[groupimp_phy_dat$taxon == 'birds',]$taxon  = 'Birds'
groupimp_phy_dat[groupimp_phy_dat$taxon == 'mammals',]$taxon  = 'Mammals'
groupimp_phy_dat[groupimp_phy_dat$taxon == 'fishes',]$taxon  = 'Fishes'

groupimp_phy_dat$variableSet = factor(groupimp_phy_dat$variableSet,
                                           levels = c("UNEXPLAINED",
                                                      "env alone",
                                                      "env intersect geo",
                                                      "geo alone"))
groupimp_phy_dat$status = factor(groupimp_phy_dat$status,
                                      levels = c("native",
                                                 "all"))
groupimp_phy_dat$taxon = factor(groupimp_phy_dat$taxon,
                                 levels = c("Plants",
                                            "Birds",
                                            "Mammals",
                                            "Fishes"
                                            ))

groupimp_phy_dat$deviance = groupimp_phy_dat$deviance/100

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_groupimp_phy = ggplot()+
  facet_wrap(vars(taxon), scales = "free", nrow = 1)+
  geom_col(aes(x = status, y = deviance, group = variableSet, fill = variableSet),
           data = groupimp_phy_dat) + 
  #geom_hline(yintercept = 0, linetype = 2, color = 'black')+
  labs(x = " ", y = "Proportion of deviance explained") +
  scale_fill_manual(values = colors_4d, 
                    labels = c('Unexplained',
                               'Enviromental only',
                               'Shared',
                               'Dispersal only'))+
  scale_color_manual(values = colors_4d, 
                     labels = c('Unexplained',
                                'Enviromental only',
                                'Shared',
                                'Dispersal only'))+
  guides(colour = guide_legend(title=' '),
         fill = guide_legend(title=' '))+
  theme(panel.grid.minor = element_blank(),
        legend.position = 'right')



relimp_phy_dat = rbind(data.table::rbindlist(imp_table_phy_GIFT) %>%
                           as_tibble() , 
                         data.table::rbindlist(imp_table_phy_TDWG) %>%
                           as_tibble(),
                       data.table::rbindlist(imp_table_phy_Basin) %>%
                         as_tibble())
relimp_phy_dat[relimp_phy_dat$taxon == 'plants',]$taxon = 'Plants'
relimp_phy_dat[relimp_phy_dat$taxon == 'birds',]$taxon  = 'Birds'
relimp_phy_dat[relimp_phy_dat$taxon == 'mammals',]$taxon  = 'Mammals'
relimp_phy_dat[relimp_phy_dat$taxon == 'fishes',]$taxon  = 'Fishes'

relimp_phy_dat$status = factor(relimp_phy_dat$status,
                                 levels = c("native",
                                            "all"))
relimp_phy_dat$taxon = factor(relimp_phy_dat$taxon,
                                levels = c("Plants",
                                           "Birds",
                                           "Mammals",
                                           "Fishes"
                                ))

relimp_phy_dat[relimp_phy_dat$predictor == 'matrix_1',]$predictor = 'Geographical linear distance'
relimp_phy_dat[relimp_phy_dat$predictor == 'matrix_2',]$predictor = 'Least-cost distance across barriers'
relimp_phy_dat[relimp_phy_dat$predictor == 'matrix_3',]$predictor = 'Least-cost distance across climatic dissimilarity surface'
relimp_phy_dat$predictor = factor(relimp_phy_dat$predictor,
                                  levels = c("Unexplained","MAT", "MAP", "TS", "PS",
                                             "AI","CECSOL",
                                             "Geographical linear distance",
                                             "Least-cost distance across barriers",
                                             "Least-cost distance across climatic dissimilarity surface"))

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_relimp_phy = ggplot()+
  facet_wrap(vars(taxon), scales = "free", nrow = 1)+
  geom_col(aes(x = status, y = deviance_percent, group = predictor, fill = predictor),
           data = relimp_phy_dat) + 
  #geom_hline(yintercept = 0, linetype = 2, color = 'black')+
  labs(x = " ", y = "Relative variable importance") +
  scale_fill_manual(values=colors_10d, 
                    labels = c('Unexplained',
                               'Mean auunal temperature',
                               'Mean auunal precipitation',
                               'Temperature seasonatity',
                               'Precipitation seasonatity',
                               'Aridity index',
                               'Cation-exchange capicity',
                               'Geographical linear distance',
                               "Least-cost distance across barriers",
                               "Least-cost distance across climatic dissimilarity surface"))+
  scale_color_manual(values=colors_10d, 
                     labels = c('Unexplained',
                                'Mean auunal temperature',
                                'Mean auunal precipitation',
                                'Temperature seasonatity',
                                'Precipitation seasonatity',
                                'Aridity index',
                                'Cation-exchange capicity',
                                'Geographical linear distance',
                                "Least-cost distance across barriers",
                                "Least-cost distance across climatic dissimilarity surface"))+
  guides(colour = guide_legend(title=' '),
         fill = guide_legend(title=' '))+
  theme(panel.grid.minor = element_blank(),
        legend.position = 'right')



#5. Export figures for paper----
library(devEMF)
library(cowplot)
library(ggpubr)

plots_gdm_phy = ggdraw() +
  # Row 1: space1 plot full width
  # Row 1: dir | indir | percent
  draw_plot(plot_relimp_phy,
            x = 0, y = 0.08,
            width = 1.01, height = 0.48) +
  draw_plot(plot_groupimp_phy,
            x = 0, y = 0.5,
            width = 0.75, height = 0.48) +
  # Optional labels
  draw_label("a", x = 0.01, y = 0.971, size = 14.5, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("b", x = 0.01, y = 0.531, size = 14.5, fontface = "bold", hjust = 0, vjust = 1)


emf('figures/plots_gdm_phy.emf',
    width = 25*0.9, height = 25*0.75, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)

plots_gdm_phy

dev.off() #turn off device and finalize file
