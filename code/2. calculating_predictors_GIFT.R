### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#0 loading required packages-------------------
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("foreach", "dplyr", 'tidyr', 'sf', 'terra', 'raster', 'scico',
                 'GIFT')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("D:/R projects/Global_ED")

# load the background maps for quantifying the predictors
df_eck4 = st_read("data/Plants/shp_glonaf_new_eck4_1.shp")
df_eck4[df_eck4$Region_id == 764,]$GID_0 = 'JPN'
df_eck4[df_eck4$Region_id == 764,]$NAME_0 = 'Japan'
df_eck4[df_eck4$Region_id == 453,]$GID_0 = 'PHL'
df_eck4[df_eck4$Region_id == 453,]$NAME_0 = 'Philippines'
df_eck4[df_eck4$Region_id == 770,]$GID_0 = 'RUS'
df_eck4[df_eck4$Region_id == 770,]$NAME_0 = 'Russian Federation'

df_country = df_eck4 %>% st_drop_geometry()
df = st_read("data/Plants/shp_glonaf_new.shp")
df = df %>% dplyr::left_join(df_country[,c('Region_id', 'GID_0', 'NAME_0')],
                      by = 'Region_id')
rm(df_eck4)
str(df)
sort(unique(df$NAME_0))
sort(unique(df$GID_0))

#1. Environemt variables----
##1.1 elevation and soil cation exchange ability----
# read the elevation dataset
ces = terra::rast("data/envir_vars/bioclim/CECSOL_M_sl2_250m_ll.tif")
elevation = terra::rast("data/envir_vars/gmted2010/GMTED2010_15n015_00625deg.nc")
# accessed at https://temis.nl/data/gmted2010/, unit = m
elevation$Elev_range = elevation$elevation_max - elevation$elevation_min

# make projection coincidence
df_for_eleva = sf::st_transform(df, crs = sf::st_crs(elevation))
df_for_ces = sf::st_transform(df, crs = sf::st_crs(ces))

# convert to SpatVector and do zonal
v_for_eleva = terra::vect(df_for_eleva)
# could be clip to high speed
elevation_crop = terra::crop(elevation, v_for_eleva, snap="out")

# zonal: by=Region_id, fun="mean"
z_for_eleva = cbind(Region_id = v_for_eleva$Region_id,
                    terra::extract(elevation_crop,
                                   v_for_eleva, fun = "mean", na.rm = TRUE))
rm(elevation_crop, v_for_eleva)

## ces data is too large, thus run this treat in Arcgis 10.8
#v_for_ces = terra::vect(df_for_ces)
# could be clip to high speed
#ces_crop = terra::crop(ces, v_for_ces, snap="out")
#z_for_ces_l = lapply(1:nrow(v_for_ces), function(i) {
#i = 1
# gc()
#terra::extract(ces_crop, v_for_ces[i,], fun = 'mean', na.rm = TRUE)
#})
#z_for_ces = do.call(rbind, results)
z_for_ces = read.table('results/primary_results/predictors_GIFT/geo_ces.txt',
                       header = T, sep = ',')
colnames(z_for_ces)[which(colnames(z_for_ces) == 'MEAN')] = 'CECSOL'

##1.2 Five main climate variables----
# accessed at https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.zhdk.cloud.switch.ch%2Fchelsav2%2F&prefix=GLOBAL%2Fclimatologies%2F1981-2010%2F
# read the bioclim multi-terra 
bioclim1 = terra::rast("data/envir_vars/bioclim/CHELSA_bio1_1981-2010_V.2.1.tif")
bioclim4 = terra::rast("data/envir_vars/bioclim/CHELSA_bio4_1981-2010_V.2.1.tif")
bioclim12 = terra::rast("data/envir_vars/bioclim/CHELSA_bio12_1981-2010_V.2.1.tif")
bioclim15 = terra::rast("data/envir_vars/bioclim/CHELSA_bio15_1981-2010_V.2.1.tif")
arid_index = terra::rast("data/envir_vars/bioclim/CHELSA_ai_1981-2010_V.2.1.tif")

# make projection coincidence
df_for_bioclim = sf::st_transform(df, crs = sf::st_crs(bioclim1))

# convert to SpatVector and do zonal
v_for_bioclim = terra::vect(df_for_bioclim)

# could be clip to high speed
bioclim1_crop = terra::crop(bioclim1, v_for_bioclim, snap="out")
bioclim4_crop = terra::crop(bioclim4, v_for_bioclim, snap="out")
bioclim12_crop = terra::crop(bioclim12, v_for_bioclim, snap="out")
bioclim15_crop = terra::crop(bioclim15, v_for_bioclim, snap="out")
arid_index_crop = terra::crop(arid_index, v_for_bioclim, snap="out")

# zonal: by=Region_id, fun="mean"
z_for_bioclim1 = terra::zonal(bioclim1_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
rm(bioclim1_crop)
gc()
z_for_bioclim4 = terra::zonal(bioclim4_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
rm(bioclim4_crop)
gc()
z_for_bioclim12 = terra::zonal(bioclim12_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
rm(bioclim12_crop)
gc()
z_for_bioclim15 = terra::zonal(bioclim15_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
rm(bioclim15_crop)
gc()
z_for_arid_index = terra::zonal(arid_index_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
rm(arid_index_crop)
gc()

### aggregate and save the 7 environmental vairables
geoentities_env = df %>% 
  left_join(z_for_eleva, by = 'Region_id') %>% 
  left_join(z_for_ces[,c('REGION_ID', 'CECSOL')],
            by = join_by('Region_id' == 'REGION_ID')) %>% 
  st_drop_geometry()

geoentities_env = cbind(geoentities_env,
                        MAT = z_for_bioclim1$`CHELSA_bio1_1981-2010_V.2.1`,
                        TS = z_for_bioclim4$`CHELSA_bio4_1981-2010_V.2.1`,
                        MAP = z_for_bioclim12$`CHELSA_bio12_1981-2010_V.2.1`,
                        PS = z_for_bioclim15$`CHELSA_bio15_1981-2010_V.2.1`,
                        AI = z_for_arid_index$`CHELSA_ai_1981-2010_V.2.1`)

summary(geoentities_env)
save(geoentities_env,
     file = 'results/primary_results/predictors_GIFT/geoentities_env.rdata')


#### calculate climate distance among different regions
load('results/primary_results/predictors_GIFT/geoentities_env.rdata')
clim_dist = as.matrix(dist(geoentities_env %>% dplyr::select(MAT, TS, MAP, PS, AI),
                 method = 'euclidean', diag = T))
rownames(clim_dist) = geoentities_env$Region_id
colnames(clim_dist) = geoentities_env$Region_id
save(clim_dist,
     file = 'results/primary_results/distances_GIFT/clim_dist.rdata')


#2. dispersal limitation ----
###(least cost distance considering physical, ecological dispersal barriers, and climate similarity)

##2.3 Calculating the distances at the regional level -----
rm(list=ls())
gc()
library(maptools)
library(rgeos)
library(geosphere)
library(rgdal)
library(dplyr)

#selected entity
geoentities_simple = readOGR(dsn = "data/Plants/shp_glonaf_new.shp")
#load equal area grid
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7.RData")

geoentities_simple@proj4string
grid_7@proj4string

geoentities_simple@proj4string = grid_7@proj4string
overlap = over(geoentities_simple, grid_7, returnList = TRUE)

###2.3.1 Pure geodistances-----
#load grid_7 data
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/geodistances_grid_7.RDATA")
# Convert vector of distances/ too large dist object into a full matrix
geodistances_full= matrix(NA, nrow(grid_7), nrow(grid_7))
rownames(geodistances_full) = grid_7$grid_7_ID
colnames(geodistances_full) = grid_7$grid_7_ID

geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = geodistances
geodistances_full = t(geodistances_full)
geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = geodistances
diag(geodistances_full) = 0

geodistances_full[1:5,1:5]

geodistances_GIFT = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)
colnames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_GIFT[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_GIFT[upper.tri(geodistances_GIFT,
                            diag=FALSE)] = t(geodistances_GIFT)[upper.tri(geodistances_GIFT,
                                                                          diag=FALSE)]
geodistances_GIFT[1:5,1:5]

geodistances_GIFT_dist = geodistances_GIFT[order(as.numeric(row.names(geodistances_GIFT))),
                                           order(as.numeric(colnames(geodistances_GIFT)))]
geodistances_GIFT_dist[1:5,1:5]
geodistances_GIFT_dist = as.dist(geodistances_GIFT_dist)

geodistances = geodistances_GIFT_dist
save(geodistances,
     file="results/primary_results/distances_GIFT/geodistances_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full,geodistances_GIFT,geodistances_GIFT_dist,i,k)

###2.3.2 water-----
#load grid_7 data
gc()
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_water_grid_7.RDATA")
# Convert vector of distances/ too large dist object into a full matrix
geodistances_full= matrix(NA, nrow(grid_7), nrow(grid_7))
rownames(geodistances_full) = grid_7$grid_7_ID
colnames(geodistances_full) = grid_7$grid_7_ID

geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_water
geodistances_full = t(geodistances_full)
geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_water
diag(geodistances_full) = 0

geodistances_full[1:5,1:5]

geodistances_GIFT = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)
colnames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_GIFT[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_GIFT[upper.tri(geodistances_GIFT,
                            diag=FALSE)] = t(geodistances_GIFT)[upper.tri(geodistances_GIFT,
                                                                          diag=FALSE)]
geodistances_GIFT[1:5,1:5]

geodistances_GIFT_dist = geodistances_GIFT[order(as.numeric(row.names(geodistances_GIFT))),
                                           order(as.numeric(colnames(geodistances_GIFT)))]
geodistances_GIFT_dist[1:5,1:5]
geodistances_GIFT_dist = as.dist(geodistances_GIFT_dist)

graphdistances_water = geodistances_GIFT_dist
save(graphdistances_water,
     file="results/primary_results/distances_GIFT/graphdistances_water_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_water, geodistances_GIFT,geodistances_GIFT_dist,i,k)


###2.3.3 elev-----
gc()
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_elev_grid_7.RDATA")
# Convert vector of distances/ too large dist object into a full matrix
geodistances_full= matrix(NA, nrow(grid_7), nrow(grid_7))
rownames(geodistances_full) = grid_7$grid_7_ID
colnames(geodistances_full) = grid_7$grid_7_ID

geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_elev
geodistances_full = t(geodistances_full)
geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_elev
diag(geodistances_full) = 0

geodistances_full[1:5,1:5]

geodistances_GIFT = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)
colnames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_GIFT[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_GIFT[upper.tri(geodistances_GIFT,
                            diag=FALSE)] = t(geodistances_GIFT)[upper.tri(geodistances_GIFT,
                                                                          diag=FALSE)]
geodistances_GIFT[1:5,1:5]

geodistances_GIFT_dist = geodistances_GIFT[order(as.numeric(row.names(geodistances_GIFT))),
                                           order(as.numeric(colnames(geodistances_GIFT)))]
geodistances_GIFT_dist[1:5,1:5]
geodistances_GIFT_dist = as.dist(geodistances_GIFT_dist)

graphdistances_elev = geodistances_GIFT_dist
save(graphdistances_elev,
     file="results/primary_results/distances_GIFT/graphdistances_elev_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_elev, geodistances_GIFT, geodistances_GIFT_dist,i,k)


###2.3.4 temp-----
gc()
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_temp_grid_7.RDATA")
# Convert vector of distances/ too large dist object into a full matrix
geodistances_full= matrix(NA, nrow(grid_7), nrow(grid_7))
rownames(geodistances_full) = grid_7$grid_7_ID
colnames(geodistances_full) = grid_7$grid_7_ID

geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_temp
geodistances_full = t(geodistances_full)
geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_temp
diag(geodistances_full) = 0

geodistances_full[1:5,1:5]

geodistances_GIFT = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)
colnames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_GIFT[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_GIFT[upper.tri(geodistances_GIFT,
                            diag=FALSE)] = t(geodistances_GIFT)[upper.tri(geodistances_GIFT,
                                                                          diag=FALSE)]
geodistances_GIFT[1:5,1:5]

geodistances_GIFT_dist = geodistances_GIFT[order(as.numeric(row.names(geodistances_GIFT))),
                                           order(as.numeric(colnames(geodistances_GIFT)))]
geodistances_GIFT_dist[1:5,1:5]
geodistances_GIFT_dist = as.dist(geodistances_GIFT_dist)

graphdistances_temp = geodistances_GIFT_dist
save(graphdistances_temp,
     file="results/primary_results/distances_GIFT/graphdistances_temp_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_temp, geodistances_GIFT, geodistances_GIFT_dist,i,k)

###2.3.5 ai-----
gc()
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_ai_grid_7.RDATA")
# Convert vector of distances/ too large dist object into a full matrix
geodistances_full= matrix(NA, nrow(grid_7), nrow(grid_7))
rownames(geodistances_full) = grid_7$grid_7_ID
colnames(geodistances_full) = grid_7$grid_7_ID

geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_ai
geodistances_full = t(geodistances_full)
geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_ai
diag(geodistances_full) = 0

geodistances_full[1:5,1:5]

geodistances_GIFT = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)
colnames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_GIFT[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_GIFT[upper.tri(geodistances_GIFT,
                            diag=FALSE)] = t(geodistances_GIFT)[upper.tri(geodistances_GIFT,
                                                                          diag=FALSE)]
geodistances_GIFT[1:5,1:5]

geodistances_GIFT_dist = geodistances_GIFT[order(as.numeric(row.names(geodistances_GIFT))),
                                           order(as.numeric(colnames(geodistances_GIFT)))]
geodistances_GIFT_dist[1:5,1:5]
geodistances_GIFT_dist = as.dist(geodistances_GIFT_dist)

graphdistances_ai = geodistances_GIFT_dist
save(graphdistances_ai,
     file="results/primary_results/distances_GIFT/graphdistances_ai_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_ai, geodistances_GIFT, geodistances_GIFT_dist,i,k)


###2.3.6 all barriers-----
gc()
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_barriers_grid_7.RDATA")
# Convert vector of distances/ too large dist object into a full matrix
geodistances_full= matrix(NA, nrow(grid_7), nrow(grid_7))
rownames(geodistances_full) = grid_7$grid_7_ID
colnames(geodistances_full) = grid_7$grid_7_ID

geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_barriers
geodistances_full = t(geodistances_full)
geodistances_full[lower.tri(geodistances_full, diag=FALSE)] = graphdistances_barriers
diag(geodistances_full) = 0

geodistances_full[1:5,1:5]

geodistances_GIFT = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)
colnames(geodistances_GIFT) = sort(geoentities_simple@data$Region_id)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_GIFT[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_GIFT[upper.tri(geodistances_GIFT,
                            diag=FALSE)] = t(geodistances_GIFT)[upper.tri(geodistances_GIFT,
                                                                          diag=FALSE)]
geodistances_GIFT[1:5,1:5]

geodistances_GIFT_dist = geodistances_GIFT[order(as.numeric(row.names(geodistances_GIFT))),
                                           order(as.numeric(colnames(geodistances_GIFT)))]
geodistances_GIFT_dist[1:5,1:5]
geodistances_GIFT_dist = as.dist(geodistances_GIFT_dist)

graphdistances_barriers = geodistances_GIFT_dist
save(graphdistances_barriers,
     file="results/primary_results/distances_GIFT/graphdistances_barriers_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_barriers, geodistances_GIFT, geodistances_GIFT_dist,i,k)



###2.3.7 climatic distances-----
gc()
#selected entity
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7_lm_isl_clipped.RData")
summary(grid_7_lm_clipped@data)

geoentities_simple@proj4string
grid_7_lm_clipped@proj4string

geoentities_simple@proj4string = grid_7_lm_clipped@proj4string

# Overlay geoentities and hexagon grid
overlap = over(geoentities_simple, grid_7_lm_clipped, returnList = TRUE)

load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids_isl_clipped/graphdistances_clim_grid_7.RDATA")
pastdistances_full = as.matrix(graphdistances_clim)
pastdistances_full[1:5,1:5]

pastdistances_gift = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))

rownames(pastdistances_gift) = geoentities_simple@data$Region_id
colnames(pastdistances_gift) = geoentities_simple@data$Region_id


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    pastdistances_gift[k,i] = min(pastdistances_full[as.character(overlap[[i]]$grid_7_isl_ID),
                                                     as.character(overlap[[k]]$grid_7_isl_ID)],
                                  na.rm=TRUE)
  }
}

pastdistances_gift[1:5,1:5]
pastdistances_gift[upper.tri(pastdistances_gift, diag=FALSE)] = t(pastdistances_gift)[upper.tri(pastdistances_gift, diag=FALSE)]
pastdistances_gift[1:5,1:5]

pastdistances_gift_dist = pastdistances_gift[order(as.numeric(row.names(pastdistances_gift))),
                                             order(as.numeric(colnames(pastdistances_gift)))]
pastdistances_gift_dist = as.dist(pastdistances_gift_dist)

graphdistances_clim = pastdistances_gift_dist
save(graphdistances_clim,
     file="results/primary_results/distances_GIFT/graphdistances_clim_grid_7.RDATA")



#3. trade flow----
gc()
library(countrycode)
library(stringr)
library(stringi)
library(tibble)
library(data.table)

trade_0 = read.csv("data/trade_flow/BACI_HS22_Y2022_V202501.csv", header = T,
                   sep = ',')
country_code = read.csv("data/trade_flow/country_codes_V202501.csv", header = T,
                        sep = ',')
#product_code = read.csv("data/trade_flow/product_codes_HS22_V202501.csv", header = T,
#                        sep = ',')
population = terra::rast("data/population/gpw_v4_population_count_rev11_2020_30_min.tif")

str(trade_0)
str(country_code)
str(df)

df = df %>% st_drop_geometry()

intersect(sort(unique(df$GID_0)), sort(unique(country_code$country_iso3)))
no_trade_id = setdiff(sort(unique(df$GID_0)), sort(unique(country_code$country_iso3)))

df_country_code = df %>% st_drop_geometry() %>% 
  left_join(country_code, by = join_by('GID_0' == 'country_iso3'))

df_no_country_code = df_country_code %>% filter(is.na(country_code))
df_no_country_code$NAME_0

## 3.1 Join country/region names of GIFT and trade data----

# modify some ISO_code to make some region in the TYDG match country_code 
df_2 = arrange(df, df$NAME_0)
df_2$my_GID_0 = df_2$GID_0
df_2$my_GID_0[df_2$NAME_0 == "United States Minor Outlying Islands"] = "PUS"
df_2$my_GID_0[df_2$NAME_0 == "Puerto Rico"] = "PUS"
df_2$my_GID_0[df_2$NAME_0 == "Taiwan"] = "CHN"

df_country_code = df_2 %>% left_join(country_code,
                                     by = join_by('my_GID_0' == 'country_iso3')) %>% 
  filter(!is.na(country_code))
table(sort(df_country_code$Region_id))
length(unique(df_country_code$Region_id))/length(df$Region_id)
# 92% of regions has trade data
df_no_country_code = df_2 %>% left_join(country_code,
                                        by = join_by('my_GID_0' == 'country_iso3')) %>% 
  filter(is.na(country_code))

Region_id_no_trade = df_no_country_code$Region_id
df_eck4 = st_read("data/Plants/shp_glonaf_new_eck4_1.shp")
df_no_country_code = df_eck4 %>% filter(Region_id %in% Region_id_no_trade)
plot(df_no_country_code$geometry)
#st_write(df_no_country_code, 
#         'results/primary_results/predictors_GIFT/df_no_country_code.shp')
#save(Region_id_no_trade,
#file = 'results/primary_results/predictors_GIFT/Region_id_no_trade.rdata')


## 3.2 Join country/region names of GIFT and population data----
# make projection coincidence
df = st_read("data/Plants/shp_glonaf_new.shp")
df_for_population = sf::st_transform(df, crs = sf::st_crs(population))

# convert to SpatVector and do zonal
v_for_population = terra::vect(df_for_population)
# could be clip to high speed
population_crop = terra::crop(population, v_for_population, snap="out")

# zonal: by=Region_id, fun="mean"
z_for_population = cbind(Region_id = v_for_population$Region_id,
                         terra::extract(population_crop,
                                        v_for_population, fun = "sum", na.rm = TRUE))
colnames(z_for_population)[which(colnames(z_for_population) == 'gpw_v4_population_count_rev11_2020_30_min')] = 'population_count_2020'

df_country_code_population = df_country_code %>% left_join(z_for_population[,c('Region_id',
                                                                               'population_count_2020')],
                                                           by = 'Region_id')
rm(population_crop, v_for_population)

## 3.3 Construct trade flow matrix----
#distance for tradeflow among regions
df_country_code_population = df_country_code_population %>% st_drop_geometry()
df_country_code_population$rel_popu = NA
df_country_code_population[is.na(df_country_code_population$population_count_2020),]$population_count_2020 = 0
list_country_code_population = split(df_country_code_population, df_country_code_population$country_code)
list_country_code_population = lapply(list_country_code_population,
                                      function(x){
                                        #x = list_country_code_population[['76']]
                                        if (nrow(x) == 1){x$rel_popu = 1
                                        } else {
                                          x$rel_popu = x$population_count_2020/sum(x$population_count_2020,
                                                                                   na.rm = T)}
                                        
                                        
                                        return(x)
                                      })
df_country_code_population = data.table::rbindlist(list_country_code_population)

df_country_code_population_direct = df_country_code_population %>% filter(rel_popu == 1)

Region_id_direct_trade = sort(df_country_code_population_direct$Region_id)
save(Region_id_direct_trade,
     file = 'results/primary_results/predictors_GIFT/Region_id_direct_trade.rdata')

### 3.3.1 Symmetrical matrix for regional ED----
tradeflow_mat_s = matrix(NA, nrow = length(unique(df_country_code_population$Region_id)),
                         ncol = length(unique(df_country_code_population$Region_id)))

rownames(tradeflow_mat_s) = sort(unique(df_country_code_population$Region_id))
colnames(tradeflow_mat_s) = sort(unique(df_country_code_population$Region_id))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_s)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_s)){
    #j = 1
    region_i = rownames(tradeflow_mat_s)[i]
    region_j = rownames(tradeflow_mat_s)[j]
    
    country_code_i = df_country_code_population %>% filter(Region_id == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(Region_id == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(Region_id == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(Region_id == region_j)%>% 
      pull(rel_popu) %>% unique()
    
    flow = (sum(tr[J(country_code_i, country_code_j), v], 
                na.rm = T) + sum(tr[J(country_code_j, country_code_i), v], 
                                 na.rm = T)) * rel_popu_i * rel_popu_j
    tradeflow_mat_s[j,i] = flow
  }
}


tradeflow_mat_s[upper.tri(tradeflow_mat_s,
                          diag=FALSE)] = t(tradeflow_mat_s)[upper.tri(tradeflow_mat_s,
                                                                      diag=FALSE)]
tradeflow_mat_s[1:5,1:5]

tradeflow_mat_s = tradeflow_mat_s[order(as.numeric(row.names(tradeflow_mat_s))),
                                  order(as.numeric(colnames(tradeflow_mat_s)))]
tradeflow_mat_s[1:10,1:10]
colSums(tradeflow_mat_s)

tradeflow_mat_s['217','147'] # trade flow value between shanghai and France
# 1163109 * 1000 USD
tradeflow_mat_s[217,147]

save(tradeflow_mat_s,
     file="results/primary_results/distances_GIFT/tradeflow_mat_s.RDATA")



### 3.3.1 Asymmetrical matrix for sp ED----
tradeflow_mat_a = matrix(NA, nrow = length(unique(df_country_code_population$Region_id)),
                         ncol = length(unique(df_country_code_population$Region_id)))

rownames(tradeflow_mat_a) = sort(unique(df_country_code_population$Region_id))
colnames(tradeflow_mat_a) = sort(unique(df_country_code_population$Region_id))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_a)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_a)){
    #j = 1
    region_i = rownames(tradeflow_mat_a)[i]
    region_j = rownames(tradeflow_mat_a)[j]
    
    country_code_i = df_country_code_population %>% filter(Region_id == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(Region_id == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(Region_id == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(Region_id == region_j)%>% 
      pull(rel_popu) %>% unique()
    
    flow_ji = sum(tr[J(country_code_i, country_code_j), v], 
                  na.rm = T)  * rel_popu_j
    tradeflow_mat_a[j,i] = flow_ji
    
    flow_ij = sum(tr[J(country_code_j, country_code_i), v], 
                  na.rm = T)  * rel_popu_j
    tradeflow_mat_a[i,j] = flow_ij
  }
}


tradeflow_mat_a[1:5,1:5]

tradeflow_mat_a = tradeflow_mat_a[order(as.numeric(row.names(tradeflow_mat_a))),
                                  order(as.numeric(colnames(tradeflow_mat_a)))]
tradeflow_mat_a[1:10,1:10]
colSums(tradeflow_mat_a)

save(tradeflow_mat_a,
     file="results/primary_results/distances_GIFT/tradeflow_mat_a.RDATA")







#4. characteristics of region's phylogenetic structure----
#library(picante)
library(PhyloMeasures)

##4.1 plant's exotics PD----
load("D:/R projects/Global_ED/data/Plants/data/phylo.fake.species.653.Rdata")
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")

colnames(df.native.650) == colnames(df.natu.650)

mat_plant_exotic = df.natu.650 %>%
  as.data.frame() %>% 
  distinct(Region_id, species) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(
    names_from  = species,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(Region_id) %>% as.matrix()

rownames(mat_plant_exotic) = mat_plant_exotic[,1]
mat_plant_exotic = mat_plant_exotic[,-1]

setdiff(colnames(mat_plant_exotic), phylo$tip.label)

pd_plant_exotic = data.frame(Region_id = rownames(mat_plant_exotic),
                              pd = pd.query(phylo, mat_plant_exotic,
                                            abundance.weights = F))

save(pd_plant_exotic,
     file="results/primary_results/predictors_GIFT/pd_plant_exotic.RDATA")
rm(pd_plant_exotic, mat_plant_exotic, df.natu.650)

##4.2 plant's natives PD----
gc()
mat_plant_native = df.native.650 %>%
  as.data.frame() %>% 
  distinct(Region_id, species) %>% 
  mutate(presence = 1) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(
    names_from  = species,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(Region_id) %>% as.matrix()

rownames(mat_plant_native) = mat_plant_native[,1]
mat_plant_native = mat_plant_native[,-1]

pd_plant_native = data.frame(Region_id = rownames(mat_plant_native),
                             pd = pd.query(phylo, mat_plant_native,
                                           abundance.weights = F))
save(pd_plant_native,
     file="results/primary_results/predictors_GIFT/pd_plant_native.RDATA")


