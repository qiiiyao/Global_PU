### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#0 loading required packages-------------------
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("foreach", "dplyr", 'tidyr', 'sf', 'terra', 'raster', 'scico')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("D:/R projects/Global_ED")

# load the background maps for quantifying the predictors
df = st_read("data/Fishes/data/Basin042017_3119.shp")
df = sf::st_make_valid(df)  

#1. Environemt variables----
##1.1 elevation and soil cation exchange ability----
# read the elevation dataset
elevation = terra::rast("data/envir_vars/gmted2010/GMTED2010_15n015_00625deg.nc")
# accessed at https://temis.nl/data/gmted2010/, unit = m
elevation$Elev_range = elevation$elevation_max - elevation$elevation_min

# make projection coincidence
df_for_eleva = sf::st_transform(df, crs = sf::st_crs(elevation))

# convert to SpatVector and do zonal
v_for_eleva = terra::vect(df_for_eleva)
# could be clip to high speed
elevation_crop = terra::crop(elevation, v_for_eleva, snap="out")

# zonal: by=BasinName, fun="mean"
z_for_eleva = cbind(BasinName = v_for_eleva$BasinName,
                    terra::extract(elevation_crop,
                             v_for_eleva, fun = "mean", na.rm = TRUE))
rm(elevation_crop, v_for_eleva)


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

# zonal: by=BasinName, fun="mean"
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
  left_join(z_for_eleva, by = 'BasinName') %>% 
  st_drop_geometry()

geoentities_env = cbind(geoentities_env,
                        MAT = z_for_bioclim1$`CHELSA_bio1_1981-2010_V.2.1`,
                        TS = z_for_bioclim4$`CHELSA_bio4_1981-2010_V.2.1`,
                        MAP = z_for_bioclim12$`CHELSA_bio12_1981-2010_V.2.1`,
                        PS = z_for_bioclim15$`CHELSA_bio15_1981-2010_V.2.1`,
                        AI = z_for_arid_index$`CHELSA_ai_1981-2010_V.2.1`)

summary(geoentities_env)
save(geoentities_env,
     file = 'results/primary_results/predictors_basins/geoentities_env.rdata')


#### calculate climate distance among different regions
load('results/primary_results/predictors_basins/geoentities_env.rdata')
clim_dist = as.matrix(dist(geoentities_env %>% dplyr::select(MAT, TS,
                                                             MAP, PS,
                                                             AI),
                           method = 'euclidean', diag = T))
rownames(clim_dist) = geoentities_env$BasinName
colnames(clim_dist) = geoentities_env$BasinName
save(clim_dist,
     file = 'results/primary_results/distances_basins/clim_dist.rdata')





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
geoentities_simple = readOGR(dsn = "data/Fishes/data/Basin042017_3119.shp")
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

geodistances_basins = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_basins) = sort(geoentities_simple@data$BasinName)
colnames(geodistances_basins) = sort(geoentities_simple@data$BasinName)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_basins[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_basins[upper.tri(geodistances_basins,
                            diag=FALSE)] = t(geodistances_basins)[upper.tri(geodistances_basins,
                                                                          diag=FALSE)]
geodistances_basins[1:5,1:5]

geodistances_basins_dist = geodistances_basins[order(as.numeric(row.names(geodistances_basins))),
                                           order(as.numeric(colnames(geodistances_basins)))]
geodistances_basins_dist[1:5,1:5]
geodistances_basins_dist = as.dist(geodistances_basins_dist)

geodistances = geodistances_basins_dist
save(geodistances,
     file="results/primary_results/distances_basins/geodistances_grid_7.RDATA")

#plot(geodistances_basins_dist, geodistances_points)
rm(geodistances_full,geodistances_basins,geodistances_basins_dist,i,k)



###2.3.2 all barriers-----
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

geodistances_basins = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_basins) = sort(geoentities_simple@data$BasinName)
colnames(geodistances_basins) = sort(geoentities_simple@data$BasinName)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_basins[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                     as.character(overlap[[k]]$grid_7_ID)],
                                   na.rm=TRUE)
  }
}

geodistances_basins[upper.tri(geodistances_basins,
                              diag=FALSE)] = t(geodistances_basins)[upper.tri(geodistances_basins,
                                                                              diag=FALSE)]
geodistances_basins[1:5,1:5]

geodistances_basins_dist = geodistances_basins[order(as.numeric(row.names(geodistances_basins))),
                                               order(as.numeric(colnames(geodistances_basins)))]
geodistances_basins_dist[1:5,1:5]
geodistances_basins_dist = as.dist(geodistances_basins_dist)

graphdistances_barriers = geodistances_basins_dist
save(graphdistances_barriers,
     file="results/primary_results/distances_basins/graphdistances_barriers_grid_7.RDATA")

#plot(geodistances_basins_dist, geodistances_points)
rm(geodistances_full, graphdistances_barriers, geodistances_basins, geodistances_basins_dist,i,k)


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

geodistances_basins = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_basins) = sort(geoentities_simple@data$BasinName)
colnames(geodistances_basins) = sort(geoentities_simple@data$BasinName)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_basins[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_basins[upper.tri(geodistances_basins,
                            diag=FALSE)] = t(geodistances_basins)[upper.tri(geodistances_basins,
                                                                          diag=FALSE)]
geodistances_basins[1:5,1:5]

geodistances_basins_dist = geodistances_basins[order(as.numeric(row.names(geodistances_basins))),
                                           order(as.numeric(colnames(geodistances_basins)))]
geodistances_basins_dist[1:5,1:5]
geodistances_basins_dist = as.dist(geodistances_basins_dist)

graphdistances_elev = geodistances_basins_dist
save(graphdistances_elev,
     file="results/primary_results/distances_basins/graphdistances_elev_grid_7.RDATA")

#plot(geodistances_basins_dist, geodistances_points)
rm(geodistances_full, graphdistances_elev, geodistances_basins, geodistances_basins_dist,i,k)


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

geodistances_basins = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_basins) = sort(geoentities_simple@data$BasinName)
colnames(geodistances_basins) = sort(geoentities_simple@data$BasinName)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_basins[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_basins[upper.tri(geodistances_basins,
                            diag=FALSE)] = t(geodistances_basins)[upper.tri(geodistances_basins,
                                                                          diag=FALSE)]
geodistances_basins[1:5,1:5]

geodistances_basins_dist = geodistances_basins[order(as.numeric(row.names(geodistances_basins))),
                                           order(as.numeric(colnames(geodistances_basins)))]
geodistances_basins_dist[1:5,1:5]
geodistances_basins_dist = as.dist(geodistances_basins_dist)

graphdistances_temp = geodistances_basins_dist
save(graphdistances_temp,
     file="results/primary_results/distances_basins/graphdistances_temp_grid_7.RDATA")

#plot(geodistances_basins_dist, geodistances_points)
rm(geodistances_full, graphdistances_temp, geodistances_basins, geodistances_basins_dist,i,k)

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

geodistances_basins = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_basins) = sort(geoentities_simple@data$BasinName)
colnames(geodistances_basins) = sort(geoentities_simple@data$BasinName)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_basins[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_basins[upper.tri(geodistances_basins,
                            diag=FALSE)] = t(geodistances_basins)[upper.tri(geodistances_basins,
                                                                          diag=FALSE)]
geodistances_basins[1:5,1:5]

geodistances_basins_dist = geodistances_basins[order(as.numeric(row.names(geodistances_basins))),
                                           order(as.numeric(colnames(geodistances_basins)))]
geodistances_basins_dist[1:5,1:5]
geodistances_basins_dist = as.dist(geodistances_basins_dist)

graphdistances_ai = geodistances_basins_dist
save(graphdistances_ai,
     file="results/primary_results/distances_basins/graphdistances_ai_grid_7.RDATA")

#plot(geodistances_basins_dist, geodistances_points)
rm(geodistances_full, graphdistances_ai, geodistances_basins, geodistances_basins_dist,i,k)


###2.3.6 climatic distances-----
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

pastdistances_basins = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))

rownames(pastdistances_basins) = geoentities_simple@data$BasinName
colnames(pastdistances_basins) = geoentities_simple@data$BasinName


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    pastdistances_basins[k,i] = min(pastdistances_full[as.character(overlap[[i]]$grid_7_isl_ID),
                                                     as.character(overlap[[k]]$grid_7_isl_ID)],
                                  na.rm=TRUE)
  }
}

pastdistances_basins[1:5,1:5]
pastdistances_basins[upper.tri(pastdistances_basins, diag=FALSE)] = t(pastdistances_basins)[upper.tri(pastdistances_basins, diag=FALSE)]
pastdistances_basins[1:5,1:5]

pastdistances_basins_dist = pastdistances_basins[order(as.character(row.names(pastdistances_basins))),
                                             order(as.character(colnames(pastdistances_basins)))]
pastdistances_basins_dist = as.dist(pastdistances_basins_dist)

graphdistances_clim = pastdistances_basins_dist
save(graphdistances_clim,
     file="results/primary_results/distances_basins/graphdistances_clim_grid_7.RDATA")





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

df = df %>%
  mutate(iso_code = countrycode(Country, origin = "country.name", destination = "iso3c"))

intersect(sort(unique(df$iso_code)), sort(unique(country_code$country_iso3)))
no_trade_id = setdiff(sort(unique(df$iso_code)), sort(unique(country_code$country_iso3)))


df_country_code = df %>% 
  left_join(country_code, by = join_by('iso_code' == 'country_iso3'))

df_no_country_code = df_country_code %>% filter(is.na(country_code))
unique(df_no_country_code$iso_code)
unique(df_no_country_code$Country)

## 3.1 Join country/region names of basins and trade data----

# modify some ISO_code to make some region in the TYDG match country_code 
df_2 = arrange(df, df$iso_code)
df_2$my_iso_code = df_2$iso_code
df_2$my_iso_code[df_2$Country == "Puerto Rico"] = "PUS"
df_2$my_iso_code[df_2$Country == "Taiwan"] = "CHN"

df_country_code = df_2 %>% left_join(country_code,
                                     by = join_by('my_iso_code' == 'country_iso3')) %>% 
  filter(!is.na(country_code))
table(sort(df_country_code$Region_id))
length(unique(df_country_code$BasinName))/length(df$BasinName)
# 0.991664 of basins has trade data
df_no_country_code = df_2 %>% left_join(country_code,
                                        by = join_by('my_iso_code' == 'country_iso3')) %>% 
  filter(is.na(country_code))
unique(df_no_country_code$Country)

basins_no_trade = sort(unique(df_no_country_code$BasinName))
df_eck4 = st_read("data/Fishes/data/Basin042017_3119_eck4.shp")
df_no_country_code = df_eck4 %>% filter(BasinName %in% basins_no_trade)
plot(df_no_country_code$geometry)
#st_write(df_no_country_code, 
#         'results/primary_results/predictors_basins/df_no_country_code.shp')
save(basins_no_trade,
file = 'results/primary_results/predictors_basins/basins_no_trade.rdata')


## 3.2 Join country/region names of basins and population data----
# make projection coincidence
df = st_read("data/Fishes/data/Basin042017_3119.shp")
df_for_population = sf::st_transform(df, crs = sf::st_crs(population))

# convert to SpatVector and do zonal
v_for_population = terra::vect(df_for_population)
# could be clip to high speed
population_crop = terra::crop(population, v_for_population, snap="out")

# zonal: by=Region_id, fun="mean"
z_for_population = cbind(BasinName = v_for_population$BasinName,
                         terra::extract(population_crop,
                                        v_for_population, fun = "sum", na.rm = TRUE))
colnames(z_for_population)[which(colnames(z_for_population) == 'gpw_v4_population_count_rev11_2020_30_min')] = 'population_count_2020'

df_country_code_population = df_country_code %>% left_join(z_for_population[,c('BasinName',
                                                                               'population_count_2020')],
                                                           by = 'BasinName')
rm(population_crop, v_for_population)

## 3.3 Construct trade flow matrix----
#distance for tradeflow among regions
df_country_code_population = df_country_code_population %>% st_drop_geometry()
df_country_code_population$rel_popu = NA

if(sum(is.na(df_country_code_population$population_count_2020)) > 0){
  df_country_code_population[is.na(df_country_code_population$population_count_2020),]$population_count_2020 = 0
}

list_country_code_population = split(df_country_code_population,
                                     df_country_code_population$country_code)
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

basin_direct_trade = sort(df_country_code_population_direct$BasinName)
save(basin_direct_trade,
     file = 'results/primary_results/predictors_basins/basin_direct_trade.rdata')

### 3.3.1 Symmetrical matrix for regional ED----
tradeflow_mat_s = matrix(NA, nrow = length(unique(df_country_code_population$BasinName)),
                         ncol = length(unique(df_country_code_population$BasinName)))

rownames(tradeflow_mat_s) = sort(unique(df_country_code_population$BasinName))
colnames(tradeflow_mat_s) = sort(unique(df_country_code_population$BasinName))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_s)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_s)){
    #j = 1
    region_i = rownames(tradeflow_mat_s)[i]
    region_j = rownames(tradeflow_mat_s)[j]
    
    country_code_i = df_country_code_population %>% filter(BasinName == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(BasinName == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
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
     file="results/primary_results/distances_basins/tradeflow_mat_s.RDATA")



### 3.3.1 Asymmetrical matrix for sp ED----
tradeflow_mat_a = matrix(NA, nrow = length(unique(df_country_code_population$BasinName)),
                         ncol = length(unique(df_country_code_population$BasinName)))

rownames(tradeflow_mat_a) = sort(unique(df_country_code_population$BasinName))
colnames(tradeflow_mat_a) = sort(unique(df_country_code_population$BasinName))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_a)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_a)){
    #j = 1
    region_i = rownames(tradeflow_mat_a)[i]
    region_j = rownames(tradeflow_mat_a)[j]
    
    country_code_i = df_country_code_population %>% filter(BasinName == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(BasinName == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
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
     file="results/primary_results/distances_basins/tradeflow_mat_a.RDATA")







#4. characteristics of basin's phylogenetic structure----
#library(picante)
library(PhyloMeasures)

##4.1 fishes' exotics PD----
load("D:/R projects/Global_ED/data/Fishes/data/my_phy.rdata")
load("D:/R projects/Global_ED/data/Fishes/data/my_data_used_final.rdata")

data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
str(data.used_final_exotics)
mat_fish_exotic = data.used_final_exotics %>%
  as.data.frame() %>% 
  distinct(X1.Basin.Name, valid_names) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(
    names_from  = valid_names,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(X1.Basin.Name) %>% as.matrix()

rownames(mat_fish_exotic) = mat_fish_exotic[,1]
mat_fish_exotic = mat_fish_exotic[,-1]

setdiff(colnames(mat_fish_exotic), phylo$tip.label)

mode(mat_fish_exotic)= "numeric"
str(mat_fish_exotic)

pd_fish_exotic = data.frame(X1.Basin.Name = rownames(mat_fish_exotic),
                            picante::pd(mat_fish_exotic, phylo,
                                 include.root = T))

save(pd_fish_exotic,
     file="results/primary_results/predictors_basins/pd_fish_exotic.RDATA")
rm(pd_fish_exotic, mat_fish_exotic, data.used_final_exotics)

##4.2 fish's natives PD----
gc()
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')

mat_fish_native = data.used_final_natives %>%
  as.data.frame() %>% 
  distinct(X1.Basin.Name, valid_names) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(
    names_from  = valid_names,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(X1.Basin.Name) %>% as.matrix()

rownames(mat_fish_native) = mat_fish_native[,1]
mat_fish_native = mat_fish_native[,-1]
setdiff(colnames(mat_fish_native), phylo$tip.label)
mode(mat_fish_native)= "numeric"
is.rooted(phylo) 

is.rooted(phylo)
is.binary(phylo)
phylo_fixed = multi2di(phylo) 

pd_fish_native = data.frame(X1.Basin.Name = rownames(mat_fish_native),
                            picante::pd(mat_fish_native, phylo_fixed,
                                        include.root = T))
save(pd_fish_native,
     file="results/primary_results/predictors_basins/pd_fish_native.RDATA")


