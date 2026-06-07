### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#0 loading required packages-------------------
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("foreach", "dplyr", 'tidyr', 'sf', 'terra', 'raster', 'scico',
                 'GIFT', 'ape')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("D:/R projects/Global_ED")

# load the background maps for quantifying the predictors
df = st_read("data/TDWG4/TDWG4_newTibet.shp")
df = st_make_valid(df)

sort(unique(df$Level_4_Na))


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

# zonal: by=regionID, fun="mean"
z_for_eleva = cbind(RegionID = v_for_eleva$RegionID,
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
z_for_ces = read.table('results/primary_results/predictors_TDWG/geo_ces.txt',
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

# zonal: by=regionID, fun="mean"
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
  left_join(z_for_eleva, by = 'RegionID') %>% 
  left_join(z_for_ces[,c('RegionID', 'CECSOL')], by = 'RegionID') %>% 
  st_drop_geometry()

geoentities_env = cbind(geoentities_env,
                        MAT = z_for_bioclim1$`CHELSA_bio1_1981-2010_V.2.1`,
                        TS = z_for_bioclim4$`CHELSA_bio4_1981-2010_V.2.1`,
                        MAP = z_for_bioclim12$`CHELSA_bio12_1981-2010_V.2.1`,
                        PS = z_for_bioclim15$`CHELSA_bio15_1981-2010_V.2.1`,
                        AI = z_for_arid_index$`CHELSA_ai_1981-2010_V.2.1`)

save(geoentities_env,
     file = 'results/primary_results/predictors_TDWG/geoentities_env.rdata')

#### calculate climate distance among different regions
load("results/primary_results/predictors_TDWG/geoentities_env.rdata")
clim_dist = as.matrix(dist(geoentities_env %>% dplyr::select(MAT, TS, MAP, PS, AI),
                           method = 'euclidean', diag = T))
rownames(clim_dist) = geoentities_env$RegionID
colnames(clim_dist) = geoentities_env$RegionID
save(clim_dist,
     file = 'results/primary_results/distances_TDWG/clim_dist.rdata')


#2. dispersal limitation ----
###(least cost distance considering physical, ecological dispersal barriers, and climate similarity)

##2.1 Prepare hexagon point shapefile and calculate point to point distances-----
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7.RData")

# calculate centroids
# point x point y
masscentroids = function(tmp_entity){
  tmp_entity_x = suppressMessages(fortify(tmp_entity))[,1] # extract all x coordinates for tmp_entity
  if (diff(range(tmp_entity_x)) > 350 & min(abs(tmp_entity_x)) > 10){  # in case entity crosses date line (only if entity does not also cross 0 degrees long)
    tmp_entity_mod = tmp_entity
    for (k in 1:length(tmp_entity_mod@polygons[[1]]@Polygons)){
      tmp_entity_mod@polygons[[1]]@Polygons[[k]]@coords[which(tmp_entity_mod@polygons[[1]]@Polygons[[k]]@coords[,1]<0),1] = tmp_entity_mod@polygons[[1]]@Polygons[[k]]@coords[which(tmp_entity_mod@polygons[[1]]@Polygons[[k]]@coords[,1]<0),1] + 360
    }
    centroid = gCentroid(tmp_entity_mod)
  } else {
    centroid = gCentroid(tmp_entity)
  }
  
  point_x = round(centroid@coords[1], 6)
  if (point_x[1] > 180){
    point_x[1] = point_x[1] - 360
  }
  point_y = round(centroid@coords[2], 6)
  return(cbind(point_x,point_y))
}

grid_7@data$point_x = NA
grid_7@data$point_y = NA

for (i in 1:nrow(grid_7)){
  grid_7@data[i,c("point_x","point_y")] = masscentroids(grid_7[i,])
}

geodistances = sapply(1:(nrow(grid_7)-1),
                      function(x) distGeo(grid_7@data[x,c("point_x","point_y")],
                                          grid_7@data[c((x+1):nrow(grid_7)),
                                                      c("point_x","point_y")])/1000)
geodistances = unlist(geodistances)
save(geodistances,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids/geodistances_grid_7.RDATA")


##2.2 Establish graph and calculate cost distances -----
###2.2.1 Calculate Barriers----
rm(list = ls())
# Add env data for calculation of barriers
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7.RData")
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/geodistances_grid_7.RDATA")

load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/env_grid_7_lm_clipped_v3.RData")
env_grid =env_grid %>%dplyr::select(grid_7_ID,mean_mn30_grd,
                                    `mean_CHELSA_bio1_1981-2010_V.2.1`,
                                    mean_ai_v3_yr)
names(env_grid)[3:4]=c("CHELSA_bio1_1981-2010_V.2.1","ai_v3_yr")

grid_7@data = left_join(grid_7@data,env_grid,by="grid_7_ID")
rm(env_grid)

# prepare measure for barriers
water = grid_7@data$mean_mn30_grd #Mean Altitude
water[!is.na(grid_7@data$mean_mn30_grd)] = 0
water[is.na(grid_7@data$mean_mn30_grd)] = 1

elev = grid_7@data$mean_mn30_grd
elev = ((elev-(min(elev, na.rm = TRUE)))/(max(elev, na.rm = TRUE)-min(elev, na.rm = TRUE)))
elev[is.na(grid_7@data$mean_mn30_grd)] = 1
summary(elev)

temp = (grid_7@data$`CHELSA_bio1_1981-2010_V.2.1`-min(grid_7@data$`CHELSA_bio1_1981-2010_V.2.1`,
                                                      na.rm=TRUE))^(1/3)#annual mean temperature
temp = (1-(temp-(min(temp, na.rm = TRUE)))/(max(temp, na.rm = TRUE)-min(temp, na.rm = TRUE)))
temp[is.na(grid_7@data$mean_mn30_grd)] = 1
all(!is.na(temp))#TRUE


ai = grid_7@data$ai_v3_yr #Global Aridity Index
ai[which(ai>1)] = 1
ai = log(ai+0.001)
ai = (1-(ai-(min(ai, na.rm = TRUE)))/(max(ai, na.rm = TRUE)-min(ai, na.rm = TRUE)))
ai[is.na(grid_7@data$mean_mn30_grd)] = 1


barriers = apply(cbind(temp, ai, elev),1, function(x) x[which.max(x)])
barriers[is.na(grid_7@data$mean_mn30_grd)] = 1

###2.2.2 Unfold distance matrix and calculate weights----
neigh_dist = data.frame(grid = paste("grid_",c(5:8), sep=""),
                        dist = c(600, 340, 210, 120)) # max distance between two neighboring gridcells 5: 600, 6: 340, 7: 210, 8: 120

system.time({
  edges = data.frame(
    from = rep(grid_7@data$grid_7_ID[1:(nrow(grid_7)-1)],
               seq(nrow(grid_7)-1,1))[geodistances < neigh_dist$dist[which(neigh_dist$grid == "grid_7")]],
    to = unlist(sapply(2:nrow(grid_7),
                       function(x) grid_7@data$grid_7_ID[x:nrow(grid_7)]))[geodistances < neigh_dist$dist[which(neigh_dist$grid == "grid_7")]],
    distance = geodistances[geodistances < neigh_dist$dist[which(neigh_dist$grid == "grid_7")]]
  )
})

range(table(c(edges$from,edges$to))) # should be 5,6
length(which(table(c(edges$from,edges$to))==5)) # should be 12

rm(grid_7, layers)


load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7_lm_isl_clipped.RData")
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7.RData")

overlap = over(grid_7_lm_clipped, grid_7, returnList = TRUE)
overlap[[2]]


for (i in 1:length(overlap)){ # For those cases where The hexagons are basically the same in grid_7 and grid_7_lm_isl_clipped we can take the one hexagon grid cell ID that matches. For smaller islands we will consider all overlapping hexagons
  if(is.na(grid_7_lm_clipped$ulm_ID[i])){
    overlap[[i]] = overlap[[i]][which(overlap[[i]]$grid_7_ID == grid_7_lm_clipped$grid_7_ID[i]),] 
  }
}

max(sapply(overlap, nrow))



#from now on separately for each cost layer:
###2.2.3 waters----
grid_7@data$value = water
edges$costs = sapply(1:nrow(edges),
                     function(x) mean(c(grid_7@data$value[which(grid_7@data$grid_7_ID==edges$to[x])],
                                        grid_7@data$value[which(grid_7@data$grid_7_ID==edges$from[x])])))

range(edges$costs)
par(mfrow=c(1,1))
hist(edges$costs)

edges$weight = edges$distance * edges$costs
nodes = grid_7@data

net = graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

graphdistances = distances(net)
graphdistances[1:5,1:5]

#distance for grid_isl
graphdistances_isl = matrix(NA, nrow = nrow(grid_7_lm_clipped), ncol = nrow(grid_7_lm_clipped))

rownames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID
colnames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID

for (i in 1:nrow(grid_7_lm_clipped)){
  if((i)%%100==0) print(paste(i, "of" , nrow(grid_7_lm_clipped)))
  for (k in i:nrow(grid_7_lm_clipped)){
    graphdistances_isl[k,i] = min(graphdistances[overlap[[i]]$grid_7_ID,overlap[[k]]$grid_7_ID] ,
                                  na.rm=TRUE)
  }
}

graphdistances_isl[1:20,1:20]


graphdistances_water = as.dist(graphdistances_isl)
save(graphdistances_water,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids_isl_clipped/graphdistances_water_grid_7.RDATA")


graphdistances_water = as.dist(graphdistances)
save(graphdistances_water,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_water_grid_7.RDATA")

rm(nodes,net,graphdistances,graphdistances_isl,i,k)


###2.2.4 all barriers----
gc()
grid_7@data$value = barriers
edges$costs = sapply(1:nrow(edges),
                     function(x) mean(c(grid_7@data$value[which(grid_7@data$grid_7_ID==edges$to[x])],
                                        grid_7@data$value[which(grid_7@data$grid_7_ID==edges$from[x])])))

range(edges$costs) 
par(mfrow=c(1,1))
hist(edges$costs)

edges$weight = edges$distance * edges$costs
nodes = grid_7@data

net = graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

graphdistances = distances(net)

graphdistances_barriers = as.dist(graphdistances)

save(graphdistances_barriers,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_barriers_grid_7.RDATA")


#distance for grid_isl
graphdistances_isl = matrix(NA, nrow = nrow(grid_7_lm_clipped), ncol = nrow(grid_7_lm_clipped))
rownames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID
colnames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID


for (i in 1:nrow(grid_7_lm_clipped)){
  if((i)%%100==0) print(paste(i, "of" , nrow(grid_7_lm_clipped)))
  for (k in i:nrow(grid_7_lm_clipped)){
    graphdistances_isl[k,i] = min(graphdistances[overlap[[i]]$grid_7_ID,overlap[[k]]$grid_7_ID],na.rm=TRUE)
  }
}

graphdistances_isl[1:5,1:5]

graphdistances_barriers = as.dist(graphdistances_isl)
save(graphdistances_barriers,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids_isl_clipped/graphdistances_barriers_grid_7.RDATA")
rm(nodes,net,graphdistances,graphdistances_isl,i,k)


###2.2.5 elev-----
gc()
grid_7@data$value = elev
edges$costs = sapply(1:nrow(edges),
                     function(x) mean(c(grid_7@data$value[which(grid_7@data$grid_7_ID==edges$to[x])],
                                        grid_7@data$value[which(grid_7@data$grid_7_ID==edges$from[x])])))

summary(edges$costs)
hist(edges$costs)

edges$weight = edges$distance * edges$costs
nodes = grid_7@data

net = graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

graphdistances = distances(net)
graphdistances_elev = as.dist(graphdistances)

save(graphdistances_elev,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_elev_grid_7.RDATA")

#distance for grid_isl
graphdistances_isl = matrix(NA, nrow = nrow(grid_7_lm_clipped), ncol = nrow(grid_7_lm_clipped))
rownames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID
colnames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID

for (i in 1:nrow(grid_7_lm_clipped)){
  if((i)%%100==0) print(paste(i, "of" , nrow(grid_7_lm_clipped)))
  for (k in i:nrow(grid_7_lm_clipped)){
    graphdistances_isl[k,i] = min(graphdistances[overlap[[i]]$grid_7_ID,overlap[[k]]$grid_7_ID],na.rm=TRUE)
  }
}


graphdistances_elev = as.dist(graphdistances_isl)
save(graphdistances_elev, file="distances_grids_isl_clipped/graphdistances_elev_grid_7.RDATA")
rm(nodes,net,graphdistances,graphdistances_isl,i,k)

###2.2.6 temp-----
gc()
grid_7@data$value = temp
head(grid_7@data$value)
edges$costs = sapply(1:nrow(edges),
                     function(x) mean(c(grid_7@data$value[which(grid_7@data$grid_7_ID==edges$to[x])],
                                        grid_7@data$value[which(grid_7@data$grid_7_ID==edges$from[x])])))

summary(edges$costs)
hist(edges$costs)

edges$weight = edges$distance * edges$costs
nodes = grid_7@data

net = graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

graphdistances = distances(net)
graphdistances_temp = as.dist(graphdistances)

save(graphdistances_temp,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_temp_grid_7.RDATA")

#distance for grid_isl
graphdistances_isl = matrix(NA, nrow = nrow(grid_7_lm_clipped), ncol = nrow(grid_7_lm_clipped))
rownames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID
colnames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID

for (i in 1:nrow(grid_7_lm_clipped)){
  if((i)%%100==0) print(paste(i, "of" , nrow(grid_7_lm_clipped)))
  for (k in i:nrow(grid_7_lm_clipped)){
    graphdistances_isl[k,i] = min(graphdistances[overlap[[i]]$grid_7_ID,overlap[[k]]$grid_7_ID],na.rm=TRUE)
  }
}


graphdistances_temp = as.dist(graphdistances_isl)
save(graphdistances_temp, file="distances_grids_isl_clipped/graphdistances_temp_grid_7.RDATA")
rm(nodes,net,graphdistances,graphdistances_isl,i,k)


###2.2.7 Aridity Index-----
gc()
grid_7@data$value = ai
head(grid_7@data$value)
edges$costs = sapply(1:nrow(edges),
                     function(x) mean(c(grid_7@data$value[which(grid_7@data$grid_7_ID==edges$to[x])],
                                        grid_7@data$value[which(grid_7@data$grid_7_ID==edges$from[x])])))

summary(edges$costs)
hist(edges$costs)

edges$weight = edges$distance * edges$costs
nodes = grid_7@data

net = graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

graphdistances = distances(net)
graphdistances_ai = as.dist(graphdistances)
save(graphdistances_ai,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids/graphdistances_ai_grid_7.RDATA")

#distance for grid_isl
graphdistances_isl = matrix(NA, nrow = nrow(grid_7_lm_clipped), ncol = nrow(grid_7_lm_clipped))
rownames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID
colnames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID

for (i in 1:nrow(grid_7_lm_clipped)){
  if((i)%%100==0) print(paste(i, "of" , nrow(grid_7_lm_clipped)))
  for (k in i:nrow(grid_7_lm_clipped)){
    graphdistances_isl[k,i] = min(graphdistances[overlap[[i]]$grid_7_ID,overlap[[k]]$grid_7_ID],na.rm=TRUE)
  }
}


graphdistances_ai = as.dist(graphdistances_isl)
save(graphdistances_ai, file="distances_grids_isl_clipped/graphdistances_ai_grid_7.RDATA")
rm(nodes,net,graphdistances,graphdistances_isl,i,k)

###2.2.8 Climate distance-----
gc()
rm(list = ls())
# PCA: use 19 climatic variables from Chelsa (https://chelsa-climate.org/)
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/env_grid_7_lm_clipped_v3.RData")
env_grid =env_grid %>%dplyr::select(grid_7_ID,contains("mean_CHELSA_bio"))
rownames(env_grid) = env_grid$grid_7_ID


nrow(env_grid[complete.cases(env_grid),])#no missing data
Hmisc::hist.data.frame(env_grid[,c(2:20)])
# transform the data to have normally distributed variables
variables_log = c("mean_CHELSA_bio12_1981-2010_V.2.1","mean_CHELSA_bio13_1981-2010_V.2.1",
                  "mean_CHELSA_bio14_1981-2010_V.2.1", "mean_CHELSA_bio15_1981-2010_V.2.1",
                  "mean_CHELSA_bio16_1981-2010_V.2.1", "mean_CHELSA_bio17_1981-2010_V.2.1",
                  "mean_CHELSA_bio18_1981-2010_V.2.1","mean_CHELSA_bio19_1981-2010_V.2.1")

for (i in 1:length(variables_log)){
  env_grid[,variables_log[i]] = log10(env_grid[,variables_log[i]]-min(env_grid[,variables_log[i]],
                                                                      na.rm=TRUE)+1)
}


# Run PCA on dataset
pca_imputed=PCA(env_grid, scale.unit=TRUE, graph=FALSE)
plot(pca_imputed,choix="var")

# Extract axis scores and eigenvalues
env_grid = cbind(env_grid, pca_imputed$ind$coord)

# make cost-surface based on climate dissimilarity of one cell to all others
# Euclidean distance
climdistances = dist(env_grid[,c("Dim.1",  "Dim.2",  "Dim.3",  "Dim.4", "Dim.5")], method = "euclidean")
summary(climdistances)
climdistances = climdistances/max(climdistances)
summary(climdistances)

climdistances = as.matrix(climdistances)
rownames(climdistances) = env_grid$grid_7_ID
colnames(climdistances) = env_grid$grid_7_ID

climdistances = as.data.frame(climdistances)
climdistances$grid_7_ID  = env_grid$grid_7_ID


#load
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7.RData")
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7_lm_isl_clipped.RData")

#overlap 
overlap = over(grid_7_lm_clipped, grid_7, returnList = TRUE)
overlap[[2]]

for (i in 1:length(overlap)){ 
  if(is.na(grid_7_lm_clipped$ulm_ID[i])){
    overlap[[i]] = overlap[[i]][which(overlap[[i]]$grid_7_ID == grid_7_lm_clipped$grid_7_ID[i]),] 
  }
}

max(sapply(overlap, nrow))

#
grid_7@data$value = NA
grid_7@data = left_join(grid_7@data,climdistances, by="grid_7_ID") 

grid_7@data[1:10,1:10]
grid_7@data[is.na(grid_7@data)] = 1
grid_7@data[1:10,1:10]



# Now individual cost distances over individual cost layers for all cells
# unfold distance matrix and calculate weights (only neighboring gridcells)
load("code/FYI/Cai_et_al_2024_NEE/data/distances_grids/geodistances_grid_7.RDATA")

neigh_dist = data.frame(grid = paste("grid_",c(5:8), sep=""), dist = c(600, 340, 210, 120)) 

system.time({
  edges = data.frame(
    from = rep(grid_7@data$grid_7_ID[1:(nrow(grid_7)-1)],seq(nrow(grid_7)-1,1))[geodistances < neigh_dist$dist[which(neigh_dist$grid == "grid_7")]],
    to = unlist(sapply(2:nrow(grid_7), function(x) grid_7@data$grid_7_ID[x:nrow(grid_7)]))[geodistances < neigh_dist$dist[which(neigh_dist$grid == "grid_7")]],
    distance = geodistances[geodistances < neigh_dist$dist[which(neigh_dist$grid == "grid_7")]]
  )
})

range(table(c(edges$from,edges$to))) # should be 5,6
length(which(table(c(edges$from,edges$to))==5)) # should be 12


# Calculated individually from and to each cell, later we can take the mean of the too triangular matrices
grid_7_data = grid_7@data
rm(list = setdiff(ls(),c("overlap","grid_7_data","edges","grid_7_lm_clipped")))

cl = makePSOCKcluster(detectCores()/2, outfile = "")
clusterExport(cl, c("graph_from_data_frame","distances","overlap","grid_7_data","edges"))
registerDoParallel(cl)

gc()
graphdistances_isl = foreach(i = 1:length(overlap), .combine = rbind) %dopar% {
  
  
  if(((i)%%100==0)|i==1) print(paste(i, "of" , length(overlap)))
  
  overlap.i = overlap[[i]]$grid_7_ID
  area.i = overlap[[i]]$area[which(as.character(overlap.i) %in% colnames(grid_7_data))]
  overlap.i = overlap.i[which(as.character(overlap.i) %in% colnames(grid_7_data))]
  
  if (length(overlap.i)>0){
    if (length(overlap.i)>1){
      grid_7_data$value = rowSums(as.matrix(grid_7_data[,as.character(overlap.i)])%*%diag(area.i))/sum(area.i)
    } else {
      grid_7_data$value = grid_7_data[,as.character(overlap.i)]
    }
    
    edges$costs = sapply(1:nrow(edges),
                         function(x) mean(c(grid_7_data$value[which(grid_7_data$grid_7_ID==edges$to[x])],
                                            grid_7_data$value[which(grid_7_data$grid_7_ID==edges$from[x])])))
    edges$weight = edges$distance * edges$costs
    nodes = grid_7_data[,c(1:3)]
    
    net = graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
    
    graphdistances = distances(net, v=which(nodes$grid_7_ID %in% overlap[[i]]$grid_7_ID))
    
    graphdistances = apply(graphdistances, 2, function(x) min(x, na.rm = TRUE))
    
    graphdistances = sapply(1:length(overlap), function(x) min(graphdistances[overlap[[x]]$grid_7_ID],
                                                               na.rm=TRUE))
  } else {
    graphdistances = rep(NA, length(overlap))
  }
  return(graphdistances)
}

stopCluster(cl)

nrow(graphdistances_isl)
graphdistances_isl[1:5,1:5]

graphdistances_isl = as.matrix(graphdistances_isl)

rownames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID
colnames(graphdistances_isl) = grid_7_lm_clipped$grid_7_isl_ID


graphdistances_isl[1:5,1:5]
graphdistances_isl = (graphdistances_isl + t(graphdistances_isl))/2
graphdistances_isl[1:5,1:5]

graphdistances_clim = as.dist(graphdistances_isl)
save(graphdistances_clim,
     file="code/FYI/Cai_et_al_2024_NEE/data/distances_grids_isl_clipped/graphdistances_clim_grid_7.RDATA")




##2.3 Calculating the distances at the regional level -----
rm(list=ls())
gc()
library(maptools)
library(rgeos)
library(geosphere)
library(rgdal)
library(dplyr)

#selected entity
geoentities_simple = readOGR(dsn = "data/TDWG4/TDWG4_newTibet.shp")
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

geodistances_tdwg = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)
colnames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_tdwg[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                    as.character(overlap[[k]]$grid_7_ID)],
                                  na.rm=TRUE)
  }
}

geodistances_tdwg[upper.tri(geodistances_tdwg,
                            diag=FALSE)] = t(geodistances_tdwg)[upper.tri(geodistances_tdwg,
                                                                           diag=FALSE)]
geodistances_tdwg[1:5,1:5]

geodistances_tdwg_dist = geodistances_tdwg[order(as.numeric(row.names(geodistances_tdwg))),
                                           order(as.numeric(colnames(geodistances_tdwg)))]
geodistances_tdwg_dist[1:5,1:5]
geodistances_tdwg_dist = as.dist(geodistances_tdwg_dist)

geodistances = geodistances_tdwg_dist
save(geodistances,
     file="results/primary_results/distances_TDWG/geodistances_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full,geodistances_tdwg,geodistances_tdwg_dist,i,k)

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

geodistances_tdwg = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)
colnames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_tdwg[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_tdwg[upper.tri(geodistances_tdwg,
                            diag=FALSE)] = t(geodistances_tdwg)[upper.tri(geodistances_tdwg,
                                                                          diag=FALSE)]
geodistances_tdwg[1:5,1:5]

geodistances_tdwg_dist = geodistances_tdwg[order(as.numeric(row.names(geodistances_tdwg))),
                                           order(as.numeric(colnames(geodistances_tdwg)))]
geodistances_tdwg_dist[1:5,1:5]
geodistances_tdwg_dist = as.dist(geodistances_tdwg_dist)

graphdistances_water = geodistances_tdwg_dist
save(graphdistances_water,
     file="results/primary_results/distances_TDWG/graphdistances_water_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_water, geodistances_tdwg,geodistances_tdwg_dist,i,k)


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

geodistances_tdwg = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)
colnames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_tdwg[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_tdwg[upper.tri(geodistances_tdwg,
                            diag=FALSE)] = t(geodistances_tdwg)[upper.tri(geodistances_tdwg,
                                                                          diag=FALSE)]
geodistances_tdwg[1:5,1:5]

geodistances_tdwg_dist = geodistances_tdwg[order(as.numeric(row.names(geodistances_tdwg))),
                                           order(as.numeric(colnames(geodistances_tdwg)))]
geodistances_tdwg_dist[1:5,1:5]
geodistances_tdwg_dist = as.dist(geodistances_tdwg_dist)

graphdistances_elev = geodistances_tdwg_dist
save(graphdistances_elev,
     file="results/primary_results/distances_TDWG/graphdistances_elev_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_elev, geodistances_tdwg, geodistances_tdwg_dist,i,k)


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

geodistances_tdwg = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)
colnames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_tdwg[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_tdwg[upper.tri(geodistances_tdwg,
                            diag=FALSE)] = t(geodistances_tdwg)[upper.tri(geodistances_tdwg,
                                                                          diag=FALSE)]
geodistances_tdwg[1:5,1:5]

geodistances_tdwg_dist = geodistances_tdwg[order(as.numeric(row.names(geodistances_tdwg))),
                                           order(as.numeric(colnames(geodistances_tdwg)))]
geodistances_tdwg_dist[1:5,1:5]
geodistances_tdwg_dist = as.dist(geodistances_tdwg_dist)

graphdistances_temp = geodistances_tdwg_dist
save(graphdistances_temp,
     file="results/primary_results/distances_TDWG/graphdistances_temp_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_temp, geodistances_tdwg, geodistances_tdwg_dist,i,k)

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

geodistances_tdwg = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)
colnames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_tdwg[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_tdwg[upper.tri(geodistances_tdwg,
                            diag=FALSE)] = t(geodistances_tdwg)[upper.tri(geodistances_tdwg,
                                                                          diag=FALSE)]
geodistances_tdwg[1:5,1:5]

geodistances_tdwg_dist = geodistances_tdwg[order(as.numeric(row.names(geodistances_tdwg))),
                                           order(as.numeric(colnames(geodistances_tdwg)))]
geodistances_tdwg_dist[1:5,1:5]
geodistances_tdwg_dist = as.dist(geodistances_tdwg_dist)

graphdistances_ai = geodistances_tdwg_dist
save(graphdistances_ai,
     file="results/primary_results/distances_TDWG/graphdistances_ai_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_ai, geodistances_tdwg, geodistances_tdwg_dist,i,k)


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

geodistances_tdwg = matrix(NA, nrow = nrow(geoentities_simple), ncol = nrow(geoentities_simple))
rownames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)
colnames(geodistances_tdwg) = sort(geoentities_simple@data$RegionID)


for (i in 1:nrow(geoentities_simple)){
  print(i)
  for (k in i:nrow(geoentities_simple)){
    geodistances_tdwg[k,i] = min(geodistances_full[as.character(overlap[[i]]$grid_7_ID),
                                                   as.character(overlap[[k]]$grid_7_ID)],
                                 na.rm=TRUE)
  }
}

geodistances_tdwg[upper.tri(geodistances_tdwg,
                            diag=FALSE)] = t(geodistances_tdwg)[upper.tri(geodistances_tdwg,
                                                                          diag=FALSE)]
geodistances_tdwg[1:5,1:5]

geodistances_tdwg_dist = geodistances_tdwg[order(as.numeric(row.names(geodistances_tdwg))),
                                           order(as.numeric(colnames(geodistances_tdwg)))]
geodistances_tdwg_dist[1:5,1:5]
geodistances_tdwg_dist = as.dist(geodistances_tdwg_dist)

graphdistances_barriers = geodistances_tdwg_dist
save(graphdistances_barriers,
     file="results/primary_results/distances_TDWG/graphdistances_barriers_grid_7.RDATA")

#plot(geodistances_gift_dist, geodistances_points)
rm(geodistances_full, graphdistances_barriers, geodistances_tdwg, geodistances_tdwg_dist,i,k)



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

rownames(pastdistances_gift) = geoentities_simple@data$RegionID
colnames(pastdistances_gift) = geoentities_simple@data$RegionID


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
     file="results/primary_results/distances_TDWG/graphdistances_clim_grid_7.RDATA")



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
plot(population$gpw_v4_population_count_rev11_2020_30_min)


intersect(sort(unique(df$ISO_Code)), sort(unique(country_code$country_iso2)))
no_trade_id = setdiff(sort(unique(df$ISO_Code)), sort(unique(country_code$country_iso2)))
setdiff(sort(unique(country_code$country_iso2)), sort(unique(df$ISO_Code)))

df_country_code = df %>% left_join(country_code,
                                   by = join_by('ISO_Code' == 'country_iso2'))

df_no_trade_id = df %>% filter(ISO_Code %in% no_trade_id)
df_trade_id = df %>% filter(!(ISO_Code %in% no_trade_id))
plot(df_trade_id$geometry)
df_no_trade_id$Level_4_Na
country_code$country_name


## 3.1 Join country/region names of TDWG and trade data----
x = unique(df_no_trade_id$Level_4_Na)
y = unique(country_code$country_name)

normalize_name = function(z){
  z %>%
    stri_trans_general("Latin-ASCII") %>%   # remove accents
    str_to_lower() %>%
    str_replace_all("&", " and ") %>%
    str_replace_all("[\\.'’`]", "") %>%
    str_replace_all("[^a-z0-9 ]+", " ") %>%
    str_squish()
}

x_df = tibble(
  name_x = x,
  name_x_n = normalize_name(x)
)

y_df = tibble(
  name_y = y,
  name_y_n = normalize_name(y)
)

x_df = x_df %>%
  mutate(iso2 = countrycode(name_x, "country.name", "iso2c", warn = FALSE))

y_df = y_df %>%
  mutate(iso2 = countrycode(name_y, "country.name", "iso2c", warn = FALSE))

sort(x_df$name_x)

alias_map = tribble(
  ~alias,                    ~canonical,
  "Taiwan",                  "China",
  "Kin-men",                 "China",
  "great britain",           "United Kingdom",
  "northern ireland",        "United Kingdom",
  "channel is",              "United Kingdom",
  "French Guiana",           "France",
  "Monaco",                  "France",  
  "east timor",              "Timor-Leste",
  "faroyar",                 "Faroe Islands",
  "virgin is",               "British Virgin Islands",
  "wake i",                  "United States Minor Outlying Islands",
  "midway is",               "United States Minor Outlying Islands",
  "johnston i",              "United States Minor Outlying Islands"
) %>%
  mutate(alias_n = normalize_name(alias))


x_df2 = x_df %>%
  left_join(alias_map, by = c("name_x_n" = "alias_n")) %>%
  mutate(
    name_for_code = ifelse(is.na(canonical), name_x, canonical),
    iso2 = ifelse(
      is.na(iso2),
      countrycode(name_for_code, "country.name", "iso2c", warn = FALSE),
      iso2
    ),
    match_type = ifelse(is.na(canonical), "direct_iso2", "alias_iso2")
  )
x_df2$my_ios2 = x_df2$iso2
sort(x_df2$name_x_n)
x_df2[x_df2$name_x_n == 'french guiana',]$my_ios2 = 'FR'
x_df2[x_df2$name_x_n == 'monaco',]$my_ios2 = 'FR'
x_df2[x_df2$name_x_n == 'taiwan',]$my_ios2 = 'CN'


matched_iso2 = x_df2 %>%
  filter(!is.na(iso2)) %>%
  left_join(
    y_df %>% filter(!is.na(iso2)) %>% distinct(iso2, name_y),
    by = join_by('my_ios2' == "iso2")
  ) %>%
  mutate(match_type = paste0(match_type, "_matched")) %>%
  dplyr::select(
    df_name = name_x,
    iso2,
    my_ios2,
    country_code_name = name_y,
    match_type
  ) %>% dplyr::filter(country_code_name != 'USA')


matched_iso2 = arrange(matched_iso2, matched_iso2$df_name)

# modify some ISO_code to make some region in the TYDG match country_code 
df_2 = arrange(df, df$Level_4_Na)
df_2$my_ISO_Code = df_2$ISO_Code
df_2[df_2$Level_4_Na %in% matched_iso2$df_name,]$my_ISO_Code = matched_iso2$my_ios2

df_country_code = df_2 %>% left_join(country_code,
                                   by = join_by('my_ISO_Code' == 'country_iso2')) %>% 
  filter(!is.na(country_code))
table(sort(df_country_code$RegionID))
length(unique(df_country_code$RegionID))/length(df$RegionID)
# 96% of regions has trade data
RegionID_no_trade = setdiff(df$RegionID, unique(df_country_code$RegionID))
#save(RegionID_no_trade,
#file = 'results/primary_results/predictors_TDWG/RegionID_no_trade.rdata')

## 3.2 Join country/region names of TDWG and population data----
# make projection coincidence
df_for_population = sf::st_transform(df, crs = sf::st_crs(population))

# convert to SpatVector and do zonal
v_for_population = terra::vect(df_for_population)
# could be clip to high speed
population_crop = terra::crop(population, v_for_population, snap="out")

# zonal: by=regionID, fun="mean"
z_for_population = cbind(RegionID = v_for_population$RegionID,
                    terra::extract(population_crop,
                                   v_for_population, fun = "sum", na.rm = TRUE))
colnames(z_for_population)[which(colnames(z_for_population) == 'gpw_v4_population_count_rev11_2020_30_min')] = 'population_count_2020'
df_country_code_population = df_country_code %>% left_join(z_for_population[,c('RegionID',
                                                                               'population_count_2020')],
                                                   by = 'RegionID')
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

RegionID_direct_trade = sort(df_country_code_population_direct$RegionID)
save(RegionID_direct_trade,
     file = 'results/primary_results/predictors_TDWG/RegionID_direct_trade.rdata')


## for Fr. South Antarctic Terr. no specific population data,
## we assume three subislands divide the total population equally
df_country_code_population[df_country_code_population$country_name == 'Fr. South Antarctic Terr.',
                           ]$rel_popu = 1/3

### 3.3.1 Symmetrical matrix for regional ED----
tradeflow_mat_s = matrix(NA, nrow = length(unique(df_country_code_population$RegionID)),
                    ncol = length(unique(df_country_code_population$RegionID)))

rownames(tradeflow_mat_s) = sort(unique(df_country_code_population$RegionID))
colnames(tradeflow_mat_s) = sort(unique(df_country_code_population$RegionID))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_s)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_s)){
    #j = 1
      region_i = rownames(tradeflow_mat_s)[i]
      region_j = rownames(tradeflow_mat_s)[j]
      
      country_code_i = df_country_code_population %>% filter(RegionID == region_i) %>% 
        pull(country_code) %>% unique()
      country_code_j = df_country_code_population %>% filter(RegionID == region_j)%>% 
        pull(country_code) %>% unique()
      
      rel_popu_i = df_country_code_population %>% filter(RegionID == region_i)%>% 
        pull(rel_popu) %>% unique()
      rel_popu_j = df_country_code_population %>% filter(RegionID == region_j)%>% 
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
     file="results/primary_results/distances_TDWG/tradeflow_mat_s.RDATA")



### 3.3.1 Asymmetrical matrix for sp ED----
tradeflow_mat_a = matrix(NA, nrow = length(unique(df_country_code_population$RegionID)),
                       ncol = length(unique(df_country_code_population$RegionID)))

rownames(tradeflow_mat_a) = sort(unique(df_country_code_population$RegionID))
colnames(tradeflow_mat_a) = sort(unique(df_country_code_population$RegionID))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_a)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_a)){
    #j = 1
    region_i = rownames(tradeflow_mat_a)[i]
    region_j = rownames(tradeflow_mat_a)[j]
    
    country_code_i = df_country_code_population %>% filter(RegionID == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(RegionID == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(RegionID == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(RegionID == region_j)%>% 
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

tradeflow_mat_a['217','147'] # trade flow value from France to shanghai 
tradeflow_mat_a['147','217'] # trade flow value from shanghai to France
# 1163109 * 1000 USD
tradeflow_mat_a[217,147]

save(tradeflow_mat_a,
     file="results/primary_results/distances_TDWG/tradeflow_mat_a.RDATA")




#4. characteristics of region's phylogenetic structure----
library(picante)
library(PhyloMeasures)

##4.1 mammal's exotics PD----

load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")

sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)

mat_mammal_exotic = sp_overlap_dat_1 %>%
  as.data.frame() %>% 
  distinct(RegionID, Binomial) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(
    names_from  = Binomial,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(RegionID) %>% as.matrix()

rownames(mat_mammal_exotic) = mat_mammal_exotic[,1]
mat_mammal_exotic = mat_mammal_exotic[,-1]

library(geiger)

setdiff(colnames(mat_mammal_exotic), spec_phy.3$tip.label)
my_pd = function(x){
  tmp.tree = treedata(spec_phy.3, x[x > 0], 
                      warnings = F)$phy
  sum(tmp.tree$edge.length)
}


#pd_mammal_exotic = data.frame(RegionID = rownames(mat_mammal_exotic),
  #                          pd = pd.query(spec_phy.3, mat_mammal_exotic,
     #                                     abundance.weights = F)) ## this method is fast,
#but do not account for the root of evolutionary history, 
# especially when some regions have only one species, which would return 0
pd_mammal_exotic = data.frame(RegionID = rownames(mat_mammal_exotic),
                          pd(mat_mammal_exotic, spec_phy.3,
                                     include.root = T))

save(pd_mammal_exotic,
     file="results/primary_results/predictors_TDWG/pd_mammal_exotic.RDATA")
rm(pd_mammal_exotic, mat_mammal_exotic, sp_overlap_dat_1)


##4.2 mammal's natives PD----
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/all_phy.rdata")

mat_mammal_native = sp_dis_5 %>%
  as.data.frame() %>% 
  distinct(Region.ID, ScientificName) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('Region.ID', 'ScientificName', 'presence')) %>% 
  pivot_wider(
    names_from  = ScientificName,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(Region.ID) %>% as.matrix()

rownames(mat_mammal_native) = mat_mammal_native[,1]
mat_mammal_native = mat_mammal_native[,-1]

#pd_mammal_native = data.frame(RegionID = rownames(mat_mammal_native),
 #                             pd = pd.query(spec_phy.3, mat_mammal_native,
   #                                         abundance.weights = F))
pd_mammal_native = data.frame(RegionID = rownames(mat_mammal_native),
                              pd(mat_mammal_native, spec_phy.3,
                                      include.root = T))

save(pd_mammal_native,
     file="results/primary_results/predictors_TDWG/pd_mammal_native.RDATA")
rm(pd_mammal_native, mat_mammal_native, spec_phy.3, sp_dis_5)

##4.3 bird's exotics PD----
all_distri_data = read.csv("data/Birds/data/Distribution_data_note.csv",
                           header = T)

all_distri_data$ScientificName = gsub(' ', '_', all_distri_data$ScientificName)

all_distri_data_c = all_distri_data %>% 
  filter(seasonal %in% c(1,2) &  # only 
           ## analysed the distribution data of birds that are resident or in breeding season
           presence %in% c(1)) %>% 
  filter(SpStatus %in% c('Native', 'alien'))# only 
## analysed the distribution data of birds that is sure they are extant

exotic_distri_data = all_distri_data_c %>% filter(SpStatus == 'alien')
phy_data = read.tree("data/Birds/data/Phylogenetic_Birds.tre")

str(exotic_distri_data)
mat_bird_exotic = exotic_distri_data %>%
  as.data.frame() %>% 
  distinct(RegionID, ScientificName) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('RegionID', 'ScientificName', 'presence')) %>% 
  pivot_wider(
    names_from  = ScientificName,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(RegionID) %>% as.matrix()

rownames(mat_bird_exotic) = mat_bird_exotic[,1]
mat_bird_exotic = mat_bird_exotic[,-1]

#pd_bird_exotic = data.frame(RegionID = rownames(mat_bird_exotic),
    #                        pd = pd.query(phy_data, mat_bird_exotic,
#                                      abundance.weights = F))

pd_bird_exotic = data.frame(RegionID = rownames(mat_bird_exotic),
                              pd(mat_bird_exotic, phy_data,
                                      include.root = T))

save(pd_bird_exotic,
     file="results/primary_results/predictors_TDWG/pd_bird_exotic.RDATA")
rm(pd_bird_exotic, mat_bird_exotic, exotic_distri_data)


##4.4 bird's natives PD----
all_distri_data = read.csv("data/Birds/data/Distribution_data_note.csv",
                           header = T)

all_distri_data$ScientificName = gsub(' ', '_', all_distri_data$ScientificName)

all_distri_data_c = all_distri_data %>% 
  filter(seasonal %in% c(1,2) &  # only 
           ## analysed the distribution data of birds that are resident or in breeding season
           presence %in% c(1)) %>% 
  filter(SpStatus %in% c('Native', 'alien'))# only 
## analysed the distribution data of birds that is sure they are extant

native_distri_data = all_distri_data_c %>% filter(SpStatus == 'Native')
phy_data = read.tree("data/Birds/data/Phylogenetic_Birds.tre")


mat_bird_native = native_distri_data %>%
  as.data.frame() %>% 
  distinct(RegionID, ScientificName) %>% 
  filter(!is.na(ScientificName)) %>% 
  mutate(presence = 1)  %>%
  dplyr::select(c('RegionID', 'ScientificName', 'presence')) %>% 
  pivot_wider(
    names_from  = ScientificName,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% dplyr::arrange(RegionID) %>% as.matrix()

rownames(mat_bird_native) = mat_bird_native[,1]
mat_bird_native = mat_bird_native[,-1]

#pd_bird_native = data.frame(RegionID = rownames(mat_bird_native),
   #                     pd = pd.query(phy_data, mat_bird_native,
           #                           abundance.weights = F))
pd_bird_native = data.frame(RegionID = rownames(mat_bird_native),
                            pd(mat_bird_native, phy_data,
                                    include.root = T))

save(pd_bird_native,
     file="results/primary_results/predictors_TDWG/pd_bird_native.RDATA")
