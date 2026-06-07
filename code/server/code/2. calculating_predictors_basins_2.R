### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#0 loading required packages-------------------
# list of packages
# Manually install the PhyloMeasures package
rm(list = ls())
requirements = c("foreach", "dplyr", 'tidyr', 'sf', 'terra', 'raster', 'scico')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("~/my_pc/Global_ED")

# load the background maps for quantifying the predictors
df = st_read("data/Fishes/data/Basin042017_3119/Basin042017_3119.shp")
df = sf::st_make_valid(df)  

#2. dispersal limitation ----
###(least cost distance considering physical, ecological dispersal barriers, and climate similarity)

##2.3 Calculating the distances at the regional level -----
rm(list=ls())
gc()

#selected entity
geoentities_simple = st_read("data/Fishes/data/Basin042017_3119/Basin042017_3119.shp", quiet = TRUE)
#load equal area grid
load("code/FYI/Cai_et_al_2024_NEE/data/Hexagons/grid_7.RData")
grid_7 = st_as_sf(grid_7)

geoentities_simple = st_transform(geoentities_simple, st_crs(grid_7))

geoentities_simple = st_make_valid(geoentities_simple)
grid_7 = st_make_valid(grid_7) 

overlap = st_intersects(geoentities_simple, grid_7)
overlap = lapply(overlap, function(ii) grid_7[ii, ])


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
rownames(geodistances_basins) = sort(geoentities_simple$BasinName)
colnames(geodistances_basins) = sort(geoentities_simple$BasinName)


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



