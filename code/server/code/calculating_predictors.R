### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("foreach", "dplyr", 'tidyr', 'sf', 'terra', 'raster', 'scico'
                 #,'GIFT'
                 )

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("~/my_pc/Global_ED")

# load the background maps for quantifying the predictors
df = st_read("data/TDWG4/TDWG4_newTibet.shp")
df = sf::st_make_valid(df)  

#### Quantifying the predictors for TDWG4 ###
### 1. Environemt variables: soil, topography and five main bioclim variables
## (1) elevation and soil cation exchange ability
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

v_for_ces = terra::vect(df_for_ces)
# could be clip to high speed
ces_crop = terra::crop(ces, v_for_ces, snap="out")

# zonal: by=regionID, fun="mean"
z_for_eleva = terra::extract(elevation_crop,
                                         v_for_eleva, fun = "mean", na.rm = TRUE)

## ces data is too large, thus using the lapply fun for each polygon
z_for_ces_l = lapply(1:nrow(v_for_ces), function(i) {
  #i = 1
  terra::extract(ces_crop, v_for_ces[i,], fun = 'mean', na.rm = TRUE)
})
z_for_ces = do.call(rbind, z_for_ces_l)


## (2) five main climate variables
# read the bioclim multi-terra 
bioclim1 = terra::rast("data/envir_vars/bioclim/CHELSA_bio1_1981-2010_V.2.1.tif")
# accessed at https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.zhdk.cloud.switch.ch%2Fchelsav2%2F&prefix=GLOBAL%2Fclimatologies%2F1981-2010%2F
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
z_for_bioclim4 = terra::zonal(bioclim4_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
z_for_bioclim12 = terra::zonal(bioclim12_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
z_for_bioclim15 = terra::zonal(bioclim15_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)
z_for_arid_index = terra::zonal(arid_index_crop, v_for_bioclim, fun = "mean", na.rm = TRUE)

### aggregate and save the 7 environmental vairables
envir_vars = list(MAT = z_for_bioclim1, TS = z_for_bioclim4, MAP = z_for_bioclim12,
               PS = z_for_bioclim15, AI = z_for_arid_index,
               CECSOL = z_for_ces,
               Elevation = z_for_eleva)

save(envir_vars, file = 'results/primary_results/envir_vars.rdata')



### 3. dispersal limitation (least cost distance considering 
### physical, ecological dispersal barriers, and climate similarity)





### 4. trade intensity


