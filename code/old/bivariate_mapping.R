### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("PhyloMeasures", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'terra', 'scico', 'ggplot2', 'gridExtra')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("D:/R projects/Global_ED")
source('code/functions/calculating_LCBD_func.R')

# load the background data for plotting the world map plot
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")


# define color gradients
scico::scico_palette_show()

colors1 = rev(scico::scico(n=8, palette = "lajolla"))
colors2 = scico::scico(n=10, palette = "vik")[1:8] 
colors3 = scico::scico(n=8, begin = 0, end = 0.4, palette = "bam")  
colors4 = scico::scico(n=8, palette = "bam")
colors5 = c(scico::scico(n=5, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=3, begin = 0.7, end = 1, direction = 1, palette = "bam"))

colors6 = c(scico::scico(n=3, begin = 0, end = 0.3, palette = "bam"), 
            scico::scico(n=5, begin = 0.6, end = 1, direction = 1, palette = "bam"))

colors7 = c(scico::scico(n=6, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=2, begin = 0.7, end = 1, direction = 1, palette = "bam"))

#plot(1:8, col = colors2, pch = 19, cex = 5)

# load supp. functions for plotting bivariate world maps
source('code/functions/plot_functions.R')
n.breaks = 10


### COASTLINES
## Shapefile downloaded from
## <https://www.ngdc.noaa.gov/mgg/shorelines/>
wrld.pol = vect("code/FYI/Baselga_et_al_2025_ELE/world-map/GSHHS_c_L1.shp", "GSHHS_c_L1")
wrld.pol.moll = project(wrld.pol, "+proj=eck4")
wrld.pol.moll = as(wrld.pol.moll, "Spatial")





#### Mammal ED & LCBD: bivariate mapping ####
load('results/primary_results/ED_mammal_native.rdata')
load('results/primary_results/ED_mammal_all.rdata')
load('results/primary_results/LCBD_mammal_native.rdata')
load('results/primary_results/LCBD_mammal_all.rdata')

load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")

#plot(df)

# spatial vectors
ext(df_trans)
empty.raster = raster(xmn=-16653225, xmx=16774959, ymn=-8460601, ymx=8375228,
                      nrows=180, ncols=360, vals=1,crs="+proj=eck4")

##### Native ######
#### ED 
ED_mammal_native_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_native$df)
ED_mammal_native_sf = df_trans %>% left_join(ED_mammal_native_1, by = 'RegionID')

# rasterize and project to Mollweide projection
#ED_mammal_native_raster = terra::rasterize(ED_mammal_native_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_mammal_native_raster = raster::rasterize(ED_mammal_native_sf,
                                              empty.raster, field = "mean_ED")

#ED_mammal_native_raster_moll = raster::projectRaster(ED_mammal_native_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_mammal_native_raster_moll_mask = terra::mask(ED_mammal_native_raster,
                                                   wrld.pol.moll)

plot(ED_mammal_native_raster_moll_mask)





#### LCBD
LCBD_mammal_native_1 = LCBD_mammal_native$LCBD_simp_geo
LCBD_mammal_native_sf = df_trans %>% left_join(LCBD_mammal_native_1, by = 'RegionID')


# spatial vectors
#LCBD_mammal_native_points = terra::vect(LCBD_mammal_native_sf)

crs(empty.raster)
crs(LCBD_mammal_native_sf) 
ext(empty.raster)
st_bbox(LCBD_mammal_native_sf)

# rasterize and project to Mollweide projection
#LCBD_mammal_native_raster = terra::rasterize(LCBD_mammal_native_points,
      #                                     empty.raster,
            #                               field="LCBD", background=NA)

LCBD_mammal_native_raster = raster::rasterize(LCBD_mammal_native_sf,
                                              empty.raster, field = "LCBD")

#LCBD_mammal_native_raster_moll = raster::projectRaster(LCBD_mammal_native_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_mammal_native_raster_moll_mask = raster::mask(LCBD_mammal_native_raster,
                                                wrld.pol.moll)

plot(LCBD_mammal_native_raster_moll_mask)


## bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_mammal_native_raster_moll_mask = rasterCM(ED_mammal_native_raster_moll_mask,
                                  LCBD_mammal_native_raster_moll_mask,
                                  n=n.breaks)



cmat.pe.pu = makeCM(n.breaks, "purple", "black", "grey95", "green4")

pdf("figures/mammal_native_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_mammal_native_raster_moll_mask,
      col=cmat.pe.pu, ylab="", xlab="", xaxt="n", yaxt="n"
      , xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
      )
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM(cm = cmat.pe.pu, main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()



##### All: native + alien #####
#### ED 
ED_mammal_all_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_all$df)
ED_mammal_all_sf = df_trans %>% left_join(ED_mammal_all_1, by = 'RegionID')

crs(empty.raster)
crs(ED_mammal_all_sf) 
ext(empty.raster)
st_bbox(ED_mammal_all_sf)

# rasterize and project to Mollweide projection
#ED_mammal_all_raster = terra::rasterize(ED_mammal_all_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_mammal_all_raster = raster::rasterize(ED_mammal_all_sf,
                                            empty.raster, field = "mean_ED")

#ED_mammal_all_raster_moll = raster::projectRaster(ED_mammal_all_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_mammal_all_raster_moll_mask = terra::mask(ED_mammal_all_raster,
                                                wrld.pol.moll)

plot(ED_mammal_all_raster_moll_mask)





#### LCBD
LCBD_mammal_all_1 = LCBD_mammal_all$LCBD_simp_geo
LCBD_mammal_all_sf = df_trans %>% left_join(LCBD_mammal_all_1, by = 'RegionID')

#LCBD_mammal_all_points = terra::vect(LCBD_mammal_all_sf)

crs(empty.raster)
crs(LCBD_mammal_all_sf) 
ext(empty.raster)
st_bbox(LCBD_mammal_all_sf)

# rasterize and project to Mollweide projection
#LCBD_mammal_all_raster = terra::rasterize(LCBD_mammal_all_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_mammal_all_raster = raster::rasterize(LCBD_mammal_all_sf,
                                              empty.raster, field = "LCBD")

#LCBD_mammal_all_raster_moll = raster::projectRaster(LCBD_mammal_all_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_mammal_all_raster_moll_mask = raster::mask(LCBD_mammal_all_raster,
                                                   wrld.pol.moll)

plot(LCBD_mammal_all_raster_moll_mask)


## bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_mammal_all_raster_moll_mask = rasterCM(ED_mammal_all_raster_moll_mask,
                                                  LCBD_mammal_all_raster_moll_mask,
                                                  n=n.breaks)

cmat.pe.pu = makeCM(n.breaks, "purple", "black", "grey95", "green4")



pdf("figures/mammal_all_ED_LCBD_map_bivar.pdf", height=10, width=15)

image(ED_LCBD_mammal_all_raster_moll_mask,
      col=cmat.pe.pu, ylab="", xlab="", xaxt="n", yaxt="n"
      , xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM(cm = cmat.pe.pu, main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()


##### Delta: (native + alien) / native #####
#### ED 
ED_mammal_all_1 = data.frame(RegionID = sort(unique(sp_dis_5$Region.ID)),
                        mean_all_ED = ED_mammal_all$df$mean_ED)
ED_mammal_native_1 = data.frame(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           mean_native_ED = ED_mammal_native$df$mean_ED)
Delta_ED_mammal = ED_mammal_native_1 %>% left_join(ED_mammal_all_1,
                                                   by = 'RegionID')

Delta_ED_mammal$delta_mean_ED = log(Delta_ED_mammal$mean_all_ED / 
                                      Delta_ED_mammal$mean_native_ED)


ED_mammal_delta_sf = df_trans %>% left_join(Delta_ED_mammal, by = 'RegionID')#%>% 
  #filter(!(delta_mean_ED == 0))


# rasterize and project to Mollweide projection
#ED_mammal_delta_raster = terra::rasterize(ED_mammal_delta_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_mammal_delta_raster = raster::rasterize(ED_mammal_delta_sf,
                                         empty.raster, field = "delta_mean_ED")

#ED_mammal_delta_raster_moll = raster::projectRaster(ED_mammal_delta_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_mammal_delta_raster_moll_mask = terra::mask(ED_mammal_delta_raster,
                                             wrld.pol.moll)

plot(ED_mammal_delta_raster_moll_mask)





#### LCBD
LCBD_mammal_all_1 = LCBD_mammal_all$LCBD_simp_geo %>% rename('all_LCBD' = 'LCBD') %>% 
  rename('all_SSI' = 'SSI')
LCBD_mammal_native_1 = LCBD_mammal_native$LCBD_simp_geo %>% rename('native_LCBD' = 'LCBD') %>% 
  rename('native_SSI' = 'SSI')
Delta_LCBD_mammal = LCBD_mammal_native_1 %>% left_join(LCBD_mammal_all_1,
                                                   by = 'RegionID')

Delta_LCBD_mammal$delta_LCBD = log(Delta_LCBD_mammal$all_LCBD / 
                                      Delta_LCBD_mammal$native_LCBD)


LCBD_mammal_delta_sf = df_trans %>% left_join(Delta_LCBD_mammal, by = 'RegionID') #%>% 
 # filter(!(delta_LCBD == 0))

# rasterize and project to Mollweide projection
#LCBD_mammal_delta_raster = terra::rasterize(LCBD_mammal_delta_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_mammal_delta_raster = raster::rasterize(LCBD_mammal_delta_sf,
                                           empty.raster, field = "delta_LCBD")

#LCBD_mammal_delta_raster_moll = raster::projectRaster(LCBD_mammal_delta_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_mammal_delta_raster_moll_mask = terra::mask(LCBD_mammal_delta_raster,
                                               wrld.pol.moll)

plot(LCBD_mammal_delta_raster_moll_mask)



#### bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_mammal_delta_raster_moll_mask = rasterCM_centered(ED_mammal_delta_raster_moll_mask,
                                               LCBD_mammal_delta_raster_moll_mask,
                                               n=2)

cmat.pe.pu = t(makeCM_bilinear(
  breaks     = 2,
  upperleft  =  "#D6ABF1",
  upperright = "#BD6432",
  lowerleft  =  "#1D6E9C",
  lowerright = "#86C486" 
))

pdf("figures/mammal_delta_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_mammal_delta_raster_moll_mask,
      col=cmat.pe.pu,
      ylab="", xlab="", xaxt="n", yaxt="n"
      , xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM_centred(cm = t(cmat.pe.pu), main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()



#### Plant ED & LCBD: bivariate mapping ####
load('results/primary_results/ED_plant_native.rdata')
load('results/primary_results/ED_plant_extant.rdata')
load('results/primary_results/LCBD_plant_native.rdata')
load('results/primary_results/LCBD_plant_extant.rdata')

load("data/Plants/data/shp.651.Rdata")
load("data/Plants/data/phylo.fake.species.653.Rdata")
load("data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)

#plot(df)

# spatial vectors
#empty.raster = rast(nrows=180, ncols=360, nlyrs=1, xmin=-16653225, xmax=16774959, 
#         ymin=-8460601, ymax=8375228, crs="+proj=eck4", vals = 1)

df.extant.650 = rbind(df.native.650, df.natu.650)

crs(shp.glonaf.trans) 
ext(shp.glonaf.trans)

empty.raster = raster(xmn=ceiling(ext(shp.glonaf.trans)[1]),
                      xmx=ceiling(ext(shp.glonaf.trans)[2]),
                      ymn=ceiling(ext(shp.glonaf.trans)[3]),
                      ymx=ceiling(ext(shp.glonaf.trans)[4]),
                      nrows=180, ncols=360, vals=1,crs="+proj=eck4")


##### Native ######
#### ED 
ED_plant_native_1 = cbind(Region_id = sort(unique(df.native.650$Region_id)),
                          ED_plant_native$df)
ED_plant_native_sf = shp.glonaf.trans %>% left_join(ED_plant_native_1, by = 'Region_id')

ED_plant_native_points = terra::vect(ED_plant_native_sf)



# rasterize and project to Mollweide projection
#ED_plant_native_raster = terra::rasterize(ED_plant_native_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_plant_native_raster = raster::rasterize(ED_plant_native_sf,
                                            empty.raster, field = "mean_ED")

#ED_plant_native_raster_moll = raster::projectRaster(ED_plant_native_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_plant_native_raster_moll_mask = terra::mask(ED_plant_native_raster,
                                                wrld.pol.moll)

plot(ED_plant_native_raster_moll_mask)





#### LCBD
LCBD_plant_native_1 = LCBD_plant_native$LCBD_simp_geo
LCBD_plant_native_sf = shp.glonaf.trans %>% left_join(LCBD_plant_native_1,
                                                      join_by('Region_id' == 'RegionID'))


# rasterize and project to Mollweide projection
#LCBD_plant_native_raster = terra::rasterize(LCBD_plant_native_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_plant_native_raster = raster::rasterize(LCBD_plant_native_sf,
                                              empty.raster, field = "LCBD")

#LCBD_plant_native_raster_moll = raster::projectRaster(LCBD_plant_native_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_plant_native_raster_moll_mask = raster::mask(LCBD_plant_native_raster,
                                                   wrld.pol.moll)

plot(LCBD_plant_native_raster_moll_mask)


## bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_plant_native_raster_moll_mask = rasterCM(ED_plant_native_raster_moll_mask,
                                                  LCBD_plant_native_raster_moll_mask,
                                                  n=n.breaks)

cmat.pe.pu = makeCM(n.breaks, "purple", "black", "grey95", "green4")

pdf("figures/plant_native_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_plant_native_raster_moll_mask,
      col=cmat.pe.pu, ylab="", xlab="", xaxt="n", yaxt="n"
      #, xlim=c(-16503000, 16454000), ylim=c(-6635000, 8375000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM(cm = cmat.pe.pu, main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()



##### All: native + alien #####
#### ED 

ED_plant_all_1 = cbind(Region_id = sort(unique(df.extant.650$Region_id)),
                          ED_plant_extant$df)
ED_plant_all_sf = shp.glonaf.trans %>% left_join(ED_plant_all_1, by = 'Region_id')

# rasterize and project to Mollweide projection
#ED_plant_all_raster = terra::rasterize(ED_plant_all_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_plant_all_raster = raster::rasterize(ED_plant_all_sf,
                                         empty.raster, field = "mean_ED")

#ED_plant_all_raster_moll = raster::projectRaster(ED_plant_all_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_plant_all_raster_moll_mask = terra::mask(ED_plant_all_raster,
                                             wrld.pol.moll)

plot(ED_plant_all_raster_moll_mask)





#### LCBD
LCBD_plant_all_1 = LCBD_plant_extant$LCBD_simp_geo
LCBD_plant_all_sf = shp.glonaf.trans %>% left_join(LCBD_plant_all_1,
                                                   join_by('Region_id' == 'RegionID'))

# rasterize and project to Mollweide projection
#LCBD_plant_all_raster = terra::rasterize(LCBD_plant_all_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_plant_all_raster = raster::rasterize(LCBD_plant_all_sf,
                                           empty.raster, field = "LCBD")

#LCBD_plant_all_raster_moll = raster::projectRaster(LCBD_plant_all_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_plant_all_raster_moll_mask = raster::mask(LCBD_plant_all_raster,
                                                wrld.pol.moll)

plot(LCBD_plant_all_raster_moll_mask)


## bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_plant_all_raster_moll_mask = rasterCM(ED_plant_all_raster_moll_mask,
                                               LCBD_plant_all_raster_moll_mask,
                                               n=n.breaks)

cmat.pe.pu = makeCM(n.breaks, "purple", "black", "grey95", "green4")

pdf("figures/plant_all_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_plant_all_raster_moll_mask,
      col=cmat.pe.pu, ylab="", xlab="", xaxt="n", yaxt="n"
      #, xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM(cm = cmat.pe.pu, main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()




##### Delta: (native + alien) / native #####
#### ED 
ED_plant_all_1 = data.frame(Region_id = sort(unique(df.extant.650$Region_id)),
                       mean_all_ED = ED_plant_extant$df$mean_ED)
ED_plant_native_1 = data.frame(Region_id = sort(unique(df.native.650$Region_id)),
                          mean_native_ED = ED_plant_native$df$mean_ED)
Delta_ED_plant = ED_plant_native_1 %>% left_join(ED_plant_all_1,
                                                   by = 'Region_id')

Delta_ED_plant$delta_mean_ED = log(Delta_ED_plant$mean_all_ED / 
                                      Delta_ED_plant$mean_native_ED)


ED_plant_delta_sf = shp.glonaf.trans %>% left_join(Delta_ED_plant, by = 'Region_id')


# rasterize and project to Mollweide projection
#ED_plant_delta_raster = terra::rasterize(ED_plant_delta_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_plant_delta_raster = raster::rasterize(ED_plant_delta_sf,
                                           empty.raster, field = "delta_mean_ED")

#ED_plant_delta_raster_moll = raster::projectRaster(ED_plant_delta_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_plant_delta_raster_moll_mask = terra::mask(ED_plant_delta_raster,
                                               wrld.pol.moll)

plot(ED_plant_delta_raster_moll_mask)





#### LCBD
LCBD_plant_all_1 = LCBD_plant_extant$LCBD_simp_geo %>% rename('all_LCBD' = 'LCBD') %>% 
  rename('all_SSI' = 'SSI')
LCBD_plant_native_1 = LCBD_plant_native$LCBD_simp_geo %>% rename('native_LCBD' = 'LCBD') %>% 
  rename('native_SSI' = 'SSI')
Delta_LCBD_plant = LCBD_plant_native_1 %>% left_join(LCBD_plant_all_1,
                                                       by = 'RegionID')

Delta_LCBD_plant$delta_LCBD = log(Delta_LCBD_plant$all_LCBD / 
                                     Delta_LCBD_plant$native_LCBD)


LCBD_plant_delta_sf = shp.glonaf.trans %>% left_join(Delta_LCBD_plant,
                                                     join_by('Region_id' == 'RegionID'))

# rasterize and project to Mollweide projection
#LCBD_plant_delta_raster = terra::rasterize(LCBD_plant_delta_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_plant_delta_raster = raster::rasterize(LCBD_plant_delta_sf,
                                             empty.raster, field = "delta_LCBD")

#LCBD_plant_delta_raster_moll = raster::projectRaster(LCBD_plant_delta_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_plant_delta_raster_moll_mask = terra::mask(LCBD_plant_delta_raster,
                                                 wrld.pol.moll)

plot(LCBD_plant_delta_raster_moll_mask)



#### bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_plant_delta_raster_moll_mask = rasterCM_centered(ED_plant_delta_raster_moll_mask,
                                                          LCBD_plant_delta_raster_moll_mask,
                                                          n=2)

cmat.pe.pu = t(makeCM_bilinear(
  breaks     = 2,
  upperleft  =  "#D6ABF1",
  upperright = "#BD6432",
  lowerleft  =  "#1D6E9C",
  lowerright = "#86C486" 
))

pdf("figures/plant_delta_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_plant_delta_raster_moll_mask,
      col=cmat.pe.pu,
      ylab="", xlab="", xaxt="n", yaxt="n"
      , xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM_centred(cm = t(cmat.pe.pu), main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()








#### Bird ED & LCBD: bivariate mapping ####
load('results/primary_results/ED_bird_native.rdata')
load('results/primary_results/ED_bird_extant.rdata')
load('results/primary_results/LCBD_bird_native.rdata')
load('results/primary_results/LCBD_bird_extant.rdata')

all_distri_data = read.csv("data/Birds/data/Distribution_data_note.csv",
                           header = T)
colnames(all_distri_data)
str(all_distri_data)
unique(all_distri_data$SpStatus)
all_distri_data$ScientificName = gsub(' ', '_', all_distri_data$ScientificName)

all_distri_data_c = all_distri_data %>% 
  filter(seasonal %in% c(1,2) &  # only 
           ## analysed the distribution data of birds that are resident or in breeding season
           presence %in% c(1)) %>% 
  filter(SpStatus %in% c('Native', 'alien'))# only 
## analysed the distribution data of birds that is sure they are extant

native_distri_data = all_distri_data_c %>% filter(SpStatus == 'Native')
exotic_distri_data = all_distri_data_c %>% filter(SpStatus == 'alien')
phy_data = read.tree("data/Birds/data/Phylogenetic_Birds.tre")


#plot(df)

# spatial vectors
ext(df_trans)
empty.raster = raster(xmn=-16653225, xmx=16774959, ymn=-8460601, ymx=8375228,
                      nrows=180, ncols=360, vals=1,crs="+proj=eck4")

##### Native ######
#### ED 
ED_bird_native_1 = cbind(RegionID = df_trans$RegionID,
                         ED_bird_native$df)
ED_bird_native_sf = df_trans %>% left_join(ED_bird_native_1, by = 'RegionID')

# rasterize and project to Mollweide projection
#ED_bird_native_raster = terra::rasterize(ED_bird_native_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_bird_native_raster = raster::rasterize(ED_bird_native_sf,
                                            empty.raster, field = "mean_ED")

#ED_bird_native_raster_moll = raster::projectRaster(ED_bird_native_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_bird_native_raster_moll_mask = terra::mask(ED_bird_native_raster,
                                                wrld.pol.moll)

plot(ED_bird_native_raster_moll_mask)





#### LCBD
LCBD_bird_native_1 = LCBD_bird_native$LCBD_simp_geo
LCBD_bird_native_sf = df_trans %>% left_join(LCBD_bird_native_1, by = 'RegionID')


# rasterize and project to Mollweide projection
#LCBD_bird_native_raster = terra::rasterize(LCBD_bird_native_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_bird_native_raster = raster::rasterize(LCBD_bird_native_sf,
                                              empty.raster, field = "LCBD")

#LCBD_bird_native_raster_moll = raster::projectRaster(LCBD_bird_native_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_bird_native_raster_moll_mask = raster::mask(LCBD_bird_native_raster,
                                                   wrld.pol.moll)

plot(LCBD_bird_native_raster_moll_mask)


## bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_bird_native_raster_moll_mask = rasterCM(ED_bird_native_raster_moll_mask,
                                                  LCBD_bird_native_raster_moll_mask,
                                                  n=n.breaks)

cmat.pe.pu = makeCM(n.breaks, "purple", "black", "grey95", "green4")

pdf("figures/bird_native_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_bird_native_raster_moll_mask,
      col=cmat.pe.pu, ylab="", xlab="", xaxt="n", yaxt="n"
      #, xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM(cm = cmat.pe.pu, main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()



##### All: native + alien #####
#### ED 
ED_bird_all_1 = cbind(RegionID = df_trans$RegionID,
                        ED_bird_extant$df)
ED_bird_all_sf = df_trans %>% left_join(ED_bird_all_1, by = 'RegionID')


# rasterize and project to Mollweide projection
#ED_bird_all_raster = terra::rasterize(ED_bird_all_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_bird_all_raster = raster::rasterize(ED_bird_all_sf,
                                         empty.raster, field = "mean_ED")

#ED_bird_all_raster_moll = raster::projectRaster(ED_bird_all_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_bird_all_raster_moll_mask = terra::mask(ED_bird_all_raster,
                                             wrld.pol.moll)

plot(ED_bird_all_raster_moll_mask)



#### LCBD
LCBD_bird_all_1 = LCBD_bird_extant$LCBD_simp_geo
LCBD_bird_all_sf = df_trans %>% left_join(LCBD_bird_all_1, by = 'RegionID')

# rasterize and project to Mollweide projection
#LCBD_bird_all_raster = terra::rasterize(LCBD_bird_all_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_bird_all_raster = raster::rasterize(LCBD_bird_all_sf,
                                           empty.raster, field = "LCBD")

#LCBD_bird_all_raster_moll = raster::projectRaster(LCBD_bird_all_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_bird_all_raster_moll_mask = raster::mask(LCBD_bird_all_raster,
                                                wrld.pol.moll)

plot(LCBD_bird_all_raster_moll_mask)


## bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_bird_all_raster_moll_mask = rasterCM(ED_bird_all_raster_moll_mask,
                                               LCBD_bird_all_raster_moll_mask,
                                               n=n.breaks)



cmat.pe.pu = makeCM(n.breaks, "purple", "black", "grey95", "green4")

pdf("figures/bird_all_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_bird_all_raster_moll_mask,
      col=cmat.pe.pu, ylab="", xlab="", xaxt="n", yaxt="n"
      , xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM(cm = cmat.pe.pu, main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()


##### Delta: (native + alien) / native #####
#### ED 
ED_bird_all_1 = data.frame(RegionID = df_trans$RegionID,
                             mean_all_ED = ED_bird_extant$df$mean_ED)
ED_bird_native_1 = data.frame(RegionID = df_trans$RegionID,
                              mean_native_ED = ED_bird_native$df$mean_ED)

Delta_ED_bird = ED_bird_native_1 %>% left_join(ED_bird_all_1,
                                                   by = 'RegionID')

Delta_ED_bird$delta_mean_ED = log(Delta_ED_bird$mean_all_ED / 
                                  Delta_ED_bird$mean_native_ED)


ED_bird_delta_sf = df_trans %>% left_join(Delta_ED_bird, by = 'RegionID')


# rasterize and project to Mollweide projection
#ED_bird_delta_raster = terra::rasterize(ED_bird_delta_points,
#                                     empty.raster,
#                               field="ED", background=NA)

ED_bird_delta_raster = raster::rasterize(ED_bird_delta_sf,
                                           empty.raster, field = "delta_mean_ED")

#ED_bird_delta_raster_moll = raster::projectRaster(ED_bird_delta_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

ED_bird_delta_raster_moll_mask = terra::mask(ED_bird_delta_raster,
                                               wrld.pol.moll)

plot(ED_bird_delta_raster_moll_mask)





#### LCBD
LCBD_bird_all_1 = LCBD_bird_extant$LCBD_simp_geo %>% rename('all_LCBD' = 'LCBD') %>% 
  rename('all_SSI' = 'SSI')
LCBD_bird_native_1 = LCBD_bird_native$LCBD_simp_geo %>% rename('native_LCBD' = 'LCBD') %>% 
  rename('native_SSI' = 'SSI')
Delta_LCBD_bird = LCBD_bird_native_1 %>% left_join(LCBD_bird_all_1,
                                                       by = 'RegionID')

Delta_LCBD_bird$delta_LCBD = log(Delta_LCBD_bird$all_LCBD / 
                                     Delta_LCBD_bird$native_LCBD)


LCBD_bird_delta_sf = df_trans %>% left_join(Delta_LCBD_bird, by = 'RegionID')

# rasterize and project to Mollweide projection
#LCBD_bird_delta_raster = terra::rasterize(LCBD_bird_delta_points,
#                                     empty.raster,
#                               field="LCBD", background=NA)

LCBD_bird_delta_raster = raster::rasterize(LCBD_bird_delta_sf,
                                             empty.raster, field = "delta_LCBD")

#LCBD_bird_delta_raster_moll = raster::projectRaster(LCBD_bird_delta_raster,
#                                       crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=km +datum=WGS84")

LCBD_bird_delta_raster_moll_mask = terra::mask(LCBD_bird_delta_raster,
                                                 wrld.pol.moll)

plot(LCBD_bird_delta_raster_moll_mask)



#### bivariate map phylo endemism vs. phylo uniqueness

# make bivariate ratster
ED_LCBD_bird_delta_raster_moll_mask = rasterCM_centered(ED_bird_delta_raster_moll_mask,
                                                          LCBD_bird_delta_raster_moll_mask,
                                                          n=2)

cmat.pe.pu = t(makeCM_bilinear(
  breaks     = 2,
  upperleft  =  "#D6ABF1",
  upperright = "#BD6432",
  lowerleft  =  "#1D6E9C",
  lowerright = "#86C486" 
))

pdf("figures/bird_delta_ED_LCBD_map_bivar.pdf", height=10, width=15)
image(ED_LCBD_bird_delta_raster_moll_mask,
      col=cmat.pe.pu,
      ylab="", xlab="", xaxt="n", yaxt="n"
      , xlim=c(-17000000, 17003000), ylim=c(-8500000, 8500000)
)
lines(wrld.pol.moll)
TeachingDemos::subplot( {
  op = par(mar = c(2, 2, 1, 1))     # enlarge margins for labels
  on.exit(par(op), add = TRUE)
  
  plotCM_centred(cm = t(cmat.pe.pu), main = "")
  
  mtext("Species Evol. distinc.",   side = 1, line = 0.5, cex = 1.6)   # x-axis title
  mtext("Regional Evol. distinc.", side = 2, line = 0.5, cex = 1.6)   # y-axis title
},
x=grconvertX(c(0.01,0.26), from='npc'),
y=grconvertY(c(0.01,0.42), from='npc'),
type='fig', pars = list(mar = c(2, 2, 0.5, 0.5)))

dev.off()



