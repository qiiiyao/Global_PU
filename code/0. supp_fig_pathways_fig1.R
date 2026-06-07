rm(list = ls())
gc()
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'scico', 'ggplot2', 'gridExtra')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("D:/R projects/Global_ED")

# load the background data for plotting the world map plot
# load the background data for plotting the world map plot
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = df
unique(df_trans$Island)
df_trans[is.na(df_trans$Island),]$Island = 0
#df_trans$area = st_area(df_trans)
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")
df_sub = st_read("data/Plants/TDWG4_Subset/TDWG4_Subset.shp")

df_trans_c = df_trans %>% st_drop_geometry()

lons = seq(-180, 180, by = 1)

# Equator line (0°)
equator_sf = st_as_sf(
  data.frame(lon = lons, lat = 0),
  coords = c("lon", "lat"), crs = 4326) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>% 
  sf::st_cast("LINESTRING")

# Tropic of Capricorn (~ -23.436°)
tropic_capricorn_sf = st_as_sf(
  data.frame(lon = lons, lat = -23.436),
  coords = c("lon", "lat"), crs = 4326) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>% 
  sf::st_cast("LINESTRING")

# Tropic of Cancer (~ +23.436°)
tropic_cancer_sf = st_as_sf(
  data.frame(lon = lons, lat = +23.436),
  coords = c("lon", "lat"), crs = 4326) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>% 
  st_cast("LINESTRING")

#1. world map ----
map_for_pathways =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "black", fill = NA) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 


library(devEMF)
library(cowplot)

emf('figures/map_for_pathways.emf',
    width = 5 * 1.2, height = 8 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways
dev.off() #turn off device and finalize file


#2. Australia map ----
library(sf)
library(terra)
library(ggplot2)
library(rnaturalearth)
library(dplyr)

cols_region = c(
  "#D7F4EF",
  "#BFD8B8",
  "#5B8B78",
  "#073B35",
  "#7A3B08")

## 1. Get Spain polygon
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")

sort(unique(df$Level_4_Na))
sort(unique(df$Level3_cod))

spain = df %>%
  filter(ISO_Code == "AU" & is.na(Island)) %>%
  st_transform(st_crs(4326))

## 2. Define source point S1
bbox_spain <- st_bbox(spain)
x_offset <- (bbox_spain["xmax"] - bbox_spain["xmin"]) * 0.08
y_offset <- (bbox_spain["ymax"] - bbox_spain["ymin"]) * 0.08

S1 <- st_as_sf(data.frame(
    id = "S1",
    lon = bbox_spain["xmin"] + x_offset,
    lat = bbox_spain["ymax"] - y_offset),
  coords = c("lon", "lat"),
  crs = 4326)

## 3. Create raster grid covering Spain
r = rast(
  ext(vect(spain)),
  resolution = 0.2,   # smaller value = smoother map
  crs = st_crs(spain)$wkt)

## 4. Convert raster cells to points
grid_pts = as.points(r) %>% 
  st_as_sf()

## 5. Keep only cells inside Spain
grid_pts = grid_pts[spain, ]

## 6. Calculate distance from each grid cell to S1
grid_pts$phy_dis = as.numeric(st_distance(grid_pts, S1)) / 1000

grid_df = grid_pts %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

spain_union = st_union(spain)

map_for_pathways_null = ggplot() +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region
  ) +
  theme_void() +
  theme(legend.position = "none")

map_for_pathways_colors1 = ggplot() +
  geom_raster(
    data = grid_df,
    aes(x = x, y = y, fill = phy_dis)
  ) +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region
  ) +
  theme_void() +
  theme(legend.position = "none")



cols_region_2 = c("#D7F4EF",
  "#BFD8B8",
  "#5B8B78",
  "#839D9A",
  "#D9BE9F")


map_for_pathways_colors2 = ggplot() +
  geom_raster(
    data = grid_df,
    aes(x = x, y = y, fill = phy_dis)
  ) +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region_2
  ) +
  theme_void() +
  theme(legend.position = "none")



cols_region_3 = c(
  "#D7F4EF",
  "#BFD8B8",
  "#5B8B78",
  "#5B8B78",
  "#D8C4B5")


map_for_pathways_colors3 = ggplot() +
  geom_raster(
    data = grid_df,
    aes(x = x, y = y, fill = phy_dis)
  ) +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region_3
  ) +
  theme_void() +
  theme(legend.position = "none")



cols_region_4 <- c(
  "#5B8B78",
  "#6F9A89",
  "#4A756A",
  "#234C45",
  "#073B35",
  "#7A3B08"
)

map_for_pathways_colors4 = ggplot() +
  geom_raster(
    data = grid_df,
    aes(x = x, y = y, fill = phy_dis)
  ) +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region_4
  ) +
  theme_void() +
  theme(legend.position = "none")



cols_region_5 = c(
    "#D8C4B5",
    "#BFD8B8",
    "#5B8B78",
    "#073B35",
    "#7A3B08")

map_for_pathways_colors5 = ggplot() +
  geom_raster(
    data = grid_df,
    aes(x = x, y = y, fill = phy_dis)
  ) +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region_5
  ) +
  theme_void() +
  theme(legend.position = "none")




cols_region_6 = c(
  "#D7F4EF",  # light cyan
  "#E8D8C5",  # very light beige
  "#D1B08B",  # light brown-yellow
  "#A8733D",  # medium brown
  "#7A3B08"   # dark brown
)

map_for_pathways_colors6 = ggplot() +
  geom_raster(
    data = grid_df,
    aes(x = x, y = y, fill = phy_dis)
  ) +
  geom_sf(
    data = spain_union,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = cols_region_6
  ) +
  theme_void() +
  theme(legend.position = "none")


library(devEMF)
library(cowplot)

emf('figures/map_for_pathways_null.emf',
    width = 5 * 1.2, height = 8 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_null
dev.off() #turn off device and finalize file


emf('figures/map_for_pathways_colors1.emf',
    width = 5 * 1.2, height = 8 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_colors1
dev.off() #turn off device and finalize file

emf('figures/map_for_pathways_colors2.emf',
    width = 5 * 1.2, height = 8 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_colors2
dev.off() #turn off device and finalize file

emf('figures/map_for_pathways_colors3.emf',
    width = 5 * 1.2, height = 8 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_colors3
dev.off() #turn off device and finalize file


emf('figures/map_for_pathways_colors4.emf',
    width = 5 * 1.2, height = 4 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_colors4
dev.off() #turn off device and finalize file

emf('figures/map_for_pathways_colors5.emf',
    width = 5 * 1.2, height = 4 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_colors5
dev.off() #turn off device and finalize file


emf('figures/map_for_pathways_colors6.emf',
    width = 5 * 1.2, height = 4 * 1.2, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
map_for_pathways_colors6
dev.off() #turn off device and finalize file


