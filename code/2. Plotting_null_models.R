### Part of the code adapted from Cai_et_al_2023_PNAS
#0. Set up R environments  ----
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
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
#df_trans$area = st_area(df_trans)
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")
df_sub = st_read("data/Plants/TDWG4_Subset/TDWG4_Subset.shp")
df_trans[is.na(df_trans$Island),]$Island = 0

df_trans_c = df_trans %>% st_drop_geometry()

lons = seq(-180, 180, by = 1)

# Equator line (0°)
equator_sf = st_as_sf(
  data.frame(lon = lons, lat = 0),
  coords = c("lon", "lat"), crs = 4326) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>% 
  sf::st_cast("LINESTRING")%>% 
  st_transform(crs = "+proj=eck4")

# Tropic of Capricorn (~ -23.436°)
tropic_capricorn_sf = st_as_sf(
  data.frame(lon = lons, lat = -23.436),
  coords = c("lon", "lat"), crs = 4326) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>% 
  sf::st_cast("LINESTRING")%>% 
  st_transform(crs = "+proj=eck4")

# Tropic of Cancer (~ +23.436°)
tropic_cancer_sf = st_as_sf(
  data.frame(lon = lons, lat = +23.436),
  coords = c("lon", "lat"), crs = 4326) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>% 
  st_cast("LINESTRING")%>% 
  st_transform(crs = "+proj=eck4")


# define color gradients
scico::scico_palette_show()

colors1 = rev(scico::scico(n=8, palette = "lajolla"))
colors2 = scico::scico(n=10, palette = "vik")[1:8] 
colors3 = scico::scico(n=8, begin = 0, end = 0.4, palette = "bam")  

colors4 = c(scico::scico(n=4, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=4, begin = 0.7, end = 1, direction = 1, palette = "bam"))

colors5 = c(scico::scico(n=5, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=3, begin = 0.7, end = 1, direction = 1, palette = "bam"))

colors6 = c(scico::scico(n=3, begin = 0, end = 0.3, palette = "bam"), 
            scico::scico(n=5, begin = 0.6, end = 1, direction = 1, palette = "bam"))

colors7 = c(scico::scico(n=6, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=2, begin = 0.7, end = 1, direction = 1, palette = "bam"))

colors8 = c(scico::scico(n=7, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=1, begin = 0.8, end = 1, direction = 1, palette = "bam"))


#plot(1:8, col = colors2, pch = 19, cex = 5)

# define a function to make legend of color gradients
legend.func = function(mycolors, mylabels) {
  group = rep("cc", 8)
  condition = letters[1:8]
  value = rep(1, 8)
  df.legend = data.frame(group, condition, value)
  mycolors.corrected = rev(mycolors)
  ggplot(df.legend, aes(fill = condition, y = value, x = group)) +
    geom_bar(position = "stack", stat = "identity", color = "white") +
    scale_fill_manual(values = mycolors.corrected) +
    theme_classic() +
    theme(
      legend.position = "none", aspect.ratio = 0.03,
      axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(size = 7, color = "black")
    ) +
    scale_y_continuous(breaks = 0:8, labels = mylabels) +
    coord_flip() +
    xlab("")
}

quantiles = c('Q1', 'Median', 'Q3')

#1. Mammal turnover: mapping q1, q2, and q3 for distribution derived from null models----
load("results/primary_results/null_models/phy_turn_mammal_delta_null_summary.rdata")

#### Delta_ED mammals mapping 
mammals_delta_PU_map_null = list()
colnames(phy_turn_mammal_delta_null_summary)[which(colnames(phy_turn_mammal_delta_null_summary) == 'median')] = 'Median'
for(i in seq_len(length(quantiles))){
  #i = 1
  quant = quantiles[i]
  PU_mammal_delta_sf = df_trans %>% left_join(phy_turn_mammal_delta_null_summary[,c('RegionID', quant)],
                                              by = 'RegionID')
  colnames(PU_mammal_delta_sf)[which(colnames(PU_mammal_delta_sf) == quant)] = 'delta_PU'
  PU_mammal_delta_sf = PU_mammal_delta_sf %>% filter(!is.na(delta_PU) & 
                                                       delta_PU != 0 
  )
  
  
  #hist(PU_mammal_delta_sf$delta_PU)
  
  breaks_delta_PU = quantile(PU_mammal_delta_sf$delta_PU,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 12)
  breaks_delta_PU[which(names(breaks_delta_PU) == '50%')] = 0
  PU_mammal_delta_sf$f_delta_PU = cut(PU_mammal_delta_sf$delta_PU,
                                      breaks = breaks_delta_PU,
                                      include.lowest = TRUE)
  
  # Ensure full coverage with a tiny buffer
  
  
  mammals_delta_PU_map =
    ggplot() +
    geom_sf(data = tropic_capricorn_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = tropic_cancer_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = equator_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = countries, color = 'gray', fill = NA) +
    geom_sf(data = bb, color = "gray", fill = NA) +
    geom_sf(data = (PU_mammal_delta_sf %>% filter(Island != 1  
                                                  #& Area < 1e4
    )),
    aes(fill = f_delta_PU),
    color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
    
    scale_fill_manual(
      values = colors4,
      na.value = 'gray90',
      name = 'PU',
      guide = guide_legend(
        title = 'PU',
        title.position = "top",
        label.position = "bottom",
        direction = "horizontal",
        nrow = 1
      )
    ) +
    #scale_fill_gradientn(
    # colors = colors4,
    # values = scales::rescale(quantile(PU_mammal_delta_sf$PU,
    #                                  probs = seq(0, 1,
    #                                            length.out = 7),
    #                                na.rm = TRUE)),
    #   na.value = 'gray90',
    #  name = 'PU'
    #  ) + 
    #ggnewscale::new_scale_fill() +
    #ggplot()+
    # Circle layer for islands, colored by PE and sized by richness
    geom_point(data = (PU_mammal_delta_sf %>% filter(Island == 1  
                                                     #&Area < 1e4
    )),
    aes(x = Lon, y = Lat, color = f_delta_PU),
    size = 3,
    shape = 21, stroke = 2, fill = NA, show.legend = F) + 
    scale_color_manual(
      values = colors4,
      na.value = 'gray90',
      name = 'PU'
    ) +
    scale_size_discrete(
      range = c(1, 5)
      #,name = "Island size indicator"
    ) + 
    coord_sf(crs = "+proj=eck4",expand = FALSE) +
    theme_void() +
    ggtitle(paste0(quant,' for Delta_PU_mammals'))+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
    ) 
  
  legend.mammals_delta_PU_map = legend.func(mycolors = colors4,
                                            mylabels = round(breaks_delta_PU,
                                                             4)) +
    ggtitle("log (PU_all / PU_native)") +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
  mammals_delta_PU_map = ggplotGrob(mammals_delta_PU_map)
  legend.mammals_delta_PU_map = ggplotGrob(legend.mammals_delta_PU_map)
  mammals_delta_PU_map_all = arrangeGrob(mammals_delta_PU_map,
                                         legend.mammals_delta_PU_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
  #plot(mammals_delta_PU_map_all)
  mammals_delta_PU_map_null[[i]] = mammals_delta_PU_map_all
  
}

plot(mammals_delta_PU_map_null[[1]])
plot(mammals_delta_PU_map_null[[2]])
plot(mammals_delta_PU_map_null[[3]])


#2. Bird turnover: mapping q1, q2, and q3 for distribution derived from null models----
load("results/primary_results/null_models/phy_turn_bird_delta_null_summary.rdata")

#### Delta_ED birds mapping 
birds_delta_PU_map_null = list()
colnames(phy_turn_bird_delta_null_summary)[which(colnames(phy_turn_bird_delta_null_summary) == 'median')] = 'Median'
for(i in seq_len(length(quantiles))){
  #i = 1
  quant = quantiles[i]
  PU_bird_delta_sf = df_trans %>% left_join(phy_turn_bird_delta_null_summary[,c('RegionID', quant)],
                                              by = 'RegionID')
  colnames(PU_bird_delta_sf)[which(colnames(PU_bird_delta_sf) == quant)] = 'delta_PU'
  PU_bird_delta_sf = PU_bird_delta_sf %>% filter(!is.na(delta_PU) & 
                                                       delta_PU != 0 
  )
  
  
  #hist(PU_bird_delta_sf$delta_PU)
  
  breaks_delta_PU = quantile(PU_bird_delta_sf$delta_PU,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 12)
  breaks_delta_PU[which(names(breaks_delta_PU) == '87.5%')] = 0
  PU_bird_delta_sf$f_delta_PU = cut(PU_bird_delta_sf$delta_PU,
                                      breaks = breaks_delta_PU,
                                      include.lowest = TRUE)
  
  # Ensure full coverage with a tiny buffer
  
  
  birds_delta_PU_map =
    ggplot() +
    geom_sf(data = tropic_capricorn_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = tropic_cancer_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = equator_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = countries, color = 'gray', fill = NA) +
    geom_sf(data = bb, color = "gray", fill = NA) +
    geom_sf(data = (PU_bird_delta_sf %>% filter(Island != 1  
                                                  #& Area < 1e4
    )),
    aes(fill = f_delta_PU),
    color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
    
    scale_fill_manual(
      values = colors8,
      na.value = 'gray90',
      name = 'PU',
      guide = guide_legend(
        title = 'PU',
        title.position = "top",
        label.position = "bottom",
        direction = "horizontal",
        nrow = 1
      )
    ) +
    #scale_fill_gradientn(
    # colors = colors8,
    # values = scales::rescale(quantile(PU_bird_delta_sf$PU,
    #                                  probs = seq(0, 1,
    #                                            length.out = 7),
    #                                na.rm = TRUE)),
    #   na.value = 'gray90',
    #  name = 'PU'
    #  ) + 
    #ggnewscale::new_scale_fill() +
    #ggplot()+
    # Circle layer for islands, colored by PE and sized by richness
    geom_point(data = (PU_bird_delta_sf %>% filter(Island == 1  
                                                     #&Area < 1e4
    )),
    aes(x = Lon, y = Lat, color = f_delta_PU),
    size = 3,
    shape = 21, stroke = 2, fill = NA, show.legend = F) + 
    scale_color_manual(
      values = colors8,
      na.value = 'gray90',
      name = 'PU'
    ) +
    scale_size_discrete(
      range = c(1, 5)
      #,name = "Island size indicator"
    ) + 
    coord_sf(crs = "+proj=eck4",expand = FALSE) +
    theme_void() +
    ggtitle(paste0(quant,' for Delta_PU_birds'))+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
    ) 
  
  legend.birds_delta_PU_map = legend.func(mycolors = colors8,
                                            mylabels = round(breaks_delta_PU,
                                                             4)) +
    ggtitle("log (PU_all / PU_native)") +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
  birds_delta_PU_map = ggplotGrob(birds_delta_PU_map)
  legend.birds_delta_PU_map = ggplotGrob(legend.birds_delta_PU_map)
  birds_delta_PU_map_all = arrangeGrob(birds_delta_PU_map,
                                         legend.birds_delta_PU_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
  #plot(birds_delta_PU_map_all)
  birds_delta_PU_map_null[[i]] = birds_delta_PU_map_all
  
}

plot(birds_delta_PU_map_null[[1]])
plot(birds_delta_PU_map_null[[2]])
plot(birds_delta_PU_map_null[[3]])



#3. Plant turnover: mapping q1, q2, and q3 for distribution derived from null models----
shp.glonaf.trans = st_read("data/Plants/shp_glonaf_new_eck4.shp")
load("results/primary_results/null_models/phy_turn_plant_delta_null_summary.rdata")

#### Delta_ED plants mapping 
plants_delta_PU_map_null = list()
colnames(phy_turn_plant_delta_null_summary)[which(colnames(phy_turn_plant_delta_null_summary) == 'median')] = 'Median'

for(i in seq_len(length(quantiles))){
  #i = 1
  quant = quantiles[i]
  PU_plant_delta_sf = shp.glonaf.trans %>% left_join(phy_turn_plant_delta_null_summary[,c('Region_id', quant)],
                                            by = 'Region_id')
  colnames(PU_plant_delta_sf)[which(colnames(PU_plant_delta_sf) == quant)] = 'delta_PU'
  PU_plant_delta_sf = PU_plant_delta_sf %>% filter(!is.na(delta_PU) & 
                                                   delta_PU != 0 
  )
  
  
  #hist(PU_plant_delta_sf$delta_PU)
  
  breaks_delta_PU = quantile(PU_plant_delta_sf$delta_PU,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 12)
  breaks_delta_PU[which(names(breaks_delta_PU) == '87.5%')] = 0
  PU_plant_delta_sf$f_delta_PU = cut(PU_plant_delta_sf$delta_PU,
                                    breaks = breaks_delta_PU,
                                    include.lowest = TRUE)
  
  # Ensure full coverage with a tiny buffer
  
  
  plants_delta_PU_map =
    ggplot() +
    geom_sf(data = tropic_capricorn_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = tropic_cancer_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = equator_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = countries, color = 'gray', fill = NA) +
    geom_sf(data = bb, color = "gray", fill = NA) +
    geom_sf(data = (PU_plant_delta_sf %>% filter(Island != 1  
                                                #& Area < 1e4
    )),
    aes(fill = f_delta_PU),
    color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
    
    scale_fill_manual(
      values = colors8,
      na.value = 'gray90',
      name = 'PU',
      guide = guide_legend(
        title = 'PU',
        title.position = "top",
        label.position = "bottom",
        direction = "horizontal",
        nrow = 1
      )
    ) +
    #scale_fill_gradientn(
    # colors = colors8,
    # values = scales::rescale(quantile(PU_plant_delta_sf$PU,
    #                                  probs = seq(0, 1,
    #                                            length.out = 7),
    #                                na.rm = TRUE)),
    #   na.value = 'gray90',
    #  name = 'PU'
    #  ) + 
    #ggnewscale::new_scale_fill() +
    #ggplot()+
    # Circle layer for islands, colored by PE and sized by richness
    geom_point(data = (PU_plant_delta_sf %>% filter(Island == 1  
                                                   #&Area < 1e4
    )),
    aes(x = Lon, y = Lat, color = f_delta_PU),
    size = 3,
    shape = 21, stroke = 2, fill = NA, show.legend = F) + 
    scale_color_manual(
      values = colors8,
      na.value = 'gray90',
      name = 'PU'
    ) +
    scale_size_discrete(
      range = c(1, 5)
      #,name = "Island size indicator"
    ) + 
    coord_sf(crs = "+proj=eck4",expand = FALSE) +
    theme_void() +
    ggtitle(paste0(quant,' for Delta_PU_plants'))+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
    ) 
  
  legend.plants_delta_PU_map = legend.func(mycolors = colors8,
                                          mylabels = round(breaks_delta_PU,
                                                           4)) +
    ggtitle("log (PU_all / PU_native)") +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
  plants_delta_PU_map = ggplotGrob(plants_delta_PU_map)
  legend.plants_delta_PU_map = ggplotGrob(legend.plants_delta_PU_map)
  plants_delta_PU_map_all = arrangeGrob(plants_delta_PU_map,
                                       legend.plants_delta_PU_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
  #plot(plants_delta_PU_map_all)
  plants_delta_PU_map_null[[i]] = plants_delta_PU_map_all
  
}

plot(plants_delta_PU_map_null[[1]])
plot(plants_delta_PU_map_null[[2]])
plot(plants_delta_PU_map_null[[3]])




#4. Fish turnover: mapping q1, q2, and q3 for distribution derived from null models----
load("results/primary_results/null_models/phy_turn_fish_delta_null_summary.rdata")
basins_trans = st_read("data/Fishes/data/Basin042017_3119_eck4.shp")
colnames(basins_trans)[which(colnames(basins_trans) == 'BasinName')] = 'X1.Basin.Name'

#### Delta_ED fishs mapping 
fishs_delta_PU_map_null = list()
colnames(phy_turn_fish_delta_null_summary)[which(colnames(phy_turn_fish_delta_null_summary) == 'median')] = 'Median'

for(i in seq_len(length(quantiles))){
  #i = 2
  quant = quantiles[i]
  PU_fish_delta_sf = basins_trans %>% left_join(phy_turn_fish_delta_null_summary[,c('X1.Basin.Name', quant)],
                                            by = 'X1.Basin.Name')
  colnames(PU_fish_delta_sf)[which(colnames(PU_fish_delta_sf) == quant)] = 'delta_PU'
  PU_fish_delta_sf = PU_fish_delta_sf %>% filter(!is.na(delta_PU) & 
                                                   delta_PU != 0 
  )
  
  
  #hist(PU_fish_delta_sf$delta_PU)
  breaks_delta_PU = quantile(PU_fish_delta_sf$delta_PU,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = TRUE,
                             digits = 12)
  
  first_positive = which(breaks_delta_PU > 0)[1]
  if (!is.na(first_positive)) {
    breaks_delta_PU[first_positive] = 0
  }
  
  if (first_positive == 6) {
    colors_my = colors5
  } else if (first_positive == 7) {
    colors_my = colors7
  } else if (first_positive == 8) {
    colors_my = colors8
  }

  
  PU_fish_delta_sf$f_delta_PU = cut(PU_fish_delta_sf$delta_PU,
                                    breaks = breaks_delta_PU,
                                    include.lowest = TRUE)
  
  # Ensure full coverage with a tiny buffer
  
  # 1. Get target CRS from polygon layer
  target_crs = st_crs(PU_fish_delta_sf)
  
  # 2. Extract attributes of island points (Island == 1) and drop polygon geometry
  island_fishes_delta = PU_fish_delta_sf %>%
    filter(Island == 1) %>%
    st_drop_geometry()   # now a plain data frame with Lon, Lat, f_turnover, etc.
  
  # 3. Create point sf from Lon/Lat (WGS84) and transform to target CRS
  island_fishes_delta_sf = island_fishes_delta %>%
    st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(crs = target_crs)
  
  # 4. Extract projected coordinates (X, Y) as numeric columns
  coords = st_coordinates(island_fishes_delta_sf)   # matrix with columns X, Y
  island_fishes_delta_sf = island_fishes_delta_sf %>%
    st_drop_geometry() %>%
    mutate(X_proj = coords[,1], Y_proj = coords[,2])
  
  
  fishes_delta_PU_map =
    ggplot() +
    geom_sf(data = tropic_capricorn_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = tropic_cancer_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = equator_sf,
            color = "gray",
            linewidth = 0.8,
            linetype = "solid") +
    geom_sf(data = countries, color = 'gray', fill = NA) +
    geom_sf(data = bb, color = "gray", fill = NA) +
    geom_sf(data = (PU_fish_delta_sf %>% filter(Island != 1  
                                                #& Area < 1e4
    )),
    aes(fill = f_delta_PU),
    color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
    
    scale_fill_manual(
      values = colors_my,
      na.value = 'gray90',
      name = 'PU',
      guide = guide_legend(
        title = 'PU',
        title.position = "top",
        label.position = "bottom",
        direction = "horizontal",
        nrow = 1
      )
    ) +
    #scale_fill_gradientn(
    # colors = colors_my,
    # values = scales::rescale(quantile(PU_fish_delta_sf$PU,
    #                                  probs = seq(0, 1,
    #                                            length.out = 7),
    #                                na.rm = TRUE)),
    #   na.value = 'gray90',
    #  name = 'PU'
    #  ) + 
    #ggnewscale::new_scale_fill() +
    #ggplot()+
    # Circle layer for islands, colored by PE and sized by richness
    geom_point(data = island_fishes_delta_sf,
               aes(x = X_proj, y = Y_proj, color = f_delta_PU),
               size = 3,
               shape = 21, stroke = 2, fill = NA, show.legend = F) + 
    scale_color_manual(
      values = colors_my,
      na.value = 'gray90',
      name = 'PU'
    ) +
    scale_size_discrete(
      range = c(1, 5)
      #,name = "Island size indicator"
    ) + 
    coord_sf(crs = "+proj=eck4",expand = FALSE) +
    theme_void() +
    ggtitle(paste0(quant,' for Delta_PU_fishes'))+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
    ) 
  
  legend.fishes_delta_PU_map = legend.func(mycolors = colors_my,
                                           mylabels = round(breaks_delta_PU,
                                                            4)) +
    ggtitle("log (PU_all / PU_native)") +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
  fishes_delta_PU_map = ggplotGrob(fishes_delta_PU_map)
  legend.fishes_delta_PU_map = ggplotGrob(legend.fishes_delta_PU_map)
  fishes_delta_PU_map_all = arrangeGrob(fishes_delta_PU_map,
                                        legend.fishes_delta_PU_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
  #plot(fishes_delta_PU_map_all)
  fishs_delta_PU_map_null[[i]] = fishes_delta_PU_map_all
  
}

plot(fishs_delta_PU_map_null[[1]])
plot(fishs_delta_PU_map_null[[2]])
plot(fishs_delta_PU_map_null[[3]])




#5. Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

##5.2 PU #####
# Compare the turnover patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.01         # small horizontal gap between columns

figs_null_PU_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_delta_PU_map_null[[1]],  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_PU_map_null[[2]],  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_PU_map_null[[3]],   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_delta_PU_map_null[[1]],   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_PU_map_null[[2]],   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_PU_map_null[[3]],    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_delta_PU_map_null[[1]], x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_PU_map_null[[2]], x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_PU_map_null[[3]],  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishs_delta_PU_map_null[[1]],   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishs_delta_PU_map_null[[2]],   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishs_delta_PU_map_null[[3]],   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
  # Labels a–l
  draw_plot_label(
    label = letters[1:12],
    size = 13,
    x = rep(c(0.01, plot_width + gap + 0.01, 2*(plot_width + gap) + 0.01), 4),
    y = rep(c(0.75 + plot_height - 0.015,
              0.5 + plot_height - 0.015,
              0.25 + plot_height - 0.015,
              0 + plot_height - 0.015), each = 3)
  )

png(filename = 'figures/figs_null_PU_nat_extant_delta_202605_1.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_null_PU_nat_extant_delta
dev.off() #turn off device and finalize file

