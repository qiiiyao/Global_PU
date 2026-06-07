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


# define color gradients
scico::scico_palette_show()

colors1 = rev(scico::scico(n=8, palette = "lajolla"))
colors2 = scico::scico(n=10, palette = "vik")[1:8] 

colors3 = scico::scico(n=8, begin = 0, end = 0.4, palette = "bam")  

colors4 = scico::scico(n=8, begin = 0.5, end = 1, palette = "bam")  

colors5 = c(scico::scico(n=4, begin = 0, end = 0.4, palette = "bam"), 
            scico::scico(n=4, begin = 0.6, end = 1, direction = 1, palette = "bam"))

colors6 = c(scico::scico(n=5, begin = 0, end = 0.3, palette = "bam"), 
            scico::scico(n=3, begin = 0.6, end = 1, direction = 1, palette = "bam"))

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


#1. Mammal sor: calculation & mapping----

##1.2 natives----
load("results/primary_results/distances_beta/phy_sor_mammal_native.rdata")

sor_mammal_native = data.frame(RegionID = rownames(phy_sor_mammal_native),
                               sor = rowMeans(phy_sor_mammal_native))
sor_mammal_native$RegionID = as.integer(sor_mammal_native$RegionID)

sor_mammal_native_sf = df_trans %>% left_join(sor_mammal_native, by = 'RegionID')
colnames(sor_mammal_native_sf)

#### native mammals mapping 
sor_mammal_native_sf = sor_mammal_native_sf %>% filter(!is.na(sor))

sor_mammal_native_sf$f_sor = cut(sor_mammal_native_sf$sor,
                                 breaks = quantile(sor_mammal_native_sf$sor,
                                                   probs = seq(0, 1,
                                                               length.out = 9),
                                                   na.rm = TRUE),
                                 include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

mammals_native_sor_map =
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
  geom_sf(data = (sor_mammal_native_sf %>% filter(Island != 1 #& 
                                                      #Area < 1e4
                                                      )),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_mammal_native_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_mammal_native_sf %>% filter(Island == 1)),
             aes(x = Lon, y = Lat, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Native_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_native_sor_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(sor_mammal_native_sf$sor,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_sor_map = ggplotGrob(mammals_native_sor_map)
legend.mammals_native_sor_map = ggplotGrob(legend.mammals_native_sor_map)
mammals_native_sor_map_all = arrangeGrob(mammals_native_sor_map,
                                         legend.mammals_native_sor_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(mammals_native_sor_map_all)



##1.3 extant regional assemlage: natives + exotics ----
load("results/primary_results/distances_beta/phy_sor_mammal_extant.rdata")

sor_mammal_all = data.frame(RegionID = rownames(phy_sor_mammal_extant),
                            sor = rowMeans(phy_sor_mammal_extant))
sor_mammal_all$RegionID = as.integer(sor_mammal_all$RegionID)

sor_mammal_all_sf = df_trans %>% left_join(sor_mammal_all, by = 'RegionID')
colnames(sor_mammal_all_sf)

#### all mammals mapping 
sor_mammal_all_sf = sor_mammal_all_sf %>% filter(!is.na(sor))

sor_mammal_all_sf$f_sor = cut(sor_mammal_all_sf$sor,
                              breaks = quantile(sor_mammal_all_sf$sor,
                                                probs = seq(0, 1,
                                                            length.out = 9),
                                                na.rm = TRUE),
                              include.lowest = TRUE)


mammals_all_sor_map =
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
  geom_sf(data =  (sor_mammal_all_sf %>% filter(Island != 1)),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_mammal_all_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_mammal_all_sf %>% filter(Island == 1)),
             aes(x = Lon, y = Lat, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Extant_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_all_sor_map = legend.func(mycolors = colors2,
                                         mylabels = round(quantile(sor_mammal_all_sf$sor,
                                                                   probs = seq(0, 1,length.out = 9),
                                                                   na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_all_sor_map = ggplotGrob(mammals_all_sor_map)
legend.mammals_all_sor_map = ggplotGrob(legend.mammals_all_sor_map)
mammals_all_sor_map_all = arrangeGrob(mammals_all_sor_map,
                                      legend.mammals_all_sor_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(mammals_all_sor_map_all)




##1.4 Delta_sor (Extant - Native) ----
load("results/primary_results/distances_beta/phy_sor_mammal_extant.rdata")
load("results/primary_results/distances_beta/phy_sor_mammal_native.rdata")

### alternative 1. detla = mean(extant - native)
sor_mammal_delta_mat = log((phy_sor_mammal_extant+0.001) /
                             (phy_sor_mammal_native+0.001))
sor_mammal_delta = data.frame(RegionID = rownames(sor_mammal_delta_mat),
                              delta_sor = rowMeans(sor_mammal_delta_mat))
sor_mammal_delta$RegionID = as.integer(sor_mammal_delta$RegionID)

sor_mammal_delta_sf = df_trans %>% left_join(sor_mammal_delta,
                                             by = 'RegionID')

colnames(sor_mammal_delta_sf)

#### Delta_ED mammals mapping 
sor_mammal_delta_sf = sor_mammal_delta_sf %>% filter(!is.na(delta_sor) & 
                                                       delta_sor != 0 
)


hist(sor_mammal_delta_sf$delta_sor)

breaks_delta_sor = quantile(sor_mammal_delta_sf$delta_sor,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_sor[which(names(breaks_delta_sor) == '87.5%')] = 0
sor_mammal_delta_sf$f_delta_sor = cut(sor_mammal_delta_sf$delta_sor,
                                      breaks = breaks_delta_sor,
                                      include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

mammals_delta_sor_map =
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
  geom_sf(data = (sor_mammal_delta_sf %>% filter(Island != 1)),
          aes(fill = f_delta_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors8,
  # values = scales::rescale(quantile(sor_mammal_delta_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_mammal_delta_sf %>% filter(Island == 1)),
             aes(x = Lon, y = Lat, color = f_delta_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Delta_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_delta_sor_map = legend.func(mycolors = colors8,
                                           mylabels = round(breaks_delta_sor,
                                                            4)) +
  ggtitle("log (sor_all / sor_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_delta_sor_map = ggplotGrob(mammals_delta_sor_map)
legend.mammals_delta_sor_map = ggplotGrob(legend.mammals_delta_sor_map)
mammals_delta_sor_map_all = arrangeGrob(mammals_delta_sor_map,
                                        legend.mammals_delta_sor_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(mammals_delta_sor_map_all)




##1.5 Delta_PU (Extant - Native) ----
load("results/primary_results/distances_beta/phy_sor_mammal_extant.rdata")
load("results/primary_results/distances_beta/phy_sor_mammal_native.rdata")

colnames(sor_mammal_all)[which(colnames(sor_mammal_all) == 'sor')] = 'all_sor'
PU_mammal_delta = sor_mammal_native %>% 
  left_join(sor_mammal_all, by = 'RegionID')
PU_mammal_delta$delta_PU = log(PU_mammal_delta$all_sor/
                                 PU_mammal_delta$sor)

head(PU_mammal_delta %>% arrange(desc(delta_PU)) %>% 
       left_join(df_trans_c[,c('RegionID', 'Level_4_Na')], by = 'RegionID'))

PU_mammal_delta_sf = df_trans %>% left_join(PU_mammal_delta,
                                            by = 'RegionID')

head(PU_mammal_delta %>% arrange(desc(delta_PU)) %>% 
       left_join(df_trans_c[,c('RegionID', 'Level_4_Na')], by = 'RegionID'))

sp_dis_5 %>% filter(Region.ID == 469) %>% pull(ScientificName) %>% unique()
sp_overlap_dat_1 %>% filter(RegionID == 469) %>% pull(Binomial) %>% unique()

PU_mammal_native %>% filter(RegionID == 469)
PU_mammal_all %>% filter(RegionID == 469)
PU_mammal_delta %>% filter(RegionID == 469)
mean(PU_mammal_delta_mat[,c('469')])
mean(phy_turn_mammal_native[,c('469')])
mean(phy_turn_mammal_extant[,c('469')])
mean(log((phy_turn_mammal_extant[,c('469')]+0.000001)/
           (phy_turn_mammal_native[,c('469')]+0.000001)))
log(0.2367636/0.2177374)

colnames(PU_mammal_delta_sf)

#### Delta_ED mammals mapping 
PU_mammal_delta_sf = PU_mammal_delta_sf %>% filter(!is.na(delta_PU) & 
                                                     delta_PU != 0 
)


hist(PU_mammal_delta_sf$delta_PU)

breaks_delta_PU = quantile(PU_mammal_delta_sf$delta_PU,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)

end_negative = sort(which(breaks_delta_PU < 0), decreasing = T)[1]
if (!is.na(end_negative)) {
  breaks_delta_PU[end_negative] = 0
}

if (end_negative == 6) {
  colors_my = colors5
} else if (end_negative == 7) {
  colors_my = colors7
} else if (end_negative == 8) {
  colors_my = colors8
}

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
  ggtitle('Delta_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_delta_PU_map = legend.func(mycolors = colors8,
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
plot(mammals_delta_PU_map_all)





## 1.6 Patitioning delta ED into 5 possible ways ----
load("results/primary_results/phy_sor_mammal_native.rdata")
colnames(sp_dis_5)
colnames(sp_overlap_dat_1)
colnames(sp_overlap_dat_1)[which(colnames(sp_overlap_dat_1) == 'Inv. Stage')] = 'Inv._stage'

sp_dis_5$Inv._stage = 'native'
colnames(sp_dis_5)[which(colnames(sp_dis_5) %in% c("ScientificName", 
                                                   "Order.1.2",
                                                   "Family.1.2",
                                                   "Region.ID"))] = c("RegionID",
                                                                      "Binomial", 
                                                                      "Order",
                                                                      "Family")

sp_overlap_dat_2 = sp_overlap_dat_1 %>% 
  dplyr::select(which(colnames(sp_overlap_dat_1) %in% colnames(sp_dis_5)))
sp_dis_5 = as_tibble(sp_dis_5)
sp_dis_6 = sp_dis_5[,colnames(sp_overlap_dat_2)]

sp_dis_6 = sp_dis_6 %>% filter(!is.na(Binomial))

sp_mammal_all = rbind(sp_dis_6,
                      sp_overlap_dat_2)

region_pairs = combn(colnames(phy_sor_mammal_native), 2)

sor_mammal_path3_5_mat = phy_sor_mammal_native
sor_mammal_path4_7_mat = phy_sor_mammal_native
sor_mammal_path6_mat = phy_sor_mammal_native

gc()
for (i in 1:ncol(region_pairs)) {
  #i = 1
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = sp_dis_6 %>% filter(RegionID == region1)
  nati_region2 = sp_dis_6 %>% filter(RegionID == region2)
  natu_region1 = sp_overlap_dat_2 %>% filter(RegionID == region1)
  natu_region2 = sp_overlap_dat_2 %>% filter(RegionID == region2)
  
  nati_region1_sps = nati_region1$Binomial
  nati_region2_sps = nati_region2$Binomial
  natu_region1_sps = setdiff(natu_region1$Binomial, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$Binomial, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    path6_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_2,
                                                              natu_2_sps_ori_1))]
    sp_mammal_path6 = rbind(nati_region1, 
                            nati_region2,
                            natu_region1 %>% filter(Binomial %in% path6_sp))
    
    comm_mammal_path6 = sp_mammal_path6 %>% 
      complete(Binomial, RegionID, fill = list(presence = 0)) %>%
      dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
      pivot_wider(names_from = Binomial,
                  values_from = presence,
                  values_fn = mean) %>% 
      dplyr::arrange('RegionID') 
    
    phy_sor_mammal_path6 = calcu_sor_simple(
      Tree = spec_phy.3,
      Comm = comm_mammal_path6[,2:ncol(comm_mammal_path6)])
    
    sor_mammal_path6_mat[region1, region2] = phy_sor_mammal_path6[1,2]
    sor_mammal_path6_mat[region2, region1] = phy_sor_mammal_path6[1,2]
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_mammal_path3_5 = rbind(nati_region1, 
                                nati_region2,
                                natu_region1 %>% filter(Binomial %in% path3_5_sp))
      
      comm_mammal_path3_5 = sp_mammal_path3_5 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_mammal_path3_5 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_mammal_path3_5[,2:ncol(comm_mammal_path3_5)])
      
      sor_mammal_path3_5_mat[region1, region2] = phy_sor_mammal_path3_5[1,2]
      sor_mammal_path3_5_mat[region2, region1] = phy_sor_mammal_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      path4_7_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_mammal_path4_7 = rbind(nati_region1, 
                                nati_region2,
                                natu_region1 %>% filter(Binomial %in% path4_7_sp))
      
      comm_mammal_path4_7 = sp_mammal_path4_7 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_mammal_path4_7 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_mammal_path4_7[,2:ncol(comm_mammal_path4_7)])
      
      sor_mammal_path4_7_mat[region1, region2] = phy_sor_mammal_path4_7[1,2]
      sor_mammal_path4_7_mat[region2, region1] = phy_sor_mammal_path4_7[1,2]
      
    }
    
  }
  
}


sor_mammal_delta_path6_mat = log((sor_mammal_path6_mat+0.001) /
                                   (phy_sor_mammal_native+0.001))
which(sor_mammal_delta_path6_mat > 0)

sor_mammal_delta_path4_7_mat = log((sor_mammal_path4_7_mat+0.001) /
                                     (sor_mammal_native_mat+0.001))
which(sor_mammal_delta_path4_7_mat > 0)
which(sor_mammal_delta_path4_7_mat < 0)

sor_mammal_delta_path3_5_mat = log((sor_mammal_path3_5_mat+0.001) /
                                     (sor_mammal_native_mat+0.001))
which(sor_mammal_delta_path3_5_mat > 0)
which(sor_mammal_delta_path3_5_mat < 0)



#2. Plant sor: calculation & mapping ----
shp.glonaf.trans = st_read("data/Plants/shp_glonaf_new_eck4.shp")
unique(shp.glonaf.trans$Island)
#phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")

colnames(df.native.650) == colnames(df.natu.650)


##2.3 natives #####
load("results/primary_results/distances_beta/phy_sor_plant_native.rdata") # load the calculated data from the server

sor_plant_native = data.frame(Region_id = rownames(phy_sor_plant_native),
                              sor = rowMeans(phy_sor_plant_native))
sor_plant_native$Region_id = as.integer(sor_plant_native$Region_id)

sor_plant_native_sf = shp.glonaf.trans %>% left_join(sor_plant_native,
                                                     by = 'Region_id')
colnames(sor_plant_native_sf)

#### native plants mapping 
sor_plant_native_sf = sor_plant_native_sf %>% filter(!is.na(sor))

sor_plant_native_sf$f_sor = cut(sor_plant_native_sf$sor,
                                breaks = quantile(sor_plant_native_sf$sor,
                                                  probs = seq(0, 1,
                                                              length.out = 9),
                                                  na.rm = TRUE),
                                include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
plants_native_sor_map =
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
  geom_sf(data = (sor_plant_native_sf %>% filter(Island != 1)),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_plant_native_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_plant_native_sf %>% filter(Island == 1)),
             aes(x = Lon, y = Lat, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Native_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_native_sor_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(sor_plant_native_sf$sor,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("Phylogenetic uniqueness (PU): Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_sor_map = ggplotGrob(plants_native_sor_map)
legend.plants_native_sor_map = ggplotGrob(legend.plants_native_sor_map)
plants_native_sor_map_all = arrangeGrob(plants_native_sor_map,
                                        legend.plants_native_sor_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(plants_native_sor_map_all)




##2.4 extant #####
load("results/primary_results/distances_beta/phy_sor_plant_extant.rdata") # load the calculated data from the server

sor_plant_extant = data.frame(Region_id = rownames(phy_sor_plant_extant),
                              sor = rowMeans(phy_sor_plant_extant))
sor_plant_extant$Region_id = as.integer(sor_plant_extant$Region_id)

sor_plant_extant_sf = shp.glonaf.trans %>% left_join(sor_plant_extant,
                                                     by = 'Region_id')
colnames(sor_plant_extant_sf)

#### extant plants mapping 
sor_plant_extant_sf = sor_plant_extant_sf %>% filter(!is.na(sor))

sor_plant_extant_sf$f_sor = cut(sor_plant_extant_sf$sor,
                                breaks = quantile(sor_plant_extant_sf$sor,
                                                  probs = seq(0, 1,
                                                              length.out = 9),
                                                  na.rm = TRUE),
                                include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

plants_extant_sor_map =
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
  geom_sf(data = (sor_plant_extant_sf %>% filter(Island != 1)),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_plant_extant_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_plant_extant_sf %>% filter(Island == 1)),
             aes(x = Lon, y = Lat, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Extant_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_extant_sor_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(sor_plant_extant_sf$sor,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_sor_map = ggplotGrob(plants_extant_sor_map)
legend.plants_extant_sor_map = ggplotGrob(legend.plants_extant_sor_map)
plants_extant_sor_map_all = arrangeGrob(plants_extant_sor_map,
                                        legend.plants_extant_sor_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(plants_extant_sor_map_all)



##2.5 Delta_sor (Extant - Native) #####
load('results/primary_results/distances_beta/phy_sor_plant_extant.rdata')
load('results/primary_results/distances_beta/phy_sor_plant_native.rdata')

### alternative 1. detla = mean(extant - native)
sor_plant_delta_mat = log((phy_sor_plant_extant+0.001) /
                            (phy_sor_plant_native+0.001))
sor_plant_delta = data.frame(Region_id = rownames(sor_plant_delta_mat),
                             delta_sor = rowMeans(sor_plant_delta_mat))
sor_plant_delta$Region_id = as.integer(sor_plant_delta$Region_id)

sor_plant_delta_sf = shp.glonaf.trans %>% left_join(sor_plant_delta,
                                                    by = 'Region_id')

colnames(sor_plant_delta_sf)


#### Delta_ED plants mapping 
sor_plant_delta_sf = sor_plant_delta_sf %>% filter(!is.na(delta_sor) & 
                                                     delta_sor != 0 
)

hist(sor_plant_delta_sf$delta_sor)

breaks_delta_sor = quantile(sor_plant_delta_sf$delta_sor,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_sor[which(names(breaks_delta_sor) == '87.5%')] = 0
sor_plant_delta_sf$f_delta_sor = cut(sor_plant_delta_sf$delta_sor,
                                     breaks = breaks_delta_sor,
                                     include.lowest = TRUE)


plants_delta_sor_map =
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
  geom_sf(data = (sor_plant_delta_sf %>% filter(!(Island == 1 & 
                                                    Area < 1e4))),
          aes(fill = f_delta_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors8,
  # values = scales::rescale(quantile(sor_plant_delta_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_plant_delta_sf %>% filter(Island == 1 & 
                                                     Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_delta_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) +
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Delta_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_delta_sor_map = legend.func(mycolors = colors8,
                                          mylabels = round(breaks_delta_sor,
                                                           4)) +
  ggtitle("log (sor_all / sor_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_sor_map = ggplotGrob(plants_delta_sor_map)
legend.plants_delta_sor_map = ggplotGrob(legend.plants_delta_sor_map)
plants_delta_sor_map_all = arrangeGrob(plants_delta_sor_map,
                                       legend.plants_delta_sor_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(plants_delta_sor_map_all)



##2.6 Delta_PU (Extant - Native) ----
load("results/primary_results/distances_beta/phy_sor_plant_extant.rdata")
load("results/primary_results/distances_beta/phy_sor_plant_native.rdata")

colnames(sor_plant_extant)[which(colnames(sor_plant_extant) == 'sor')] = 'all_sor'
PU_plant_delta = sor_plant_native %>% 
  left_join(sor_plant_extant, by = 'Region_id')
PU_plant_delta$delta_PU = log(PU_plant_delta$all_sor/
                                PU_plant_delta$sor)

head(PU_plant_delta %>% arrange(desc(delta_PU)) %>% 
       left_join(shp.glonaf.trans[,c('Region_id', 'Area')], by = 'Region_id'))

PU_plant_delta_sf = shp.glonaf.trans %>% left_join(PU_plant_delta,
                                                   by = 'Region_id')


colnames(PU_plant_delta_sf)

#### Delta_ED plants mapping 
PU_plant_delta_sf = PU_plant_delta_sf %>% filter(!is.na(delta_PU) & 
                                                   delta_PU != 0 
)


hist(PU_plant_delta_sf$delta_PU)

breaks_delta_PU = quantile(PU_plant_delta_sf$delta_PU,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)

end_negative = sort(which(breaks_delta_PU < 0), decreasing = T)[1]
if (!is.na(end_negative)) {
  breaks_delta_PU[end_negative] = 0
}

if (end_negative == 6) {
  colors_my = colors5
} else if (end_negative == 7) {
  colors_my = colors7
} else if (end_negative == 8) {
  colors_my = colors8
}

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
  ggtitle('Delta_plants')+
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
plot(plants_delta_PU_map_all)


##2.6 Patitioning delta ED into 5 possible ways ----
colnames(sp_dis_5)
colnames(sp_overlap_dat_1)
colnames(sp_overlap_dat_1)[which(colnames(sp_overlap_dat_1) == 'Inv. Stage')] = 'Inv._stage'

sp_dis_5$Inv._stage = 'native'
colnames(sp_dis_5)[which(colnames(sp_dis_5) %in% c("ScientificName", 
                                                   "Order.1.2",
                                                   "Family.1.2",
                                                   "Region.ID"))] = c("RegionID",
                                                                      "Binomial", 
                                                                      "Order",
                                                                      "Family")

sp_overlap_dat_2 = sp_overlap_dat_1 %>% 
  dplyr::select(which(colnames(sp_overlap_dat_1) %in% colnames(sp_dis_5)))
sp_dis_5 = as_tibble(sp_dis_5)
sp_dis_6 = sp_dis_5[,colnames(sp_overlap_dat_2)]

sp_dis_6 = sp_dis_6 %>% filter(!is.na(Binomial))

sp_plant_all = rbind(sp_dis_6,
                     sp_overlap_dat_2)

region_pairs = combn(colnames(sor_plant_delta_mat), 2)

sor_plant_path3_5_mat = sor_plant_native_mat
sor_plant_path4_7_mat = sor_plant_native_mat
sor_plant_path6_mat = sor_plant_native_mat

for (i in 1:ncol(region_pairs)) {
  #i = 1
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = sp_dis_6 %>% filter(RegionID == region1)
  nati_region2 = sp_dis_6 %>% filter(RegionID == region2)
  natu_region1 = sp_overlap_dat_2 %>% filter(RegionID == region1)
  natu_region2 = sp_overlap_dat_2 %>% filter(RegionID == region2)
  
  nati_region1_sps = nati_region1$Binomial
  nati_region2_sps = nati_region2$Binomial
  natu_region1_sps = setdiff(natu_region1$Binomial, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$Binomial, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    path6_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_2,
                                                              natu_2_sps_ori_1))]
    sp_plant_path6 = rbind(nati_region1, 
                           nati_region2,
                           natu_region1 %>% filter(Binomial %in% path6_sp))
    
    comm_plant_path6 = sp_plant_path6 %>% 
      complete(Binomial, RegionID, fill = list(presence = 0)) %>%
      dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
      pivot_wider(names_from = Binomial,
                  values_from = presence,
                  values_fn = mean) %>% 
      dplyr::arrange('RegionID') 
    
    phy_sor_plant_path6 = calcu_sor_simple(
      Tree = spec_phy.3,
      Comm = comm_plant_path6[,2:ncol(comm_plant_path6)])
    
    sor_plant_path6_mat[region1, region2] = phy_sor_plant_path6[1,2]
    sor_plant_path6_mat[region2, region1] = phy_sor_plant_path6[1,2]
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_plant_path3_5 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(Binomial %in% path3_5_sp))
      
      comm_plant_path3_5 = sp_plant_path3_5 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_plant_path3_5 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_plant_path3_5[,2:ncol(comm_plant_path3_5)])
      
      sor_plant_path3_5_mat[region1, region2] = phy_sor_plant_path3_5[1,2]
      sor_plant_path3_5_mat[region2, region1] = phy_sor_plant_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      path4_7_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_plant_path4_7 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(Binomial %in% path4_7_sp))
      
      comm_plant_path4_7 = sp_plant_path4_7 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_plant_path4_7 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_plant_path4_7[,2:ncol(comm_plant_path4_7)])
      
      sor_plant_path4_7_mat[region1, region2] = phy_sor_plant_path4_7[1,2]
      sor_plant_path4_7_mat[region2, region1] = phy_sor_plant_path4_7[1,2]
      
    }
    
  }
  
}


sor_plant_delta_path6_mat = log((sor_plant_path6_mat+0.001) /
                                  (sor_plant_native_mat+0.001))

sor_plant_delta_path4_7_mat = log((sor_plant_path4_7_mat+0.001) /
                                    (sor_plant_native_mat+0.001))
which(sor_plant_delta_path4_7_mat > 0)
which(sor_plant_delta_path4_7_mat < 0)

sor_plant_delta_path3_5_mat = log((sor_plant_path3_5_mat+0.001) /
                                    (sor_plant_native_mat+0.001))
which(sor_plant_delta_path3_5_mat > 0)
which(sor_plant_delta_path3_5_mat < 0)




#3. Bird sor: calculation & mapping ####
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
#phy_data = read.tree("data/Birds/data/Phylogenetic_Birds.tre")
#is.rooted(phy_data)


##3.2 natives #####
load("results/primary_results/distances_beta/phy_sor_bird_native.rdata") # load the calculated data from the server

sor_bird_native = data.frame(RegionID = rownames(phy_sor_bird_native),
                             sor = rowMeans(phy_sor_bird_native))
sor_bird_native$RegionID = as.integer(sor_bird_native$RegionID)

sor_bird_native_sf = df_trans %>% left_join(sor_bird_native, by = 'RegionID')
colnames(sor_bird_native_sf)

#### native birds mapping 
sor_bird_native_sf = sor_bird_native_sf %>% filter(!is.na(sor))

sor_bird_native_sf$f_sor = cut(sor_bird_native_sf$sor,
                               breaks = quantile(sor_bird_native_sf$sor,
                                                 probs = seq(0, 1,
                                                             length.out = 9),
                                                 na.rm = TRUE),
                               include.lowest = TRUE)

sor_bird_native_sf = sor_bird_native_sf %>% filter(RegionID != 37)
birds_native_sor_map =
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
  geom_sf(data = (sor_bird_native_sf %>% filter(!(Island == 1 & 
                                                    Area < 1e4))),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_bird_native_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_bird_native_sf %>% filter(Island == 1 & 
                                                     Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Native_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_native_sor_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(sor_bird_native_sf$sor,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_sor_map = ggplotGrob(birds_native_sor_map)
legend.birds_native_sor_map = ggplotGrob(legend.birds_native_sor_map)
birds_native_sor_map_all = arrangeGrob(birds_native_sor_map,
                                       legend.birds_native_sor_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(birds_native_sor_map_all)




##3.3 extant #####
load("results/primary_results/distances_beta/phy_sor_bird_extant.rdata") # load the calculated data from the server

sor_bird_extant = data.frame(RegionID = rownames(phy_sor_bird_extant),
                             sor = rowMeans(phy_sor_bird_extant))
sor_bird_extant$RegionID = as.integer(sor_bird_extant$RegionID)

sor_bird_extant_sf = df_trans %>% left_join(sor_bird_extant, by = 'RegionID')
colnames(sor_bird_extant_sf)

#### extant birds mapping 
sor_bird_extant_sf = sor_bird_extant_sf %>% filter(!is.na(sor))

sor_bird_extant_sf$f_sor = cut(sor_bird_extant_sf$sor,
                               breaks = quantile(sor_bird_extant_sf$sor,
                                                 probs = seq(0, 1,
                                                             length.out = 9),
                                                 na.rm = TRUE),
                               include.lowest = TRUE)

sor_bird_extant_sf = sor_bird_extant_sf %>% filter(RegionID != 37)
birds_extant_sor_map =
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
  geom_sf(data = (sor_bird_extant_sf %>% filter(!(Island == 1 & 
                                                    Area < 1e4))),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_bird_extant_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_bird_extant_sf %>% filter(Island == 1 & 
                                                     Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Extant_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_extant_sor_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(sor_bird_extant_sf$sor,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_extant_sor_map = ggplotGrob(birds_extant_sor_map)
legend.birds_extant_sor_map = ggplotGrob(legend.birds_extant_sor_map)
birds_extant_sor_map_all = arrangeGrob(birds_extant_sor_map,
                                       legend.birds_extant_sor_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(birds_extant_sor_map_all)





##3.4 Delta_sor (Extant - Native) #####
load('results/primary_results/distances_beta/phy_sor_bird_extant.rdata')
load('results/primary_results/distances_beta/phy_sor_bird_native.rdata')

### alternative 1. detla = mean(extant - native)
sor_bird_delta_mat = log((phy_sor_bird_extant+0.001) /
                           (phy_sor_bird_native+0.001))
sor_bird_delta = data.frame(RegionID = rownames(sor_bird_delta_mat),
                            delta_sor = rowMeans(sor_bird_delta_mat))
sor_bird_delta$RegionID = as.integer(sor_bird_delta$RegionID)

sor_bird_delta_sf = df_trans %>% left_join(sor_bird_delta,
                                           by = 'RegionID')

colnames(sor_bird_delta_sf)


#### Delta_ED birds mapping 
sor_bird_delta_sf = sor_bird_delta_sf %>% filter(!is.na(delta_sor) & 
                                                   delta_sor != 0 
)

hist(sor_bird_delta_sf$delta_sor)

breaks_delta_sor = quantile(sor_bird_delta_sf$delta_sor,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_sor[which(names(breaks_delta_sor) == '87.5%')] = 0
sor_bird_delta_sf$f_delta_sor = cut(sor_bird_delta_sf$delta_sor,
                                    breaks = breaks_delta_sor,
                                    include.lowest = TRUE)


birds_delta_sor_map =
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
  geom_sf(data = (sor_bird_delta_sf %>% filter(!(Island == 1 & 
                                                   Area < 1e4))), aes(fill = f_delta_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors8,
  # values = scales::rescale(quantile(sor_bird_delta_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_bird_delta_sf %>% filter(Island == 1 & 
                                                    Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_delta_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Delta_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_delta_sor_map = legend.func(mycolors = colors8,
                                         mylabels = round(breaks_delta_sor,
                                                          4)) +
  ggtitle("log (sor_all / sor_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_delta_sor_map = ggplotGrob(birds_delta_sor_map)
legend.birds_delta_sor_map = ggplotGrob(legend.birds_delta_sor_map)
birds_delta_sor_map_all = arrangeGrob(birds_delta_sor_map,
                                      legend.birds_delta_sor_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(birds_delta_sor_map_all)




##3.5 Delta_PU (Extant - Native) ----
load("results/primary_results/distances_beta/phy_turn_bird_extant.rdata")
load("results/primary_results/distances_beta/phy_turn_bird_native.rdata")

colnames(sor_bird_extant)[which(colnames(sor_bird_extant) == 'sor')] = 'all_sor'
PU_bird_delta = sor_bird_native %>% 
  left_join(sor_bird_extant, by = 'RegionID')
PU_bird_delta$delta_PU = log(PU_bird_delta$all_sor/
                               PU_bird_delta$sor)

PU_bird_delta_sf = df_trans %>% left_join(PU_bird_delta,
                                          by = 'RegionID')

head(PU_bird_delta %>% arrange(desc(delta_PU)) %>% 
       left_join(df_trans_c[,c('RegionID', 'Level_4_Na')], by = 'RegionID'))


#### Delta_ED birds mapping 
PU_bird_delta_sf = PU_bird_delta_sf %>% filter(!is.na(delta_PU) & 
                                                 delta_PU != 0 
)


hist(PU_bird_delta_sf$delta_PU)

breaks_delta_PU = quantile(PU_bird_delta_sf$delta_PU,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)


end_negative = sort(which(breaks_delta_PU < 0), decreasing = T)[1]
if (!is.na(end_negative)) {
  breaks_delta_PU[end_negative] = 0
}

if (end_negative == 6) {
  colors_my = colors5
} else if (end_negative == 7) {
  colors_my = colors7
} else if (end_negative == 8) {
  colors_my = colors8
}

PU_bird_delta_sf$f_delta_PU = cut(PU_bird_delta_sf$delta_PU,
                                  breaks = breaks_delta_PU,
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
PU_bird_delta_sf = PU_bird_delta_sf %>% filter(RegionID != 37)

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
  ggtitle('Delta_birds')+
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
plot(birds_delta_PU_map_all)







## 3.6 Patitioning delta ED into 5 possible ways ----
colnames(sp_dis_5)
colnames(sp_overlap_dat_1)
colnames(sp_overlap_dat_1)[which(colnames(sp_overlap_dat_1) == 'Inv. Stage')] = 'Inv._stage'

sp_dis_5$Inv._stage = 'native'
colnames(sp_dis_5)[which(colnames(sp_dis_5) %in% c("ScientificName", 
                                                   "Order.1.2",
                                                   "Family.1.2",
                                                   "Region.ID"))] = c("RegionID",
                                                                      "Binomial", 
                                                                      "Order",
                                                                      "Family")

sp_overlap_dat_2 = sp_overlap_dat_1 %>% 
  dplyr::select(which(colnames(sp_overlap_dat_1) %in% colnames(sp_dis_5)))
sp_dis_5 = as_tibble(sp_dis_5)
sp_dis_6 = sp_dis_5[,colnames(sp_overlap_dat_2)]

sp_dis_6 = sp_dis_6 %>% filter(!is.na(Binomial))

sp_bird_all = rbind(sp_dis_6,
                    sp_overlap_dat_2)

region_pairs = combn(colnames(sor_bird_delta_mat), 2)

sor_bird_path3_5_mat = sor_bird_native_mat
sor_bird_path4_7_mat = sor_bird_native_mat
sor_bird_path6_mat = sor_bird_native_mat

for (i in 1:ncol(region_pairs)) {
  #i = 1
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = sp_dis_6 %>% filter(RegionID == region1)
  nati_region2 = sp_dis_6 %>% filter(RegionID == region2)
  natu_region1 = sp_overlap_dat_2 %>% filter(RegionID == region1)
  natu_region2 = sp_overlap_dat_2 %>% filter(RegionID == region2)
  
  nati_region1_sps = nati_region1$Binomial
  nati_region2_sps = nati_region2$Binomial
  natu_region1_sps = setdiff(natu_region1$Binomial, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$Binomial, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    path6_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_2,
                                                              natu_2_sps_ori_1))]
    sp_bird_path6 = rbind(nati_region1, 
                          nati_region2,
                          natu_region1 %>% filter(Binomial %in% path6_sp))
    
    comm_bird_path6 = sp_bird_path6 %>% 
      complete(Binomial, RegionID, fill = list(presence = 0)) %>%
      dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
      pivot_wider(names_from = Binomial,
                  values_from = presence,
                  values_fn = mean) %>% 
      dplyr::arrange('RegionID') 
    
    phy_sor_bird_path6 = calcu_sor_simple(
      Tree = spec_phy.3,
      Comm = comm_bird_path6[,2:ncol(comm_bird_path6)])
    
    sor_bird_path6_mat[region1, region2] = phy_sor_bird_path6[1,2]
    sor_bird_path6_mat[region2, region1] = phy_sor_bird_path6[1,2]
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_bird_path3_5 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(Binomial %in% path3_5_sp))
      
      comm_bird_path3_5 = sp_bird_path3_5 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_bird_path3_5 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_bird_path3_5[,2:ncol(comm_bird_path3_5)])
      
      sor_bird_path3_5_mat[region1, region2] = phy_sor_bird_path3_5[1,2]
      sor_bird_path3_5_mat[region2, region1] = phy_sor_bird_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      path4_7_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_bird_path4_7 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(Binomial %in% path4_7_sp))
      
      comm_bird_path4_7 = sp_bird_path4_7 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_bird_path4_7 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_bird_path4_7[,2:ncol(comm_bird_path4_7)])
      
      sor_bird_path4_7_mat[region1, region2] = phy_sor_bird_path4_7[1,2]
      sor_bird_path4_7_mat[region2, region1] = phy_sor_bird_path4_7[1,2]
      
    }
    
  }
  
}


sor_bird_delta_path6_mat = log((sor_bird_path6_mat+0.001) /
                                 (sor_bird_native_mat+0.001))

sor_bird_delta_path4_7_mat = log((sor_bird_path4_7_mat+0.001) /
                                   (sor_bird_native_mat+0.001))
which(sor_bird_delta_path4_7_mat > 0)
which(sor_bird_delta_path4_7_mat < 0)

sor_bird_delta_path3_5_mat = log((sor_bird_path3_5_mat+0.001) /
                                   (sor_bird_native_mat+0.001))
which(sor_bird_delta_path3_5_mat > 0)
which(sor_bird_delta_path3_5_mat < 0)





#4. Fish sor: calculation & mapping ####
#load("D:/R projects/Global_ED/data/Fishes/data/my_phy.rdata")
#is.rooted(phylo)
load("D:/R projects/Global_ED/data/Fishes/data/my_data_used_final.rdata")
#basin = st_read("data/Fishes/data/Basin042017_3119.shp")
basins_trans = st_read("data/Fishes/data/Basin042017_3119_eck4.shp")
#save(df, file = "data/Fishes/data/Basin042017_3119.rdata")
colnames(basins_trans)[which(colnames(basins_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')

##4.3 natives #####
load('results/primary_results/distances_beta/phy_sor_fish_native.rdata')

sor_fish_native = data.frame(X1.Basin.Name = rownames(phy_sor_fish_native),
                             sor = rowMeans(phy_sor_fish_native, 
                                            na.rm = T))

sor_fish_native_sf = basins_trans %>% left_join(sor_fish_native,
                                                by = 'X1.Basin.Name')
colnames(sor_fish_native_sf)

#### native fishes mapping 
sor_fish_native_sf = sor_fish_native_sf %>% filter(!is.na(sor))

sor_fish_native_sf$f_sor = cut(sor_fish_native_sf$sor,
                               breaks = quantile(sor_fish_native_sf$sor,
                                                 probs = seq(0, 1,
                                                             length.out = 9),
                                                 na.rm = TRUE),
                               include.lowest = TRUE)


# 1. Get target CRS from polygon layer
target_crs = st_crs(sor_fish_native_sf)

# 2. Extract attributes of island points (Island == 1) and drop polygon geometry
island_fishes_native = sor_fish_native_sf %>%
  filter(Island == 1) %>%
  st_drop_geometry()   # now a plain data frame with Lon, Lat, f_sor, etc.

# 3. Create point sf from Lon/Lat (WGS84) and transform to target CRS
island_fishes_native_sf = island_fishes_native %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
  st_transform(crs = target_crs)

# 4. Extract projected coordinates (X, Y) as numeric columns
coords = st_coordinates(island_fishes_native_sf)   # matrix with columns X, Y
island_fishes_native_sf = island_fishes_native_sf %>%
  st_drop_geometry() %>%
  mutate(X_proj = coords[,1], Y_proj = coords[,2])

fishes_native_sor_map =
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
  geom_sf(data = (sor_fish_native_sf %>% filter(!(Island == 1 & 
                                                    Surf_area < 1e4))),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_fish_native_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = island_fishes_native_sf,
             aes(x = X_proj, y = Y_proj, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Native_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_native_sor_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(sor_fish_native_sf$sor,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_native_sor_map = ggplotGrob(fishes_native_sor_map)
legend.fishes_native_sor_map = ggplotGrob(legend.fishes_native_sor_map)
fishes_native_sor_map_all = arrangeGrob(fishes_native_sor_map,
                                        legend.fishes_native_sor_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(fishes_native_sor_map_all)



##4.4 extant #####
load('results/primary_results/distances_beta/phy_sor_fish_extant.rdata')

sor_fish_extant = data.frame(X1.Basin.Name = rownames(phy_sor_fish_extant),
                             sor = rowMeans(phy_sor_fish_extant, 
                                            na.rm = T))

sor_fish_extant_sf = basins_trans %>% left_join(sor_fish_extant,
                                                by = 'X1.Basin.Name')
colnames(sor_fish_extant_sf)

#### extant fishes mapping 
sor_fish_extant_sf = sor_fish_extant_sf %>% filter(!is.na(sor))

sor_fish_extant_sf$f_sor = cut(sor_fish_extant_sf$sor,
                               breaks = quantile(sor_fish_extant_sf$sor,
                                                 probs = seq(0, 1,
                                                             length.out = 9),
                                                 na.rm = TRUE),
                               include.lowest = TRUE)

# 2. Extract attributes of island points (Island == 1) and drop polygon geometry
island_fishes_extant = sor_fish_extant_sf %>%
  filter(Island == 1) %>%
  st_drop_geometry()   # now a plain data frame with Lon, Lat, f_sor, etc.

# 3. Create point sf from Lon/Lat (WGS84) and transform to target CRS
island_fishes_extant_sf = island_fishes_extant %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
  st_transform(crs = target_crs)

# 4. Extract projected coordinates (X, Y) as numeric columns
coords = st_coordinates(island_fishes_extant_sf)   # matrix with columns X, Y
island_fishes_extant_sf = island_fishes_extant_sf %>%
  st_drop_geometry() %>%
  mutate(X_proj = coords[,1], Y_proj = coords[,2])

fishes_extant_sor_map =
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
  geom_sf(data = (sor_fish_extant_sf %>% filter(!(Island == 1 & 
                                                    Surf_area < 1e4))),
          aes(fill = f_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(sor_fish_extant_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = island_fishes_extant_sf,
             aes(x = X_proj, y = Y_proj, color = f_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Extant_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_extant_sor_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(sor_fish_extant_sf$sor,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("PU: Sorensen's Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_extant_sor_map = ggplotGrob(fishes_extant_sor_map)
legend.fishes_extant_sor_map = ggplotGrob(legend.fishes_extant_sor_map)
fishes_extant_sor_map_all = arrangeGrob(fishes_extant_sor_map,
                                        legend.fishes_extant_sor_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(fishes_extant_sor_map_all)




##4.5 Delta_sor (Extant - Native) #####
load('results/primary_results/distances_beta/phy_sor_fish_extant.rdata')
load('results/primary_results/distances_beta/phy_sor_fish_native.rdata')

### alternative 1. detla = mean(extant - native)

common_basins = which(colnames(phy_sor_fish_extant) %in% intersect(colnames(phy_sor_fish_extant),
                                                                   colnames(phy_sor_fish_native)))
phy_sor_fish_extant = phy_sor_fish_extant[common_basins,common_basins]
sor_fish_delta_mat = log((phy_sor_fish_extant+0.001) /
                           (phy_sor_fish_native+0.001))
sor_fish_delta = data.frame(X1.Basin.Name = rownames(sor_fish_delta_mat),
                            delta_sor = rowMeans(sor_fish_delta_mat, 
                                                 na.rm = T))

sor_fish_delta_sf = basins_trans %>% left_join(sor_fish_delta,
                                               by = 'X1.Basin.Name')

colnames(sor_fish_delta_sf)



#### Delta_ED fishes mapping 
sor_fish_delta_sf = sor_fish_delta_sf %>% filter(!is.na(delta_sor) & 
                                                   delta_sor != 0 
)

hist(sor_fish_delta_sf$delta_sor)

breaks_delta_sor = quantile(sor_fish_delta_sf$delta_sor,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_sor[which(names(breaks_delta_sor) == '87.5%')] = 0
sor_fish_delta_sf$f_delta_sor = cut(sor_fish_delta_sf$delta_sor,
                                    breaks = breaks_delta_sor,
                                    include.lowest = TRUE)


fishes_delta_sor_map =
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
  geom_sf(data = (sor_fish_delta_sf %>% filter(!(Island == 1 & 
                                                   Surf_area < 1e4))),
          aes(fill = f_delta_sor),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor',
    guide = guide_legend(
      title = 'sor',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors8,
  # values = scales::rescale(quantile(sor_fish_delta_sf$sor,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'sor'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (sor_fish_delta_sf %>% filter(Island == 1 & 
                                                    Surf_area < 1e4)),
             aes(x = Lon, y = Lat, color = f_delta_sor),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'sor'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Delta_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_delta_sor_map = legend.func(mycolors = colors8,
                                          mylabels = round(breaks_delta_sor,
                                                           4)) +
  ggtitle("log (sor_all / sor_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_delta_sor_map = ggplotGrob(fishes_delta_sor_map)
legend.fishes_delta_sor_map = ggplotGrob(legend.fishes_delta_sor_map)
fishes_delta_sor_map_all = arrangeGrob(fishes_delta_sor_map,
                                       legend.fishes_delta_sor_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(fishes_delta_sor_map_all)



##4.6 Delta_PU (Extant - Native) ----
load("results/primary_results/distances_beta/phy_turn_fish_extant.rdata")
load("results/primary_results/distances_beta/phy_turn_fish_native.rdata")

colnames(sor_fish_extant)[which(colnames(sor_fish_extant) == 'sor')] = 'all_sor'
PU_fish_delta = sor_fish_native %>% 
  left_join(sor_fish_extant, by = 'X1.Basin.Name')
PU_fish_delta$delta_PU = log(PU_fish_delta$all_sor/
                               PU_fish_delta$sor)

PU_fish_delta_sf = basins_trans %>% left_join(PU_fish_delta,
                                              by = 'X1.Basin.Name')
colnames(PU_fish_delta_sf)
head(PU_fish_delta %>% arrange(desc(delta_PU)) %>% 
       left_join(basins_trans[,c('X1.Basin.Name', 'Country')], by = 'X1.Basin.Name'))


#### Delta_ED fishes mapping 
PU_fish_delta_sf = PU_fish_delta_sf %>% filter(!is.na(delta_PU) & 
                                                 delta_PU != 0 
)


hist(PU_fish_delta_sf$delta_PU)

breaks_delta_PU = quantile(PU_fish_delta_sf$delta_PU,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)

end_negative = sort(which(breaks_delta_PU < 0), decreasing = T)[1]
if (!is.na(end_negative)) {
  breaks_delta_PU[end_negative] = 0
}

if (end_negative == 6) {
  colors_my = colors5
} else if (end_negative == 7) {
  colors_my = colors7
} else if (end_negative == 8) {
  colors_my = colors8
}

PU_fish_delta_sf$f_delta_PU = cut(PU_fish_delta_sf$delta_PU,
                                  breaks = breaks_delta_PU,
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer


# 2. Extract attributes of island points (Island == 1) and drop polygon geometry
island_fishes_delta = PU_fish_delta_sf %>%
  filter(Island == 1) %>%
  st_drop_geometry()   # now a plain data frame with Lon, Lat, f_sor, etc.

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
  ggtitle('Delta_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_delta_PU_map = legend.func(mycolors = colors8,
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
plot(fishes_delta_PU_map_all)



## 4.7 Patitioning delta ED into 5 possible ways ----
colnames(sp_dis_5)
colnames(sp_overlap_dat_1)
colnames(sp_overlap_dat_1)[which(colnames(sp_overlap_dat_1) == 'Inv. Stage')] = 'Inv._stage'

sp_dis_5$Inv._stage = 'native'
colnames(sp_dis_5)[which(colnames(sp_dis_5) %in% c("ScientificName", 
                                                   "Order.1.2",
                                                   "Family.1.2",
                                                   "Region.ID"))] = c("RegionID",
                                                                      "Binomial", 
                                                                      "Order",
                                                                      "Family")

sp_overlap_dat_2 = sp_overlap_dat_1 %>% 
  dplyr::select(which(colnames(sp_overlap_dat_1) %in% colnames(sp_dis_5)))
sp_dis_5 = as_tibble(sp_dis_5)
sp_dis_6 = sp_dis_5[,colnames(sp_overlap_dat_2)]

sp_dis_6 = sp_dis_6 %>% filter(!is.na(Binomial))

sp_fish_all = rbind(sp_dis_6,
                    sp_overlap_dat_2)

region_pairs = combn(colnames(sor_fish_delta_mat), 2)

sor_fish_path3_5_mat = sor_fish_native_mat
sor_fish_path4_7_mat = sor_fish_native_mat
sor_fish_path6_mat = sor_fish_native_mat

for (i in 1:ncol(region_pairs)) {
  #i = 1
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = sp_dis_6 %>% filter(RegionID == region1)
  nati_region2 = sp_dis_6 %>% filter(RegionID == region2)
  natu_region1 = sp_overlap_dat_2 %>% filter(RegionID == region1)
  natu_region2 = sp_overlap_dat_2 %>% filter(RegionID == region2)
  
  nati_region1_sps = nati_region1$Binomial
  nati_region2_sps = nati_region2$Binomial
  natu_region1_sps = setdiff(natu_region1$Binomial, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$Binomial, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    path6_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_2,
                                                              natu_2_sps_ori_1))]
    sp_fish_path6 = rbind(nati_region1, 
                          nati_region2,
                          natu_region1 %>% filter(Binomial %in% path6_sp))
    
    comm_fish_path6 = sp_fish_path6 %>% 
      complete(Binomial, RegionID, fill = list(presence = 0)) %>%
      dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
      pivot_wider(names_from = Binomial,
                  values_from = presence,
                  values_fn = mean) %>% 
      dplyr::arrange('RegionID') 
    
    phy_sor_fish_path6 = calcu_sor_simple(
      Tree = spec_phy.3,
      Comm = comm_fish_path6[,2:ncol(comm_fish_path6)])
    
    sor_fish_path6_mat[region1, region2] = phy_sor_fish_path6[1,2]
    sor_fish_path6_mat[region2, region1] = phy_sor_fish_path6[1,2]
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_fish_path3_5 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(Binomial %in% path3_5_sp))
      
      comm_fish_path3_5 = sp_fish_path3_5 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_fish_path3_5 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_fish_path3_5[,2:ncol(comm_fish_path3_5)])
      
      sor_fish_path3_5_mat[region1, region2] = phy_sor_fish_path3_5[1,2]
      sor_fish_path3_5_mat[region2, region1] = phy_sor_fish_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      path4_7_sp = natu_region1_sps[which(natu_region1_sps %in% c(natu_1_sps_ori_3,
                                                                  natu_2_sps_ori_3))]
      sp_fish_path4_7 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(Binomial %in% path4_7_sp))
      
      comm_fish_path4_7 = sp_fish_path4_7 %>% 
        complete(Binomial, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
        pivot_wider(names_from = Binomial,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_sor_fish_path4_7 = calcu_sor_simple(
        Tree = spec_phy.3,
        Comm = comm_fish_path4_7[,2:ncol(comm_fish_path4_7)])
      
      sor_fish_path4_7_mat[region1, region2] = phy_sor_fish_path4_7[1,2]
      sor_fish_path4_7_mat[region2, region1] = phy_sor_fish_path4_7[1,2]
      
    }
    
  }
  
}


sor_fish_delta_path6_mat = log((sor_fish_path6_mat+0.001) /
                                 (sor_fish_native_mat+0.001))

sor_fish_delta_path4_7_mat = log((sor_fish_path4_7_mat+0.001) /
                                   (sor_fish_native_mat+0.001))
which(sor_fish_delta_path4_7_mat > 0)
which(sor_fish_delta_path4_7_mat < 0)

sor_fish_delta_path3_5_mat = log((sor_fish_path3_5_mat+0.001) /
                                   (sor_fish_native_mat+0.001))
which(sor_fish_delta_path3_5_mat > 0)
which(sor_fish_delta_path3_5_mat < 0)




#5. Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

##5.1 sor #####
# Compare the sor patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.01         # small horizontal gap between columns

figs_sor_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_sor_map_all,  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_sor_map_all,  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_sor_map_all,   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_sor_map_all,   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_extant_sor_map_all,   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_sor_map_all,    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_sor_map_all, x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_all_sor_map_all, x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_sor_map_all,  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_sor_map_all,   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_sor_map_all,   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_sor_map_all,   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
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

png(filename = 'figures/figs_sor_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_sor_nat_extant_delta
dev.off() #turn off device and finalize file


### export for PPT

figs_sor_delta = ggarrange(plants_delta_sor_map_all, 
                           birds_delta_sor_map_all,
                           mammals_delta_sor_map_all,
                           fishes_delta_sor_map_all,
                           nrow = 2, ncol = 2, 
                           labels = c('a', 'b', 'c', 'd'))


png(filename = 'figures/figs_sor_delta.png',
    height=15, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_sor_delta
dev.off() #turn off device and finalize file


##5.2 PU #####
# Compare the sor patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.01         # small horizontal gap between columns

figs_PU_sor_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_sor_map_all,  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_sor_map_all,  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_PU_map_all,   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_sor_map_all,   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_extant_sor_map_all,   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_PU_map_all,    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_sor_map_all, x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_all_sor_map_all, x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_PU_map_all,  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_sor_map_all,   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_sor_map_all,   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_PU_map_all,   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
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

png(filename = 'figures/figs_PU_sor_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_PU_sor_nat_extant_delta
dev.off() #turn off device and finalize file


