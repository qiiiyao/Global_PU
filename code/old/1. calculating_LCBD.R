### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'scico', 'ggplot2', 'gridExtra')

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


#### Mammal LCBD: calculation & mapping ####
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")
is.rooted(spec_phy.3)

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)
sp_overlap_dat_1$presence = 1

comm_mammal_exotic = sp_overlap_dat_1 %>% 
  #complete(Binomial, RegionID, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(names_from = Binomial,
              values_from = presence,
              values_fn = mean) %>% 
  #complete(RegionID = unique(df$RegionID)) %>% ## assume that 
  ## the absence regions have no naturalized aliens
  dplyr::mutate(
    across(everything(), \(x) {
      if (is.factor(x)) {
        fct_expand(x, "0") %>%  tidyr::replace_na("0")
      } else if (is.character(x)) {
        tidyr::replace_na(x, "0")
      } else if (inherits(x, "units")) {
        unit = units(x)
        tidyr::replace_na(x, set_units(0, unit, mode = "standard"))
      } else {
        tidyr::replace_na(x, 0)
      }
    })
  )  %>% left_join(df_trans[,c('RegionID', 'area')],
                   by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  dplyr::arrange(RegionID)

LCBD_mammal_exotic = calcu_LCBD_parallel(
  Tree = spec_phy.3,
  Comm = comm_mammal_exotic,
  Region_posi = which(colnames(comm_mammal_exotic) == 'RegionID'))

save(LCBD_mammal_exotic,
     file = 'results/primary_results/LCBD_mammal_exotic.rdata')

load('results/primary_results/LCBD_mammal_exotic.rdata')

LCBD_mammal_exotic_sf = df_trans %>% left_join(LCBD_mammal_exotic$LCBD_simp_geo,
                                               by = 'RegionID')
colnames(LCBD_mammal_exotic_sf)

#### exotic mammals mapping 
LCBD_mammal_exotic_sf = LCBD_mammal_exotic_sf %>% filter(!is.na(LCBD))

LCBD_mammal_exotic_sf$f_LCBD = cut(LCBD_mammal_exotic_sf$LCBD,
                                   breaks = quantile(LCBD_mammal_exotic_sf$LCBD,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_mammal_exotic_sf_2 = LCBD_mammal_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


SSI_mammal_exotic_sf = LCBD_mammal_exotic_sf %>% filter(!is.na(SSI))

SSI_mammal_exotic_sf$f_delta_SSI = cut(SSI_mammal_exotic_sf$SSI,
                                 breaks = quantile(SSI_mammal_exotic_sf$SSI,
                                                   probs = seq(0, 1,
                                                               length.out = 9),
                                                   na.rm = TRUE),
                                 include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_mammal_exotic_sf_2 = SSI_mammal_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


mammals_exotic_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_mammal_exotic_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_mammal_exotic_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_mammal_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_exotic_LCBD_map = legend.func(mycolors = colors2,
                                             mylabels = round(quantile(LCBD_mammal_exotic_sf$LCBD,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_exotic_LCBD_map = ggplotGrob(mammals_exotic_LCBD_map)
legend.mammals_exotic_LCBD_map = ggplotGrob(legend.mammals_exotic_LCBD_map)
mammals_exotic_LCBD_map_all = arrangeGrob(mammals_exotic_LCBD_map,
                                          legend.mammals_exotic_LCBD_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(mammals_exotic_LCBD_map_all)




mammals_exotic_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_mammal_exotic_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_mammal_exotic_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_mammal_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_exotic_SSI_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(SSI_mammal_exotic_sf$SSI,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_exotic_SSI_map = ggplotGrob(mammals_exotic_SSI_map)
legend.mammals_exotic_SSI_map = ggplotGrob(legend.mammals_exotic_SSI_map)
mammals_exotic_SSI_map_all = arrangeGrob(mammals_exotic_SSI_map,
                                         legend.mammals_exotic_SSI_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(mammals_exotic_SSI_map_all)




##### natives #####
colnames(sp_dis_5)
sp_dis_5$presence = 1
comm_mammal_native = sp_dis_5 %>% 
  complete(Region.ID, ScientificName, fill = list(presence = 0)) %>%
  dplyr::select(c('Region.ID', 'ScientificName', 'presence')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = presence) %>% 
  left_join(df_trans[,c('RegionID', 'area')],
            by = join_by('Region.ID' == 'RegionID')) %>% 
  relocate(area, .after = Region.ID) %>% 
  dplyr::arrange('Region.ID') 

LCBD_mammal_native = calcu_LCBD_parallel(Tree = spec_phy.3,
                                         Comm = comm_mammal_native,
                                         Region_posi = which(colnames(comm_mammal_native) == 'Region.ID'))

save(LCBD_mammal_native,
     file = 'results/primary_results/LCBD_mammal_native.rdata')

load('results/primary_results/LCBD_mammal_native.rdata')

LCBD_mammal_native_sf = df_trans %>% left_join(LCBD_mammal_native$LCBD_simp_geo, by = 'RegionID')
colnames(LCBD_mammal_native_sf)

#### native mammals mapping 
LCBD_mammal_native_sf = LCBD_mammal_native_sf %>% filter(!is.na(LCBD))

LCBD_mammal_native_sf$f_LCBD = cut(LCBD_mammal_native_sf$LCBD,
                                    breaks = quantile(LCBD_mammal_native_sf$LCBD,
                                                      probs = seq(0, 1,
                                                                  length.out = 9),
                                                      na.rm = TRUE),
                                    include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_mammal_native_sf_2 = LCBD_mammal_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


SSI_mammal_native_sf = LCBD_mammal_native_sf %>% filter(!is.na(SSI))

SSI_mammal_native_sf$f_delta_SSI = cut(SSI_mammal_native_sf$SSI,
                                   breaks = quantile(SSI_mammal_native_sf$SSI,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_mammal_native_sf_2 = SSI_mammal_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


mammals_native_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_mammal_native_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_mammal_native_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_mammal_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.mammals_native_LCBD_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(LCBD_mammal_native_sf$LCBD,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_LCBD_map = ggplotGrob(mammals_native_LCBD_map)
legend.mammals_native_LCBD_map = ggplotGrob(legend.mammals_native_LCBD_map)
mammals_native_LCBD_map_all = arrangeGrob(mammals_native_LCBD_map,
                                        legend.mammals_native_LCBD_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(mammals_native_LCBD_map_all)




mammals_native_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_mammal_native_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_mammal_native_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_mammal_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.mammals_native_SSI_map = legend.func(mycolors = colors2,
                                             mylabels = round(quantile(SSI_mammal_native_sf$SSI,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_SSI_map = ggplotGrob(mammals_native_SSI_map)
legend.mammals_native_SSI_map = ggplotGrob(legend.mammals_native_SSI_map)
mammals_native_SSI_map_all = arrangeGrob(mammals_native_SSI_map,
                                          legend.mammals_native_SSI_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(mammals_native_SSI_map_all)


##### extant regional assemlage: natives + exotics #####
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

sp_mammal_all = rbind(sp_dis_6,
                      sp_overlap_dat_2)

unique(sp_mammal_all$Inv._stage)
colnames(sp_mammal_all)

#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_mammal_all$Binomial = gsub(' ', '_', sp_mammal_all$Binomial)
sp_mammal_all$presence = 1

comm_mammal_all = sp_mammal_all %>% 
  complete(Binomial, RegionID, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(names_from = Binomial,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_trans[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  dplyr::arrange('RegionID') 

LCBD_mammal_all = calcu_LCBD_parallel(
  Tree = spec_phy.3,
  Comm = comm_mammal_all,
  Region_posi = which(colnames(comm_mammal_all) == 'RegionID'))

save(LCBD_mammal_all,
     file = 'results/primary_results/LCBD_mammal_all.rdata')

load('results/primary_results/LCBD_mammal_all.rdata')

LCBD_mammal_extant_sf = df_trans %>% left_join(LCBD_mammal_all$LCBD_simp_geo, by = 'RegionID')
colnames(LCBD_mammal_extant_sf)

#### all mammals mapping 
LCBD_mammal_extant_sf = LCBD_mammal_extant_sf %>% filter(!is.na(LCBD))

LCBD_mammal_extant_sf$f_LCBD = cut(LCBD_mammal_extant_sf$LCBD,
                                   breaks = quantile(LCBD_mammal_extant_sf$LCBD,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_mammal_extant_sf_2 = LCBD_mammal_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



SSI_mammal_extant_sf = LCBD_mammal_extant_sf %>% filter(!is.na(SSI))

SSI_mammal_extant_sf$f_delta_SSI = cut(SSI_mammal_extant_sf$SSI,
                                   breaks = quantile(SSI_mammal_extant_sf$SSI,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_mammal_extant_sf_2 = SSI_mammal_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



mammals_extant_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_mammal_extant_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_mammal_extant_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_mammal_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.mammals_extant_LCBD_map = legend.func(mycolors = colors2,
                                             mylabels = round(quantile(LCBD_mammal_extant_sf$LCBD,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_extant_LCBD_map = ggplotGrob(mammals_extant_LCBD_map)
legend.mammals_extant_LCBD_map = ggplotGrob(legend.mammals_extant_LCBD_map)
mammals_extant_LCBD_map_all = arrangeGrob(mammals_extant_LCBD_map,
                                          legend.mammals_extant_LCBD_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(mammals_extant_LCBD_map_all)



mammals_extant_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_mammal_extant_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_mammal_extant_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_mammal_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.mammals_extant_SSI_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(SSI_mammal_extant_sf$SSI,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_extant_SSI_map = ggplotGrob(mammals_extant_SSI_map)
legend.mammals_extant_SSI_map = ggplotGrob(legend.mammals_extant_SSI_map)
mammals_extant_SSI_map_all = arrangeGrob(mammals_extant_SSI_map,
                                          legend.mammals_extant_SSI_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(mammals_extant_SSI_map_all)




##### Delta_LCBD (Extant - Native) #####
load('results/primary_results/LCBD_mammal_all.rdata')
load('results/primary_results/LCBD_mammal_native.rdata')

LCBD_mammal_extant_1 = LCBD_mammal_all$LCBD_simp_geo
colnames(LCBD_mammal_extant_1)[2:3] = c("all_LCBD", "all_SSI")

LCBD_mammal_native_1 = LCBD_mammal_native$LCBD_simp_geo
colnames(LCBD_mammal_native_1)[2:3] = c("native_LCBD", "native_SSI")

LCBD_mammal_delta = LCBD_mammal_native_1 %>% left_join(LCBD_mammal_extant_1,
                                                   by = 'RegionID')

LCBD_mammal_delta$delta_LCBD = log(LCBD_mammal_delta$all_LCBD / 
                                      LCBD_mammal_delta$native_LCBD)

LCBD_mammal_delta$delta_SSI = log(LCBD_mammal_delta$all_SSI / 
                                       LCBD_mammal_delta$native_SSI)

LCBD_mammal_delta_sf = df_trans %>% left_join(LCBD_mammal_delta, by = 'RegionID')

colnames(LCBD_mammal_delta_sf)


#### Delta_ED mammals mapping 
LCBD_mammal_delta_sf = LCBD_mammal_delta_sf %>% filter(!is.na(delta_LCBD) & 
                                                         delta_LCBD != 0 
)

LCBD_mammal_delta_sf = LCBD_mammal_delta_sf %>% filter(!is.na(delta_SSI) & 
                                                         delta_SSI != 0 
)

hist(LCBD_mammal_delta_sf$delta_LCBD)

breaks_delta_LCBD = quantile(LCBD_mammal_delta_sf$delta_LCBD,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)
breaks_delta_LCBD[which(names(breaks_delta_LCBD) == '50%')] = 0
LCBD_mammal_delta_sf$f_delta_LCBD = cut(LCBD_mammal_delta_sf$delta_LCBD,
                                         breaks = breaks_delta_LCBD,
                                         include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks_delta_SSI = quantile(LCBD_mammal_delta_sf$delta_SSI,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 4)
breaks_delta_SSI[which(names(breaks_delta_SSI) == '62.5%')] = 0
LCBD_mammal_delta_sf$f_delta_SSI = cut(LCBD_mammal_delta_sf$delta_SSI,
                                        breaks = breaks_delta_SSI,
                                        include.lowest = TRUE)

hist(LCBD_mammal_delta_sf$delta_SSI)

LCBD_mammal_delta_sf_2 = LCBD_mammal_delta_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])




mammals_delta_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_mammal_delta_sf, aes(fill = f_delta_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors4,
  # values = scales::rescale(quantile(LCBD_mammal_delta_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_mammal_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_LCBD, color = f_delta_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.mammals_delta_LCBD_map = legend.func(mycolors = colors4,
                                             mylabels = round(breaks_delta_LCBD,
                                                              4)) +
  ggtitle("log (LCBD_all / LCBD_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_delta_LCBD_map = ggplotGrob(mammals_delta_LCBD_map)
legend.mammals_delta_LCBD_map = ggplotGrob(legend.mammals_delta_LCBD_map)
mammals_delta_LCBD_map_all = arrangeGrob(mammals_delta_LCBD_map,
                                          legend.mammals_delta_LCBD_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(mammals_delta_LCBD_map_all)



mammals_delta_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_mammal_delta_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors5,
  # values = scales::rescale(quantile(SSI_mammal_delta_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_mammal_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'Delta_SSI'
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

legend.mammals_delta_SSI_map = legend.func(mycolors = colors5,
                                            mylabels =  round(breaks_delta_SSI,
                                                              4)) +
  ggtitle("log (SSI_all / SSI_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_delta_SSI_map = ggplotGrob(mammals_delta_SSI_map)
legend.mammals_delta_SSI_map = ggplotGrob(legend.mammals_delta_SSI_map)
mammals_delta_SSI_map_all = arrangeGrob(mammals_delta_SSI_map,
                                         legend.mammals_delta_SSI_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(mammals_delta_SSI_map_all)





#### Plant LCBD: calculation & mapping ####
load("D:/R projects/Global_ED/data/Plants/data/shp.651.Rdata")
phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)

colnames(df.native.650) == colnames(df.natu.650)



##### 1.3 Remove apomictic species #####
# Download latest version of the Apomixis Database:
# (i) click on "Search" in query box. (ii) scroll down in species list. 
# (iii) click "Export Report as CSV" at bottom of species list. 
# https://uni-goettingen.de/en/433689.html

Apomixis = read.table("data/Plants/data/Apomixis_data.txt",
                      header = T, sep = "")

Apomixis$Genus = as.character(Apomixis$Genus)
Apomixis = Apomixis%>%filter(Apomixis.Yes.Uncertain.=="Y")
#Genus with parenthesis
Apomixis$Genus1 =  sapply(Apomixis$Genus, function(x)strsplit(x, split = ' [()]')[[1]][1])
Apomixis$Genus2 =  sapply(Apomixis$Genus, function(x)strsplit(x, split = '[()]')[[1]][2])
#get all genera in one vector
Apomixis_Genus =  Apomixis$Genus[-grep("\\s*\\([^\\)]+\\)",Apomixis$Genus)]
Apomixis_Genus = unique(c(Apomixis_Genus,Apomixis$Genus1[which(!is.na(Apomixis$Genus1))]))
Apomixis_Genus = unique(c(Apomixis_Genus,Apomixis$Genus2[which(!is.na(Apomixis$Genus2))]))
Apomixis$Genus2[which(Apomixis$Genus2%in%Apomixis$Genus)]#3 genera

#remove genus containing apomictic species
length(unique(df.native.650$species))
df.native.650 = df.native.650 %>% filter(!(genus%in%Apomixis_Genus))
length(unique(df.native.650$species))

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
colnames(df.natu.650)
df.natu.650$presence = 1
phylo_plant_exotic = drop.tip(phylo_big, setdiff(phylo_big$tip.label, df.natu.650$species))

comm_plant_exotic = df.natu.650 %>% 
  complete(species, Region_id, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  dplyr::arrange('Region_id') 

LCBD_plant_exotic = calcu_LCBD_parallel(Tree = phylo_plant_exotic,
                                    Comm = comm_plant_exotic,
                                    Region_posi = which(colnames(comm_plant_exotic) == 'Region_id'))

save(LCBD_plant_exotic,
     file = 'results/primary_results/LCBD_plant_exotic.rdata')

load("results/primary_results/LCBD_plant_exotic.rdata")

LCBD_plant_exotic_sf = shp.glonaf.trans %>% 
  left_join(LCBD_plant_exotic$LCBD_simp_geo,
            by = join_by('Region_id' == 'RegionID'))

colnames(LCBD_plant_exotic_sf)

#### exotic plants mapping 
LCBD_plant_exotic_sf = LCBD_plant_exotic_sf %>% filter(!is.na(LCBD))

LCBD_plant_exotic_sf$f_LCBD = cut(LCBD_plant_exotic_sf$LCBD,
                                   breaks = quantile(LCBD_plant_exotic_sf$LCBD,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_plant_exotic_sf_2 = LCBD_plant_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


SSI_plant_exotic_sf = LCBD_plant_exotic_sf %>% filter(!is.na(SSI))

SSI_plant_exotic_sf$f_delta_SSI = cut(SSI_plant_exotic_sf$SSI,
                                       breaks = quantile(SSI_plant_exotic_sf$SSI,
                                                         probs = seq(0, 1,
                                                                     length.out = 9),
                                                         na.rm = TRUE),
                                       include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_plant_exotic_sf_2 = SSI_plant_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


plants_exotic_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_plant_exotic_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_plant_exotic_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_plant_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_exotic_LCBD_map = legend.func(mycolors = colors2,
                                             mylabels = round(quantile(LCBD_plant_exotic_sf$LCBD,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_LCBD_map = ggplotGrob(plants_exotic_LCBD_map)
legend.plants_exotic_LCBD_map = ggplotGrob(legend.plants_exotic_LCBD_map)
plants_exotic_LCBD_map_all = arrangeGrob(plants_exotic_LCBD_map,
                                          legend.plants_exotic_LCBD_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(plants_exotic_LCBD_map_all)


plants_exotic_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_plant_exotic_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_plant_exotic_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_plant_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_exotic_SSI_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(SSI_plant_exotic_sf$SSI,
                                                           probs = seq(0, 1,length.out = 9),
                                                           na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_SSI_map = ggplotGrob(plants_exotic_SSI_map)
legend.plants_exotic_SSI_map = ggplotGrob(legend.plants_exotic_SSI_map)
plants_exotic_SSI_map_all = arrangeGrob(plants_exotic_SSI_map,
                                         legend.plants_exotic_SSI_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(plants_exotic_SSI_map_all)




##### natives #####
colnames(df.native.650)
df.native.650$presence = 1
phylo_plant_native = drop.tip(phylo_big, setdiff(phylo_big$tip.label, df.native.650$species))
comm_plant_native = df.native.650 %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

LCBD_plant_native = calcu_LCBD_parallel(Tree = phylo_plant_native,
                                        Comm = comm_plant_native,
                                        Region_posi = which(colnames(comm_plant_native) == 'Region_id'))

#save(LCBD_plant_native,
#     file = 'results/primary_results/LCBD_plant_native.rdata') # run in the server

load("results/primary_results/LCBD_plant_native.rdata") # load the calculated data from the server

LCBD_plant_native_sf = shp.glonaf.trans %>% left_join(LCBD_plant_native$LCBD_simp_geo,
                                                      by = join_by('Region_id' == 'RegionID'))
colnames(LCBD_plant_native_sf)

#### native plants mapping 
LCBD_plant_native_sf = LCBD_plant_native_sf %>% filter(!is.na(LCBD))

LCBD_plant_native_sf$f_LCBD = cut(LCBD_plant_native_sf$LCBD,
                                   breaks = quantile(LCBD_plant_native_sf$LCBD,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_plant_native_sf_2 = LCBD_plant_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


SSI_plant_native_sf = LCBD_plant_native_sf %>% filter(!is.na(SSI))

SSI_plant_native_sf$f_delta_SSI = cut(SSI_plant_native_sf$SSI,
                                       breaks = quantile(SSI_plant_native_sf$SSI,
                                                         probs = seq(0, 1,
                                                                     length.out = 9),
                                                         na.rm = TRUE),
                                       include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_plant_native_sf_2 = SSI_plant_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


plants_native_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_plant_native_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_plant_native_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_plant_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.plants_native_LCBD_map = legend.func(mycolors = colors2,
                                             mylabels = round(quantile(LCBD_plant_native_sf$LCBD,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_LCBD_map = ggplotGrob(plants_native_LCBD_map)
legend.plants_native_LCBD_map = ggplotGrob(legend.plants_native_LCBD_map)
plants_native_LCBD_map_all = arrangeGrob(plants_native_LCBD_map,
                                          legend.plants_native_LCBD_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(plants_native_LCBD_map_all)




plants_native_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_plant_native_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_plant_native_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_plant_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.plants_native_SSI_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(SSI_plant_native_sf$SSI,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_SSI_map = ggplotGrob(plants_native_SSI_map)
legend.plants_native_SSI_map = ggplotGrob(legend.plants_native_SSI_map)
plants_native_SSI_map_all = arrangeGrob(plants_native_SSI_map,
                                         legend.plants_native_SSI_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(plants_native_SSI_map_all)




##### extant #####
df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1
comm_plant_extant = df.extant.650 %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

LCBD_plant_extant = calcu_LCBD_parallel(Tree = phylo,
                                        Comm = comm_plant_extant,
                                        Region_posi = which(colnames(comm_plant_extant) == 'Region_id'))
save(LCBD_plant_extant,
     file = 'results/primary_results/LCBD_plant_extant.rdata')

load("results/primary_results/LCBD_plant_extant.rdata") # load the calculated data from the server

LCBD_plant_extant_sf = shp.glonaf.trans %>% left_join(LCBD_plant_extant$LCBD_simp_geo,
                                                      by = join_by('Region_id' == 'RegionID'))
colnames(LCBD_plant_extant_sf)
hist(LCBD_plant_extant_sf$LCBD)

#### extant plants mapping 
LCBD_plant_extant_sf = LCBD_plant_extant_sf %>% filter(!is.na(LCBD))

LCBD_plant_extant_sf$f_LCBD = cut(LCBD_plant_extant_sf$LCBD,
                                   breaks = quantile(LCBD_plant_extant_sf$LCBD,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_plant_extant_sf_2 = LCBD_plant_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



SSI_plant_extant_sf = LCBD_plant_extant_sf %>% filter(!is.na(SSI))

SSI_plant_extant_sf$f_delta_SSI = cut(SSI_plant_extant_sf$SSI,
                                       breaks = quantile(SSI_plant_extant_sf$SSI,
                                                         probs = seq(0, 1,
                                                                     length.out = 9),
                                                         na.rm = TRUE),
                                       include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_plant_extant_sf_2 = SSI_plant_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



plants_extant_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_plant_extant_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_plant_extant_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_plant_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.plants_extant_LCBD_map = legend.func(mycolors = colors2,
                                             mylabels = round(quantile(LCBD_plant_extant_sf$LCBD,
                                                                       probs = seq(0, 1,length.out = 9),
                                                                       na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_LCBD_map = ggplotGrob(plants_extant_LCBD_map)
legend.plants_extant_LCBD_map = ggplotGrob(legend.plants_extant_LCBD_map)
plants_extant_LCBD_map_all = arrangeGrob(plants_extant_LCBD_map,
                                          legend.plants_extant_LCBD_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10),
                                                                c(NA, rep(2, 8), NA)))
plot(plants_extant_LCBD_map_all)



plants_extant_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_plant_extant_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_plant_extant_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_plant_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.plants_extant_SSI_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(SSI_plant_extant_sf$SSI,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_SSI_map = ggplotGrob(plants_extant_SSI_map)
legend.plants_extant_SSI_map = ggplotGrob(legend.plants_extant_SSI_map)
plants_extant_SSI_map_all = arrangeGrob(plants_extant_SSI_map,
                                         legend.plants_extant_SSI_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(plants_extant_SSI_map_all)


##### Delta_ED (Extant - Native) #####
load('results/primary_results/LCBD_plant_extant.rdata')
load('results/primary_results/LCBD_plant_native.rdata')

LCBD_plant_extant_1 = LCBD_plant_extant$LCBD_simp_geo
colnames(LCBD_plant_extant_1)[2:3] = c("all_LCBD", "all_SSI")

LCBD_plant_native_1 = LCBD_plant_native$LCBD_simp_geo
colnames(LCBD_plant_native_1)[2:3] = c("native_LCBD", "native_SSI")

LCBD_plant_delta = LCBD_plant_native_1 %>% left_join(LCBD_plant_extant_1,
                                                       by = 'RegionID')

LCBD_plant_delta$delta_LCBD = log(LCBD_plant_delta$all_LCBD / 
                                     LCBD_plant_delta$native_LCBD)

LCBD_plant_delta$delta_SSI = log(LCBD_plant_delta$all_SSI / 
                                    LCBD_plant_delta$native_SSI)

LCBD_plant_delta_sf = shp.glonaf.trans %>% left_join(LCBD_plant_delta,
                                                     by = join_by('Region_id' == 'RegionID'))

colnames(LCBD_plant_delta_sf)


#### Delta_ED plants mapping 
LCBD_plant_delta_sf = LCBD_plant_delta_sf %>% filter(!is.na(delta_LCBD) & 
                                                         delta_LCBD != 0 
)

LCBD_plant_delta_sf = LCBD_plant_delta_sf %>% filter(!is.na(delta_SSI) & 
                                                         delta_SSI != 0 
)

hist(LCBD_plant_delta_sf$delta_LCBD)

breaks_delta_LCBD = quantile(LCBD_plant_delta_sf$delta_LCBD,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 12)
breaks_delta_LCBD[which(names(breaks_delta_LCBD) == '50%')] = 0
LCBD_plant_delta_sf$f_delta_LCBD = cut(LCBD_plant_delta_sf$delta_LCBD,
                                        breaks = breaks_delta_LCBD,
                                        include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks_delta_SSI = quantile(LCBD_plant_delta_sf$delta_SSI,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 4)
#breaks_delta_SSI[which(names(breaks_delta_SSI) == '62.5%')] = 0
LCBD_plant_delta_sf$f_delta_SSI = cut(LCBD_plant_delta_sf$delta_SSI,
                                       breaks = breaks_delta_SSI,
                                       include.lowest = TRUE)

hist(LCBD_plant_delta_sf$delta_SSI)

LCBD_plant_delta_sf_2 = LCBD_plant_delta_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



plants_delta_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_plant_delta_sf, aes(fill = f_delta_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors4,
  # values = scales::rescale(quantile(LCBD_plant_delta_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_plant_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_LCBD, color = f_delta_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.plants_delta_LCBD_map = legend.func(mycolors = colors4,
                                            mylabels = round(breaks_delta_LCBD,
                                                             4)) +
  ggtitle("log (LCBD_all / LCBD_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_LCBD_map = ggplotGrob(plants_delta_LCBD_map)
legend.plants_delta_LCBD_map = ggplotGrob(legend.plants_delta_LCBD_map)
plants_delta_LCBD_map_all = arrangeGrob(plants_delta_LCBD_map,
                                         legend.plants_delta_LCBD_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(plants_delta_LCBD_map_all)



plants_delta_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_plant_delta_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors3,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors3,
  # values = scales::rescale(quantile(SSI_plant_delta_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_plant_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors3,
    na.value = 'gray90',
    name = 'Delta_SSI'
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

legend.plants_delta_SSI_map = legend.func(mycolors = colors3,
                                           mylabels =  round(breaks_delta_SSI,
                                                             4)) +
  ggtitle("log (SSI_all / SSI_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_SSI_map = ggplotGrob(plants_delta_SSI_map)
legend.plants_delta_SSI_map = ggplotGrob(legend.plants_delta_SSI_map)
plants_delta_SSI_map_all = arrangeGrob(plants_delta_SSI_map,
                                        legend.plants_delta_SSI_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(plants_delta_SSI_map_all)








#### Bird LCBD: calculation & mapping ####
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
is.rooted(phy_data)

df.trans = st_transform(df, crs = "+proj=eck4") 
df.trans$area = st_area(df.trans)
sum(as.numeric(df.trans$area) == 0)

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
exotic_distri_data$present = 1
comm_bird_exotic = exotic_distri_data %>% 
  complete(ScientificName, RegionID, fill = list(present = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'present')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = present,
              values_fn = mean) %>% 
  complete(RegionID = unique(df$RegionID)) %>% ## assume that 
  ##  the absence regions have no naturalized aliens
  mutate(
    across(everything(), \(x) {
      if (is.factor(x)) {
        fct_expand(x, "0") %>%  tidyr::replace_na("0")
      } else if (is.character(x)) {
        tidyr::replace_na(x, "0")
      } else if (inherits(x, "units")) {
        unit = units(x)
        tidyr::replace_na(x, set_units(0, unit, mode = "standard"))
      } else {
        tidyr::replace_na(x, 0)
      }
    })
  )  %>% left_join(df.trans[,c('RegionID', 'area')],
                   by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID') 

LCBD_bird_exotic = calcu_LCBD_parallel(Tree = phy_data,
                                        Comm = comm_bird_exotic[,unique(exotic_distri_data$ScientificName)],
                                        Area = comm_bird_exotic$area)
save(LCBD_bird_exotic,
     file = 'results/primary_results/LCBD_bird_exotic.rdata')
load('results/primary_results/LCBD_bird_exotic.rdata')
LCBD_bird_exotic$df$delta_LCBDR = LCBD_bird_exotic$df$delta_LCBDR * 1e+06 * 100 
LCBD_bird_exotic_1 = cbind(RegionID = df$RegionID,
                         LCBD_bird_exotic$df)
LCBD_bird_exotic_sf = df.trans %>% left_join(LCBD_bird_exotic_1, by = 'RegionID')
colnames(LCBD_bird_exotic_sf)

#### exotic birds mapping 
LCBD_bird_exotic_sf = LCBD_bird_exotic_sf %>% filter(!is.na(delta_LCBD))
LCBD_bird_exotic_sf$log_EDR = log10(LCBD_bird_exotic_sf$delta_LCBDR)

LCBD_bird_exotic_sf$f_delta_LCBD = cut(LCBD_bird_exotic_sf$delta_LCBD,
                                  breaks = quantile(LCBD_bird_exotic_sf$delta_LCBD,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(LCBD_bird_exotic_sf$delta_LCBDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(LCBD_bird_exotic_sf$delta_LCBDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(LCBD_bird_exotic_sf$delta_LCBDR, na.rm = TRUE) * 1.01

LCBD_bird_exotic_sf$f_log10_delta_LCBDR = cut(LCBD_bird_exotic_sf$delta_LCBDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# delta_LCBDR = 2.352844e-07
LCBD_bird_exotic_sf_2 = LCBD_bird_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

LCBD_bird_exotic_sf_2 %>% filter(is.na(f_log10_delta_LCBDR))

birds_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_bird_exotic_sf, aes(fill = f_delta_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_bird_exotic_sf$delta_LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_LCBD, color = f_delta_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_exotic_ED_map = legend.func(mycolors = colors2,
                                         mylabels = round(quantile(LCBD_bird_exotic_sf$delta_LCBD,
                                                                   probs = seq(0, 1,length.out = 9),
                                                                   na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_exotic_ED_map = ggplotGrob(birds_exotic_ED_map)
legend.birds_exotic_ED_map = ggplotGrob(legend.birds_exotic_ED_map)
birds_exotic_ED_map_all = arrangeGrob(birds_exotic_ED_map,
                                      legend.birds_exotic_ED_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(birds_exotic_ED_map_all)

birds_exotic_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_bird_exotic_sf, aes(fill = f_log10_delta_LCBDR),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'EDR',
    guide = guide_legend(
      title = "EDR",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorEDR by PE and sizEDR by richness
  geom_point(data = subset(LCBD_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_delta_LCBDR, color = f_log10_delta_LCBDR),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'EDR'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_exotic_EDR_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 5)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_exotic_EDR_map = ggplotGrob(birds_exotic_EDR_map)
legend.birds_exotic_EDR_map = ggplotGrob(legend.birds_exotic_EDR_map)
birds_exotic_EDR_map_all = arrangeGrob(birds_exotic_EDR_map,
                                       legend.birds_exotic_EDR_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_exotic_EDR_map_all)




##### natives #####
colnames(native_distri_data)
native_distri_data$present = 1
comm_bird_native = native_distri_data %>% 
  complete(RegionID, ScientificName, fill = list(present = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'present')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = present,
              values_fn = mean) %>% 
  left_join(df.trans[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID') 

LCBD_bird_native = calcu_LCBD_parallel(Tree = phy_data,
                                       Comm = comm_bird_native,
                                       Region_posi = which(colnames(comm_bird_native) == 'RegionID'))

save(LCBD_bird_native, file = 'results/primary_results/LCBD_bird_native.rdata')

load("results/primary_results/LCBD_bird_native.rdata") # load the calculated data from the server

LCBD_bird_native_sf = df.trans %>% left_join(LCBD_bird_native$LCBD_simp_geo, by = 'RegionID')
colnames(LCBD_bird_native_sf)

#### native birds mapping 
LCBD_bird_native_sf = LCBD_bird_native_sf %>% filter(!is.na(LCBD))

LCBD_bird_native_sf$f_LCBD = cut(LCBD_bird_native_sf$LCBD,
                                  breaks = quantile(LCBD_bird_native_sf$LCBD,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_bird_native_sf_2 = LCBD_bird_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


SSI_bird_native_sf = LCBD_bird_native_sf %>% filter(!is.na(SSI))

SSI_bird_native_sf$f_delta_SSI = cut(SSI_bird_native_sf$SSI,
                                      breaks = quantile(SSI_bird_native_sf$SSI,
                                                        probs = seq(0, 1,
                                                                    length.out = 9),
                                                        na.rm = TRUE),
                                      include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_bird_native_sf_2 = SSI_bird_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


birds_native_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_bird_native_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_bird_native_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_bird_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.birds_native_LCBD_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(LCBD_bird_native_sf$LCBD,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_LCBD_map = ggplotGrob(birds_native_LCBD_map)
legend.birds_native_LCBD_map = ggplotGrob(legend.birds_native_LCBD_map)
birds_native_LCBD_map_all = arrangeGrob(birds_native_LCBD_map,
                                         legend.birds_native_LCBD_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(birds_native_LCBD_map_all)




birds_native_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_bird_native_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_bird_native_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_bird_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.birds_native_SSI_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(SSI_bird_native_sf$SSI,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_SSI_map = ggplotGrob(birds_native_SSI_map)
legend.birds_native_SSI_map = ggplotGrob(legend.birds_native_SSI_map)
birds_native_SSI_map_all = arrangeGrob(birds_native_SSI_map,
                                        legend.birds_native_SSI_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(birds_native_SSI_map_all)




##### extant #####
colnames(all_distri_data_c)
all_distri_data_c$present = 1
comm_bird_extant = all_distri_data_c %>% 
  complete(RegionID, ScientificName, fill = list(present = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'present')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = present,
              values_fn = mean) %>% 
  left_join(df.trans[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID') 

LCBD_bird_extant = calcu_LCBD_parallel(Tree = phy_data,
                                      Comm = comm_bird_extant,
                                      Region_posi = which(colnames(comm_bird_extant) == 'RegionID'))

save(LCBD_bird_extant, file = 'results/primary_results/LCBD_bird_extant.rdata')

load("results/primary_results/LCBD_bird_extant.rdata") # load the calculated data from the server

LCBD_bird_extant_sf = df.trans %>% left_join(LCBD_bird_extant$LCBD_simp_geo, by = 'RegionID')
colnames(LCBD_bird_extant_sf)

#### extant birds mapping 
LCBD_bird_extant_sf = LCBD_bird_extant_sf %>% filter(!is.na(LCBD))

LCBD_bird_extant_sf$f_LCBD = cut(LCBD_bird_extant_sf$LCBD,
                                  breaks = quantile(LCBD_bird_extant_sf$LCBD,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_bird_extant_sf_2 = LCBD_bird_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



SSI_bird_extant_sf = LCBD_bird_extant_sf %>% filter(!is.na(SSI))

SSI_bird_extant_sf$f_delta_SSI = cut(SSI_bird_extant_sf$SSI,
                                      breaks = quantile(SSI_bird_extant_sf$SSI,
                                                        probs = seq(0, 1,
                                                                    length.out = 9),
                                                        na.rm = TRUE),
                                      include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_bird_extant_sf_2 = SSI_bird_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



birds_extant_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_bird_extant_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_bird_extant_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_bird_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.birds_extant_LCBD_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(LCBD_bird_extant_sf$LCBD,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_extant_LCBD_map = ggplotGrob(birds_extant_LCBD_map)
legend.birds_extant_LCBD_map = ggplotGrob(legend.birds_extant_LCBD_map)
birds_extant_LCBD_map_all = arrangeGrob(birds_extant_LCBD_map,
                                         legend.birds_extant_LCBD_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(birds_extant_LCBD_map_all)



birds_extant_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_bird_extant_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_bird_extant_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_bird_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.birds_extant_SSI_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(SSI_bird_extant_sf$SSI,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_extant_SSI_map = ggplotGrob(birds_extant_SSI_map)
legend.birds_extant_SSI_map = ggplotGrob(legend.birds_extant_SSI_map)
birds_extant_SSI_map_all = arrangeGrob(birds_extant_SSI_map,
                                        legend.birds_extant_SSI_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(birds_extant_SSI_map_all)





##### Delta_ED (Extant - Native) #####
load('results/primary_results/LCBD_bird_extant.rdata')
load('results/primary_results/LCBD_bird_native.rdata')

LCBD_bird_extant_1 = LCBD_bird_extant$LCBD_simp_geo
colnames(LCBD_bird_extant_1)[2:3] = c("all_LCBD", "all_SSI")

LCBD_bird_native_1 = LCBD_bird_native$LCBD_simp_geo
colnames(LCBD_bird_native_1)[2:3] = c("native_LCBD", "native_SSI")

LCBD_bird_delta = LCBD_bird_native_1 %>% left_join(LCBD_bird_extant_1,
                                                       by = 'RegionID')

LCBD_bird_delta$delta_LCBD = log(LCBD_bird_delta$all_LCBD / 
                                     LCBD_bird_delta$native_LCBD)

LCBD_bird_delta$delta_SSI = log(LCBD_bird_delta$all_SSI / 
                                    LCBD_bird_delta$native_SSI)

LCBD_bird_delta_sf = df_trans %>% left_join(LCBD_bird_delta, by = 'RegionID')

colnames(LCBD_bird_delta_sf)


#### Delta_ED birds mapping 
LCBD_bird_delta_sf = LCBD_bird_delta_sf %>% filter(!is.na(delta_LCBD) & 
                                                         delta_LCBD != 0 
)

LCBD_bird_delta_sf = LCBD_bird_delta_sf %>% filter(!is.na(delta_SSI) & 
                                                         delta_SSI != 0 
)

hist(LCBD_bird_delta_sf$delta_LCBD)

breaks_delta_LCBD = quantile(LCBD_bird_delta_sf$delta_LCBD,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 12)
breaks_delta_LCBD[which(names(breaks_delta_LCBD) == '37.5%')] = 0
LCBD_bird_delta_sf$f_delta_LCBD = cut(LCBD_bird_delta_sf$delta_LCBD,
                                        breaks = breaks_delta_LCBD,
                                        include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks_delta_SSI = quantile(LCBD_bird_delta_sf$delta_SSI,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 4)
breaks_delta_SSI[which(names(breaks_delta_SSI) == '75%')] = 0
LCBD_bird_delta_sf$f_delta_SSI = cut(LCBD_bird_delta_sf$delta_SSI,
                                       breaks = breaks_delta_SSI,
                                       include.lowest = TRUE)

hist(LCBD_bird_delta_sf$delta_SSI)

LCBD_bird_delta_sf_2 = LCBD_bird_delta_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])




birds_delta_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_bird_delta_sf, aes(fill = f_delta_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors6,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors6,
  # values = scales::rescale(quantile(LCBD_bird_delta_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_bird_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_LCBD, color = f_delta_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors6,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.birds_delta_LCBD_map = legend.func(mycolors = colors6,
                                            mylabels = round(breaks_delta_LCBD,
                                                             4)) +
  ggtitle("log (LCBD_all / LCBD_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_delta_LCBD_map = ggplotGrob(birds_delta_LCBD_map)
legend.birds_delta_LCBD_map = ggplotGrob(legend.birds_delta_LCBD_map)
birds_delta_LCBD_map_all = arrangeGrob(birds_delta_LCBD_map,
                                         legend.birds_delta_LCBD_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(birds_delta_LCBD_map_all)



birds_delta_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_bird_delta_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors7,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors7,
  # values = scales::rescale(quantile(SSI_bird_delta_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_bird_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors7,
    na.value = 'gray90',
    name = 'Delta_SSI'
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

legend.birds_delta_SSI_map = legend.func(mycolors = colors7,
                                           mylabels =  round(breaks_delta_SSI,
                                                             4)) +
  ggtitle("log (SSI_all / SSI_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_delta_SSI_map = ggplotGrob(birds_delta_SSI_map)
legend.birds_delta_SSI_map = ggplotGrob(legend.birds_delta_SSI_map)
birds_delta_SSI_map_all = arrangeGrob(birds_delta_SSI_map,
                                        legend.birds_delta_SSI_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(birds_delta_SSI_map_all)







#### Fish LCBD: calculation & mapping ####
load("D:/R projects/Global_ED/data/Fishes/data/my_phy.rdata")
is.rooted(phylo)
load("D:/R projects/Global_ED/data/Fishes/data/my_data_used_final.rdata")
df = st_read("data/Fishes/data/Basin042017_3119.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
#save(df, file = "data/Fishes/data/Basin042017_3119.rdata")
colnames(df_trans)[which(colnames(df_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
data.used_final_exotics$presence = 1

comm_fish_exotic = data.used_final_exotics %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  complete(X1.Basin.Name = unique(df$BasinName)) %>% ## assume that 
  ##  the absence regions have no naturalized aliens
  mutate(
    across(everything(), \(x) {
      if (is.factor(x)) {
        fct_expand(x, "0") %>%  tidyr::replace_na("0")
      } else if (is.character(x)) {
        tidyr::replace_na(x, "0")
      } else if (inherits(x, "units")) {
        unit = units(x)
        tidyr::replace_na(x, set_units(0, unit, mode = "standard"))
      } else {
        tidyr::replace_na(x, 0)
      }
    })
  )  %>% left_join(df_trans[,c('X1.Basin.Name', 'Surf_area')], 
                   #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
                   by = 'X1.Basin.Name') %>% 
  relocate(Surf_area, .after = X1.Basin.Name) %>% 
  arrange('X1.Basin.Name') 

LCBD_fish_exotic = calcu_LCBD_parallel(
  Tree = phylo,
  Comm = comm_fish_exotic[,unique(data.used_final_exotics$valid_names)],
  Area = comm_fish_exotic$Surf_area)


save(LCBD_fish_exotic_sf, file = 'results/primary_results/LCBD_fish_exotic.rdata')
load('results/primary_results/LCBD_fish_exotic.rdata')
LCBD_fish_exotic$df$delta_LCBDR = LCBD_fish_exotic$df$delta_LCBDR * 1e+6
LCBD_fish_exotic_1 = cbind(Basin.Name = sort(unique(df_trans$X1.Basin.Name)),
                         LCBD_fish_exotic$df)
LCBD_fish_exotic_sf = df_trans %>% left_join(LCBD_fish_exotic_1,
                                           by = join_by('X1.Basin.Name' == 'Basin.Name'))
colnames(LCBD_fish_exotic_sf)

#### exotic fishes mapping 
LCBD_fish_exotic_sf = LCBD_fish_exotic_sf %>% filter(!is.na(delta_LCBD))
LCBD_fish_exotic_sf$log_EDR = log10(LCBD_fish_exotic_sf$delta_LCBDR)

LCBD_fish_exotic_sf$f_delta_LCBD = cut(LCBD_fish_exotic_sf$delta_LCBD,
                                  breaks = quantile(LCBD_fish_exotic_sf$delta_LCBD,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(LCBD_fish_exotic_sf$delta_LCBDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(LCBD_fish_exotic_sf$delta_LCBDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(LCBD_fish_exotic_sf$delta_LCBDR, na.rm = TRUE) * 1.01

LCBD_fish_exotic_sf$f_log10_delta_LCBDR = cut(LCBD_fish_exotic_sf$delta_LCBDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# delta_LCBDR = 2.352844e-07
LCBD_fish_exotic_sf_2 = LCBD_fish_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

LCBD_fish_exotic_sf_2 %>% filter(is.na(f_log10_delta_LCBDR))

fishes_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_fish_exotic_sf, aes(fill = f_delta_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_fish_exotic_sf$delta_LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_fish_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_LCBD, color = f_delta_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_exotic_ED_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(LCBD_fish_exotic_sf$delta_LCBD,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_exotic_ED_map = ggplotGrob(fishes_exotic_ED_map)
legend.fishes_exotic_ED_map = ggplotGrob(legend.fishes_exotic_ED_map)
fishes_exotic_ED_map_all = arrangeGrob(fishes_exotic_ED_map,
                                       legend.fishes_exotic_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(fishes_exotic_ED_map_all)

fishes_exotic_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_fish_exotic_sf, aes(fill = f_log10_delta_LCBDR),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'EDR',
    guide = guide_legend(
      title = "EDR",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorEDR by PE and sizEDR by richness
  geom_point(data = subset(LCBD_fish_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_delta_LCBDR, color = f_log10_delta_LCBDR),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'EDR'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Exotic_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_exotic_EDR_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("MEDR (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_exotic_EDR_map = ggplotGrob(fishes_exotic_EDR_map)
legend.fishes_exotic_EDR_map = ggplotGrob(legend.fishes_exotic_EDR_map)
fishes_exotic_EDR_map_all = arrangeGrob(fishes_exotic_EDR_map,
                                        legend.fishes_exotic_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishes_exotic_EDR_map_all)




##### natives #####
data.used_final_natives$presence = 1

comm_fish_native = data.used_final_natives %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_trans[,c('X1.Basin.Name', 'Surf_area')], 
            #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
            by = 'X1.Basin.Name') %>% 
  relocate(Surf_area, .after = X1.Basin.Name)%>% 
  arrange('X1.Basin.Name') 

LCBD_fish_native = calcu_LCBD_parallel(
  Tree = phylo,
  Comm = comm_fish_native,
  Region_posi = which(colnames(comm_fish_native) == 'X1.Basin.Name'))

save(LCBD_fish_native, file = 'results/primary_results/LCBD_fish_native.rdata')

load('results/primary_results/LCBD_fish_native.rdata')

LCBD_fish_native_sf = df_trans %>% left_join(LCBD_fish_native$LCBD_simp_geo,
                                           by = join_by('X1.Basin.Name' == 'RegionID'))
colnames(LCBD_fish_native_sf)

#### native fishes mapping 
LCBD_fish_native_sf = LCBD_fish_native_sf %>% filter(!is.na(LCBD))

LCBD_fish_native_sf$f_LCBD = cut(LCBD_fish_native_sf$LCBD,
                                  breaks = quantile(LCBD_fish_native_sf$LCBD,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_fish_native_sf_2 = LCBD_fish_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


SSI_fish_native_sf = LCBD_fish_native_sf %>% filter(!is.na(SSI))

SSI_fish_native_sf$f_delta_SSI = cut(SSI_fish_native_sf$SSI,
                                      breaks = quantile(SSI_fish_native_sf$SSI,
                                                        probs = seq(0, 1,
                                                                    length.out = 9),
                                                        na.rm = TRUE),
                                      include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_fish_native_sf_2 = SSI_fish_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


fishes_native_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_fish_native_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_fish_native_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_fish_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.fishes_native_LCBD_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(LCBD_fish_native_sf$LCBD,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_native_LCBD_map = ggplotGrob(fishes_native_LCBD_map)
legend.fishes_native_LCBD_map = ggplotGrob(legend.fishes_native_LCBD_map)
fishes_native_LCBD_map_all = arrangeGrob(fishes_native_LCBD_map,
                                         legend.fishes_native_LCBD_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(fishes_native_LCBD_map_all)




fishes_native_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_fish_native_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_fish_native_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_fish_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.fishes_native_SSI_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(SSI_fish_native_sf$SSI,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_native_SSI_map = ggplotGrob(fishes_native_SSI_map)
legend.fishes_native_SSI_map = ggplotGrob(legend.fishes_native_SSI_map)
fishes_native_SSI_map_all = arrangeGrob(fishes_native_SSI_map,
                                        legend.fishes_native_SSI_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(fishes_native_SSI_map_all)



##### extant #####
data.used_final$presence = 1

comm_fish_extant = data.used_final %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_trans[,c('X1.Basin.Name', 'Surf_area')], 
            #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
            by = 'X1.Basin.Name') %>% 
  relocate(Surf_area, .after = X1.Basin.Name) %>% 
  arrange('X1.Basin.Name') 

LCBD_fish_extant = calcu_LCBD_parallel(
  Tree = phylo,
  Comm = comm_fish_extant,
  Region_posi = which(colnames(comm_fish_extant) == 'X1.Basin.Name'))

save(LCBD_fish_extant, file = 'results/primary_results/LCBD_fish_extant.rdata')

load('results/primary_results/LCBD_fish_extant.rdata')

LCBD_fish_extant_sf = df_trans %>% left_join(LCBD_fish_extant$LCBD_simp_geo,
                                             by = join_by('X1.Basin.Name' == 'RegionID'))
colnames(LCBD_fish_extant_sf)


#### extant fishes mapping 
LCBD_fish_extant_sf = LCBD_fish_extant_sf %>% filter(!is.na(LCBD))

LCBD_fish_extant_sf$f_LCBD = cut(LCBD_fish_extant_sf$LCBD,
                                  breaks = quantile(LCBD_fish_extant_sf$LCBD,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

LCBD_fish_extant_sf_2 = LCBD_fish_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



SSI_fish_extant_sf = LCBD_fish_extant_sf %>% filter(!is.na(SSI))

SSI_fish_extant_sf$f_delta_SSI = cut(SSI_fish_extant_sf$SSI,
                                      breaks = quantile(SSI_fish_extant_sf$SSI,
                                                        probs = seq(0, 1,
                                                                    length.out = 9),
                                                        na.rm = TRUE),
                                      include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

SSI_fish_extant_sf_2 = SSI_fish_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



fishes_extant_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_fish_extant_sf, aes(fill = f_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(LCBD_fish_extant_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_fish_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_LCBD, color = f_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.fishes_extant_LCBD_map = legend.func(mycolors = colors2,
                                            mylabels = round(quantile(LCBD_fish_extant_sf$LCBD,
                                                                      probs = seq(0, 1,length.out = 9),
                                                                      na.rm = TRUE), 5)) +
  ggtitle("LCBD") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_extant_LCBD_map = ggplotGrob(fishes_extant_LCBD_map)
legend.fishes_extant_LCBD_map = ggplotGrob(legend.fishes_extant_LCBD_map)
fishes_extant_LCBD_map_all = arrangeGrob(fishes_extant_LCBD_map,
                                         legend.fishes_extant_LCBD_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10),
                                                               c(NA, rep(2, 8), NA)))
plot(fishes_extant_LCBD_map_all)



fishes_extant_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = SSI_fish_extant_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(SSI_fish_extant_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(SSI_fish_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'SSI'
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

legend.fishes_extant_SSI_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(SSI_fish_extant_sf$SSI,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 5)) +
  ggtitle("SSI") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_extant_SSI_map = ggplotGrob(fishes_extant_SSI_map)
legend.fishes_extant_SSI_map = ggplotGrob(legend.fishes_extant_SSI_map)
fishes_extant_SSI_map_all = arrangeGrob(fishes_extant_SSI_map,
                                        legend.fishes_extant_SSI_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(fishes_extant_SSI_map_all)






##### Delta_ED (Extant - Native) #####
load('results/primary_results/LCBD_fish_extant.rdata')
load('results/primary_results/LCBD_fish_native.rdata')

LCBD_fish_extant_1 = LCBD_fish_extant$LCBD_simp_geo
colnames(LCBD_fish_extant_1)[2:3] = c("all_LCBD", "all_SSI")

LCBD_fish_native_1 = LCBD_fish_native$LCBD_simp_geo
colnames(LCBD_fish_native_1)[2:3] = c("native_LCBD", "native_SSI")

LCBD_fish_delta = LCBD_fish_native_1 %>% left_join(LCBD_fish_extant_1,
                                                     by = 'RegionID')

LCBD_fish_delta$delta_LCBD = log(LCBD_fish_delta$all_LCBD / 
                                    LCBD_fish_delta$native_LCBD)

LCBD_fish_delta$delta_SSI = log(LCBD_fish_delta$all_SSI / 
                                   LCBD_fish_delta$native_SSI)

LCBD_fish_delta_sf = df_trans %>% 
  left_join(LCBD_fish_delta, by = join_by('X1.Basin.Name' == 'RegionID'))

colnames(LCBD_fish_delta_sf)


#### Delta_ED fishes mapping 
LCBD_fish_delta_sf = LCBD_fish_delta_sf %>% filter(!is.na(delta_LCBD) & 
                                                       delta_LCBD != 0 
)

LCBD_fish_delta_sf = LCBD_fish_delta_sf %>% filter(!is.na(delta_SSI) & 
                                                       delta_SSI != 0 
)

hist(LCBD_fish_delta_sf$delta_LCBD)

breaks_delta_LCBD = quantile(LCBD_fish_delta_sf$delta_LCBD,
                             probs = seq(0, 1,
                                         length.out = 9),
                             na.rm = TRUE,
                             names = T,
                             digits = 12)
breaks_delta_LCBD[which(names(breaks_delta_LCBD) == '50%')] = 0
LCBD_fish_delta_sf$f_delta_LCBD = cut(LCBD_fish_delta_sf$delta_LCBD,
                                       breaks = breaks_delta_LCBD,
                                       include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks_delta_SSI = quantile(LCBD_fish_delta_sf$delta_SSI,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 4)
breaks_delta_SSI[which(names(breaks_delta_SSI) == '62.5%')] = 0
LCBD_fish_delta_sf$f_delta_SSI = cut(LCBD_fish_delta_sf$delta_SSI,
                                      breaks = breaks_delta_SSI,
                                      include.lowest = TRUE)

hist(LCBD_fish_delta_sf$delta_SSI)

LCBD_fish_delta_sf_2 = LCBD_fish_delta_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])



fishes_delta_LCBD_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_fish_delta_sf, aes(fill = f_delta_LCBD),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'LCBD',
    guide = guide_legend(
      title = 'LCBD',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors4,
  # values = scales::rescale(quantile(LCBD_fish_delta_sf$LCBD,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'LCBD'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_fish_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_LCBD, color = f_delta_LCBD),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'LCBD'
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

legend.fishes_delta_LCBD_map = legend.func(mycolors = colors4,
                                           mylabels = round(breaks_delta_LCBD,
                                                            4)) +
  ggtitle("log (LCBD_all / LCBD_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_delta_LCBD_map = ggplotGrob(fishes_delta_LCBD_map)
legend.fishes_delta_LCBD_map = ggplotGrob(legend.fishes_delta_LCBD_map)
fishes_delta_LCBD_map_all = arrangeGrob(fishes_delta_LCBD_map,
                                        legend.fishes_delta_LCBD_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(fishes_delta_LCBD_map_all)



fishes_delta_SSI_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = LCBD_fish_delta_sf, aes(fill = f_delta_SSI),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'SSI',
    guide = guide_legend(
      title = 'SSI',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors5,
  # values = scales::rescale(quantile(SSI_fish_delta_sf$SSI,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'SSI'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(LCBD_fish_delta_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_SSI, color = f_delta_SSI),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'Delta_SSI'
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

legend.fishes_delta_SSI_map = legend.func(mycolors = colors5,
                                          mylabels =  round(breaks_delta_SSI,
                                                            4)) +
  ggtitle("log (SSI_all / SSI_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_delta_SSI_map = ggplotGrob(fishes_delta_SSI_map)
legend.fishes_delta_SSI_map = ggplotGrob(legend.fishes_delta_SSI_map)
fishes_delta_SSI_map_all = arrangeGrob(fishes_delta_SSI_map,
                                       legend.fishes_delta_SSI_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(fishes_delta_SSI_map_all)








#### Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

##### LCBD #####
# Compare the LCBD patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.015         # small horizontal gap between columns

figs_LCBD_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_LCBD_map_all,  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_LCBD_map_all,  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_LCBD_map_all,   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_LCBD_map_all,   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_extant_LCBD_map_all,   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_LCBD_map_all,    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_LCBD_map_all, x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_LCBD_map_all, x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_LCBD_map_all,  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_LCBD_map_all,   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_LCBD_map_all,   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_LCBD_map_all,   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
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

png(filename = 'figures/figs_LCBD_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_LCBD_nat_extant_delta
dev.off() #turn off device and finalize file




plot_width = 0.5
plot_height = 0.20

figs_LCBD_nat_exo = ggdraw() +
  draw_plot(plants_native_LCBD_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_LCBD_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_LCBD_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_LCBD_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_LCBD_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_LCBD_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_LCBD_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_exotic_LCBD_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_LCBD_nat_exo.png',
    #plot = figs_LCBD_nat_exo,
    height=40, width=25, # setting emfPlusFontToPath=TRUE to 
    res = 300,
    # ensure text looks correct on the viewing system
    units = 'cm')
figs_LCBD_nat_exo
dev.off() #turn off device and finalize file


#emf('figures/figs_LCBD_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_LCBD_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_LCBD_nat_extant = ggdraw() +
  draw_plot(plants_native_LCBD_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_LCBD_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_LCBD_map_all, x = 0, y = 0.53, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_LCBD_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_LCBD_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_LCBD_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_LCBD_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_LCBD_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_LCBD_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_LCBD_nat_extant
dev.off() #turn off device and finalize file




##### SSI #####
# Compare the SSI patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.015         # small horizontal gap between columns

figs_SSI_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_SSI_map_all,  x = 0,                 y = 0.75,
            width = plot_width, height = plot_height) +
  draw_plot(plants_extant_SSI_map_all,  x = plot_width + gap,  y = 0.75,
            width = plot_width, height = plot_height) +
  draw_plot(plants_delta_SSI_map_all,   x = 2*(plot_width + gap), y = 0.75,
            width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_SSI_map_all,   x = 0,                 y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(birds_extant_SSI_map_all,   x = plot_width + gap,  y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(birds_delta_SSI_map_all,    x = 2*(plot_width + gap), y = 0.5,
            width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_SSI_map_all, x = 0,                 y = 0.25,
            width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_SSI_map_all, x = plot_width + gap,  y = 0.25,
            width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_SSI_map_all,  x = 2*(plot_width + gap), y = 0.25,
            width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_SSI_map_all,   x = 0,                 y = 0,
            width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_SSI_map_all,   x = plot_width + gap,  y = 0,
            width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_SSI_map_all,   x = 2*(plot_width + gap), y = 0,
            width = plot_width, height = plot_height) +
  
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

png(filename = 'figures/figs_SSI_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_SSI_nat_extant_delta
dev.off() #turn off device and finalize file



plot_width = 0.5
plot_height = 0.20

figs_SSI_nat_exo = ggdraw() +
  draw_plot(plants_native_SSI_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_SSI_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_SSI_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_SSI_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_SSI_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_SSI_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_SSI_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_exotic_SSI_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_SSI_nat_exo.png',
    #plot = figs_SSI_nat_exo,
    height=40, width=25, # setting emfPlusFontToPath=TRUE to 
    res = 300,
    # ensure text looks correct on the viewing system
    units = 'cm')
figs_SSI_nat_exo
dev.off() #turn off device and finalize file


#emf('figures/figs_SSI_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_SSI_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_SSI_nat_extant = ggdraw() +
  draw_plot(plants_native_SSI_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_SSI_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_SSI_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_SSI_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_SSI_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_SSI_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_SSI_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_SSI_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_SSI_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_SSI_nat_extant
dev.off() #turn off device and finalize file

