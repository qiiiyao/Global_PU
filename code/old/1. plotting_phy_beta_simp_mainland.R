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
source('code/functions/calculating_phy_turnover_func_2.R')

# load the background data for plotting the world map plot
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")

# load the background data for plotting the world map plot
df_trans = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_mainland = df_trans %>% filter(!(Island == 1 & 
                           Area < 5e3))
plot(df_mainland$geometry)
tdwg_mainlands = as.character(sort(unique(df_mainland$RegionID)))
#df_trans$area = st_area(df_trans)
#load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")
df_sub = st_read("data/Plants/TDWG4_Subset/TDWG4_Subset.shp")


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
colors4 = scico::scico(n=8, palette = "bam")
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


#1. Mammal turnover: calculation & mapping----
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")
is.rooted(spec_phy.3)

##1.1 exotics-----
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)
sp_overlap_dat_1$presence = 1

comm_mammal_exotic = sp_overlap_dat_1 %>% 
  complete(Binomial, RegionID, fill = list(presence = 0)) %>%
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
  )  %>% left_join(df_trans[,c('RegionID', 'Area')],
                   by = 'RegionID') %>% 
  relocate(Area, .after = RegionID) %>% 
  dplyr::arrange(RegionID)

phy_turn_mammal_exotic = calcu_phy_turn_simple(
  tree = spec_phy.3,
  x = comm_mammal_exotic,
  Region_posi = which(colnames(comm_mammal_exotic) == 'RegionID'))

save(phy_turn_mammal_exotic,
     file = 'results/primary_results/distances_beta/phy_turn_mammal_exotic.rdata')

load('results/primary_results/distances_beta/phy_turn_mammal_exotic.rdata')
which(phy_turn_mammal_exotic < 0)

tdwg_mainlands_exo_mammal = intersect(colnames(phy_turn_mammal_exotic), tdwg_mainlands)
phy_turn_mammal_exotic = phy_turn_mammal_exotic[tdwg_mainlands_exo_mammal, tdwg_mainlands_exo_mammal]

phy_turn_mammal_exotic_df = data.frame(RegionID = as.integer(tdwg_mainlands_exo_mammal),
                                       phy_turn = colMeans(phy_turn_mammal_exotic,
                                                           na.rm = T))
phy_turn_mammal_exotic_sf = df_trans %>% left_join(phy_turn_mammal_exotic_df,
                                                   by = 'RegionID')
colnames(phy_turn_mammal_exotic_sf)

#### exotic mammals mapping 
phy_turn_mammal_exotic_sf = phy_turn_mammal_exotic_sf %>% filter(!is.na(phy_turn))

phy_turn_mammal_exotic_sf$f_phy_turn = cut(phy_turn_mammal_exotic_sf$phy_turn,
                                           breaks = quantile(phy_turn_mammal_exotic_sf$phy_turn,
                                                             probs = seq(0, 1,
                                                                         length.out = 9),
                                                             na.rm = TRUE),
                                           include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

phy_turn_mammal_exotic_sf_2 = phy_turn_mammal_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])


mammals_exotic_phy_turn_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = phy_turn_mammal_exotic_sf, aes(fill = f_phy_turn),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'phy_turn',
    guide = guide_legend(
      title = 'phy_turn',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(phy_turn_mammal_exotic_sf$phy_turn,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'phy_turn'
  #  ) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'phy_turn'
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

legend.mammals_exotic_phy_turn_map = legend.func(mycolors = colors2,
                                                 mylabels = round(quantile(phy_turn_mammal_exotic_sf$phy_turn,
                                                                           probs = seq(0, 1,length.out = 9),
                                                                           na.rm = TRUE), 5)) +
  ggtitle("phy_turn") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_exotic_phy_turn_map = ggplotGrob(mammals_exotic_phy_turn_map)
legend.mammals_exotic_phy_turn_map = ggplotGrob(legend.mammals_exotic_phy_turn_map)
mammals_exotic_phy_turn_map_all = arrangeGrob(mammals_exotic_phy_turn_map,
                                              legend.mammals_exotic_phy_turn_map,
                                              ncol = 1,
                                              layout_matrix = rbind(matrix(1, 4, 10),
                                                                    c(NA, rep(2, 8), NA)))
plot(mammals_exotic_phy_turn_map_all)



##1.2 natives----
load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')

tdwg_mainlands_nati_mammal = intersect(colnames(phy_turn_mammal_native),
                                      tdwg_mainlands)
turnover_mammal_native_mat = phy_turn_mammal_native[tdwg_mainlands_nati_mammal,
                                                        tdwg_mainlands_nati_mammal]

turnover_mammal_native_df = data.frame(RegionID = as.integer(colnames(turnover_mammal_native_mat)),
                                    turnover = colMeans(turnover_mammal_native_mat,
                                                        na.rm = T))

turnover_mammal_native_sf = df_trans %>% left_join(turnover_mammal_native_df,
                                                   by = 'RegionID')
colnames(turnover_mammal_native_sf)

#### native mammals mapping 
turnover_mammal_native_sf = turnover_mammal_native_sf %>% filter(!is.na(turnover))

turnover_mammal_native_sf$f_turnover = cut(turnover_mammal_native_sf$turnover,
                                           breaks = quantile(turnover_mammal_native_sf$turnover,
                                                             probs = seq(0, 1,
                                                                         length.out = 9),
                                                             na.rm = TRUE),
                                           include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

mammals_native_turnover_map =
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
  geom_sf(data = (turnover_mammal_native_sf %>% filter(!(Island == 1 & 
                                                           Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_mammal_native_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.mammals_native_turnover_map = legend.func(mycolors = colors2,
                                                 mylabels = round(quantile(turnover_mammal_native_sf$turnover,
                                                                           probs = seq(0, 1,length.out = 9),
                                                                           na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_turnover_map = ggplotGrob(mammals_native_turnover_map)
legend.mammals_native_turnover_map = ggplotGrob(legend.mammals_native_turnover_map)
mammals_native_turnover_map_all = arrangeGrob(mammals_native_turnover_map,
                                              legend.mammals_native_turnover_map,
                                              ncol = 1,
                                              layout_matrix = rbind(matrix(1, 4, 10),
                                                                    c(NA, rep(2, 8), NA)))
plot(mammals_native_turnover_map_all)



##1.3 extant regional assemlage: natives + exotics ----
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

unique(sp_mammal_all$Inv._stage)
colnames(sp_mammal_all)

#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_mammal_all$Binomial = gsub(' ', '_', sp_mammal_all$Binomial)
sp_mammal_all$presence = 1

region1_2 = sp_mammal_all %>% filter(RegionID %in% c('1', '2'))

comm_mammal_all = sp_mammal_all %>% 
  complete(Binomial, RegionID, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(names_from = Binomial,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_trans[,c('RegionID', 'Area')],
            by = 'RegionID') %>% 
  relocate(Area, .after = RegionID) %>% 
  dplyr::arrange('RegionID') 

phy_turn_mammal_extant = calcu_phy_turn_simple(
  Tree = spec_phy.3,
  Comm = comm_mammal_all,
  Region_posi = which(colnames(comm_mammal_all) == 'RegionID'))

save(phy_turn_mammal_extant,
     file = 'results/primary_results/distances_beta/phy_turn_mammal_extant.rdata')

load('results/primary_results/distances_beta/phy_turn_mammal_extant.rdata')

tdwg_mainlands_all_mammal = intersect(colnames(phy_turn_mammal_extant),
                                       tdwg_mainlands)
turnover_mammal_all_mat = phy_turn_mammal_extant[tdwg_mainlands_all_mammal,
                                              tdwg_mainlands_all_mammal]

turnover_mammal_all_df = data.frame(RegionID = as.integer(colnames(turnover_mammal_all_mat)),
                                       turnover = colMeans(turnover_mammal_all_mat,
                                                           na.rm = T))

turnover_mammal_all_sf = df_trans %>% left_join(turnover_mammal_all_df, by = 'RegionID')
colnames(turnover_mammal_all_sf)

#### all mammals mapping 
turnover_mammal_all_sf = turnover_mammal_all_sf %>% filter(!is.na(turnover))

turnover_mammal_all_sf$f_turnover = cut(turnover_mammal_all_sf$turnover,
                                        breaks = quantile(turnover_mammal_all_sf$turnover,
                                                          probs = seq(0, 1,
                                                                      length.out = 9),
                                                          na.rm = TRUE),
                                        include.lowest = TRUE)


mammals_all_turnover_map =
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
  geom_sf(data =  (turnover_mammal_all_sf %>% filter(!(Island == 1 & 
                                                         Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_mammal_all_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.mammals_all_turnover_map = legend.func(mycolors = colors2,
                                              mylabels = round(quantile(turnover_mammal_all_sf$turnover,
                                                                        probs = seq(0, 1,length.out = 9),
                                                                        na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_all_turnover_map = ggplotGrob(mammals_all_turnover_map)
legend.mammals_all_turnover_map = ggplotGrob(legend.mammals_all_turnover_map)
mammals_all_turnover_map_all = arrangeGrob(mammals_all_turnover_map,
                                           legend.mammals_all_turnover_map,
                                           ncol = 1,
                                           layout_matrix = rbind(matrix(1, 4, 10),
                                                                 c(NA, rep(2, 8), NA)))
plot(mammals_all_turnover_map_all)




##1.4 Delta_turnover (Extant - Native) ----
load("results/primary_results/rstudio-export/phy_turn_mammal_all.rdata")
load("results/primary_results/rstudio-export/phy_turn_mammal_native.rdata")

### alternative 1. detla = mean(extant - native)
turnover_mammal_extant_mat = phy_turn_mammal_extant
turnover_mammal_native_mat = phy_turn_mammal_native

turnover_mammal_delta_mat = log((turnover_mammal_extant_mat+0.001) /
                                  (turnover_mammal_native_mat+0.001))

tdwg_mainlands_delta_mammal = intersect(colnames(turnover_mammal_delta_mat),
                                      tdwg_mainlands)
turnover_mammal_delta_mat = turnover_mammal_delta_mat[tdwg_mainlands_delta_mammal,
                                              tdwg_mainlands_delta_mammal]

turnover_mammal_delta = data.frame(RegionID = as.integer(colnames(turnover_mammal_delta_mat)),
                                   delta_turnover = colMeans(turnover_mammal_delta_mat,
                                                             na.rm = T))

turnover_mammal_delta_sf = df_trans %>% left_join(turnover_mammal_delta,
                                                  by = 'RegionID')

colnames(turnover_mammal_delta_sf)

#### Delta_ED mammals mapping 
turnover_mammal_delta_sf = turnover_mammal_delta_sf %>% filter(!is.na(delta_turnover))


hist(turnover_mammal_delta_sf$delta_turnover)

breaks_delta_turnover = quantile(turnover_mammal_delta_sf$delta_turnover,
                                 probs = seq(0, 1,
                                             length.out = 9),
                                 na.rm = TRUE,
                                 names = T,
                                 digits = 12)
breaks_delta_turnover[which(names(breaks_delta_turnover) == '50%')] = 0
turnover_mammal_delta_sf$f_delta_turnover = cut(turnover_mammal_delta_sf$delta_turnover,
                                                breaks = breaks_delta_turnover,
                                                include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

mammals_delta_turnover_map =
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
  geom_sf(data = (turnover_mammal_delta_sf %>% filter(!(Island == 1 & 
                                                          Area < 1e4))),
          aes(fill = f_delta_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors4,
  # values = scales::rescale(quantile(turnover_mammal_delta_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
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

legend.mammals_delta_turnover_map = legend.func(mycolors = colors4,
                                                mylabels = round(breaks_delta_turnover,
                                                                 4)) +
  ggtitle("log (Turnover_all / Turnover_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_delta_turnover_map = ggplotGrob(mammals_delta_turnover_map)
legend.mammals_delta_turnover_map = ggplotGrob(legend.mammals_delta_turnover_map)
mammals_delta_turnover_map_all = arrangeGrob(mammals_delta_turnover_map,
                                             legend.mammals_delta_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(mammals_delta_turnover_map_all)


## 1.5 Patitioning delta ED into 5 possible ways ----
load("results/primary_results/phy_turn_mammal_native.rdata")
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

region_pairs = combn(colnames(phy_turn_mammal_native), 2)

turnover_mammal_path3_5_mat = phy_turn_mammal_native
turnover_mammal_path4_7_mat = phy_turn_mammal_native
turnover_mammal_path6_mat = phy_turn_mammal_native

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
    
    phy_turn_mammal_path6 = calcu_turnover_simple(
      Tree = spec_phy.3,
      Comm = comm_mammal_path6[,2:ncol(comm_mammal_path6)])
    
    turnover_mammal_path6_mat[region1, region2] = phy_turn_mammal_path6[1,2]
    turnover_mammal_path6_mat[region2, region1] = phy_turn_mammal_path6[1,2]
    
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
      
      phy_turn_mammal_path3_5 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_mammal_path3_5[,2:ncol(comm_mammal_path3_5)])
      
      turnover_mammal_path3_5_mat[region1, region2] = phy_turn_mammal_path3_5[1,2]
      turnover_mammal_path3_5_mat[region2, region1] = phy_turn_mammal_path3_5[1,2]
      
      
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
      
      phy_turn_mammal_path4_7 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_mammal_path4_7[,2:ncol(comm_mammal_path4_7)])
      
      turnover_mammal_path4_7_mat[region1, region2] = phy_turn_mammal_path4_7[1,2]
      turnover_mammal_path4_7_mat[region2, region1] = phy_turn_mammal_path4_7[1,2]
      
    }
    
  }
  
}


turnover_mammal_delta_path6_mat = log((turnover_mammal_path6_mat+0.001) /
                                        (phy_turn_mammal_native+0.001))
which(turnover_mammal_delta_path6_mat > 0)

turnover_mammal_delta_path4_7_mat = log((turnover_mammal_path4_7_mat+0.001) /
                                          (turnover_mammal_native_mat+0.001))
which(turnover_mammal_delta_path4_7_mat > 0)
which(turnover_mammal_delta_path4_7_mat < 0)

turnover_mammal_delta_path3_5_mat = log((turnover_mammal_path3_5_mat+0.001) /
                                          (turnover_mammal_native_mat+0.001))
which(turnover_mammal_delta_path3_5_mat > 0)
which(turnover_mammal_delta_path3_5_mat < 0)

#2. Plant turnover: calculation & mapping ----
shp.glonaf.new = st_read("data/Plants/shp_glonaf_new_eck4.shp")
shp.glonaf.mainland = shp.glonaf.new %>% filter(!(Island == 1 & 
                                                    Area < 5e3)) 
plot(shp.glonaf.mainland$geometry)
glonaf_mainlands = as.character(sort(unique(shp.glonaf.mainland$Region_id)))

phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")

load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
colnames(df.native.650) == colnames(df.natu.650)

## 2.1 Remove apomictic species ----
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

##2.2 exotics #####
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

phy_turn_plant_exotic = calcu_phy_turn_simple(Tree = phylo_plant_exotic,
                                                Comm = comm_plant_exotic,
                                                Region_posi = which(colnames(comm_plant_exotic) == 'Region_id'))

save(phy_turn_plant_exotic,
     file = 'results/primary_results/distances_beta/phy_turn_plant_exotic.rdata')

load('results/primary_results/distances_beta/phy_turn_plant_exotic.rdata')

turnover_plant_exotic_mat = phy_turn_plant_exotic
rm(phy_turn_plant_exotic)

glonaf_mainlands_exo_plant = intersect(colnames(turnover_plant_exotic_mat),
                                       glonaf_mainlands)
turnover_plant_exotic_mat = turnover_plant_exotic_mat[glonaf_mainlands_exo_plant,
                                                      glonaf_mainlands_exo_plant]

turnover_plant_exotic_df = data.frame(Region_id = as.integer(colnames(turnover_plant_exotic_mat)),
                                          turnover = colMeans(turnover_plant_exotic_mat,
                                                           na.rm = T))

turnover_plant_exotic_sf = shp.glonaf.mainland %>% left_join(turnover_plant_exotic_df,
                                                          by = 'Region_id')
colnames(turnover_plant_exotic_sf)

##2.3 exotic plants mapping 
turnover_plant_exotic_sf = turnover_plant_exotic_sf %>% filter(!is.na(turnover))

turnover_plant_exotic_sf$f_turnover = cut(turnover_plant_exotic_sf$turnover,
                                          breaks = quantile(turnover_plant_exotic_sf$turnover,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE),
                                          include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

turnover_plant_exotic_sf_2 = turnover_plant_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

plants_exotic_turnover_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = turnover_plant_exotic_sf, aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('exotic_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_exotic_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_plant_exotic_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_turnover_map = ggplotGrob(plants_exotic_turnover_map)
legend.plants_exotic_turnover_map = ggplotGrob(legend.plants_exotic_turnover_map)
plants_exotic_turnover_map_all = arrangeGrob(plants_exotic_turnover_map,
                                             legend.plants_exotic_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(plants_exotic_turnover_map_all)



##2.3 natives #####
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

phy_turn_plant_native = calcu_phy_turn_simple(Tree = phylo_plant_native,
                                                Comm = comm_plant_native,
                                                Region_posi = which(colnames(comm_plant_native) == 'Region_id'))

#save(phy_turn_plant_native,
#     file = 'results/primary_results/distances_beta/phy_turn_plant_native.rdata') # run in the server

load("results/primary_results/phy_turn_plant_native.rdata") # load the calculated data from the server
load("results/primary_results/distances_beta/phy_turn_plant_native.rdata") # load the calculated data from the server

turnover_plant_native_mat = phy_turn_plant_native
rm(phy_turn_plant_native)

glonaf_mainlands_nati_plant = intersect(colnames(turnover_plant_native_mat),
                                       glonaf_mainlands)
turnover_plant_native_mat = turnover_plant_native_mat[glonaf_mainlands_nati_plant,
                                                      glonaf_mainlands_nati_plant]

turnover_plant_native_df = data.frame(Region_id = as.integer(colnames(turnover_plant_native_mat)),
                                          turnover = colMeans(turnover_plant_native_mat,
                                                              na.rm = T))

turnover_plant_native_sf = shp.glonaf.mainland %>% left_join(turnover_plant_native_df,
                                                          by = 'Region_id')
colnames(turnover_plant_native_sf)

#### native plants mapping 
turnover_plant_native_sf = turnover_plant_native_sf %>% filter(!is.na(turnover))

turnover_plant_native_sf$f_turnover = cut(turnover_plant_native_sf$turnover,
                                          breaks = quantile(turnover_plant_native_sf$turnover,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE),
                                          include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
plants_native_turnover_map =
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
  geom_sf(data = (turnover_plant_native_sf %>% filter(!(Island == 1 & 
                                                          Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.plants_native_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_plant_native_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_turnover_map = ggplotGrob(plants_native_turnover_map)
legend.plants_native_turnover_map = ggplotGrob(legend.plants_native_turnover_map)
plants_native_turnover_map_all = arrangeGrob(plants_native_turnover_map,
                                             legend.plants_native_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(plants_native_turnover_map_all)




##2.4 extant #####
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

phy_turn_plant_extant = calcu_phy_turn_simple(Tree = phylo,
                                                Comm = comm_plant_extant,
                                                Region_posi = which(colnames(comm_plant_extant) == 'Region_id'))
#save(phy_turn_plant_extant,
#     file = 'results/primary_results/distances_beta/phy_turn_plant_extant.rdata')

load("results/primary_results/distances_beta/phy_turn_plant_extant.rdata") # load the calculated data from the server

turnover_plant_extant_mat = phy_turn_plant_extant
rm(phy_turn_plant_extant)

glonaf_mainlands_all_plant = intersect(colnames(turnover_plant_extant_mat),
                                        glonaf_mainlands)
turnover_plant_extant_mat = turnover_plant_extant_mat[glonaf_mainlands_all_plant,
                                                      glonaf_mainlands_all_plant]

turnover_plant_extant_df = data.frame(Region_id = as.integer(colnames(turnover_plant_extant_mat)),
                                      turnover = colMeans(turnover_plant_extant_mat,
                                                          na.rm = T))

turnover_plant_extant_sf = shp.glonaf.mainland %>% left_join(turnover_plant_extant_df,
                                                          by = 'Region_id')
colnames(turnover_plant_extant_sf)

#### extant plants mapping 
turnover_plant_extant_sf = turnover_plant_extant_sf %>% filter(!is.na(turnover))

turnover_plant_extant_sf$f_turnover = cut(turnover_plant_extant_sf$turnover,
                                          breaks = quantile(turnover_plant_extant_sf$turnover,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE),
                                          include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
plants_extant_turnover_map =
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
  geom_sf(data = (turnover_plant_extant_sf %>% filter(!(Island == 1 & 
                                                          Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.plants_extant_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_plant_extant_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_turnover_map = ggplotGrob(plants_extant_turnover_map)
legend.plants_extant_turnover_map = ggplotGrob(legend.plants_extant_turnover_map)
plants_extant_turnover_map_all = arrangeGrob(plants_extant_turnover_map,
                                             legend.plants_extant_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(plants_extant_turnover_map_all)



##2.5 Delta_turnover (Extant - Native) #####
load('results/primary_results/distances_beta/phy_turn_plant_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_plant_native.rdata')

### 1. detla = mean(extant) - mean(native)
#turnover_plant_extant_mat = phy_turn_plant_extant$beta_mat
#turnover_plant_extant = data.frame(Region_id = colnames(turnover_plant_extant_mat),
#                                   extant_turnover = colMeans(turnover_plant_extant_mat))
#turnover_plant_extant$Region_id = as.integer(turnover_plant_extant$Region_id)

#turnover_plant_native_mat = phy_turn_plant_native$beta_mat
#turnover_plant_native = data.frame(Region_id = colnames(turnover_plant_native_mat),
##                                   native_turnover = colMeans(turnover_plant_native_mat))
#turnover_plant_native$Region_id = as.integer(turnover_plant_native$Region_id)

#turnover_plant_delta = turnover_plant_native %>% left_join(turnover_plant_extant,
#                                                   by = 'Region_id')

#turnover_plant_delta$delta_turnover = log(turnover_plant_delta$extant_turnover / 
#                                   turnover_plant_delta$native_turnover)

### alternative 1. detla = mean(extant - native)
turnover_plant_delta_mat = log((turnover_plant_extant_mat+0.001) /
                                 (turnover_plant_native_mat+0.001))
rm(phy_turn_plant_extant)
rm(phy_turn_plant_native)

glonaf_mainlands_delta_plant = intersect(colnames(turnover_plant_delta_mat),
                                        glonaf_mainlands)
turnover_plant_delta_mat = turnover_plant_delta_mat[glonaf_mainlands_delta_plant,
                                                      glonaf_mainlands_delta_plant]

turnover_plant_delta = data.frame(Region_id = as.integer(colnames(turnover_plant_delta_mat)),
                                   delta_turnover = colMeans(turnover_plant_delta_mat,
                                                             na.rm = T))

turnover_plant_delta_sf = shp.glonaf.mainland %>% left_join(turnover_plant_delta,
                                                  by = 'Region_id')

colnames(turnover_plant_delta_sf)


#### Delta_ED plants mapping 
turnover_plant_delta_sf = turnover_plant_delta_sf %>% filter(!is.na(delta_turnover) & 
                                                               delta_turnover != 0 
)

hist(turnover_plant_delta_sf$delta_turnover)

breaks_delta_turnover = quantile(turnover_plant_delta_sf$delta_turnover,
                                 probs = seq(0, 1,
                                             length.out = 9),
                                 na.rm = TRUE,
                                 names = T,
                                 digits = 12)
breaks_delta_turnover[which(names(breaks_delta_turnover) == '87.5%')] = 0
turnover_plant_delta_sf$f_delta_turnover = cut(turnover_plant_delta_sf$delta_turnover,
                                               breaks = breaks_delta_turnover,
                                               include.lowest = TRUE)


plants_delta_turnover_map =
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
  geom_sf(data = (turnover_plant_delta_sf %>% filter(!(Island == 1 & 
                                                         Area < 1e4))),
          aes(fill = f_delta_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Delta_plants')+
  theme(legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5)) 

legend.plants_delta_turnover_map = legend.func(mycolors = colors8,
                                               mylabels = round(breaks_delta_turnover,
                                                                4)) +
  ggtitle("log (Turnover_all / Turnover_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_turnover_map = ggplotGrob(plants_delta_turnover_map)
legend.plants_delta_turnover_map = ggplotGrob(legend.plants_delta_turnover_map)
plants_delta_turnover_map_all = arrangeGrob(plants_delta_turnover_map,
                                            legend.plants_delta_turnover_map,
                                            ncol = 1,
                                            layout_matrix = rbind(matrix(1, 4, 10),
                                                                  c(NA, rep(2, 8), NA)))
plot(plants_delta_turnover_map_all)


## 2.6 Patitioning delta ED into 5 possible ways ----
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

region_pairs = combn(colnames(turnover_plant_delta_mat), 2)

turnover_plant_path3_5_mat = turnover_plant_native_mat
turnover_plant_path4_7_mat = turnover_plant_native_mat
turnover_plant_path6_mat = turnover_plant_native_mat

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
    
    phy_turn_plant_path6 = calcu_turnover_simple(
      Tree = spec_phy.3,
      Comm = comm_plant_path6[,2:ncol(comm_plant_path6)])
    
    turnover_plant_path6_mat[region1, region2] = phy_turn_plant_path6[1,2]
    turnover_plant_path6_mat[region2, region1] = phy_turn_plant_path6[1,2]
    
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
      
      phy_turn_plant_path3_5 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_plant_path3_5[,2:ncol(comm_plant_path3_5)])
      
      turnover_plant_path3_5_mat[region1, region2] = phy_turn_plant_path3_5[1,2]
      turnover_plant_path3_5_mat[region2, region1] = phy_turn_plant_path3_5[1,2]
      
      
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
      
      phy_turn_plant_path4_7 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_plant_path4_7[,2:ncol(comm_plant_path4_7)])
      
      turnover_plant_path4_7_mat[region1, region2] = phy_turn_plant_path4_7[1,2]
      turnover_plant_path4_7_mat[region2, region1] = phy_turn_plant_path4_7[1,2]
      
    }
    
  }
  
}


turnover_plant_delta_path6_mat = log((turnover_plant_path6_mat+0.001) /
                                       (turnover_plant_native_mat+0.001))

turnover_plant_delta_path4_7_mat = log((turnover_plant_path4_7_mat+0.001) /
                                         (turnover_plant_native_mat+0.001))
which(turnover_plant_delta_path4_7_mat > 0)
which(turnover_plant_delta_path4_7_mat < 0)

turnover_plant_delta_path3_5_mat = log((turnover_plant_path3_5_mat+0.001) /
                                         (turnover_plant_native_mat+0.001))
which(turnover_plant_delta_path3_5_mat > 0)
which(turnover_plant_delta_path3_5_mat < 0)




#3. Plant turnover_TDWG4: calculation & mapping ####
phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative_TDWG.Rdata")
df.native.650 = df.native.natu.species.650.nonative_TDWG[[1]]
df.natu.650 = df.native.natu.species.650.nonative_TDWG[[2]]
colnames(df.native.650) == colnames(df.natu.650)
RegionID = df_sub$RegionID


## 3.1 Remove apomictic species #####
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

##3.2 exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
colnames(df.natu.650)
df.natu.650$presence = 1
phylo_plant_exotic = drop.tip(phylo_big, setdiff(phylo_big$tip.label, df.natu.650$species))

comm_plant_exotic = df.natu.650 %>% 
  complete(species, RegionID, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_sub[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  dplyr::arrange('RegionID') 

phy_turn_plant_exotic = calcu_phy_turn_simple(Tree = phylo_plant_exotic,
                                                Comm = comm_plant_exotic,
                                                Region_posi = which(colnames(comm_plant_exotic) == 'RegionID'))

save(phy_turn_plant_exotic,
     file = 'results/primary_results/phy_turn_plant_exotic.rdata')

load('results/primary_results/phy_turn_plant_exotic.rdata')

turnover_plant_exotic_mat = phy_turn_plant_exotic$beta_mat
turnover_plant_exotic = data.frame(RegionID = colnames(turnover_plant_exotic_mat),
                                   turnover = colMeans(turnover_plant_exotic_mat, 
                                                       na.rm = T))
turnover_plant_exotic$RegionID = as.integer(turnover_plant_exotic$RegionID)

turnover_plant_exotic_sf = df_sub %>% left_join(turnover_plant_exotic,
                                                by = 'RegionID')
colnames(turnover_plant_exotic_sf)

#### exotic plants mapping 
turnover_plant_exotic_sf = turnover_plant_exotic_sf %>% filter(!is.na(turnover))

turnover_plant_exotic_sf$f_turnover = cut(turnover_plant_exotic_sf$turnover,
                                          breaks = quantile(turnover_plant_exotic_sf$turnover,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE),
                                          include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

turnover_plant_exotic_sf_2 = turnover_plant_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

plants_exotic_turnover_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = turnover_plant_exotic_sf, aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_plant_exotic_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(turnover_plant_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_turnover, color = f_turnover),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('exotic_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_exotic_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_plant_exotic_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_turnover_map = ggplotGrob(plants_exotic_turnover_map)
legend.plants_exotic_turnover_map = ggplotGrob(legend.plants_exotic_turnover_map)
plants_exotic_turnover_map_all = arrangeGrob(plants_exotic_turnover_map,
                                             legend.plants_exotic_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(plants_exotic_turnover_map_all)





##3.3 natives #####
colnames(df.native.650)
df.native.650$presence = 1
phylo_plant_native = drop.tip(phylo_big, setdiff(phylo_big$tip.label, df.native.650$species))
comm_plant_native = df.native.650 %>% 
  ungroup() %>% 
  tidyr::complete(RegionID, species, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(df_sub[,c('RegionID', 'Area')],
            by = 'RegionID') %>% 
  relocate(Area, .after = RegionID) %>% 
  filter(!is.na(RegionID))

comm_plant_native = arrange(comm_plant_native, comm_plant_native$RegionID)

phy_turn_plant_native = calcu_phy_turnover_parallel(Tree = phylo_plant_native,
                                                    Comm = comm_plant_native,
                                                    Region_posi = which(colnames(comm_plant_native) == 'RegionID'))

save(phy_turn_plant_native,
     file = 'results/primary_results/phy_turn_plant_native.rdata') # run in the server

load("results/primary_results/phy_turn_plant_native.rdata") # load the calculated data from the server

turnover_plant_native_mat = phy_turn_plant_native
turnover_plant_native = data.frame(RegionID = colnames(turnover_plant_native_mat),
                                   turnover = colMeans(turnover_plant_native_mat,
                                                       na.rm = T))
turnover_plant_native$RegionID = as.integer(turnover_plant_native$RegionID)

turnover_plant_native_sf = df_sub %>% left_join(turnover_plant_native,
                                                by = 'RegionID') %>% 
  left_join((df[,c("RegionID", "Lon", "Lat")] %>% 
               st_drop_geometry),by = 'RegionID')
colnames(turnover_plant_native_sf)

#### native plants mapping 
turnover_plant_native_sf = turnover_plant_native_sf %>% filter(!is.na(turnover))

turnover_plant_native_sf$f_turnover = cut(turnover_plant_native_sf$turnover,
                                          breaks = quantile(turnover_plant_native_sf$turnover,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE),
                                          include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
plants_native_turnover_map =
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
  geom_sf(data = (turnover_plant_native_sf %>% filter(!(Island == 1 & 
                                                          Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_plant_native_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (turnover_plant_native_sf %>% filter(Island == 1 & 
                                                           Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_turnover),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.plants_native_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_plant_native_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_turnover_map = ggplotGrob(plants_native_turnover_map)
legend.plants_native_turnover_map = ggplotGrob(legend.plants_native_turnover_map)
plants_native_turnover_map_all = arrangeGrob(plants_native_turnover_map,
                                             legend.plants_native_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(plants_native_turnover_map_all)




##3.4 extant #####
df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1
comm_plant_extant = df.extant.650 %>% 
  complete(RegionID, species, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(df_sub[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID') 

phy_turn_plant_extant = calcu_phy_turn_simple(Tree = phylo,
                                                Comm = comm_plant_extant,
                                                Region_posi = which(colnames(comm_plant_extant) == 'RegionID'))
save(phy_turn_plant_extant,
     file = 'results/primary_results/phy_turn_plant_extant.rdata')

load("results/primary_results/phy_turn_plant_extant.rdata") # load the calculated data from the server

turnover_plant_extant_mat = phy_turn_plant_extant$beta_mat
turnover_plant_extant = data.frame(RegionID = colnames(turnover_plant_extant_mat),
                                   turnover = colMeans(turnover_plant_extant_mat))
turnover_plant_extant$RegionID = as.integer(turnover_plant_extant$RegionID)

turnover_plant_extant_sf = df_sub %>% left_join(turnover_plant_extant,
                                                by = 'RegionID')
colnames(turnover_plant_extant_sf)

#### extant plants mapping 
turnover_plant_extant_sf = turnover_plant_extant_sf %>% filter(!is.na(turnover))

turnover_plant_extant_sf$f_turnover = cut(turnover_plant_extant_sf$turnover,
                                          breaks = quantile(turnover_plant_extant_sf$turnover,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE),
                                          include.lowest = TRUE)

# Ensure full coverage with a tiny buffer

plants_extant_turnover_map =
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
  geom_sf(data = (turnover_plant_extant_sf %>% filter(!(Island == 1 & 
                                                          Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_plant_extant_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (turnover_plant_extant_sf %>% filter(Island == 1 & 
                                                           Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_turnover),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.plants_extant_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_plant_extant_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_turnover_map = ggplotGrob(plants_extant_turnover_map)
legend.plants_extant_turnover_map = ggplotGrob(legend.plants_extant_turnover_map)
plants_extant_turnover_map_all = arrangeGrob(plants_extant_turnover_map,
                                             legend.plants_extant_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(plants_extant_turnover_map_all)



##3.5 Delta_turnover (Extant - Native) #####
load('results/primary_results/phy_turn_plant_extant.rdata')
load('results/primary_results/phy_turn_plant_native.rdata')

### 1. detla = mean(extant) - mean(native)
#turnover_plant_extant_mat = phy_turn_plant_extant$beta_mat
#turnover_plant_extant = data.frame(RegionID = colnames(turnover_plant_extant_mat),
#                                   extant_turnover = colMeans(turnover_plant_extant_mat))
#turnover_plant_extant$RegionID = as.integer(turnover_plant_extant$RegionID)

#turnover_plant_native_mat = phy_turn_plant_native$beta_mat
#turnover_plant_native = data.frame(RegionID = colnames(turnover_plant_native_mat),
##                                   native_turnover = colMeans(turnover_plant_native_mat))
#turnover_plant_native$RegionID = as.integer(turnover_plant_native$RegionID)

#turnover_plant_delta = turnover_plant_native %>% left_join(turnover_plant_extant,
#                                                   by = 'RegionID')

#turnover_plant_delta$delta_turnover = log(turnover_plant_delta$extant_turnover / 
#                                   turnover_plant_delta$native_turnover)

### alternative 1. detla = mean(extant - native)
turnover_plant_extant_mat = phy_turn_plant_extant$beta_mat
turnover_plant_native_mat = phy_turn_plant_native$beta_mat
turnover_plant_delta_mat = log((turnover_plant_extant_mat+0.001) /
                                 (turnover_plant_native_mat+0.001))
turnover_plant_delta = data.frame(RegionID = colnames(turnover_plant_delta_mat),
                                  delta_turnover = colMeans(turnover_plant_delta_mat))
turnover_plant_delta$RegionID = as.integer(turnover_plant_delta$RegionID)

turnover_plant_delta_sf = df_sub %>% left_join(turnover_plant_delta,
                                               by = 'RegionID')

colnames(turnover_plant_delta_sf)


#### Delta_ED plants mapping 
turnover_plant_delta_sf = turnover_plant_delta_sf %>% filter(!is.na(delta_turnover) & 
                                                               delta_turnover != 0 
)

hist(turnover_plant_delta_sf$delta_turnover)

breaks_delta_turnover = quantile(turnover_plant_delta_sf$delta_turnover,
                                 probs = seq(0, 1,
                                             length.out = 9),
                                 na.rm = TRUE,
                                 names = T,
                                 digits = 12)
breaks_delta_turnover[which(names(breaks_delta_turnover) == '87.5%')] = 0
turnover_plant_delta_sf$f_delta_turnover = cut(turnover_plant_delta_sf$delta_turnover,
                                               breaks = breaks_delta_turnover,
                                               include.lowest = TRUE)


plants_delta_turnover_map =
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
  geom_sf(data = (turnover_plant_delta_sf %>% filter(!(Island == 1 & 
                                                         Area < 1e4))),
          aes(fill = f_delta_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors8,
  # values = scales::rescale(quantile(turnover_plant_delta_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = (turnover_plant_delta_sf %>% filter(Island == 1 & 
                                                          Area < 1e4)),
             aes(x = Lon, y = Lat, color = f_delta_turnover),
             size = 3,
             shape = 21, stroke = 2, fill = NA, show.legend = F) +
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover'
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

legend.plants_delta_turnover_map = legend.func(mycolors = colors8,
                                               mylabels = round(breaks_delta_turnover,
                                                                4)) +
  ggtitle("log (Turnover_all / Turnover_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_turnover_map = ggplotGrob(plants_delta_turnover_map)
legend.plants_delta_turnover_map = ggplotGrob(legend.plants_delta_turnover_map)
plants_delta_turnover_map_all = arrangeGrob(plants_delta_turnover_map,
                                            legend.plants_delta_turnover_map,
                                            ncol = 1,
                                            layout_matrix = rbind(matrix(1, 4, 10),
                                                                  c(NA, rep(2, 8), NA)))
plot(plants_delta_turnover_map_all)



#4. Bird turnover: calculation & mapping ####
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

##4.1 exotics #####
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

phy_turn_bird_exotic = calcu_phy_turn_simple(Tree = phy_data,
                                               Comm = comm_bird_exotic[,unique(exotic_distri_data$ScientificName)],
                                               Area = comm_bird_exotic$area)
save(phy_turn_bird_exotic,
     file = 'results/primary_results/phy_turn_bird_exotic.rdata')

load('results/primary_results/phy_turn_bird_exotic.rdata')


phy_turn_bird_exotic$df$delta_phy_turnR = phy_turn_bird_exotic$df$delta_phy_turnR * 1e+06 * 100 
phy_turn_bird_exotic_1 = cbind(RegionID = df$RegionID,
                               phy_turn_bird_exotic$df)
phy_turn_bird_exotic_sf = df.trans %>% left_join(phy_turn_bird_exotic_1, by = 'RegionID')
colnames(phy_turn_bird_exotic_sf)

#### exotic birds mapping 
phy_turn_bird_exotic_sf = phy_turn_bird_exotic_sf %>% filter(!is.na(delta_phy_turn))
phy_turn_bird_exotic_sf$log_EDR = log10(phy_turn_bird_exotic_sf$delta_phy_turnR)

phy_turn_bird_exotic_sf$f_delta_phy_turn = cut(phy_turn_bird_exotic_sf$delta_phy_turn,
                                               breaks = quantile(phy_turn_bird_exotic_sf$delta_phy_turn,
                                                                 probs = seq(0, 1,
                                                                             length.out = 9),
                                                                 na.rm = TRUE),
                                               include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(phy_turn_bird_exotic_sf$delta_phy_turnR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(phy_turn_bird_exotic_sf$delta_phy_turnR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(phy_turn_bird_exotic_sf$delta_phy_turnR, na.rm = TRUE) * 1.01

phy_turn_bird_exotic_sf$f_log10_delta_phy_turnR = cut(phy_turn_bird_exotic_sf$delta_phy_turnR,
                                                      breaks = breaks,
                                                      include.lowest = TRUE,
                                                      include.highest = T,
                                                      right = T,
                                                      dig.lab = 2)

# delta_phy_turnR = 2.352844e-07
phy_turn_bird_exotic_sf_2 = phy_turn_bird_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

phy_turn_bird_exotic_sf_2 %>% filter(is.na(f_log10_delta_phy_turnR))

birds_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = phy_turn_bird_exotic_sf, aes(fill = f_delta_phy_turn),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'phy_turn',
    guide = guide_legend(
      title = 'phy_turn',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(phy_turn_bird_exotic_sf$delta_phy_turn,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'phy_turn'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(phy_turn_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_phy_turn, color = f_delta_phy_turn),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'phy_turn'
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
                                         mylabels = round(quantile(phy_turn_bird_exotic_sf$delta_phy_turn,
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
  geom_sf(data = phy_turn_bird_exotic_sf, aes(fill = f_log10_delta_phy_turnR),
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
  geom_point(data = subset(phy_turn_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_delta_phy_turnR, color = f_log10_delta_phy_turnR),
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




##4.2 natives #####
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

phy_turn_bird_native = calcu_phy_turn_simple(Tree = phy_data,
                                               Comm = comm_bird_native,
                                               Region_posi = which(colnames(comm_bird_native) == 'RegionID'))

save(phy_turn_bird_native,
     file = 'results/primary_results/distances_beta/phy_turn_bird_native.rdata')

load("results/primary_results/distances_beta/phy_turn_bird_native.rdata") # load the calculated data from the server

turnover_bird_native_mat = phy_turn_bird_native
rm(phy_turn_bird_native)

tdwg_mainlands_nati_bird = intersect(colnames(turnover_bird_native_mat),
                                       tdwg_mainlands)
turnover_bird_native_mat = turnover_bird_native_mat[tdwg_mainlands_nati_bird,
                                                tdwg_mainlands_nati_bird]

turnover_bird_native_df = data.frame(RegionID = as.integer(colnames(turnover_bird_native_mat)),
                                       turnover = colMeans(turnover_bird_native_mat,
                                                           na.rm = T))

turnover_bird_native_sf = df_mainland %>% left_join(turnover_bird_native_df,
                                                   by = 'RegionID')
colnames(turnover_bird_native_sf)

#### native birds mapping 
turnover_bird_native_sf = turnover_bird_native_sf %>% filter(!is.na(turnover))

turnover_bird_native_sf$f_turnover = cut(turnover_bird_native_sf$turnover,
                                         breaks = quantile(turnover_bird_native_sf$turnover,
                                                           probs = seq(0, 1,
                                                                       length.out = 9),
                                                           na.rm = TRUE),
                                         include.lowest = TRUE)


birds_native_turnover_map =
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
  geom_sf(data = (turnover_bird_native_sf %>% filter(!(Island == 1 & 
                                                         Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.birds_native_turnover_map = legend.func(mycolors = colors2,
                                               mylabels = round(quantile(turnover_bird_native_sf$turnover,
                                                                         probs = seq(0, 1,length.out = 9),
                                                                         na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_turnover_map = ggplotGrob(birds_native_turnover_map)
legend.birds_native_turnover_map = ggplotGrob(legend.birds_native_turnover_map)
birds_native_turnover_map_all = arrangeGrob(birds_native_turnover_map,
                                            legend.birds_native_turnover_map,
                                            ncol = 1,
                                            layout_matrix = rbind(matrix(1, 4, 10),
                                                                  c(NA, rep(2, 8), NA)))
plot(birds_native_turnover_map_all)




##4.3 extant #####
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

phy_turn_bird_extant = calcu_phy_turn_simple(Tree = phy_data,
                                               Comm = comm_bird_extant,
                                               Region_posi = which(colnames(comm_bird_extant) == 'RegionID'))

save(phy_turn_bird_extant,
     file = 'results/primary_results/distances_beta/phy_turn_bird_extant.rdata')

load("results/primary_results/distances_beta/phy_turn_bird_extant.rdata") # load the calculated data from the server

turnover_bird_extant_mat = phy_turn_bird_extant
rm(phy_turn_bird_extant)

tdwg_mainlands_all_bird = intersect(colnames(turnover_bird_extant_mat),
                                     tdwg_mainlands)
turnover_bird_extant_mat = turnover_bird_extant_mat[tdwg_mainlands_all_bird,
                                                    tdwg_mainlands_all_bird]

turnover_bird_extant_df = data.frame(RegionID = as.integer(colnames(turnover_bird_extant_mat)),
                                     turnover = colMeans(turnover_bird_extant_mat,
                                                         na.rm = T))

turnover_bird_extant_sf = df_mainland %>% left_join(turnover_bird_extant_df,
                                                    by = 'RegionID')
colnames(turnover_bird_extant_sf)

#### extant birds mapping 
turnover_bird_extant_sf = turnover_bird_extant_sf %>% filter(!is.na(turnover))

turnover_bird_extant_sf$f_turnover = cut(turnover_bird_extant_sf$turnover,
                                         breaks = quantile(turnover_bird_extant_sf$turnover,
                                                           probs = seq(0, 1,
                                                                       length.out = 9),
                                                           na.rm = TRUE),
                                         include.lowest = TRUE)


birds_extant_turnover_map =
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
  geom_sf(data = (turnover_bird_extant_sf %>% filter(!(Island == 1 & 
                                                         Area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover'
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

legend.birds_extant_turnover_map = legend.func(mycolors = colors2,
                                               mylabels = round(quantile(turnover_bird_extant_sf$turnover,
                                                                         probs = seq(0, 1,length.out = 9),
                                                                         na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_extant_turnover_map = ggplotGrob(birds_extant_turnover_map)
legend.birds_extant_turnover_map = ggplotGrob(legend.birds_extant_turnover_map)
birds_extant_turnover_map_all = arrangeGrob(birds_extant_turnover_map,
                                            legend.birds_extant_turnover_map,
                                            ncol = 1,
                                            layout_matrix = rbind(matrix(1, 4, 10),
                                                                  c(NA, rep(2, 8), NA)))
plot(birds_extant_turnover_map_all)





##4.4 Delta_turnover (Extant - Native) #####
load('results/primary_results/distances_beta/phy_turn_bird_extant.rdata')
load('results/primary_results/distances_beta/phy_turn_bird_native.rdata')

### alternative 1. detla = mean(extant - native)
turnover_bird_extant_mat = phy_turn_bird_extant
turnover_bird_native_mat = phy_turn_bird_native

turnover_bird_delta_mat = log((turnover_bird_extant_mat+0.001) /
                                  (turnover_bird_native_mat+0.001))

tdwg_mainlands_delta_bird = intersect(colnames(turnover_bird_delta_mat),
                                        tdwg_mainlands)
turnover_bird_delta_mat = turnover_bird_delta_mat[tdwg_mainlands_delta_bird,
                                                      tdwg_mainlands_delta_bird]

turnover_bird_delta = data.frame(RegionID = as.integer(colnames(turnover_bird_delta_mat)),
                                   delta_turnover = colMeans(turnover_bird_delta_mat,
                                                             na.rm = T))

turnover_bird_delta_sf = df_mainland %>% left_join(turnover_bird_delta,
                                                  by = 'RegionID')

colnames(turnover_bird_delta_sf)


#### Delta_ED birds mapping 
turnover_bird_delta_sf = turnover_bird_delta_sf %>% filter(!is.na(delta_turnover))

hist(turnover_bird_delta_sf$delta_turnover)

breaks_delta_turnover = quantile(turnover_bird_delta_sf$delta_turnover,
                                 probs = seq(0, 1,
                                             length.out = 9),
                                 na.rm = TRUE,
                                 names = T,
                                 digits = 12)
breaks_delta_turnover[which(names(breaks_delta_turnover) == '87.5%')] = 0
turnover_bird_delta_sf$f_delta_turnover = cut(turnover_bird_delta_sf$delta_turnover,
                                              breaks = breaks_delta_turnover,
                                              include.lowest = TRUE)


birds_delta_turnover_map =
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
  geom_sf(data = (turnover_bird_delta_sf %>% filter(!(Island == 1 & 
                                                        Area < 1e4))), aes(fill = f_delta_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  scale_color_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover'
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

legend.birds_delta_turnover_map = legend.func(mycolors = colors8,
                                              mylabels = round(breaks_delta_turnover,
                                                               4)) +
  ggtitle("log (turnover_all / turnover_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_delta_turnover_map = ggplotGrob(birds_delta_turnover_map)
legend.birds_delta_turnover_map = ggplotGrob(legend.birds_delta_turnover_map)
birds_delta_turnover_map_all = arrangeGrob(birds_delta_turnover_map,
                                           legend.birds_delta_turnover_map,
                                           ncol = 1,
                                           layout_matrix = rbind(matrix(1, 4, 10),
                                                                 c(NA, rep(2, 8), NA)))
plot(birds_delta_turnover_map_all)






## 4.5 Patitioning delta ED into 5 possible ways ----
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

region_pairs = combn(colnames(turnover_bird_delta_mat), 2)

turnover_bird_path3_5_mat = turnover_bird_native_mat
turnover_bird_path4_7_mat = turnover_bird_native_mat
turnover_bird_path6_mat = turnover_bird_native_mat

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
    
    phy_turn_bird_path6 = calcu_turnover_simple(
      Tree = spec_phy.3,
      Comm = comm_bird_path6[,2:ncol(comm_bird_path6)])
    
    turnover_bird_path6_mat[region1, region2] = phy_turn_bird_path6[1,2]
    turnover_bird_path6_mat[region2, region1] = phy_turn_bird_path6[1,2]
    
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
      
      phy_turn_bird_path3_5 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_bird_path3_5[,2:ncol(comm_bird_path3_5)])
      
      turnover_bird_path3_5_mat[region1, region2] = phy_turn_bird_path3_5[1,2]
      turnover_bird_path3_5_mat[region2, region1] = phy_turn_bird_path3_5[1,2]
      
      
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
      
      phy_turn_bird_path4_7 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_bird_path4_7[,2:ncol(comm_bird_path4_7)])
      
      turnover_bird_path4_7_mat[region1, region2] = phy_turn_bird_path4_7[1,2]
      turnover_bird_path4_7_mat[region2, region1] = phy_turn_bird_path4_7[1,2]
      
    }
    
  }
  
}


turnover_bird_delta_path6_mat = log((turnover_bird_path6_mat+0.001) /
                                      (turnover_bird_native_mat+0.001))

turnover_bird_delta_path4_7_mat = log((turnover_bird_path4_7_mat+0.001) /
                                        (turnover_bird_native_mat+0.001))
which(turnover_bird_delta_path4_7_mat > 0)
which(turnover_bird_delta_path4_7_mat < 0)

turnover_bird_delta_path3_5_mat = log((turnover_bird_path3_5_mat+0.001) /
                                        (turnover_bird_native_mat+0.001))
which(turnover_bird_delta_path3_5_mat > 0)
which(turnover_bird_delta_path3_5_mat < 0)





#5. Fish turnover: calculation & mapping ####
load("D:/R projects/Global_ED/data/Fishes/data/my_phy.rdata")
is.rooted(phylo)
load("D:/R projects/Global_ED/data/Fishes/data/my_data_used_final.rdata")

colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')


basin_trans = st_read("data/Fishes/data/Basin042017_3119_eck4.shp")
colnames(basin_trans)[which(colnames(basin_trans) == 'BasinName')] = 'X1.Basin.Name'
basin_mainland = basin_trans %>% filter(!(Island == 1 & 
                                Area < 5e3))

plot(basin_mainland$geometry)
plot(basin_trans$geometry)
basin_mainlands = as.character(sort(unique(basin_mainland$RegionID)))


##5.1 exotics #####
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

phy_turn_fish_exotic = calcu_phy_turn_simple(
  Tree = phylo,
  Comm = comm_fish_exotic[,unique(data.used_final_exotics$valid_names)],
  Area = comm_fish_exotic$Surf_area)


save(phy_turn_fish_exotic_sf, file = 'results/primary_results/phy_turn_fish_exotic.rdata')
load('results/primary_results/phy_turn_fish_exotic.rdata')
phy_turn_fish_exotic$df$delta_phy_turnR = phy_turn_fish_exotic$df$delta_phy_turnR * 1e+6
phy_turn_fish_exotic_1 = cbind(Basin.Name = sort(unique(df_trans$X1.Basin.Name)),
                               phy_turn_fish_exotic$df)
phy_turn_fish_exotic_sf = df_trans %>% left_join(phy_turn_fish_exotic_1,
                                                 by = join_by('X1.Basin.Name' == 'Basin.Name'))
colnames(phy_turn_fish_exotic_sf)

#### exotic fishes mapping 
phy_turn_fish_exotic_sf = phy_turn_fish_exotic_sf %>% filter(!is.na(delta_phy_turn))
phy_turn_fish_exotic_sf$log_EDR = log10(phy_turn_fish_exotic_sf$delta_phy_turnR)

phy_turn_fish_exotic_sf$f_delta_phy_turn = cut(phy_turn_fish_exotic_sf$delta_phy_turn,
                                               breaks = quantile(phy_turn_fish_exotic_sf$delta_phy_turn,
                                                                 probs = seq(0, 1,
                                                                             length.out = 9),
                                                                 na.rm = TRUE),
                                               include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(phy_turn_fish_exotic_sf$delta_phy_turnR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(phy_turn_fish_exotic_sf$delta_phy_turnR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(phy_turn_fish_exotic_sf$delta_phy_turnR, na.rm = TRUE) * 1.01

phy_turn_fish_exotic_sf$f_log10_delta_phy_turnR = cut(phy_turn_fish_exotic_sf$delta_phy_turnR,
                                                      breaks = breaks,
                                                      include.lowest = TRUE,
                                                      include.highest = T,
                                                      right = T,
                                                      dig.lab = 2)

# delta_phy_turnR = 2.352844e-07
phy_turn_fish_exotic_sf_2 = phy_turn_fish_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

phy_turn_fish_exotic_sf_2 %>% filter(is.na(f_log10_delta_phy_turnR))

fishes_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = phy_turn_fish_exotic_sf, aes(fill = f_delta_phy_turn),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'phy_turn',
    guide = guide_legend(
      title = 'phy_turn',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(phy_turn_fish_exotic_sf$delta_phy_turn,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'phy_turn'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(phy_turn_fish_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_phy_turn, color = f_delta_phy_turn),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'phy_turn'
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
                                          mylabels = round(quantile(phy_turn_fish_exotic_sf$delta_phy_turn,
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
  geom_sf(data = phy_turn_fish_exotic_sf, aes(fill = f_log10_delta_phy_turnR),
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
  geom_point(data = subset(phy_turn_fish_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_delta_phy_turnR, color = f_log10_delta_phy_turnR),
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





##5.2 natives #####
data.used_final_natives$presence = 1

data.used_final_natives_sub = data.used_final_natives %>% 
  filter(X1.Basin.Name %in% c("Aa",
                              "Abant.lake",
                              "Abashiri"))

comm_fish_native = data.used_final_natives_sub %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(basin_trans[,c('X1.Basin.Name', 'Surf_area')], 
            #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
            by = 'X1.Basin.Name') %>% 
  relocate(Surf_area, .after = X1.Basin.Name)%>% 
  arrange('X1.Basin.Name') 

phy_turn_fish_native = calcu_turnover_simple(Tree = phylo,
  Comm = comm_fish_native)

save(phy_turn_fish_native, file = 'results/primary_results/phy_turn_fish_native.rdata')

load('results/primary_results/phy_turn_fish_native.rdata')

turnover_fish_native_mat = phy_turn_fish_native
turnover_fish_native = data.frame(X1.Basin.Name = colnames(turnover_fish_native_mat),
                                  turnover = colMeans(turnover_fish_native_mat, 
                                                      na.rm = T))

turnover_fish_native_sf = basin_trans %>% left_join(turnover_fish_native,
                                                 by = join_by('BasinName' == 'X1.Basin.Name'))
colnames(turnover_fish_native_sf)

#### native fishes mapping 
turnover_fish_native_sf = turnover_fish_native_sf %>% filter(!is.na(turnover))

turnover_fish_native_sf$f_turnover = cut(turnover_fish_native_sf$turnover,
                                         breaks = quantile(turnover_fish_native_sf$turnover,
                                                           probs = seq(0, 1,
                                                                       length.out = 9),
                                                           na.rm = TRUE),
                                         include.lowest = TRUE)


fishes_native_turnover_map =
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
  geom_sf(data = (turnover_fish_native_sf %>% filter(!(Island == 1 & 
                                                         Surf_area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_fish_native_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Native_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_native_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_fish_native_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_native_turnover_map = ggplotGrob(fishes_native_turnover_map)
legend.fishes_native_turnover_map = ggplotGrob(legend.fishes_native_turnover_map)
fishes_native_turnover_map_all = arrangeGrob(fishes_native_turnover_map,
                                             legend.fishes_native_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(fishes_native_turnover_map_all)



##5.3 extant #####
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

phy_turn_fish_extant = calcu_phy_turn_simple(
  Tree = phylo,
  Comm = comm_fish_extant,
  Region_posi = which(colnames(comm_fish_extant) == 'X1.Basin.Name'))

save(phy_turn_fish_extant, file = 'results/primary_results/phy_turn_fish_extant.rdata')

load('results/primary_results/phy_turn_fish_extant.rdata')

turnover_fish_extant_mat = phy_turn_fish_extant
turnover_fish_extant = data.frame(X1.Basin.Name = colnames(turnover_fish_extant_mat),
                                  turnover = colMeans(turnover_fish_extant_mat, 
                                                      na.rm = T))

turnover_fish_extant_sf = basin_trans %>% left_join(turnover_fish_extant,
                                                 by = join_by('BasinName' == 'X1.Basin.Name'))
colnames(turnover_fish_extant_sf)

#### extant fishes mapping 
turnover_fish_extant_sf = turnover_fish_extant_sf %>% filter(!is.na(turnover))

turnover_fish_extant_sf$f_turnover = cut(turnover_fish_extant_sf$turnover,
                                         breaks = quantile(turnover_fish_extant_sf$turnover,
                                                           probs = seq(0, 1,
                                                                       length.out = 9),
                                                           na.rm = TRUE),
                                         include.lowest = TRUE)

fishes_extant_turnover_map =
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
  geom_sf(data = (turnover_fish_extant_sf %>% filter(!(Island == 1 & 
                                                         Surf_area < 1e4))),
          aes(fill = f_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(turnover_fish_extant_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Extant_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_extant_turnover_map = legend.func(mycolors = colors2,
                                                mylabels = round(quantile(turnover_fish_extant_sf$turnover,
                                                                          probs = seq(0, 1,length.out = 9),
                                                                          na.rm = TRUE), 5)) +
  ggtitle("Turnover") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_extant_turnover_map = ggplotGrob(fishes_extant_turnover_map)
legend.fishes_extant_turnover_map = ggplotGrob(legend.fishes_extant_turnover_map)
fishes_extant_turnover_map_all = arrangeGrob(fishes_extant_turnover_map,
                                             legend.fishes_extant_turnover_map,
                                             ncol = 1,
                                             layout_matrix = rbind(matrix(1, 4, 10),
                                                                   c(NA, rep(2, 8), NA)))
plot(fishes_extant_turnover_map_all)




##5.4 Delta_turnover (Extant - Native) #####
load('results/primary_results/phy_turn_fish_extant.rdata')
load('results/primary_results/phy_turn_fish_native.rdata')

### alternative 1. detla = mean(extant - native)

turnover_fish_extant_mat = phy_turn_fish_extant
turnover_fish_native_mat = phy_turn_fish_native
common_basins = which(colnames(turnover_fish_extant_mat) %in% intersect(colnames(turnover_fish_extant_mat),
                                                                        colnames(turnover_fish_native_mat)))
turnover_fish_extant_mat = turnover_fish_extant_mat[common_basins,common_basins]
turnover_fish_delta_mat = log((turnover_fish_extant_mat+0.001) /
                                (turnover_fish_native_mat+0.001))
turnover_fish_delta = data.frame(X1.Basin.Name = colnames(turnover_fish_delta_mat),
                                 delta_turnover = colMeans(turnover_fish_delta_mat, 
                                                           na.rm = T))

turnover_fish_delta_sf = basin_trans %>% left_join(turnover_fish_delta,
                                                   by = join_by('BasinName' == 'X1.Basin.Name'))
colnames(turnover_fish_delta_sf)


#### Delta_ED fishes mapping 
turnover_fish_delta_sf = turnover_fish_delta_sf %>% filter(!is.na(delta_turnover) & 
                                                             delta_turnover != 0 
)

hist(turnover_fish_delta_sf$delta_turnover)

breaks_delta_turnover = quantile(turnover_fish_delta_sf$delta_turnover,
                                 probs = seq(0, 1,
                                             length.out = 9),
                                 na.rm = TRUE,
                                 names = T,
                                 digits = 12)
breaks_delta_turnover[which(names(breaks_delta_turnover) == '87.5%')] = 0
turnover_fish_delta_sf$f_delta_turnover = cut(turnover_fish_delta_sf$delta_turnover,
                                              breaks = breaks_delta_turnover,
                                              include.lowest = TRUE)


fishes_delta_turnover_map =
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
  geom_sf(data = (turnover_fish_delta_sf %>% filter(!(Island == 1 & 
                                                        Surf_area < 1e4))),
          aes(fill = f_delta_turnover),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors8,
    na.value = 'gray90',
    name = 'turnover',
    guide = guide_legend(
      title = 'turnover',
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors8,
  # values = scales::rescale(quantile(turnover_fish_delta_sf$turnover,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'turnover'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Delta_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_delta_turnover_map = legend.func(mycolors = colors8,
                                               mylabels = round(breaks_delta_turnover,
                                                                4)) +
  ggtitle("log (Turnover_all / Turnover_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_delta_turnover_map = ggplotGrob(fishes_delta_turnover_map)
legend.fishes_delta_turnover_map = ggplotGrob(legend.fishes_delta_turnover_map)
fishes_delta_turnover_map_all = arrangeGrob(fishes_delta_turnover_map,
                                            legend.fishes_delta_turnover_map,
                                            ncol = 1,
                                            layout_matrix = rbind(matrix(1, 4, 10),
                                                                  c(NA, rep(2, 8), NA)))
plot(fishes_delta_turnover_map_all)





## 5.5 Patitioning delta ED into 5 possible ways ----
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

region_pairs = combn(colnames(turnover_fish_delta_mat), 2)

turnover_fish_path3_5_mat = turnover_fish_native_mat
turnover_fish_path4_7_mat = turnover_fish_native_mat
turnover_fish_path6_mat = turnover_fish_native_mat

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
    
    phy_turn_fish_path6 = calcu_turnover_simple(
      Tree = spec_phy.3,
      Comm = comm_fish_path6[,2:ncol(comm_fish_path6)])
    
    turnover_fish_path6_mat[region1, region2] = phy_turn_fish_path6[1,2]
    turnover_fish_path6_mat[region2, region1] = phy_turn_fish_path6[1,2]
    
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
      
      phy_turn_fish_path3_5 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_fish_path3_5[,2:ncol(comm_fish_path3_5)])
      
      turnover_fish_path3_5_mat[region1, region2] = phy_turn_fish_path3_5[1,2]
      turnover_fish_path3_5_mat[region2, region1] = phy_turn_fish_path3_5[1,2]
      
      
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
      
      phy_turn_fish_path4_7 = calcu_turnover_simple(
        Tree = spec_phy.3,
        Comm = comm_fish_path4_7[,2:ncol(comm_fish_path4_7)])
      
      turnover_fish_path4_7_mat[region1, region2] = phy_turn_fish_path4_7[1,2]
      turnover_fish_path4_7_mat[region2, region1] = phy_turn_fish_path4_7[1,2]
      
    }
    
  }
  
}


turnover_fish_delta_path6_mat = log((turnover_fish_path6_mat+0.001) /
                                      (turnover_fish_native_mat+0.001))

turnover_fish_delta_path4_7_mat = log((turnover_fish_path4_7_mat+0.001) /
                                        (turnover_fish_native_mat+0.001))
which(turnover_fish_delta_path4_7_mat > 0)
which(turnover_fish_delta_path4_7_mat < 0)

turnover_fish_delta_path3_5_mat = log((turnover_fish_path3_5_mat+0.001) /
                                        (turnover_fish_native_mat+0.001))
which(turnover_fish_delta_path3_5_mat > 0)
which(turnover_fish_delta_path3_5_mat < 0)




#6. Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

##6.1 turnover #####
# Compare the turnover patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.01         # small horizontal gap between columns

figs_turnover_nat_extant_delta_mainland = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_turnover_map_all,  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_turnover_map_all,  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_turnover_map_all,   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_turnover_map_all,   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_extant_turnover_map_all,   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_turnover_map_all,    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_turnover_map_all, x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_all_turnover_map_all, x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_turnover_map_all,  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_turnover_map_all,   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_turnover_map_all,   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_turnover_map_all,   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
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

png(filename = 'figures/figs_turnover_nat_extant_delta_mainland.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_turnover_nat_extant_delta_mainland
dev.off() #turn off device and finalize file


### export for PPT

figs_turnover_delta = ggarrange(plants_delta_turnover_map_all, 
                                birds_delta_turnover_map_all,
                                mammals_delta_turnover_map_all,
                                fishes_delta_turnover_map_all,
                                nrow = 2, ncol = 2, 
                                labels = c('a', 'b', 'c', 'd'))


png(filename = 'figures/figs_turnover_delta.png',
    height=15, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_turnover_delta
dev.off() #turn off device and finalize file



plot_width = 0.5
plot_height = 0.20

figs_turnover_nat_exo = ggdraw() +
  draw_plot(plants_native_turnover_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_turnover_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_turnover_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_turnover_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_turnover_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_turnover_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_turnover_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_exotic_turnover_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_turnover_nat_exo.png',
    #plot = figs_turnover_nat_exo,
    height=40, width=25, # setting emfPlusFontToPath=TRUE to 
    res = 300,
    # ensure text looks correct on the viewing system
    units = 'cm')
figs_turnover_nat_exo
dev.off() #turn off device and finalize file


#emf('figures/figs_turnover_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_turnover_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_turnover_nat_extant = ggdraw() +
  draw_plot(plants_native_turnover_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_turnover_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_turnover_map_all, x = 0, y = 0.53, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_turnover_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_turnover_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_turnover_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_turnover_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_turnover_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_turnover_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_turnover_nat_extant
dev.off() #turn off device and finalize file


