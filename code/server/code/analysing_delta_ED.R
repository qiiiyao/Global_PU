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
source('code/calculating_ED_func.R')

# load the background data for plotting the world map plot
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
load("code/FYI/Codes_and_Data_Fan_et_al_2023/data.for.shp.plots.4.Rdata")


# define color gradients
scico::scico_palette_show()

colors1 = rev(scico::scico(n=8, palette = "lajolla"))
colors2 = scico::scico(n=10, palette = "vik")[1:8] 
colors3 = scico::scico(n=8, begin = 0, end = 0.4, palette = "bam")
colors4 = scico::scico(n=8, palette = "bam")
colors5 = scico::scico(n=7, begin = 0, end = 0.4, palette = "bam")
colors5[length(colors5)+1] = colors4[length(colors4)-2]

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


#### Mammal ED: calculation & mapping ####
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")
phy_mammal = spec_phy.3
initial_div_t_mammal = max(phyloregion::evol_distinct(tree = phy_mammal,
                                                      type = "fair.proportion"))

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)
sp_overlap_dat_1$presence = 1

comm_mammal_exotic = sp_overlap_dat_1 %>% 
  complete(Binomial, RegionID, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(names_from = Binomial,
              values_from = presence,
              values_fn = mean) %>% 
  complete(RegionID = unique(df$RegionID)) %>% ## assume that 
  ## the absence regions have no naturalized aliens
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
  )  %>% left_join(df_trans[,c('RegionID', 'area')],
                   by = 'RegionID') %>% 
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID')

ED_mammal_exotic = calcu_mean_ED_parallel(
  Tree = spec_phy.3,
  Comm = comm_mammal_exotic[,unique(sp_overlap_dat_1$Binomial)],
  Area = comm_mammal_exotic$area)

save(ED_mammal_exotic,
     file = 'results/primary_results/ED_mammal_exotic.rdata')

load('results/primary_results/ED_mammal_exotic.rdata')
# convert the unit of EDR to Mya 10^4 KM^(-2)
ED_mammal_exotic$df$mean_EDR = ED_mammal_exotic$df$mean_EDR * 1e+06 * 10000 
ED_mammal_exotic_1 = cbind(RegionID = df$RegionID,
                           ED_mammal_exotic$df)
ED_mammal_exotic_sf = df_trans %>% left_join(ED_mammal_exotic_1, by = 'RegionID')
colnames(ED_mammal_exotic_sf)

#### exotic mammals mapping 
ED_mammal_exotic_sf = ED_mammal_exotic_sf %>% filter(!is.na(mean_ED))
ED_mammal_exotic_sf$log_EDR = log10(ED_mammal_exotic_sf$mean_EDR)

ED_mammal_exotic_sf$f_mean_ED = cut(ED_mammal_exotic_sf$mean_ED,
                                    breaks = quantile(ED_mammal_exotic_sf$mean_ED,
                                                      probs = seq(0, 1,
                                                                  length.out = 9),
                                                      na.rm = TRUE),
                                    include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_mammal_exotic_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_mammal_exotic_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_mammal_exotic_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_mammal_exotic_sf$f_log10_mean_EDR = cut(ED_mammal_exotic_sf$mean_EDR,
                                           breaks = breaks,
                                           include.lowest = TRUE,
                                           include.highest = T,
                                           right = T,
                                           dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_mammal_exotic_sf_2 = ED_mammal_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_mammal_exotic_sf_2 %>% filter(is.na(f_log10_mean_EDR))

mammals_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_mammal_exotic_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_mammal_exotic_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_mammal_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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
    plot.margin = margin(0, 0, 0, 0)
  ) 

legend.mammals_exotic_ED_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(ED_mammal_exotic_sf$mean_ED,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = .5, vjust = -22, size = 9))
mammals_exotic_ED_map = ggplotGrob(mammals_exotic_ED_map)
legend.mammals_exotic_ED_map = ggplotGrob(legend.mammals_exotic_ED_map)
mammals_exotic_ED_map_all = arrangeGrob(mammals_exotic_ED_map,
                                        legend.mammals_exotic_ED_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA))
)
plot(mammals_exotic_ED_map_all)

mammals_exotic_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_mammal_exotic_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_mammal_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Exotic_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_exotic_EDR_map = legend.func(mycolors = colors2,
                                            mylabels = round(breaks, 6)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_exotic_EDR_map = ggplotGrob(mammals_exotic_EDR_map)
legend.mammals_exotic_EDR_map = ggplotGrob(legend.mammals_exotic_EDR_map)
mammals_exotic_EDR_map_all = arrangeGrob(mammals_exotic_EDR_map,
                                         legend.mammals_exotic_EDR_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_exotic_EDR_map_all)



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
  arrange('Region.ID') 

ED_mammal_native = calcu_mean_ED_parallel(Tree = spec_phy.3,
                                          Comm = as.matrix(comm_mammal_native[,unique(sp_dis_5$ScientificName)]),
                                          Area = comm_mammal_native$area)

save(ED_mammal_native,
     file = 'results/primary_results/ED_mammal_native.rdata')

load('results/primary_results/ED_mammal_native.rdata')
ED_mammal_native$df$mean_EDR = ED_mammal_native$df$mean_EDR * 1e+06 * 10000 
ED_mammal_native_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_native$df)
ED_mammal_native_sf = df_trans %>% left_join(ED_mammal_native_1, by = 'RegionID')
colnames(ED_mammal_native_sf)

#### native mammals mapping 
ED_mammal_native_sf = ED_mammal_native_sf %>% filter(!is.na(mean_ED))
ED_mammal_native_sf$log_EDR = log10(ED_mammal_native_sf$mean_EDR)

ED_mammal_native_sf$f_mean_ED = cut(ED_mammal_native_sf$mean_ED,
                                    breaks = quantile(ED_mammal_native_sf$mean_ED,
                                                      probs = seq(0, 1,
                                                                  length.out = 9),
                                                      na.rm = TRUE),
                                    include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_mammal_native_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_mammal_native_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_mammal_native_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_mammal_native_sf$f_log10_mean_EDR = cut(ED_mammal_native_sf$mean_EDR,
                                           breaks = breaks,
                                           include.lowest = TRUE,
                                           include.highest = T,
                                           right = T,
                                           dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_mammal_native_sf_2 = ED_mammal_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_mammal_native_sf_2 %>% filter(is.na(f_log10_mean_EDR))

mammals_native_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_mammal_native_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_mammal_native_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_mammal_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.mammals_native_ED_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(ED_mammal_native_sf$mean_ED,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_ED_map = ggplotGrob(mammals_native_ED_map)
legend.mammals_native_ED_map = ggplotGrob(legend.mammals_native_ED_map)
mammals_native_ED_map_all = arrangeGrob(mammals_native_ED_map,
                                        legend.mammals_native_ED_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(mammals_native_ED_map_all)

mammals_native_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_mammal_native_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_mammal_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Native_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_native_EDR_map = legend.func(mycolors = colors2,
                                            mylabels = round(breaks, 4)) +
  ggtitle("MEDR (Mya 10^4 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_EDR_map = ggplotGrob(mammals_native_EDR_map)
legend.mammals_native_EDR_map = ggplotGrob(legend.mammals_native_EDR_map)
mammals_native_EDR_map_all = arrangeGrob(mammals_native_EDR_map,
                                         legend.mammals_native_EDR_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_native_EDR_map_all)


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
  arrange('RegionID') 

ED_mammal_all = calcu_mean_ED_parallel(
  Tree = spec_phy.3,
  Comm = comm_mammal_all[,unique(sp_mammal_all$Binomial)],
  Area = comm_mammal_all$area)

save(ED_mammal_all,
     file = 'results/primary_results/ED_mammal_all.rdata')

load('results/primary_results/ED_mammal_all.rdata')
ED_mammal_all$df$mean_EDR = ED_mammal_all$df$mean_EDR * 1e+06 * 10000 
ED_mammal_extant_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_all$df)
ED_mammal_extant_sf = df_trans %>% left_join(ED_mammal_extant_1, by = 'RegionID')
colnames(ED_mammal_extant_sf)

#### all mammals mapping 
ED_mammal_extant_sf = ED_mammal_extant_sf %>% filter(!is.na(mean_ED))
ED_mammal_extant_sf$log_EDR = log10(ED_mammal_extant_sf$mean_EDR)

ED_mammal_extant_sf$f_mean_ED = cut(ED_mammal_extant_sf$mean_ED,
                                    breaks = quantile(ED_mammal_extant_sf$mean_ED,
                                                      probs = seq(0, 1,
                                                                  length.out = 9),
                                                      na.rm = TRUE),
                                    include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_mammal_extant_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_mammal_extant_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_mammal_extant_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_mammal_extant_sf$f_log10_mean_EDR = cut(ED_mammal_extant_sf$mean_EDR,
                                           breaks = breaks,
                                           include.lowest = TRUE,
                                           include.highest = T,
                                           right = T,
                                           dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_mammal_extant_sf_2 = ED_mammal_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_mammal_extant_sf_2 %>% filter(is.na(f_log10_mean_EDR))

mammals_extant_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_mammal_extant_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_mammal_extant_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_mammal_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.mammals_extant_ED_map = legend.func(mycolors = colors2,
                                           mylabels = round(quantile(ED_mammal_extant_sf$mean_ED,
                                                                     probs = seq(0, 1,length.out = 9),
                                                                     na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_extant_ED_map = ggplotGrob(mammals_extant_ED_map)
legend.mammals_extant_ED_map = ggplotGrob(legend.mammals_extant_ED_map)
mammals_extant_ED_map_all = arrangeGrob(mammals_extant_ED_map,
                                        legend.mammals_extant_ED_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10),
                                                              c(NA, rep(2, 8), NA)))
plot(mammals_extant_ED_map_all)

mammals_extant_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_mammal_extant_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_mammal_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Extant_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_extant_EDR_map = legend.func(mycolors = colors2,
                                            mylabels = round(breaks, 4)) +
  ggtitle("MEDR (Mya 10^4 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_extant_EDR_map = ggplotGrob(mammals_extant_EDR_map)
legend.mammals_extant_EDR_map = ggplotGrob(legend.mammals_extant_EDR_map)
mammals_extant_EDR_map_all = arrangeGrob(mammals_extant_EDR_map,
                                         legend.mammals_extant_EDR_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_extant_EDR_map_all)



##### Delta_ED (Extant - Native) #####
load('results/primary_results/ED_mammal_all.rdata')
load('results/primary_results/ED_mammal_native.rdata')

ED_mammal_all$df$mean_EDR = ED_mammal_all$df$mean_EDR * 1e+06 * 100 
ED_mammal_extant_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_all$df)
colnames(ED_mammal_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_mammal_native$df$mean_EDR = ED_mammal_native$df$mean_EDR * 1e+06 * 10000 
ED_mammal_native_1 = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                           ED_mammal_native$df)
colnames(ED_mammal_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_mammal = ED_mammal_native_1 %>% left_join(ED_mammal_extant_1,
                                                   by = 'RegionID')


Delta_ED_mammal_sf = df_trans %>% left_join(Delta_ED_mammal, by = 'RegionID')
colnames(Delta_ED_mammal_sf)
Delta_ED_mammal_sf$delta_mean_ED = log(Delta_ED_mammal_sf$mean_all_ED / 
                                         Delta_ED_mammal_sf$mean_native_ED)

Delta_ED_mammal_sf$delta_mean_EDR = log(Delta_ED_mammal_sf$mean_all_EDR / 
                                          Delta_ED_mammal_sf$mean_native_EDR)

#### Delta_ED mammals mapping 
Delta_ED_mammal_sf = Delta_ED_mammal_sf %>% filter(!is.na(delta_mean_ED) & 
                                                     delta_mean_ED != 0 
)

hist(Delta_ED_mammal_sf$delta_mean_ED)

breaks_delta_ED = quantile(Delta_ED_mammal_sf$delta_mean_ED,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)
breaks_delta_ED[which(names(breaks_delta_ED) == '50%')] = 0
Delta_ED_mammal_sf$f_delta_mean_ED = cut(Delta_ED_mammal_sf$delta_mean_ED,
                                         breaks = breaks_delta_ED,
                                         include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
hist(Delta_ED_mammal_sf$delta_mean_EDR)
Delta_ED_mammal_sf$f_delta_mean_EDR = cut(Delta_ED_mammal_sf$delta_mean_EDR,
                                          breaks = quantile(Delta_ED_mammal_sf$delta_mean_EDR,
                                                            probs = seq(0, 1,
                                                                        length.out = 9),
                                                            na.rm = TRUE,
                                                            names = T,
                                                            digits = 12),
                                          include.lowest = TRUE,
                                          include.highest = T,
                                          right = T,
                                          dig.lab = 2)

# delta_mean_EDR = 2.352844e-07
Delta_ED_mammal_sf_2 = Delta_ED_mammal_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

Delta_ED_mammal_sf_2 = Delta_ED_mammal_sf_2 %>% filter(!is.na(delta_mean_EDR) & 
                                                         delta_mean_EDR != 0 
)

mammals_delta_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_mammal_sf, aes(fill = f_delta_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors4,
  # values = scales::rescale(quantile(Delta_ED_mammal_sf$delta_mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(Delta_ED_mammal_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_ED, color = f_delta_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors4,
    na.value = 'gray90',
    name = 'ED'
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

legend.mammals_delta_ED_map = legend.func(mycolors = colors4,
                                          mylabels = round(breaks_delta_ED, 4)) +
  ggtitle("log (MED_all / MED_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_delta_ED_map = ggplotGrob(mammals_delta_ED_map)
legend.mammals_delta_ED_map = ggplotGrob(legend.mammals_delta_ED_map)
mammals_delta_ED_map_all = arrangeGrob(mammals_delta_ED_map,
                                       legend.mammals_delta_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(mammals_delta_ED_map_all)

mammals_delta_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_mammal_sf, aes(fill = f_delta_mean_EDR),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors3,
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
  geom_point(data = subset(Delta_ED_mammal_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_EDR, color = f_delta_mean_EDR),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors3,
    na.value = 'gray90',
    name = 'EDR'
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

legend.mammals_delta_EDR_map = legend.func(mycolors = colors3,
                                           mylabels = round(quantile(Delta_ED_mammal_sf$delta_mean_EDR,
                                                                     probs = seq(0, 1,
                                                                                 length.out = 9),
                                                                     na.rm = TRUE,
                                                                     names = T,
                                                                     digits = 12), 2)) +
  ggtitle("log (MEDR_all / MEDR_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_delta_EDR_map = ggplotGrob(mammals_delta_EDR_map)
legend.mammals_delta_EDR_map = ggplotGrob(legend.mammals_delta_EDR_map)
mammals_delta_EDR_map_all = arrangeGrob(mammals_delta_EDR_map,
                                        legend.mammals_delta_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_delta_EDR_map_all)




#### Plant ED: calculation & mapping ####
load("D:/R projects/Global_ED/data/Plants/data/shp.651.Rdata")
load("D:/R projects/Global_ED/data/Plants/data/phylo.fake.species.653.Rdata")
phy_plant = phylo
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)

initial_div_t_plant = max(phyloregion::evol_distinct(tree = phy_plant,
                                                     type = "fair.proportion"))

colnames(df.native.650) == colnames(df.natu.650)

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
colnames(df.natu.650)
df.natu.650$presence = 1
comm_plant_exotic = df.natu.650 %>% 
  complete(species, Region_id, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

ED_plant_exotic = calcu_mean_ED_parallel(Tree = phylo,
                                         Comm = comm_plant_exotic[,unique(df.natu.650$species)],
                                         Area = comm_plant_exotic$area)
save(ED_plant_exotic,
     file = 'results/primary_results/ED_plant_exotic.rdata')

load("results/primary_results/ED_plant_exotic.rdata")
ED_plant_exotic$df$mean_EDR = ED_plant_exotic$df$mean_EDR * 1e+06 * 100 
ED_plant_exotic_1 = cbind(Region_id = sort(unique(df.natu.650$Region_id)),
                          ED_plant_exotic$df)
ED_plant_exotic_sf = shp.glonaf.trans %>% left_join(ED_plant_exotic_1, by = 'Region_id')
colnames(ED_plant_exotic_sf)

#### exotic plants mapping 
ED_plant_exotic_sf = ED_plant_exotic_sf %>% filter(!is.na(mean_ED))
ED_plant_exotic_sf$log_EDR = log10(ED_plant_exotic_sf$mean_EDR)

ED_plant_exotic_sf$f_mean_ED = cut(ED_plant_exotic_sf$mean_ED,
                                   breaks = quantile(ED_plant_exotic_sf$mean_ED,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_plant_exotic_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_plant_exotic_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_plant_exotic_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_plant_exotic_sf$f_log10_mean_EDR = cut(ED_plant_exotic_sf$mean_EDR,
                                          breaks = breaks,
                                          include.lowest = TRUE,
                                          include.highest = T,
                                          right = T,
                                          dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_plant_exotic_sf_2 = ED_plant_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_plant_exotic_sf_2 %>% filter(is.na(f_log10_mean_EDR))

plants_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_plant_exotic_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_plant_exotic_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_plant_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.plants_exotic_ED_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(ED_plant_exotic_sf$mean_ED,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_ED_map = ggplotGrob(plants_exotic_ED_map)
legend.plants_exotic_ED_map = ggplotGrob(legend.plants_exotic_ED_map)
plants_exotic_ED_map_all = arrangeGrob(plants_exotic_ED_map,
                                       legend.plants_exotic_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(plants_exotic_ED_map_all)

plants_exotic_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_plant_exotic_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_plant_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Exotic_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_exotic_EDR_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 3)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_EDR_map = ggplotGrob(plants_exotic_EDR_map)
legend.plants_exotic_EDR_map = ggplotGrob(legend.plants_exotic_EDR_map)
plants_exotic_EDR_map_all = arrangeGrob(plants_exotic_EDR_map,
                                        legend.plants_exotic_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_exotic_EDR_map_all)


##### natives #####
colnames(df.native.650)
df.native.650$presence = 1
comm_plant_native = df.native.650 %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

ED_plant_native = calcu_mean_ED_parallel(Tree = phylo,
                                         Comm = as.matrix(comm_plant_native[,unique(df.native.650$species)]),
                                         Area = comm_plant_native$area)

save(ED_plant_native,
     file = 'results/primary_results/ED_plant_native.rdata')

load("results/primary_results/ED_plant_native.rdata") # load the calculated data from the server
ED_plant_native$df$mean_EDR = ED_plant_native$df$mean_EDR * 1e+06 * 100 
ED_plant_native_1 = cbind(Region_id = sort(unique(df.native.650$Region_id)),
                          ED_plant_native$df)
ED_plant_native_sf = shp.glonaf.trans %>% left_join(ED_plant_native_1, by = 'Region_id')
colnames(ED_plant_native_sf)

#### native plants mapping 
ED_plant_native_sf = ED_plant_native_sf %>% filter(!is.na(mean_ED))
ED_plant_native_sf$log_EDR = log10(ED_plant_native_sf$mean_EDR)

ED_plant_native_sf$f_mean_ED = cut(ED_plant_native_sf$mean_ED,
                                   breaks = quantile(ED_plant_native_sf$mean_ED,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_plant_native_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_plant_native_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_plant_native_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_plant_native_sf$f_log10_mean_EDR = cut(ED_plant_native_sf$mean_EDR,
                                          breaks = breaks,
                                          include.lowest = TRUE,
                                          include.highest = T,
                                          right = T,
                                          dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_plant_native_sf_2 = ED_plant_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_plant_native_sf_2 %>% filter(is.na(f_log10_mean_EDR))

plants_native_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_plant_native_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_plant_native_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_plant_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.plants_native_ED_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(ED_plant_native_sf$mean_ED,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_ED_map = ggplotGrob(plants_native_ED_map)
legend.plants_native_ED_map = ggplotGrob(legend.plants_native_ED_map)
plants_native_ED_map_all = arrangeGrob(plants_native_ED_map,
                                       legend.plants_native_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(plants_native_ED_map_all)

plants_native_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_plant_native_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_plant_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Native_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_native_EDR_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 4)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_EDR_map = ggplotGrob(plants_native_EDR_map)
legend.plants_native_EDR_map = ggplotGrob(legend.plants_native_EDR_map)
plants_native_EDR_map_all = arrangeGrob(plants_native_EDR_map,
                                        legend.plants_native_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_native_EDR_map_all)


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

ED_plant_extant = calcu_mean_ED_parallel(Tree = phylo,
                                         Comm = as.matrix(comm_plant_extant[,unique(df.extant.650$species)]),
                                         Area = comm_plant_extant$area)
save(ED_plant_extant,
     file = 'results/primary_results/ED_plant_extant.rdata')

load("results/primary_results/ED_plant_extant.rdata") # load the calculated data from the server
ED_plant_extant$df$mean_EDR = ED_plant_extant$df$mean_EDR * 1e+06 * 100 
ED_plant_extant_1 = cbind(Region_id = sort(unique(df.extant.650$Region_id)),
                          ED_plant_extant$df)
ED_plant_extant_sf = shp.glonaf.trans %>% left_join(ED_plant_extant_1, by = 'Region_id')
colnames(ED_plant_extant_sf)

#### extant plants mapping 
ED_plant_extant_sf = ED_plant_extant_sf %>% filter(!is.na(mean_ED))
ED_plant_extant_sf$log_EDR = log10(ED_plant_extant_sf$mean_EDR)

ED_plant_extant_sf$f_mean_ED = cut(ED_plant_extant_sf$mean_ED,
                                   breaks = quantile(ED_plant_extant_sf$mean_ED,
                                                     probs = seq(0, 1,
                                                                 length.out = 9),
                                                     na.rm = TRUE),
                                   include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_plant_extant_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_plant_extant_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_plant_extant_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_plant_extant_sf$f_log10_mean_EDR = cut(ED_plant_extant_sf$mean_EDR,
                                          breaks = breaks,
                                          include.lowest = TRUE,
                                          include.highest = T,
                                          right = T,
                                          dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_plant_extant_sf_2 = ED_plant_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_plant_extant_sf_2 %>% filter(is.na(f_log10_mean_EDR))

plants_extant_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_plant_extant_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_plant_extant_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_plant_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.plants_extant_ED_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(ED_plant_extant_sf$mean_ED,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_ED_map = ggplotGrob(plants_extant_ED_map)
legend.plants_extant_ED_map = ggplotGrob(legend.plants_extant_ED_map)
plants_extant_ED_map_all = arrangeGrob(plants_extant_ED_map,
                                       legend.plants_extant_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(plants_extant_ED_map_all)

plants_extant_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_plant_extant_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_plant_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Extant_plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.plants_extant_EDR_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 5)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_EDR_map = ggplotGrob(plants_extant_EDR_map)
legend.plants_extant_EDR_map = ggplotGrob(legend.plants_extant_EDR_map)
plants_extant_EDR_map_all = arrangeGrob(plants_extant_EDR_map,
                                        legend.plants_extant_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_extant_EDR_map_all)



##### Delta_ED (Extant - Native) #####
load('results/primary_results/ED_plant_extant.rdata')
load('results/primary_results/ED_plant_native.rdata')
df.extant.650 = rbind(df.native.650, df.natu.650)

ED_plant_extant$df$mean_EDR = ED_plant_extant$df$mean_EDR * 1e+06 * 100 
ED_plant_extant_1 = cbind(Region_id = sort(unique(df.extant.650$Region_id)),
                          ED_plant_extant$df)
colnames(ED_plant_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")


ED_plant_native$df$mean_EDR = ED_plant_native$df$mean_EDR * 1e+06 * 100 
ED_plant_native_1 = cbind(Region_id = sort(unique(df.native.650$Region_id)),
                          ED_plant_native$df)
colnames(ED_plant_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_plant = ED_plant_native_1 %>% left_join(ED_plant_extant_1,
                                                 by = 'Region_id')


Delta_ED_plant_sf = shp.glonaf.trans %>% left_join(Delta_ED_plant, by = 'Region_id')
colnames(Delta_ED_plant_sf)
Delta_ED_plant_sf$delta_mean_ED = log(Delta_ED_plant_sf$mean_all_ED / 
                                        Delta_ED_plant_sf$mean_native_ED)

Delta_ED_plant_sf$delta_mean_EDR = log(Delta_ED_plant_sf$mean_all_EDR / 
                                         Delta_ED_plant_sf$mean_native_EDR)

#### Delta_ED plants mapping 
Delta_ED_plant_sf = Delta_ED_plant_sf %>% filter(!is.na(delta_mean_ED) & 
                                                   delta_mean_ED != 0 
)

hist(Delta_ED_plant_sf$delta_mean_ED)

breaks_delta_ED = quantile(Delta_ED_plant_sf$delta_mean_ED,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)
breaks_delta_ED[which(names(breaks_delta_ED) == '87.5%')] = 0
Delta_ED_plant_sf$f_delta_mean_ED = cut(Delta_ED_plant_sf$delta_mean_ED,
                                        breaks = breaks_delta_ED,
                                        include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
hist(Delta_ED_plant_sf$delta_mean_EDR)

breaks_delta_EDR = quantile(Delta_ED_plant_sf$delta_mean_EDR,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_EDR[which(names(breaks_delta_EDR) == '87.5%')] = 0
Delta_ED_plant_sf$f_delta_mean_EDR = cut(Delta_ED_plant_sf$delta_mean_EDR,
                                         breaks = breaks_delta_EDR,
                                         include.lowest = T,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# delta_mean_EDR = 2.352844e-07
Delta_ED_plant_sf_2 = Delta_ED_plant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

Delta_ED_plant_sf_2 = Delta_ED_plant_sf_2 %>% filter(!is.na(delta_mean_EDR) & 
                                                       delta_mean_EDR != 0 
)

plants_delta_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_plant_sf, aes(fill = f_delta_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors5,
  # values = scales::rescale(quantile(Delta_ED_plant_sf$delta_mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(Delta_ED_plant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_ED, color = f_delta_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'ED'
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

legend.plants_delta_ED_map = legend.func(mycolors = colors5,
                                         mylabels = round(breaks_delta_ED, 4)) +
  ggtitle("log (MED_all / MED_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_ED_map = ggplotGrob(plants_delta_ED_map)
legend.plants_delta_ED_map = ggplotGrob(legend.plants_delta_ED_map)
plants_delta_ED_map_all = arrangeGrob(plants_delta_ED_map,
                                      legend.plants_delta_ED_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(plants_delta_ED_map_all)

plants_delta_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_plant_sf, aes(fill = f_delta_mean_EDR),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
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
  geom_point(data = subset(Delta_ED_plant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_EDR, color = f_delta_mean_EDR),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'EDR'
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

legend.plants_delta_EDR_map = legend.func(mycolors = colors5,
                                          mylabels = round(breaks_delta_EDR, 2)) +
  ggtitle("log (MEDR_all / MEDR_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_delta_EDR_map = ggplotGrob(plants_delta_EDR_map)
legend.plants_delta_EDR_map = ggplotGrob(legend.plants_delta_EDR_map)
plants_delta_EDR_map_all = arrangeGrob(plants_delta_EDR_map,
                                       legend.plants_delta_EDR_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_delta_EDR_map_all)





#### Bird ED: calculation & mapping ####
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
phy_bird = read.tree("data/Birds/data/Phylogenetic_Birds.tre")

initial_div_t_bird = max(phyloregion::evol_distinct(tree = phy_bird,
                                                    type = "fair.proportion"))

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

ED_bird_exotic = calcu_mean_ED_parallel(Tree = phy_data,
                                        Comm = comm_bird_exotic[,unique(exotic_distri_data$ScientificName)],
                                        Area = comm_bird_exotic$area)
save(ED_bird_exotic,
     file = 'results/primary_results/ED_bird_exotic.rdata')
load('results/primary_results/ED_bird_exotic.rdata')
ED_bird_exotic$df$mean_EDR = ED_bird_exotic$df$mean_EDR * 1e+06 * 100 
ED_bird_exotic_1 = cbind(RegionID = df$RegionID,
                         ED_bird_exotic$df)
ED_bird_exotic_sf = df.trans %>% left_join(ED_bird_exotic_1, by = 'RegionID')
colnames(ED_bird_exotic_sf)

#### exotic birds mapping 
ED_bird_exotic_sf = ED_bird_exotic_sf %>% filter(!is.na(mean_ED))
ED_bird_exotic_sf$log_EDR = log10(ED_bird_exotic_sf$mean_EDR)

ED_bird_exotic_sf$f_mean_ED = cut(ED_bird_exotic_sf$mean_ED,
                                  breaks = quantile(ED_bird_exotic_sf$mean_ED,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_bird_exotic_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_bird_exotic_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_bird_exotic_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_bird_exotic_sf$f_log10_mean_EDR = cut(ED_bird_exotic_sf$mean_EDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_bird_exotic_sf_2 = ED_bird_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_bird_exotic_sf_2 %>% filter(is.na(f_log10_mean_EDR))

birds_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_bird_exotic_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_bird_exotic_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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
                                         mylabels = round(quantile(ED_bird_exotic_sf$mean_ED,
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
  geom_sf(data = ED_bird_exotic_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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

ED_bird_native = calcu_mean_ED_parallel(Tree = phy_data,
                                        Comm = as.matrix(comm_bird_native[,unique(native_distri_data$ScientificName)]),
                                        Area = comm_bird_native$area)

save(ED_bird_native, file = 'results/primary_results/ED_bird_native.rdata')

load("results/primary_results/ED_bird_native.rdata") # load the calculated data from the server
ED_bird_native$df$mean_EDR = ED_bird_native$df$mean_EDR * 1e+06 * 100 
ED_bird_native_1 = cbind(RegionID = df.trans$RegionID,
                         ED_bird_native$df)
ED_bird_native_sf = df.trans %>% left_join(ED_bird_native_1, by = 'RegionID')
colnames(ED_bird_native_sf)

#### native birds mapping 
ED_bird_native_sf = ED_bird_native_sf %>% filter(!is.na(mean_ED))
ED_bird_native_sf$log_EDR = log10(ED_bird_native_sf$mean_EDR)

ED_bird_native_sf$f_mean_ED = cut(ED_bird_native_sf$mean_ED,
                                  breaks = quantile(ED_bird_native_sf$mean_ED,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_bird_native_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_bird_native_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_bird_native_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_bird_native_sf$f_log10_mean_EDR = cut(ED_bird_native_sf$mean_EDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_bird_native_sf_2 = ED_bird_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_bird_native_sf_2 %>% filter(is.na(f_log10_mean_EDR))

birds_native_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_bird_native_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_bird_native_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_bird_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.birds_native_ED_map = legend.func(mycolors = colors2,
                                         mylabels = round(quantile(ED_bird_native_sf$mean_ED,
                                                                   probs = seq(0, 1,length.out = 9),
                                                                   na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_ED_map = ggplotGrob(birds_native_ED_map)
legend.birds_native_ED_map = ggplotGrob(legend.birds_native_ED_map)
birds_native_ED_map_all = arrangeGrob(birds_native_ED_map,
                                      legend.birds_native_ED_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(birds_native_ED_map_all)

birds_native_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_bird_native_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_bird_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Native_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_native_EDR_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 5)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_EDR_map = ggplotGrob(birds_native_EDR_map)
legend.birds_native_EDR_map = ggplotGrob(legend.birds_native_EDR_map)
birds_native_EDR_map_all = arrangeGrob(birds_native_EDR_map,
                                       legend.birds_native_EDR_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_native_EDR_map_all)





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

ED_bird_extant = calcu_mean_ED_parallel(Tree = phy_data,
                                        Comm = as.matrix(comm_bird_extant[,unique(all_distri_data_c$ScientificName)]),
                                        Area = comm_bird_extant$area)

save(ED_bird_extant, file = 'results/primary_results/ED_bird_extant.rdata')

load("results/primary_results/ED_bird_extant.rdata") # load the calculated data from the server
ED_bird_extant$df$mean_EDR = ED_bird_extant$df$mean_EDR * 1e+06 * 100 
ED_bird_extant_1 = cbind(RegionID = df.trans$RegionID,
                         ED_bird_extant$df)
ED_bird_extant_sf = df.trans %>% left_join(ED_bird_extant_1, by = 'RegionID')
colnames(ED_bird_extant_sf)

#### extant birds mapping 
ED_bird_extant_sf = ED_bird_extant_sf %>% filter(!is.na(mean_ED))
ED_bird_extant_sf$log_EDR = log10(ED_bird_extant_sf$mean_EDR)

ED_bird_extant_sf$f_mean_ED = cut(ED_bird_extant_sf$mean_ED,
                                  breaks = quantile(ED_bird_extant_sf$mean_ED,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_bird_extant_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_bird_extant_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_bird_extant_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_bird_extant_sf$f_log10_mean_EDR = cut(ED_bird_extant_sf$mean_EDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_bird_extant_sf_2 = ED_bird_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_bird_extant_sf_2 %>% filter(is.na(f_log10_mean_EDR))

birds_extant_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_bird_extant_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_bird_extant_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_bird_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.birds_extant_ED_map = legend.func(mycolors = colors2,
                                         mylabels = round(quantile(ED_bird_extant_sf$mean_ED,
                                                                   probs = seq(0, 1,length.out = 9),
                                                                   na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_extant_ED_map = ggplotGrob(birds_extant_ED_map)
legend.birds_extant_ED_map = ggplotGrob(legend.birds_extant_ED_map)
birds_extant_ED_map_all = arrangeGrob(birds_extant_ED_map,
                                      legend.birds_extant_ED_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(birds_extant_ED_map_all)

birds_extant_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_bird_extant_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_bird_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Extant_birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.birds_extant_EDR_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 5)) +
  ggtitle("MEDR (Mya 10^2 KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_extant_EDR_map = ggplotGrob(birds_extant_EDR_map)
legend.birds_extant_EDR_map = ggplotGrob(legend.birds_extant_EDR_map)
birds_extant_EDR_map_all = arrangeGrob(birds_extant_EDR_map,
                                       legend.birds_extant_EDR_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_extant_EDR_map_all)



##### Delta_ED (Extant - Native) #####
load('results/primary_results/ED_bird_extant.rdata')
load('results/primary_results/ED_bird_native.rdata')

ED_bird_extant$df$mean_EDR = ED_bird_extant$df$mean_EDR * 1e+06 * 100 
ED_bird_extant_1 = cbind(RegionID = df.trans$RegionID,
                         ED_bird_extant$df)
colnames(ED_bird_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_bird_native$df$mean_EDR = ED_bird_native$df$mean_EDR * 1e+06 * 100 
ED_bird_native_1 = cbind(RegionID = df.trans$RegionID,
                         ED_bird_native$df)
colnames(ED_bird_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_bird = ED_bird_native_1 %>% left_join(ED_bird_extant_1,
                                               by = 'RegionID')


Delta_ED_bird_sf = df_trans %>% left_join(Delta_ED_bird, by = 'RegionID')
colnames(Delta_ED_bird_sf)
Delta_ED_bird_sf$delta_mean_ED = log(Delta_ED_bird_sf$mean_all_ED / 
                                       Delta_ED_bird_sf$mean_native_ED)

Delta_ED_bird_sf$delta_mean_EDR = log(Delta_ED_bird_sf$mean_all_EDR / 
                                        Delta_ED_bird_sf$mean_native_EDR)

#### Delta_ED birds mapping 
Delta_ED_bird_sf = Delta_ED_bird_sf %>% filter(!is.na(delta_mean_ED) & 
                                                 delta_mean_ED != 0 
)

hist(Delta_ED_bird_sf$delta_mean_ED)

breaks_delta_ED = quantile(Delta_ED_bird_sf$delta_mean_ED,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)
breaks_delta_ED[which(names(breaks_delta_ED) == '87.5%')] = 0
Delta_ED_bird_sf$f_delta_mean_ED = cut(Delta_ED_bird_sf$delta_mean_ED,
                                       breaks = breaks_delta_ED,
                                       include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
hist(Delta_ED_bird_sf$delta_mean_EDR)
breaks_delta_EDR = quantile(Delta_ED_bird_sf$delta_mean_EDR,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_EDR[which(names(breaks_delta_EDR) == '87.5%')] = 0
Delta_ED_bird_sf$f_delta_mean_EDR = cut(Delta_ED_bird_sf$delta_mean_EDR,
                                        breaks = breaks_delta_EDR,
                                        include.lowest = TRUE,
                                        include.highest = T,
                                        right = T,
                                        dig.lab = 2)

# delta_mean_EDR = 2.352844e-07
Delta_ED_bird_sf_2 = Delta_ED_bird_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

Delta_ED_bird_sf_2 = Delta_ED_bird_sf_2 %>% filter(!is.na(delta_mean_EDR) & 
                                                     delta_mean_EDR != 0 
)

birds_delta_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_bird_sf, aes(fill = f_delta_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors5,
  # values = scales::rescale(quantile(Delta_ED_bird_sf$delta_mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(Delta_ED_bird_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_ED, color = f_delta_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'ED'
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

legend.birds_delta_ED_map = legend.func(mycolors = colors5,
                                        mylabels = round(breaks_delta_ED, 4)) +
  ggtitle("log (MED_all / MED_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_delta_ED_map = ggplotGrob(birds_delta_ED_map)
legend.birds_delta_ED_map = ggplotGrob(legend.birds_delta_ED_map)
birds_delta_ED_map_all = arrangeGrob(birds_delta_ED_map,
                                     legend.birds_delta_ED_map,
                                     ncol = 1,
                                     layout_matrix = rbind(matrix(1, 4, 10),
                                                           c(NA, rep(2, 8), NA)))
plot(birds_delta_ED_map_all)

birds_delta_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_bird_sf, aes(fill = f_delta_mean_EDR),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
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
  geom_point(data = subset(Delta_ED_bird_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_EDR, color = f_delta_mean_EDR),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'EDR'
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

legend.birds_delta_EDR_map = legend.func(mycolors = colors5,
                                         mylabels = round(breaks_delta_EDR, 4)) +
  ggtitle("log (MEDR_all / MEDR_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_delta_EDR_map = ggplotGrob(birds_delta_EDR_map)
legend.birds_delta_EDR_map = ggplotGrob(legend.birds_delta_EDR_map)
birds_delta_EDR_map_all = arrangeGrob(birds_delta_EDR_map,
                                      legend.birds_delta_EDR_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_delta_EDR_map_all)





#### Fish ED: calculation & mapping ####
load("D:/R projects/Global_ED/data/Fishes/data/my_phy.rdata")
phy_fish = phylo
initial_div_t_fish = max(phyloregion::evol_distinct(tree = phy_fish,
                                                    type = "fair.proportion"))

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

ED_fish_exotic = calcu_mean_ED_parallel(
  Tree = phylo,
  Comm = comm_fish_exotic[,unique(data.used_final_exotics$valid_names)],
  Area = comm_fish_exotic$Surf_area)


save(ED_fish_exotic_sf, file = 'results/primary_results/ED_fish_exotic.rdata')
load('results/primary_results/ED_fish_exotic.rdata')
ED_fish_exotic$df$mean_EDR = ED_fish_exotic$df$mean_EDR * 1e+6
ED_fish_exotic_1 = cbind(Basin.Name = sort(unique(df_trans$X1.Basin.Name)),
                         ED_fish_exotic$df)
ED_fish_exotic_sf = df_trans %>% left_join(ED_fish_exotic_1,
                                           by = join_by('X1.Basin.Name' == 'Basin.Name'))
colnames(ED_fish_exotic_sf)

#### exotic fishes mapping 
ED_fish_exotic_sf = ED_fish_exotic_sf %>% filter(!is.na(mean_ED))
ED_fish_exotic_sf$log_EDR = log10(ED_fish_exotic_sf$mean_EDR)

ED_fish_exotic_sf$f_mean_ED = cut(ED_fish_exotic_sf$mean_ED,
                                  breaks = quantile(ED_fish_exotic_sf$mean_ED,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_fish_exotic_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_fish_exotic_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_fish_exotic_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_fish_exotic_sf$f_log10_mean_EDR = cut(ED_fish_exotic_sf$mean_EDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_fish_exotic_sf_2 = ED_fish_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_fish_exotic_sf_2 %>% filter(is.na(f_log10_mean_EDR))

fishes_exotic_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_fish_exotic_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_fish_exotic_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_fish_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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
                                          mylabels = round(quantile(ED_fish_exotic_sf$mean_ED,
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
  geom_sf(data = ED_fish_exotic_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_fish_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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

ED_fish_native = calcu_mean_ED_parallel(
  Tree = phylo,
  Comm = comm_fish_native[,unique(data.used_final_natives$valid_names)],
  Area = comm_fish_native$Surf_area)

save(ED_fish_native, file = 'results/primary_results/ED_fish_native.rdata')

load('results/primary_results/ED_fish_native.rdata')
ED_fish_native$df$mean_EDR = ED_fish_native$df$mean_EDR * 1e+06
ED_fish_native_1 = cbind(Basin.Name = sort(unique(data.used_final_natives$X1.Basin.Name)),
                         ED_fish_native$df)
ED_fish_native_sf = df_trans %>% left_join(ED_fish_native_1,
                                           by = join_by('X1.Basin.Name' == 'Basin.Name'))
colnames(ED_fish_native_sf)

#### native fishes mapping 
ED_fish_native_sf = ED_fish_native_sf %>% filter(!is.na(mean_ED))
ED_fish_native_sf$log_EDR = log10(ED_fish_native_sf$mean_EDR)

ED_fish_native_sf$f_mean_ED = cut(ED_fish_native_sf$mean_ED,
                                  breaks = quantile(ED_fish_native_sf$mean_ED,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_fish_native_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_fish_native_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_fish_native_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_fish_native_sf$f_log10_mean_EDR = cut(ED_fish_native_sf$mean_EDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_fish_native_sf_2 = ED_fish_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_fish_native_sf_2 %>% filter(is.na(f_log10_mean_EDR))

fishes_native_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_fish_native_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_fish_native_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_fish_native_sf_2, area < units::set_units(8e6, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.fishes_native_ED_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(ED_fish_native_sf$mean_ED,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_native_ED_map = ggplotGrob(fishes_native_ED_map)
legend.fishes_native_ED_map = ggplotGrob(legend.fishes_native_ED_map)
fishes_native_ED_map_all = arrangeGrob(fishes_native_ED_map,
                                       legend.fishes_native_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(fishes_native_ED_map_all)

fishes_native_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_fish_native_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_fish_native_sf_2, area < units::set_units(8e6, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Native_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_native_EDR_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("MEDR (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_native_EDR_map = ggplotGrob(fishes_native_EDR_map)
legend.fishes_native_EDR_map = ggplotGrob(legend.fishes_native_EDR_map)
fishes_native_EDR_map_all = arrangeGrob(fishes_native_EDR_map,
                                        legend.fishes_native_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishes_native_EDR_map_all)


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

ED_fish_extant = calcu_mean_ED_parallel(
  Tree = phylo,
  Comm = comm_fish_extant[,unique(data.used_final$valid_names)],
  Area = comm_fish_extant$Surf_area)

save(ED_fish_extant, file = 'results/primary_results/ED_fish_extant.rdata')

load('results/primary_results/ED_fish_extant.rdata')
ED_fish_extant$df$mean_EDR = ED_fish_extant$df$mean_EDR * 1e+06
ED_fish_extant_1 = cbind(Basin.Name = sort(unique(data.used_final$X1.Basin.Name)),
                         ED_fish_extant$df)
ED_fish_extant_sf = df_trans %>% left_join(ED_fish_extant_1,
                                           by = join_by('X1.Basin.Name' == 'Basin.Name'))
colnames(ED_fish_extant_sf)

#### extant fishes mapping 
ED_fish_extant_sf = ED_fish_extant_sf %>% filter(!is.na(mean_ED))
ED_fish_extant_sf$log_EDR = log10(ED_fish_extant_sf$mean_EDR)

ED_fish_extant_sf$f_mean_ED = cut(ED_fish_extant_sf$mean_ED,
                                  breaks = quantile(ED_fish_extant_sf$mean_ED,
                                                    probs = seq(0, 1,
                                                                length.out = 9),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
breaks = 10^quantile(log10(ED_fish_extant_sf$mean_EDR),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(ED_fish_extant_sf$mean_EDR, na.rm = TRUE)) * 0.99
breaks[length(breaks)] = max(ED_fish_extant_sf$mean_EDR, na.rm = TRUE) * 1.01

ED_fish_extant_sf$f_log10_mean_EDR = cut(ED_fish_extant_sf$mean_EDR,
                                         breaks = breaks,
                                         include.lowest = TRUE,
                                         include.highest = T,
                                         right = T,
                                         dig.lab = 2)

# mean_EDR = 2.352844e-07
ED_fish_extant_sf_2 = ED_fish_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

ED_fish_extant_sf_2 %>% filter(is.na(f_log10_mean_EDR))

fishes_extant_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_fish_extant_sf, aes(fill = f_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors2,
  # values = scales::rescale(quantile(ED_fish_extant_sf$mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(ED_fish_extant_sf_2, area < units::set_units(8e6, "m^2")),
             aes(x = lon, y = lat, size = f_mean_ED, color = f_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'ED'
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

legend.fishes_extant_ED_map = legend.func(mycolors = colors2,
                                          mylabels = round(quantile(ED_fish_extant_sf$mean_ED,
                                                                    probs = seq(0, 1,length.out = 9),
                                                                    na.rm = TRUE), 2)) +
  ggtitle("MED (Mya)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_extant_ED_map = ggplotGrob(fishes_extant_ED_map)
legend.fishes_extant_ED_map = ggplotGrob(legend.fishes_extant_ED_map)
fishes_extant_ED_map_all = arrangeGrob(fishes_extant_ED_map,
                                       legend.fishes_extant_ED_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10),
                                                             c(NA, rep(2, 8), NA)))
plot(fishes_extant_ED_map_all)

fishes_extant_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = ED_fish_extant_sf, aes(fill = f_log10_mean_EDR),
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
  geom_point(data = subset(ED_fish_extant_sf_2, area < units::set_units(8e6, "m^2")),
             aes(x = lon, y = lat, size = f_log10_mean_EDR, color = f_log10_mean_EDR),
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
  ggtitle('Extant_fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.fishes_extant_EDR_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("MEDR (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_extant_EDR_map = ggplotGrob(fishes_extant_EDR_map)
legend.fishes_extant_EDR_map = ggplotGrob(legend.fishes_extant_EDR_map)
fishes_extant_EDR_map_all = arrangeGrob(fishes_extant_EDR_map,
                                        legend.fishes_extant_EDR_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishes_extant_EDR_map_all)





##### Delta_ED (Extant - Native) #####
load('results/primary_results/ED_fish_extant.rdata')
load('results/primary_results/ED_fish_native.rdata')

ED_fish_extant$df$mean_EDR = ED_fish_extant$df$mean_EDR * 1e+06
ED_fish_extant_1 = cbind(Basin.Name = sort(unique(data.used_final$X1.Basin.Name)),
                         ED_fish_extant$df)
colnames(ED_fish_extant_1)[2:3] = c("mean_all_ED", "mean_all_EDR")

ED_fish_native$df$mean_EDR = ED_fish_native$df$mean_EDR * 1e+06
ED_fish_native_1 = cbind(Basin.Name = sort(unique(data.used_final_natives$X1.Basin.Name)),
                         ED_fish_native$df)
colnames(ED_fish_native_1)[2:3] = c("mean_native_ED", "mean_native_EDR")

Delta_ED_fish = ED_fish_native_1 %>% left_join(ED_fish_extant_1,
                                               by = 'Basin.Name')


Delta_ED_fish_sf = df_trans %>% left_join(Delta_ED_fish,
                                          join_by('X1.Basin.Name' == 'Basin.Name'))
colnames(Delta_ED_fish_sf)
Delta_ED_fish_sf$delta_mean_ED = log(Delta_ED_fish_sf$mean_all_ED / 
                                       Delta_ED_fish_sf$mean_native_ED)

Delta_ED_fish_sf$delta_mean_EDR = log(Delta_ED_fish_sf$mean_all_EDR / 
                                        Delta_ED_fish_sf$mean_native_EDR)

#### Delta_ED fishes mapping 
Delta_ED_fish_sf = Delta_ED_fish_sf %>% filter(!is.na(delta_mean_ED) & 
                                                 delta_mean_ED != 0 
)

hist(Delta_ED_fish_sf$delta_mean_ED)

breaks_delta_ED = quantile(Delta_ED_fish_sf$delta_mean_ED,
                           probs = seq(0, 1,
                                       length.out = 9),
                           na.rm = TRUE,
                           names = T,
                           digits = 12)
breaks_delta_ED[which(names(breaks_delta_ED) == '87.5%')] = 0
Delta_ED_fish_sf$f_delta_mean_ED = cut(Delta_ED_fish_sf$delta_mean_ED,
                                       breaks = breaks_delta_ED,
                                       include.lowest = TRUE)

# Ensure full coverage with a tiny buffer
hist(Delta_ED_fish_sf$delta_mean_EDR)

breaks_delta_EDR = quantile(Delta_ED_fish_sf$delta_mean_EDR,
                            probs = seq(0, 1,
                                        length.out = 9),
                            na.rm = TRUE,
                            names = T,
                            digits = 12)
breaks_delta_EDR[which(names(breaks_delta_EDR) == '87.5%')] = 0
Delta_ED_fish_sf$f_delta_mean_EDR = cut(Delta_ED_fish_sf$delta_mean_EDR,
                                        breaks = breaks_delta_EDR,
                                        include.lowest = T,
                                        include.highest = T,
                                        right = T,
                                        dig.lab = 2)

# delta_mean_EDR = 2.352844e-07
Delta_ED_fish_sf_2 = Delta_ED_fish_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

Delta_ED_fish_sf_2 = Delta_ED_fish_sf_2 %>% filter(!is.na(delta_mean_EDR) & 
                                                     delta_mean_EDR != 0 
)

fishes_delta_ED_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_fish_sf, aes(fill = f_delta_mean_ED),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'ED',
    guide = guide_legend(
      title = "ED",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  #scale_fill_gradientn(
  # colors = colors5,
  # values = scales::rescale(quantile(Delta_ED_fish_sf$delta_mean_ED,
  #                                  probs = seq(0, 1,
  #                                            length.out = 7),
  #                                na.rm = TRUE)),
  #   na.value = 'gray90',
  #  name = 'ED'
  #  ) + 
  #ggnewscale::new_scale_fill() +
  #ggplot()+
  # Circle layer for islands, colored by PE and sized by richness
  geom_point(data = subset(Delta_ED_fish_sf_2, area < units::set_units(8e6, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_ED, color = f_delta_mean_ED),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'ED'
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

legend.fishes_delta_ED_map = legend.func(mycolors = colors5,
                                         mylabels = round(breaks_delta_ED, 4)) +
  ggtitle("log (MED_all / MED_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_delta_ED_map = ggplotGrob(fishes_delta_ED_map)
legend.fishes_delta_ED_map = ggplotGrob(legend.fishes_delta_ED_map)
fishes_delta_ED_map_all = arrangeGrob(fishes_delta_ED_map,
                                      legend.fishes_delta_ED_map,
                                      ncol = 1,
                                      layout_matrix = rbind(matrix(1, 4, 10),
                                                            c(NA, rep(2, 8), NA)))
plot(fishes_delta_ED_map_all)

fishes_delta_EDR_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = Delta_ED_fish_sf, aes(fill = f_delta_mean_EDR),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors5,
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
  geom_point(data = subset(Delta_ED_fish_sf_2, area < units::set_units(8e6, "m^2")),
             aes(x = lon, y = lat, size = f_delta_mean_EDR, color = f_delta_mean_EDR),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors5,
    na.value = 'gray90',
    name = 'EDR'
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

legend.fishes_delta_EDR_map = legend.func(mycolors = colors5,
                                          mylabels = round(breaks_delta_EDR, 2)) +
  ggtitle("log (MEDR_all / MEDR_native)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishes_delta_EDR_map = ggplotGrob(fishes_delta_EDR_map)
legend.fishes_delta_EDR_map = ggplotGrob(legend.fishes_delta_EDR_map)
fishes_delta_EDR_map_all = arrangeGrob(fishes_delta_EDR_map,
                                       legend.fishes_delta_EDR_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishes_delta_EDR_map_all)








#### Which factors affect Delta_ED across 4 taxonomic groups ####
Delta_ED_mammal$ini_div_t = initial_div_t_mammal
Delta_ED_plant$ini_div_t = initial_div_t_plant
Delta_ED_bird$ini_div_t = initial_div_t_bird
Delta_ED_fish$ini_div_t = initial_div_t_fish

Delta_ED_mammal$taxa = 'mammal'
Delta_ED_plant$taxa = 'plant'
Delta_ED_bird$taxa = 'bird'
Delta_ED_fish$taxa = 'fish'

colnames(Delta_ED_plant)[which(colnames(Delta_ED_plant) == "Region_id")] = "RegionID"
colnames(Delta_ED_fish)[which(colnames(Delta_ED_fish) == "Basin.Name")] = "RegionID"

Delta_ED_all_taxa = rbind(Delta_ED_mammal, Delta_ED_plant, 
                          Delta_ED_bird, Delta_ED_fish)

Delta_ED_all_taxa$delta_mean_ED = log(Delta_ED_all_taxa$mean_all_ED / 
                                        Delta_ED_all_taxa$mean_native_ED)

Delta_ED_all_taxa$delta_mean_EDR = log(Delta_ED_all_taxa$mean_all_EDR / 
                                         Delta_ED_all_taxa$mean_native_EDR)

Delta_ED_all_taxa$ini_div_t = factor(Delta_ED_all_taxa$ini_div_t)
Delta_ED_all_taxa$taxa = factor(Delta_ED_all_taxa$taxa)

plot(Delta_ED_all_taxa$ini_div_t, Delta_ED_all_taxa$delta_mean_ED)
plot(Delta_ED_all_taxa$taxa, Delta_ED_all_taxa$delta_mean_ED)


#### Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

##### Mean ED #####
plot_width = 0.5
plot_height = 0.20

figs_ED_nat_exo = ggdraw() +
  draw_plot(plants_native_ED_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_ED_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_ED_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_ED_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_ED_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_ED_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_ED_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_exotic_ED_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_ED_nat_exo.png',
    #plot = figs_ED_nat_exo,
    height=40, width=25, # setting emfPlusFontToPath=TRUE to 
    res = 300,
    # ensure text looks correct on the viewing system
    units = 'cm')
figs_ED_nat_exo
dev.off() #turn off device and finalize file


#emf('figures/figs_ED_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_ED_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_ED_nat_extant = ggdraw() +
  draw_plot(plants_native_ED_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_ED_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_ED_map_all, x = 0, y = 0.53, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_ED_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_ED_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_ED_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_ED_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_ED_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_ED_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_ED_nat_extant
dev.off() #turn off device and finalize file


# Compare the PE patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.015         # small horizontal gap between columns

figs_ED_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_ED_map_all,  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_ED_map_all,  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_ED_map_all,   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_ED_map_all,   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_extant_ED_map_all,   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_ED_map_all,    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_ED_map_all, x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_ED_map_all, x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_ED_map_all,  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_ED_map_all,   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_ED_map_all,   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_ED_map_all,   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
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

png(filename = 'figures/figs_ED_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_ED_nat_extant_delta
dev.off() #turn off device and finalize file



##### Mean EDR #####
plot_width = 0.5
plot_height = 0.20

figs_EDR_nat_exo = ggdraw() +
  draw_plot(plants_native_EDR_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_EDR_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_EDR_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_EDR_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_EDR_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_EDR_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_EDR_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_exotic_EDR_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_EDR_nat_exo.png',
    #plot = figs_EDR_nat_exo,
    height=40, width=25, # setting emfPlusFontToPath=TRUE to 
    res = 300,
    # ensure text looks correct on the viewing system
    units = 'cm')
figs_EDR_nat_exo
dev.off() #turn off device and finalize file


#emf('figures/figs_EDR_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_EDR_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_EDR_nat_extant = ggdraw() +
  draw_plot(plants_native_EDR_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_EDR_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_EDR_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_EDR_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_EDR_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_EDR_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_EDR_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_EDR_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_EDR_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_EDR_nat_extant
dev.off() #turn off device and finalize file


# Compare the PE patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.015         # small horizontal gap between columns

figs_EDR_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_EDR_map_all,  x = 0,                 y = 0.75,
            width = plot_width, height = plot_height) +
  draw_plot(plants_extant_EDR_map_all,  x = plot_width + gap,  y = 0.75,
            width = plot_width, height = plot_height) +
  draw_plot(plants_delta_EDR_map_all,   x = 2*(plot_width + gap), y = 0.75,
            width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_EDR_map_all,   x = 0,                 y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(birds_extant_EDR_map_all,   x = plot_width + gap,  y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(birds_delta_EDR_map_all,    x = 2*(plot_width + gap), y = 0.5,
            width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_EDR_map_all, x = 0,                 y = 0.25,
            width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_EDR_map_all, x = plot_width + gap,  y = 0.25,
            width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_EDR_map_all,  x = 2*(plot_width + gap), y = 0.25,
            width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_EDR_map_all,   x = 0,                 y = 0,
            width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_EDR_map_all,   x = plot_width + gap,  y = 0,
            width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_EDR_map_all,   x = 2*(plot_width + gap), y = 0,
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

png(filename = 'figures/figs_EDR_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)

figs_EDR_nat_extant_delta
dev.off() #turn off device and finalize file

