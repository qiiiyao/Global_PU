### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Globa_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'scico', 'ggplot', 'gridExtra')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})

setwd("D:/R projects/Global_ED")
source('code/calculating_PE_func.R')
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)


# define color gradients
colors1 = rev(scico::scico(n=8, palette = "lajolla"))
colors2 = scico::scico(n=10, palette = "vik")[1:8] 
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

# load the background data for plotting the world map plot
load("code/FYI/Codes_and_Data_Fan_et_al_2023/data.for.shp.plots.4.Rdata")


#### Mammal PE: calculation & mapping ####
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)

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
  relocate(area, .after = RegionID)

pe_mammal_exotic = calcu_PE_PE_parallel(
                      Tree = spec_phy.3 %>% 
                        keep.tip(unique(sp_overlap_dat_1$Binomial)),
                      Comm = comm_mammal_exotic[,unique(sp_overlap_dat_1$Binomial)],
                      Area = comm_mammal_exotic$area)

pe_mammal_exotic_1 = cbind(RegionID = comm_mammal_exotic$RegionID,
                           pe_mammal_exotic)
pe_mammal_exotic_sf = df %>% left_join(pe_mammal_exotic_1, by = 'RegionID')
colnames(pe_mammal_exotic_sf)

save(pe_mammal_exotic_sf,
     file = 'results/primary_results/pe_mammal_exotic_sf.rdata')

#### exotic mammals mapping 
load('results/primary_results/pe_mammal_exotic_sf.rdata')

pe_mammal_exotic_sf = df_trans %>% left_join((pe_mammal_exotic_sf %>% 
                                               st_drop_geometry() %>% 
                                               dplyr::select(c('PE',
                                                      'PE.numerator',
                                                      'PE.denominator',
                                                      'RPE',
                                                      'RegionID'))), by = 'RegionID')

pe_mammal_exotic_sf = pe_mammal_exotic_sf %>% filter(PE !=0 )
pe_mammal_exotic_sf$log_PE = log10(pe_mammal_exotic_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_mammal_exotic_sf$PE),
                     probs = seq(0, 1,
                                 length.out = 9),
                     na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_mammal_exotic_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_mammal_exotic_sf$log_PE, na.rm = TRUE) * 0.99

pe_mammal_exotic_sf$f_log10_PE = cut(pe_mammal_exotic_sf$log_PE,
                                           breaks = breaks,
                                           include.lowest = TRUE,
                                           include.highest = T,
                                           right = T,
                                           dig.lab = 2)

# PE = 2.352844e-07
pe_mammal_exotic_sf_2 = pe_mammal_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_mammal_exotic_sf_2 %>% filter(is.na(f_log10_PE))

mammals_exotic_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_mammal_exotic_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_mammal_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.mammals_exotic_PE_map = legend.func(mycolors = colors2,
                                            mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_exotic_PE_map = ggplotGrob(mammals_exotic_PE_map)
legend.mammals_exotic_PE_map = ggplotGrob(legend.mammals_exotic_PE_map)
mammals_exotic_PE_map_all = arrangeGrob(mammals_exotic_PE_map,
                                         legend.mammals_exotic_PE_map,
                                         ncol = 1,
                                         layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_exotic_PE_map_all)

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
  relocate(area, .after = Region.ID)

pe_mammal_native = calcu_PE_PE_parallel(Tree = spec_phy.3 %>% keep.tip(unique(sp_dis_5$ScientificName)),
                                  Comm = as.matrix(comm_mammal_native[,unique(sp_dis_5$ScientificName)]),
                                  Area = comm_mammal_native$area)
pe_mammal_native_1 = cbind(RegionID = comm_mammal_native$Region.ID,
                    pe_mammal_native)
pe_mammal_native_sf = df %>% left_join(pe_mammal_native_1, by = 'RegionID')
colnames(pe_mammal_native_sf)
save(pe_mammal_native_sf,
     file = 'results/primary_results/pe_mammal_native_sf.rdata')

#### native mammals mapping 
load('results/primary_results/pe_mammal_native_sf.rdata')

pe_mammal_native_sf = df_trans %>% left_join((pe_mammal_native_sf %>% 
                                                st_drop_geometry() %>% 
                                                dplyr::select(c('PE',
                                                                'PE.numerator',
                                                                'PE.denominator',
                                                                'RPE',
                                                                'RegionID'))), by = 'RegionID')

pe_mammal_native_sf = pe_mammal_native_sf %>% filter(PE !=0 )
pe_mammal_native_sf$log_PE = log10(pe_mammal_native_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_mammal_native_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_mammal_native_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_mammal_native_sf$log_PE, na.rm = TRUE) * 0.99

pe_mammal_native_sf$f_log10_PE = cut(pe_mammal_native_sf$log_PE,
                                     breaks = breaks,
                                     include.lowest = TRUE,
                                     include.highest = T,
                                     right = T,
                                     dig.lab = 2)

# PE = 2.352844e-07
pe_mammal_native_sf_2 = pe_mammal_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_mammal_native_sf_2 %>% filter(is.na(f_log10_PE))

mammals_native_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_mammal_native_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_mammal_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.mammals_native_PE_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_native_PE_map = ggplotGrob(mammals_native_PE_map)
legend.mammals_native_PE_map = ggplotGrob(legend.mammals_native_PE_map)
mammals_native_PE_map_all = arrangeGrob(mammals_native_PE_map,
                                        legend.mammals_native_PE_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_native_PE_map_all)

##### extant regional assemlage: natives + exotics #####
colnames(sp_dis_5)
colnames(sp_overlap_dat_1)

colnames(sp_overlap_dat_1)[which(colnames(sp_overlap_dat_1) == 'Inv. Stage')] = 'Inv._stage'
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)

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
  relocate(area, .after = RegionID)

pe_mammal_all = calcu_PE_PE_parallel(
  Tree = spec_phy.3 %>% 
    keep.tip(unique(sp_mammal_all$Binomial)),
  Comm = comm_mammal_all[,unique(sp_mammal_all$Binomial)],
  Area = comm_mammal_all$area)

pe_mammal_all_1 = cbind(RegionID = comm_mammal_all$RegionID,
                        pe_mammal_all)
pe_mammal_all_sf = df %>% left_join(pe_mammal_all_1, by = 'RegionID')
colnames(pe_mammal_all_sf)

save(pe_mammal_all_sf,
     file = 'results/primary_results/pe_mammal_all_sf.rdata')

#### all mammals mapping 
load('results/primary_results/pe_mammal_all_sf.rdata')

pe_mammal_all_sf = df_trans %>% left_join((pe_mammal_all_sf %>% 
                                                st_drop_geometry() %>% 
                                                dplyr::select(c('PE',
                                                                'PE.numerator',
                                                                'PE.denominator',
                                                                'RPE',
                                                                'RegionID'))), by = 'RegionID')

pe_mammal_all_sf = pe_mammal_all_sf %>% filter(PE !=0 )
pe_mammal_all_sf$log_PE = log10(pe_mammal_all_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_mammal_all_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_mammal_all_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_mammal_all_sf$log_PE, na.rm = TRUE) * 0.99

pe_mammal_all_sf$f_log10_PE = cut(pe_mammal_all_sf$log_PE,
                                     breaks = breaks,
                                     include.lowest = TRUE,
                                     include.highest = T,
                                     right = T,
                                     dig.lab = 2)

# PE = 2.352844e-07
pe_mammal_all_sf_2 = pe_mammal_all_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_mammal_all_sf_2 %>% filter(is.na(f_log10_PE))

mammals_all_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_mammal_all_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_mammal_all_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
  ) +
  scale_size_discrete(
    range = c(1, 5)
    #,name = "Island size indicator"
  ) + 
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('all_mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) 

legend.mammals_all_PE_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
mammals_all_PE_map = ggplotGrob(mammals_all_PE_map)
legend.mammals_all_PE_map = ggplotGrob(legend.mammals_all_PE_map)
mammals_all_PE_map_all = arrangeGrob(mammals_all_PE_map,
                                        legend.mammals_all_PE_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(mammals_all_PE_map_all)


#### Plant PE: calculation & mapping ####
load("D:/R projects/Globa_ED/data/Plants/data/shp.651.Rdata")
load("D:/R projects/Globa_ED/data/Plants/data/phylo.fake.species.653.Rdata")
load("D:/R projects/Globa_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)

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
  relocate(area, .after = Region_id)

pe_plant_exotic = calcu_PE_PE_parallel(Tree = phylo %>% keep.tip(unique(df.natu.650$species)),
                                  Comm = comm_plant_exotic[,unique(df.natu.650$species)],
                                  Area = comm_plant_exotic$area)
pe_plant_exotic_1 = cbind(Region_id = comm_plant_exotic$Region_id,
                          pe_plant_exotic)
pe_plant_exotic_sf = shp.glonaf.trans %>% left_join(pe_plant_exotic_1, by = 'Region_id')
colnames(pe_plant_exotic_sf)

#### exotic plants mapping 
load("results/primary_results/pe_plant_exotic_sf.rdata") # load the calculated data from the server

pe_plant_exotic_sf = pe_plant_exotic_sf %>% filter(PE !=0 )
pe_plant_exotic_sf$log_PE = log10(pe_plant_exotic_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_plant_exotic_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_plant_exotic_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_plant_exotic_sf$log_PE, na.rm = TRUE) * 0.99

pe_plant_exotic_sf$f_log10_PE = cut(pe_plant_exotic_sf$log_PE,
                                     breaks = breaks,
                                     include.lowest = TRUE,
                                     include.highest = T,
                                     right = T,
                                     dig.lab = 2)

# PE = 2.352844e-07
pe_plant_exotic_sf_2 = pe_plant_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_plant_exotic_sf_2 %>% filter(is.na(f_log10_PE))

plants_exotic_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_plant_exotic_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_plant_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.plants_exotic_PE_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_exotic_PE_map = ggplotGrob(plants_exotic_PE_map)
legend.plants_exotic_PE_map = ggplotGrob(legend.plants_exotic_PE_map)
plants_exotic_PE_map_all = arrangeGrob(plants_exotic_PE_map,
                                        legend.plants_exotic_PE_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_exotic_PE_map_all)

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
  relocate(area, .after = Region_id)

pe_plant_native = calcu_PE_PE_parallel(Tree = phylo %>% keep.tip(unique(df.native.650$species)),
                                  Comm = as.matrix(comm_plant_native[,unique(df.native.650$species)]),
                                  Area = comm_plant_native$area)
pe_plant_native_1 = cbind(Region_id = comm_plant_native$Region_id,
                    pe_plant_native)
pe_plant_native_sf = shp.glonaf.trans %>% left_join(pe_plant_native_1, by = 'Region_id')
colnames(pe_plant_native_sf)

#### native plants mapping 
load("results/primary_results/pe_plant_native_sf.rdata") # load the calculated data from the server

pe_plant_native_sf = pe_plant_native_sf %>% filter(PE !=0 )
pe_plant_native_sf$log_PE = log10(pe_plant_native_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_plant_native_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_plant_native_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_plant_native_sf$log_PE, na.rm = TRUE) * 0.99

pe_plant_native_sf$f_log10_PE = cut(pe_plant_native_sf$log_PE,
                                    breaks = breaks,
                                    include.lowest = TRUE,
                                    include.highest = T,
                                    right = T,
                                    dig.lab = 2)

# PE = 2.352844e-07
pe_plant_native_sf_2 = pe_plant_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_plant_native_sf_2 %>% filter(is.na(f_log10_PE))

plants_native_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_plant_native_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_plant_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.plants_native_PE_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_native_PE_map = ggplotGrob(plants_native_PE_map)
legend.plants_native_PE_map = ggplotGrob(legend.plants_native_PE_map)
plants_native_PE_map_all = arrangeGrob(plants_native_PE_map,
                                       legend.plants_native_PE_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_native_PE_map_all)


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
  relocate(area, .after = Region_id)

pe_plant_extant = calcu_PE_PE_parallel(Tree = phylo %>% keep.tip(unique(df.extant.650$species)),
                                        Comm = as.matrix(comm_plant_extant[,unique(df.extant.650$species)]),
                                        Area = comm_plant_extant$area)
pe_plant_extant_1 = cbind(Region_id = comm_plant_extant$Region_id,
                          pe_plant_extant)
pe_plant_extant_sf = shp.glonaf.trans %>% left_join(pe_plant_extant_1, by = 'Region_id')
colnames(pe_plant_extant_sf)

save(pe_plant_extant_sf, file = 'results/primary_results/pe_plant_all_sf.rdata')

#### extant plants mapping 
load("results/primary_results/pe_plant_all_sf.rdata") # load the calculated data from the server

pe_plant_extant_sf = pe_plant_extant_sf %>% filter(PE !=0 )
pe_plant_extant_sf$log_PE = log10(pe_plant_extant_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_plant_extant_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_plant_extant_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_plant_extant_sf$log_PE, na.rm = TRUE) * 0.99

pe_plant_extant_sf$f_log10_PE = cut(pe_plant_extant_sf$log_PE,
                                    breaks = breaks,
                                    include.lowest = TRUE,
                                    include.highest = T,
                                    right = T,
                                    dig.lab = 2)

# PE = 2.352844e-07
pe_plant_extant_sf_2 = pe_plant_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_plant_extant_sf_2 %>% filter(is.na(f_log10_PE))

plants_extant_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_plant_extant_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_plant_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.plants_extant_PE_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
plants_extant_PE_map = ggplotGrob(plants_extant_PE_map)
legend.plants_extant_PE_map = ggplotGrob(legend.plants_extant_PE_map)
plants_extant_PE_map_extant = arrangeGrob(plants_extant_PE_map,
                                       legend.plants_extant_PE_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(plants_extant_PE_map_extant)


#### Bird PE: calculation & mapping ####
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
phy_data = ape::read.tree("data/Birds/data/Phylogenetic_Birds.tre")

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
  relocate(area, .after = RegionID)

pe_bird_exotic = calcu_PE_PE_parallel(Tree = phy_data %>% keep.tip(unique(exotic_distri_data$ScientificName)),
                                        Comm = comm_bird_exotic[,unique(exotic_distri_data$ScientificName)],
                                        Area = comm_bird_exotic$area)
pe_bird_exotic_1 = cbind(RegionID = comm_bird_exotic$RegionID,
                          pe_bird_exotic)
pe_bird_exotic_sf = df.trans %>% left_join(pe_bird_exotic_1, by = 'RegionID')
colnames(pe_bird_exotic_sf)
save(pe_bird_exotic_sf,
     file = 'results/primary_results/pe_bird_exotic_sf.rdata')

#### exotic birds mapping 
load('results/primary_results/pe_bird_exotic_sf.rdata')

pe_bird_exotic_sf = df_trans %>% left_join((pe_bird_exotic_sf %>% 
                                                st_drop_geometry() %>% 
                                                dplyr::select(c('PE',
                                                                'PE.numerator',
                                                                'PE.denominator',
                                                                'RPE',
                                                                'RegionID'))), by = 'RegionID')

pe_bird_exotic_sf = pe_bird_exotic_sf %>% filter(PE !=0 )
pe_bird_exotic_sf$log_PE = log10(pe_bird_exotic_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_bird_exotic_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_bird_exotic_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_bird_exotic_sf$log_PE, na.rm = TRUE) * 0.99

pe_bird_exotic_sf$f_log10_PE = cut(pe_bird_exotic_sf$log_PE,
                                     breaks = breaks,
                                     include.lowest = TRUE,
                                     include.highest = T,
                                     right = T,
                                     dig.lab = 2)

# PE = 2.352844e-07
pe_bird_exotic_sf_2 = pe_bird_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_bird_exotic_sf_2 %>% filter(is.na(f_log10_PE))

birds_exotic_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_bird_exotic_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_bird_exotic_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.birds_exotic_PE_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_exotic_PE_map = ggplotGrob(birds_exotic_PE_map)
legend.birds_exotic_PE_map = ggplotGrob(legend.birds_exotic_PE_map)
birds_exotic_PE_map_all = arrangeGrob(birds_exotic_PE_map,
                                        legend.birds_exotic_PE_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_exotic_PE_map_all)


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
  relocate(area, .after = RegionID)

pe_bird_native = calcu_PE_PE_parallel(Tree = phy_data %>% keep.tip(unique(native_distri_data$ScientificName)),
                                        Comm = as.matrix(comm_bird_native[,unique(native_distri_data$ScientificName)]),
                                        Area = comm_bird_native$area)
pe_bird_native_1 = cbind(RegionID = comm_bird_native$RegionID,
                          pe_bird_native)
pe_bird_native_sf = df.trans %>% left_join(pe_bird_native_1, by = 'RegionID')
colnames(pe_bird_native_sf)
save(pe_bird_native_sf, file = 'results/primary_results/pe_bird_native_sf.rdata')

#### native birds mapping 
load("results/primary_results/pe_bird_native_sf.rdata") # load the calculated data from the server

pe_bird_native_sf = df_trans %>% left_join((pe_bird_native_sf %>% 
                                                st_drop_geometry() %>% 
                                                dplyr::select(c('PE',
                                                                'PE.numerator',
                                                                'PE.denominator',
                                                                'RPE',
                                                                'RegionID'))), by = 'RegionID')

pe_bird_native_sf = pe_bird_native_sf %>% filter(PE !=0 )
pe_bird_native_sf$log_PE = log10(pe_bird_native_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_bird_native_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_bird_native_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_bird_native_sf$log_PE, na.rm = TRUE) * 0.99

pe_bird_native_sf$f_log10_PE = cut(pe_bird_native_sf$log_PE,
                                     breaks = breaks,
                                     include.lowest = TRUE,
                                     include.highest = T,
                                     right = T,
                                     dig.lab = 2)

# PE = 2.352844e-07
pe_bird_native_sf_2 = pe_bird_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_bird_native_sf_2 %>% filter(is.na(f_log10_PE))

birds_native_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_bird_native_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_bird_native_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.birds_native_PE_map = legend.func(mycolors = colors2,
                                           mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_native_PE_map = ggplotGrob(birds_native_PE_map)
legend.birds_native_PE_map = ggplotGrob(legend.birds_native_PE_map)
birds_native_PE_map_all = arrangeGrob(birds_native_PE_map,
                                        legend.birds_native_PE_map,
                                        ncol = 1,
                                        layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_native_PE_map_all)


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
  relocate(area, .after = RegionID)

pe_bird_extant = calcu_PE_PE_parallel(Tree = phy_data %>% keep.tip(unique(all_distri_data_c$ScientificName)),
                                       Comm = as.matrix(comm_bird_extant[,unique(all_distri_data_c$ScientificName)]),
                                       Area = comm_bird_extant$area)
pe_bird_extant_1 = cbind(RegionID = comm_bird_extant$RegionID,
                         pe_bird_extant)
pe_bird_extant_sf = df.trans %>% left_join(pe_bird_extant_1, by = 'RegionID')
colnames(pe_bird_extant_sf)
save(pe_bird_extant_sf, file = 'results/primary_results/pe_bird_all_sf.rdata')

#### extant birds mapping 
load("results/primary_results/pe_bird_all_sf.rdata") # load the calculated data from the server

pe_bird_extant_sf = df_trans %>% left_join((pe_bird_extant_sf %>% 
                                             st_drop_geometry() %>% 
                                             dplyr::select(c('PE',
                                                             'PE.numerator',
                                                             'PE.denominator',
                                                             'RPE',
                                                             'RegionID'))), by = 'RegionID')

pe_bird_extant_sf = pe_bird_extant_sf %>% filter(PE !=0 )
pe_bird_extant_sf$log_PE = log10(pe_bird_extant_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_bird_extant_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_bird_extant_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_bird_extant_sf$log_PE, na.rm = TRUE) * 0.99

pe_bird_extant_sf$f_log10_PE = cut(pe_bird_extant_sf$log_PE,
                                  breaks = breaks,
                                  include.lowest = TRUE,
                                  include.highest = T,
                                  right = T,
                                  dig.lab = 2)

# PE = 2.352844e-07
pe_bird_extant_sf_2 = pe_bird_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_bird_extant_sf_2 %>% filter(is.na(f_log10_PE))

birds_all_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_bird_extant_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_bird_extant_sf_2, area < units::set_units(1e8, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.birds_all_PE_map = legend.func(mycolors = colors2,
                                        mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
birds_all_PE_map = ggplotGrob(birds_all_PE_map)
legend.birds_all_PE_map = ggplotGrob(legend.birds_all_PE_map)
birds_all_PE_map_all = arrangeGrob(birds_all_PE_map,
                                     legend.birds_all_PE_map,
                                     ncol = 1,
                                     layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(birds_all_PE_map_all)


#### Fish PE: calculation & mapping ####
load("D:/R projects/Globa_ED/data/Fishes/data/my_phy.rdata")
load("D:/R projects/Globa_ED/data/Fishes/data/my_data_used_final.rdata")
df = st_read("data/Fishes/data/Basin042017_3119.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
save(df, file = "data/Fishes/data/Basin042017_3119.rdata")
colnames(df_trans)[which(colnames(df_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% dplyr::filter(X3.Native.Exotic.Status == 
                                                       'exotic')
data.used_final_natives = data.used_final %>% dplyr::filter(X3.Native.Exotic.Status == 
                                                       'native')

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
  relocate(Surf_area, .after = X1.Basin.Name)

pe_fish_exotic = calcu_PE_PE_parallel(
  Tree = phylo %>% 
    keep.tip(unique(data.used_final_exotics$valid_names)),
  Comm = comm_fish_exotic[,unique(data.used_final_exotics$valid_names)],
  Area = comm_fish_exotic$Surf_area)

pe_fish_exotic_1 = cbind(Surf_area = comm_fish_exotic$Surf_area,
                    pe_fish_exotic)
pe_fish_exotic_sf = df_trans %>% left_join(pe_fish_exotic_1, by = 'Surf_area')
colnames(pe_fish_exotic_sf)
save(pe_fish_exotic_sf, file = 'results/primary_results/pe_fish_exotic_sf.rdata')

#### exotic fishs mapping 
load('results/primary_results/pe_fish_exotic_sf.rdata')

pe_fish_exotic_sf = pe_fish_exotic_sf %>% filter(PE !=0 )
pe_fish_exotic_sf$log_PE = log10(pe_fish_exotic_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_fish_exotic_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_fish_exotic_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_fish_exotic_sf$log_PE, na.rm = TRUE) * 0.99

pe_fish_exotic_sf$f_log10_PE = cut(pe_fish_exotic_sf$log_PE,
                                    breaks = breaks,
                                    include.lowest = TRUE,
                                    include.highest = T,
                                    right = T,
                                    dig.lab = 2)

# PE = 2.352844e-07
pe_fish_exotic_sf_2 = pe_fish_exotic_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_fish_exotic_sf_2 %>% filter(is.na(f_log10_PE))

fishs_exotic_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_fish_exotic_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_fish_exotic_sf_2, area < units::set_units(1e7, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.fishs_exotic_PE_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishs_exotic_PE_map = ggplotGrob(fishs_exotic_PE_map)
legend.fishs_exotic_PE_map = ggplotGrob(legend.fishs_exotic_PE_map)
fishs_exotic_PE_map_all = arrangeGrob(fishs_exotic_PE_map,
                                       legend.fishs_exotic_PE_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishs_exotic_PE_map_all)


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
  relocate(Surf_area, .after = X1.Basin.Name)

pe_fish_native = calcu_PE_PE_parallel(
  Tree = phylo %>% 
    keep.tip(unique(data.used_final_natives$valid_names)),
  Comm = comm_fish_native[,unique(data.used_final_natives$valid_names)],
  Area = comm_fish_native$Surf_area)

pe_fish_native_1 = cbind(Surf_area = comm_fish_native$Surf_area,
                    pe_fish_native)
pe_fish_native_sf = df_trans %>% left_join(pe_fish_native_1, by = 'Surf_area')
colnames(pe_fish_native_sf)
save(pe_fish_native_sf, file = 'results/primary_results/pe_fish_native_sf.rdata')

#### native fishs mapping 
load('results/primary_results/pe_fish_native_sf.rdata')

pe_fish_native_sf = pe_fish_native_sf %>% filter(PE !=0 )
pe_fish_native_sf$log_PE = log10(pe_fish_native_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_fish_native_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_fish_native_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_fish_native_sf$log_PE, na.rm = TRUE) * 0.99

pe_fish_native_sf$f_log10_PE = cut(pe_fish_native_sf$log_PE,
                                    breaks = breaks,
                                    include.lowest = TRUE,
                                    include.highest = T,
                                    right = T,
                                    dig.lab = 2)

# PE = 2.352844e-07
pe_fish_native_sf_2 = pe_fish_native_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_fish_native_sf_2 %>% filter(is.na(f_log10_PE))

fishs_native_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_fish_native_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_fish_native_sf_2, area < units::set_units(1e7, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.fishs_native_PE_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishs_native_PE_map = ggplotGrob(fishs_native_PE_map)
legend.fishs_native_PE_map = ggplotGrob(legend.fishs_native_PE_map)
fishs_native_PE_map_all = arrangeGrob(fishs_native_PE_map,
                                       legend.fishs_native_PE_map,
                                       ncol = 1,
                                       layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishs_native_PE_map_all)


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
  relocate(Surf_area, .after = X1.Basin.Name)

pe_fish_extant = calcu_PE_PE_parallel(
  Tree = phylo %>% 
    keep.tip(unique(data.used_final$valid_names)),
  Comm = comm_fish_extant[,unique(data.used_final$valid_names)],
  Area = comm_fish_extant$Surf_area)

pe_fish_extant_1 = cbind(Surf_area = comm_fish_extant$Surf_area,
                         pe_fish_extant)
pe_fish_extant_sf = df_trans %>% left_join(pe_fish_extant_1, by = 'Surf_area')
colnames(pe_fish_extant_sf)
save(pe_fish_extant_sf, file = 'results/primary_results/pe_fish_all_sf.rdata')

#### extant fishs mapping 
load('results/primary_results/pe_fish_all_sf.rdata')

pe_fish_extant_sf = pe_fish_extant_sf %>% filter(PE !=0 )
pe_fish_extant_sf$log_PE = log10(pe_fish_extant_sf$PE)
# Ensure full coverage with a tiny buffer
breaks = quantile(log10(pe_fish_extant_sf$PE),
                  probs = seq(0, 1,
                              length.out = 9),
                  na.rm = TRUE)
breaks[1] = min(breaks[1], min(pe_fish_extant_sf$log_PE, na.rm = TRUE)) * 1.01
breaks[length(breaks)] = max(pe_fish_extant_sf$log_PE, na.rm = TRUE) * 0.99

pe_fish_extant_sf$f_log10_PE = cut(pe_fish_extant_sf$log_PE,
                                    breaks = breaks,
                                    include.lowest = TRUE,
                                    include.highest = T,
                                    right = T,
                                    dig.lab = 2)

# PE = 2.352844e-07
pe_fish_extant_sf_2 = pe_fish_extant_sf %>%
  st_centroid() %>% 
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2])

pe_fish_extant_sf_2 %>% filter(is.na(f_log10_PE))

fishs_extant_PE_map =
  ggplot() +
  geom_sf(data = countries, color = 'gray', fill = NA) +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = pe_fish_extant_sf, aes(fill = f_log10_PE),
          color = NA, dTolerance = 2, linewidth = 0.2, show.legend = F) +
  
  scale_fill_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE',
    guide = guide_legend(
      title = "PE",
      title.position = "top",
      label.position = "bottom",
      direction = "horizontal",
      nrow = 1
    )
  ) +
  # Circle layer for islands, colorPE by PE and sizPE by richness
  geom_point(data = subset(pe_fish_extant_sf_2, area < units::set_units(1e7, "m^2")),
             aes(x = lon, y = lat, size = f_log10_PE, color = f_log10_PE),
             shape = 21, stroke = 2, fill = NA, show.legend = F) + 
  scale_color_manual(
    values = colors2,
    na.value = 'gray90',
    name = 'PE'
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

legend.fishs_extant_PE_map = legend.func(mycolors = colors2,
                                          mylabels = round(breaks, 2)) +
  ggtitle("PE (Mya KM^-2)") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
fishs_extant_PE_map = ggplotGrob(fishs_extant_PE_map)
legend.fishs_extant_PE_map = ggplotGrob(legend.fishs_extant_PE_map)
fishs_extant_PE_map_extant = arrangeGrob(fishs_extant_PE_map,
                                          legend.fishs_extant_PE_map,
                                          ncol = 1,
                                          layout_matrix = rbind(matrix(1, 4, 10), c(NA, rep(2, 8), NA)))
plot(fishs_extant_PE_map_extant)



#### Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Dimensions for 2 columns and 4 rows
plot_width = 0.5
plot_height = 0.20

figs_pe_nat_exo = ggdraw() +
  draw_plot(plants_native_PE_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_PE_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_PE_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_PE_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_PE_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_PE_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishs_native_PE_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishs_exotic_PE_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_pe_nat_exo.png',
       #plot = figs_pe_nat_exo,
       height=40, width=25, # setting emfPlusFontToPath=TRUE to 
       res = 300,
       # ensure text looks correct on the viewing system
       units = 'cm')
figs_pe_nat_exo
dev.off() #turn off device and finalize file

#emf('figures/figs_pe_nat_exo.emf',
   # height=40, width=25, coordDPI = 1, 
   # emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
   # units = 'cm')
#figs_pe_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_pe_nat_extant = ggdraw() +
  draw_plot(plants_native_PE_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_PE_map_extant, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_PE_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_all_PE_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_PE_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_all_PE_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishs_native_PE_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishs_extant_PE_map_extant, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_pe_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_pe_nat_extant
dev.off() #turn off device and finalize file
