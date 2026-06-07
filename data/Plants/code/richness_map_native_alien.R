### Clear work-space
rm(list = ls())

## Package loading
require(sf)
require(raster)
require(ggplot2)
require(dplyr)
require(tidyr)
require(data.table)

## Set up work place
setwd("D:/R projects/Global_ED")

load("data/Plants/data/df.native.natu.species.650.nonative_TDWG.Rdata")
df = st_read("data/Plants/TDWG4_Subset/TDWG4_Subset.shp")
df_trans = df

df.native.650 = df.native.natu.species.650.nonative_TDWG[[1]]
df.natu.650 = df.native.natu.species.650.nonative_TDWG[[2]]
colnames(df.native.650) == colnames(df.natu.650)

df_natu_native = rbind(df.native.650, df.natu.650)
str(df_natu_native)
richness_dat = df_natu_native %>% group_by(RegionID, status) %>% 
  summarise_at(vars(species), c(richness = function(x){length(unique(x))})) %>% 
  pivot_wider(names_from = status, values_from = richness)

colnames(richness_dat)[2:3] = c('native_richenss', 'exotic_richness')

nrow(richness_dat %>% filter(!is.na(exotic_richness)))

df_joined_merge = df_trans %>% left_join(richness_dat,
                                         by = 'RegionID') %>% 
  mutate(richness_all = exotic_richness + native_richenss)

df_joined_merge[is.na(df_joined_merge$exotic_richness),]$exotic_richness = 0
df_joined_merge[is.na(df_joined_merge$native_richenss),]$native_richenss = 0
df_joined_merge[is.na(df_joined_merge$richness_all),]$richness_all = 0


exotic_richness_map =
  ggplot() +
  #geom_sf(data = countries, color = NA, fill = "gray") +
  #geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = df_joined_merge, aes(fill = exotic_richness), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12)
  ) + 
  scale_fill_scico(palette = 'vik', direction = 1)
exotic_richness_map


save(df_joined_merge, file = 'Plants/Results_data/richness_map_data.Rdata')

### mapping the exotic and native richness of fishes for multiple basins ####