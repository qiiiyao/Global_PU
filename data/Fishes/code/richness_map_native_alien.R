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
setwd("D:/Globe taxa distribution/Multi_taxa")
df = st_read("Fishes/3119 basins.shp/3119 basins.shp")

richness_dat = read.csv('Fishes/Basin_richness_summary.csv',
                        header = T)
nrow(richness_dat %>% filter(!is.na(Exotic_richness)))

colnames(richness_dat)[which(colnames(richness_dat) == 'X1.Basin.Name')] = 'BasinName'

df_joined_merge = df %>% left_join(richness_dat,
                                   by = 'BasinName')

save(df_joined_merge, file = 'Fishes/Results_data/richness_map_data.Rdata')


### mapping the exotic and native richness of fishes for multiple basins ####