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
df = st_read("../TDWG4/TDWG4_newTibet.shp")

distribution_dat = read.csv('Birds/Distribution_Data.csv',
                        header = T)
str(distribution_dat)
num_native = length(unique((distribution_dat %>% filter(SpStatus == 'Native'))$BinomialName))
richness_dat = distribution_dat %>% group_by(RegionID, SpStatus) %>% 
  summarise_at(vars(BinomialName), c(richness = function(x){length(unique(x))})) %>% 
  pivot_wider(names_from = SpStatus, values_from = richness)

colnames(richness_dat)[2:3] = c('native_richness', 'exotic_richness')

nrow(richness_dat %>% filter(!is.na(exotic_richness)))

richness_dat[is.na(richness_dat$exotic_richness),]$exotic_richness = 0

df_joined_merge = df %>% left_join(richness_dat,
                                   by = 'RegionID') %>% 
  mutate(richness_all = native_richness + exotic_richness)

save(df_joined_merge, file = 'Birds/Results_data/richness_map_data.Rdata')


### mapping the exotic and native richness of fishes for multiple basins ####