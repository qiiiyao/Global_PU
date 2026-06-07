### Clear work-space
rm(list = ls())

## Package loading
require(sf)
require(raster)
require(ggplot2)
require(dplyr)
require(tidyr)
require(data.table)
library(scico)

## Set up work place
setwd("D:/R projects/Multi_taxa")


# define color gradients
colors1 = scico::scico(n=8, palette = "lajolla") 
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
scale_manual_richness = list(
  scale_fill_manual(values = colors1, drop = FALSE),
  scale_color_manual(values = colors1, drop = FALSE)
)

load("data/Mammals/Sample_code/data.for.shp.plots.4.Rdata")

# Import the primary result
df = st_read("data/TDWG4/TDWG4_newTibet.shp")


#### mammals mapping 
load("data/Mammals/Results_data/richness_map_data.Rdata")
mammals_richness_map =
  ggplot() +
  geom_sf(data = countries, color = NA, fill = "gray") +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = df_joined_merge, aes(fill = richness_all), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Mammals')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5),
  ) + 
  scale_fill_scico(palette = 'vik', direction = 1, 
                   name = 'Richness')
mammals_richness_map

#### plants mapping
load("data/Plants/Results_data/richness_map_data.Rdata")
plants_richness_map =
  ggplot() +
  #geom_sf(data = countries, color = NA, fill = "gray") +
  #geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = df_joined_merge, aes(fill = richness_all), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Plants')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5)
  ) + 
  scale_fill_scico(palette = 'vik', direction = 1, 
                   name = 'Richness', 
                   breaks = seq(min(df_joined_merge$richness_all),
                                  max(df_joined_merge$richness_all),
                                length.out = 3 ))
plants_richness_map


#### birds mapping
load("Birds/Results_data/richness_map_data.Rdata")
birds_richness_map =
  ggplot() +
  geom_sf(data = countries, color = NA, fill = "gray") +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = df_joined_merge, aes(fill = richness_all), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Birds')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5)
  ) + 
  scale_fill_scico(palette = 'vik', direction = 1, 
                   name = 'Richness')
birds_richness_map


#### fishes mapping
load("Fishes/Results_data/richness_map_data.Rdata")

fishes_richness_map =
  ggplot() +
  geom_sf(data = countries, color = NA, fill = "gray") +
  geom_sf(data = bb, color = "gray", fill = NA) +
  geom_sf(data = df_joined_merge, aes(fill = Total_richness), color = NA, dTolerance = 2) +
  coord_sf(crs = "+proj=eck4",expand = FALSE) +
  theme_void() +
  ggtitle('Fishes')+
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, 
                              hjust = 0.5)
  ) + 
  scale_fill_scico(palette = 'vik', direction = 1, 
                   name = 'Richness', 
                   breaks = seq(min(df_joined_merge$Total_richness),
                                max(df_joined_merge$Total_richness),
                                length.out = 3 ))
fishes_richness_map

#### Merge all plots
library(ggpubr)
library(export)
library(devEMF)
multi_taxa = ggarrange(plants_richness_map,
                       mammals_richness_map,
                       birds_richness_map,
                       fishes_richness_map,
                       labels = c('a', 'b',
                                  'c', 'd'),
                       nrow = 2, ncol = 2)

ggsave('figures/richness_multi_taxa.jpeg',
       plot = multi_taxa,
       height=9*1.8, width=13*1.8,
       units = 'cm')

emf('figures/richness_multi_taxa.emf',
    height=9*1.8, width=13*1.8, coordDPI = 300, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
multi_taxa
dev.off() #turn off device and finalize file



