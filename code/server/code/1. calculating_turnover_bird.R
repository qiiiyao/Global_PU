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

setwd("~/my_pc/Global_ED")
source('code/functions/calculating_phy_turnover_func.R')

# load the background data for plotting the world map plot
# load the background data for plotting the world map plot
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = df
#df_trans$area = st_area(df_trans)
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")

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

## 4.5 Patitioning delta ED into 5 possible ways ----
load("results/primary_results/distances_beta/LCBD_bird_native.rdata")

exotic_distri_data$present = 1
colnames(native_distri_data)
native_distri_data$present = 1

turnover_bird_native_mat = LCBD_bird_native$beta_mat
region_pairs = combn(colnames(turnover_bird_native_mat), 2)

turnover_bird_path3_5_mat = turnover_bird_native_mat
turnover_bird_path4_7_mat = turnover_bird_native_mat
turnover_bird_path6_mat = turnover_bird_native_mat

for (i in 1:ncol(region_pairs)) {
  #i = 1
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = native_distri_data %>% filter(RegionID == region1)
  nati_region2 = native_distri_data %>% filter(RegionID == region2)
  natu_region1 = exotic_distri_data %>% filter(RegionID == region1)
  natu_region2 = exotic_distri_data %>% filter(RegionID == region2)
  
  nati_region1_sps = nati_region1$ScientificName
  nati_region2_sps = nati_region2$ScientificName
  natu_region1_sps = setdiff(natu_region1$ScientificName, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$ScientificName, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    
    if(turnover_bird_native_mat[region1, region2] == 0) { ## if regional pairs have no different native ScientificName,
      #then whatever two regions spread ScientificName, which did not change their phy_turnover!
      turnover_bird_path6_mat[region1, region2] = 0
      turnover_bird_path6_mat[region2, region1] = 0
      
    } else {
      path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                          natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
      
      sp_bird_path6 = rbind(nati_region1, 
                            nati_region2,
                            natu_region1 %>% filter(ScientificName %in% path6_sp),
                            natu_region2 %>% filter(ScientificName %in% path6_sp))
      
      comm_bird_path6 = sp_bird_path6 %>% 
        complete(RegionID, ScientificName, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'ScientificName', 'presence')) %>% 
        pivot_wider(names_from = ScientificName,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_turn_bird_path6 = calcu_turnover_simple(
        Tree = phy_data,
        Comm = comm_bird_path6[,2:ncol(comm_bird_path6)])
      
      turnover_bird_path6_mat[region1, region2] = phy_turn_bird_path6[1,2]
      turnover_bird_path6_mat[region2, region1] = phy_turn_bird_path6[1,2]
    }
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_bird_path3_5 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(ScientificName %in% path3_5_sp),
                              natu_region2 %>% filter(ScientificName %in% path3_5_sp))
      
      comm_bird_path3_5 = sp_bird_path3_5 %>% 
        complete(ScientificName, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'ScientificName', 'presence')) %>% 
        pivot_wider(names_from = ScientificName,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_turn_bird_path3_5 = calcu_turnover_simple(
        Tree = phy_data,
        Comm = comm_bird_path3_5[,2:ncol(comm_bird_path3_5)])
      
      turnover_bird_path3_5_mat[region1, region2] = phy_turn_bird_path3_5[1,2]
      turnover_bird_path3_5_mat[region2, region1] = phy_turn_bird_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_bird_path4_7 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(ScientificName %in% path4_7_sp),
                              natu_region2 %>% filter(ScientificName %in% path4_7_sp))
      
      comm_bird_path4_7 = sp_bird_path4_7 %>% 
        complete(ScientificName, RegionID, fill = list(presence = 0)) %>%
        dplyr::select(c('RegionID', 'ScientificName', 'presence')) %>% 
        pivot_wider(names_from = ScientificName,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('RegionID') 
      
      phy_turn_bird_path4_7 = calcu_turnover_simple(
        Tree = phy_data,
        Comm = comm_bird_path4_7[,2:ncol(comm_bird_path4_7)])
      
      turnover_bird_path4_7_mat[region1, region2] = phy_turn_bird_path4_7[1,2]
      turnover_bird_path4_7_mat[region2, region1] = phy_turn_bird_path4_7[1,2]
      
    }
    
  }
  
}


save(turnover_bird_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/birds/turnover_bird_path3_5_mat.rdata')

save(turnover_bird_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/birds/turnover_bird_path6_mat.rdata')

save(turnover_bird_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/birds/turnover_bird_path4_7_mat.rdata')


