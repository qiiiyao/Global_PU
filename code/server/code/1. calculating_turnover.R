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




#1. Mammal turnover: calculation & mapping----
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")
is.rooted(spec_phy.3)

## 1.5 Patitioning delta ED into 5 possible ways ----
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)
sp_overlap_dat_1$presence = 1

sp_dis_5$presence = 1

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
  #i = 226
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
    
    if(phy_turn_mammal_native[region1, region2] == 0) { ## if regional pairs have no different native species,
      #then whatever two regions spread species, which did not change their phy_turnover!
      turnover_mammal_path6_mat[region1, region2] = 0
      turnover_mammal_path6_mat[region2, region1] = 0
      
    } else {
      path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                          natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
      sp_mammal_path6 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(Binomial %in% path6_sp),
                              natu_region2 %>% filter(Binomial %in% path6_sp))
      
      comm_mammal_path6 = sp_mammal_path6 %>% 
        complete(RegionID, Binomial, fill = list(presence = 0)) %>%
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
      
    }
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_mammal_path3_5 = rbind(nati_region1, 
                                nati_region2,
                                natu_region1 %>% filter(Binomial %in% path3_5_sp),
                                natu_region2 %>% filter(Binomial %in% path3_5_sp))
      
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
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_mammal_path4_7 = rbind(nati_region1, 
                                nati_region2,
                                natu_region1 %>% filter(Binomial %in% path4_7_sp),
                                natu_region2 %>% filter(Binomial %in% path4_7_sp))
      
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

save(turnover_mammal_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path3_5_mat.rdata')

save(turnover_mammal_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path6_mat.rdata')

save(turnover_mammal_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path4_7_mat.rdata')




#2. Plant turnover: calculation & mapping ----
shp.glonaf.new = st_read("data/Plants/shp_glonaf_new_eck4.shp")
phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")
load("data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = shp.glonaf.new 

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


## 2.6 Patitioning delta ED into 5 possible ways ----
colnames(df.natu.650)
df.natu.650$presence = 1

colnames(df.native.650)
df.native.650$presence = 1


df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1

load("results/primary_results/distances_beta/LCBD_plant_native.rdata")

turnover_plant_native_mat = LCBD_plant_native$beta_mat
region_pairs = combn(colnames(turnover_plant_native_mat), 2)

turnover_plant_path3_5_mat = turnover_plant_native_mat
turnover_plant_path4_7_mat = turnover_plant_native_mat
turnover_plant_path6_mat = turnover_plant_native_mat

for (i in 1:ncol(region_pairs)) {
  #i = 1
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = df.native.650 %>% filter(Region_id == region1)
  nati_region2 = df.native.650 %>% filter(Region_id == region2)
  natu_region1 = df.natu.650 %>% filter(Region_id == region1)
  natu_region2 = df.natu.650 %>% filter(Region_id == region2)
  
  nati_region1_sps = nati_region1$species
  nati_region2_sps = nati_region2$species
  natu_region1_sps = setdiff(natu_region1$species, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$species, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    
    if(turnover_plant_native_mat[region1, region2] == 0) { ## if regional pairs have no different native species,
      #then whatever two regions spread species, which did not change their phy_turnover!
      turnover_plant_path6_mat[region1, region2] = 0
      turnover_plant_path6_mat[region2, region1] = 0
      
    } else {
      
      path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                          natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
      
      sp_plant_path6 = rbind(nati_region1, 
                             nati_region2,
                             natu_region1 %>% filter(species %in% path6_sp),
                             natu_region2 %>% filter(species %in% path6_sp))
      
      comm_plant_path6 = sp_plant_path6 %>% 
        complete(Region_id, species, fill = list(presence = 0)) %>%
        dplyr::select(c('Region_id', 'species', 'presence')) %>% 
        pivot_wider(names_from = species,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('Region_id') 
      
      phylo_plant_path6 = drop.tip(phylo_big, setdiff(phylo_big$tip.label, sp_plant_path6$species))
      
      phy_turn_plant_path6 = calcu_turnover_simple(
        Tree = phylo_plant_path6,
        Comm = comm_plant_path6[,2:ncol(comm_plant_path6)])
      
      turnover_plant_path6_mat[region1, region2] = phy_turn_plant_path6[1,2]
      turnover_plant_path6_mat[region2, region1] = phy_turn_plant_path6[1,2]
    }
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_plant_path3_5 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(species %in% path3_5_sp),
                               natu_region2 %>% filter(species %in% path3_5_sp))
      
      comm_plant_path3_5 = sp_plant_path3_5 %>% 
        complete(species, Region_id, fill = list(presence = 0)) %>%
        dplyr::select(c('Region_id', 'species', 'presence')) %>% 
        pivot_wider(names_from = species,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('Region_id') 
      
      phylo_plant_path3_5 = drop.tip(phylo_big, setdiff(phylo_big$tip.label, sp_plant_path3_5$species))
      
      phy_turn_plant_path3_5 = calcu_turnover_simple(
        Tree = phylo_plant_path3_5,
        Comm = comm_plant_path3_5[,2:ncol(comm_plant_path3_5)])
      
      turnover_plant_path3_5_mat[region1, region2] = phy_turn_plant_path3_5[1,2]
      turnover_plant_path3_5_mat[region2, region1] = phy_turn_plant_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_plant_path4_7 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(species %in% path4_7_sp),
                               natu_region2 %>% filter(species %in% path4_7_sp))
      
      comm_plant_path4_7 = sp_plant_path4_7 %>% 
        complete(species, Region_id, fill = list(presence = 0)) %>%
        dplyr::select(c('Region_id', 'species', 'presence')) %>% 
        pivot_wider(names_from = species,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('Region_id') 
      
      phylo_plant_path4_7 = drop.tip(phylo_big, setdiff(phylo_big$tip.label, sp_plant_path4_7$species))
      
      phy_turn_plant_path4_7 = calcu_turnover_simple(
        Tree = phylo_plant_path4_7,
        Comm = comm_plant_path4_7[,2:ncol(comm_plant_path4_7)])
      
      turnover_plant_path4_7_mat[region1, region2] = phy_turn_plant_path4_7[1,2]
      turnover_plant_path4_7_mat[region2, region1] = phy_turn_plant_path4_7[1,2]
      
    }
    
  }
  
}


save(turnover_plant_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/plants/turnover_plant_path3_5_mat.rdata')

save(turnover_plant_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/plants/turnover_plant_path6_mat.rdata')

save(turnover_plant_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/plants/turnover_plant_path4_7_mat.rdata')


#3. Plant turnover_TDWG4: calculation & mapping ####
phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative_TDWG.Rdata")
df.native.650 = df.native.natu.species.650.nonative_TDWG[[1]]
df.natu.650 = df.native.natu.species.650.nonative_TDWG[[2]]
colnames(df.native.650) == colnames(df.natu.650)
Region_id = df_sub$Region_id


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
  complete(species, Region_id, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_sub[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  dplyr::arrange('Region_id') 

phy_turn_plant_exotic = calcu_phy_turn_parallel(Tree = phylo_plant_exotic,
                                                Comm = comm_plant_exotic,
                                                Region_posi = which(colnames(comm_plant_exotic) == 'Region_id'))

save(phy_turn_plant_exotic,
     file = 'results/primary_results/phy_turn_plant_exotic.rdata')

load('results/primary_results/phy_turn_plant_exotic.rdata')

turnover_plant_exotic_mat = phy_turn_plant_exotic$beta_mat
turnover_plant_exotic = data.frame(Region_id = colnames(turnover_plant_exotic_mat),
                                   turnover = colMeans(turnover_plant_exotic_mat, 
                                                       na.rm = T))
turnover_plant_exotic$Region_id = as.integer(turnover_plant_exotic$Region_id)

turnover_plant_exotic_sf = df_sub %>% left_join(turnover_plant_exotic,
                                                by = 'Region_id')
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
  tidyr::complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(df_sub[,c('Region_id', 'Area')],
            by = 'Region_id') %>% 
  relocate(Area, .after = Region_id) %>% 
  filter(!is.na(Region_id))

comm_plant_native = arrange(comm_plant_native, comm_plant_native$Region_id)

phy_turn_plant_native = calcu_phy_turnover_parallel(Tree = phylo_plant_native,
                                                    Comm = comm_plant_native,
                                                    Region_posi = which(colnames(comm_plant_native) == 'Region_id'))

save(phy_turn_plant_native,
     file = 'results/primary_results/phy_turn_plant_native.rdata') # run in the server

load("results/primary_results/phy_turn_plant_native.rdata") # load the calculated data from the server

turnover_plant_native_mat = phy_turn_plant_native
turnover_plant_native = data.frame(Region_id = colnames(turnover_plant_native_mat),
                                   turnover = colMeans(turnover_plant_native_mat,
                                                       na.rm = T))
turnover_plant_native$Region_id = as.integer(turnover_plant_native$Region_id)

turnover_plant_native_sf = df_sub %>% left_join(turnover_plant_native,
                                                by = 'Region_id') %>% 
  left_join((df[,c("Region_id", "Lon", "Lat")] %>% 
               st_drop_geometry),by = 'Region_id')
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
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(df_sub[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

phy_turn_plant_extant = calcu_phy_turn_parallel(Tree = phylo,
                                                Comm = comm_plant_extant,
                                                Region_posi = which(colnames(comm_plant_extant) == 'Region_id'))
save(phy_turn_plant_extant,
     file = 'results/primary_results/phy_turn_plant_extant.rdata')

load("results/primary_results/phy_turn_plant_extant.rdata") # load the calculated data from the server

turnover_plant_extant_mat = phy_turn_plant_extant$beta_mat
turnover_plant_extant = data.frame(Region_id = colnames(turnover_plant_extant_mat),
                                   turnover = colMeans(turnover_plant_extant_mat))
turnover_plant_extant$Region_id = as.integer(turnover_plant_extant$Region_id)

turnover_plant_extant_sf = df_sub %>% left_join(turnover_plant_extant,
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
turnover_plant_extant_mat = phy_turn_plant_extant$beta_mat
turnover_plant_native_mat = phy_turn_plant_native$beta_mat
turnover_plant_delta_mat = log((turnover_plant_extant_mat+0.001) /
                                 (turnover_plant_native_mat+0.001))
turnover_plant_delta = data.frame(Region_id = colnames(turnover_plant_delta_mat),
                                  delta_turnover = colMeans(turnover_plant_delta_mat))
turnover_plant_delta$Region_id = as.integer(turnover_plant_delta$Region_id)

turnover_plant_delta_sf = df_sub %>% left_join(turnover_plant_delta,
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





#5. Fish turnover: calculation & mapping ####
load("data/Fishes/data/my_phy.rdata")
is.rooted(phylo)
load("data/Fishes/data/my_data_used_final.rdata")
df_trans = st_read("data/Fishes/data/Basin042017_3119_eck4/Basin042017_3119_eck4.shp")
#save(df, file = "data/Fishes/data/Basin042017_3119.rdata")
colnames(df_trans)[which(colnames(df_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')


## 5.5 Patitioning delta ED into 5 possible ways ----
load("results/primary_results/distances_beta/LCBD_fish_native.rdata")

data.used_final_exotics$presence = 1
data.used_final_natives$presence = 1
colnames(data.used_final_exotics)

turnover_fish_native_mat = LCBD_fish_native$beta_mat
region_pairs = combn(colnames(turnover_fish_native_mat), 2)

turnover_fish_native_mat[is.na(turnover_fish_native_mat)] = 1
sum(is.na(turnover_fish_native_mat))

turnover_fish_path3_5_mat = turnover_fish_native_mat
turnover_fish_path4_7_mat = turnover_fish_native_mat
turnover_fish_path6_mat = turnover_fish_native_mat

for (i in 1:ncol(region_pairs)) {
  #i = 4
  region1 = region_pairs[1,i]
  region2 = region_pairs[2,i]
  
  nati_region1 = data.used_final_natives %>% filter(X1.Basin.Name == region1)
  nati_region2 = data.used_final_natives %>% filter(X1.Basin.Name == region2)
  natu_region1 = data.used_final_exotics %>% filter(X1.Basin.Name == region1)
  natu_region2 = data.used_final_exotics %>% filter(X1.Basin.Name == region2)
  
  nati_region1_sps = nati_region1$valid_names
  nati_region2_sps = nati_region2$valid_names
  natu_region1_sps = setdiff(natu_region1$valid_names, nati_region1_sps)
  natu_region2_sps = setdiff(natu_region2$valid_names, nati_region2_sps)
  
  natu_1_sps_ori_2 = natu_region1_sps[which(natu_region1_sps %in% nati_region2_sps)]
  natu_1_sps_ori_3 = natu_region1_sps[which(!(natu_region1_sps %in% nati_region2_sps))]
  
  natu_2_sps_ori_1 = natu_region2_sps[which(natu_region2_sps %in% nati_region1_sps)]
  natu_2_sps_ori_3 = natu_region2_sps[which(!(natu_region2_sps %in% nati_region1_sps))]
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else if (length(natu_1_sps_ori_2) > 0 |
             length(natu_2_sps_ori_1) > 0) {
    
    if(turnover_fish_native_mat[region1, region2] == 0 |
       !is.na(turnover_fish_native_mat[region1, region2] == 0)) { ## if regional pairs have no different native valid_names,
      #then whatever two regions spread valid_names, which did not change their phy_turnover!
      turnover_fish_path6_mat[region1, region2] = 0
      turnover_fish_path6_mat[region2, region1] = 0
      
    } else {
      path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                          natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
      
      sp_fish_path6 = rbind(nati_region1, 
                            nati_region2#,
                            #natu_region1 %>% filter(valid_names %in% path6_sp),
                            #natu_region2 %>% filter(valid_names %in% path6_sp)
                            )
      
      comm_fish_path6 = sp_fish_path6 %>% 
        complete(X1.Basin.Name, valid_names, fill = list(presence = 0)) %>%
        dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
        pivot_wider(names_from = valid_names,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('X1.Basin.Name') 
      
      phy_turn_fish_path6 = calcu_turnover_simple(
        Tree = phylo,
        Comm = comm_fish_path6[,2:ncol(comm_fish_path6)])
      
      turnover_fish_path6_mat[region1, region2] = phy_turn_fish_path6[1,2]
      turnover_fish_path6_mat[region2, region1] = phy_turn_fish_path6[1,2]
    }
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_fish_path3_5 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(valid_names %in% path3_5_sp),
                              natu_region2 %>% filter(valid_names %in% path3_5_sp))
      
      comm_fish_path3_5 = sp_fish_path3_5 %>% 
        complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
        dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
        pivot_wider(names_from = valid_names,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('X1.Basin.Name') 
      
      phy_turn_fish_path3_5 = calcu_turnover_simple(
        Tree = phylo,
        Comm = comm_fish_path3_5[,2:ncol(comm_fish_path3_5)])
      
      turnover_fish_path3_5_mat[region1, region2] = phy_turn_fish_path3_5[1,2]
      turnover_fish_path3_5_mat[region2, region1] = phy_turn_fish_path3_5[1,2]
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_fish_path4_7 = rbind(nati_region1, 
                              nati_region2,
                              natu_region1 %>% filter(valid_names %in% path4_7_sp),
                              natu_region2 %>% filter(valid_names %in% path4_7_sp))
      
      comm_fish_path4_7 = sp_fish_path4_7 %>% 
        complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
        dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
        pivot_wider(names_from = valid_names,
                    values_from = presence,
                    values_fn = mean) %>% 
        dplyr::arrange('X1.Basin.Name') 
      
      phy_turn_fish_path4_7 = calcu_turnover_simple(
        Tree = phylo,
        Comm = comm_fish_path4_7[,2:ncol(comm_fish_path4_7)])
      
      turnover_fish_path4_7_mat[region1, region2] = phy_turn_fish_path4_7[1,2]
      turnover_fish_path4_7_mat[region2, region1] = phy_turn_fish_path4_7[1,2]
      
    }
    
  }
  
}


save(turnover_fish_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/fishs/turnover_fish_path3_5_mat.rdata')

save(turnover_fish_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/fishs/turnover_fish_path6_mat.rdata')

save(turnover_fish_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/fishs/turnover_fish_path4_7_mat.rdata')




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

figs_turnover_nat_extant_delta = ggdraw() +
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

png(filename = 'figures/figs_turnover_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_turnover_nat_extant_delta
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


