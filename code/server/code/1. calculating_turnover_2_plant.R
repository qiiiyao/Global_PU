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
source('code/functions/calculating_phy_turnover_func_2.R')

# load the background data for plotting the world map plot
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")

# load the background data for plotting the world map plot
shp.glonaf.trans = st_read("data/Plants/shp_glonaf_new_eck4.shp")

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



#2. Plant turnover: calculation & mapping ----
shp.glonaf.trans = st_read("data/Plants/shp_glonaf_new_eck4.shp")
shp.glonaf.mainland = shp.glonaf.trans %>% filter(!(Island == 1 & 
                                                      Area < 5e3)) 
plot(shp.glonaf.mainland$geometry)
glonaf_mainlands = as.character(sort(unique(shp.glonaf.mainland$Region_id)))

phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")

load("data/Plants/data/df.native.natu.species.650.nonative.Rdata")
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

##2.4 extant #####
df.native.650$presence = 1
df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1

df.extant.650 = arrange(df.extant.650, df.extant.650$Region_id)
comm_plant_extant = df.extant.650 %>% 
  #filter(Region_id %in% region_ids[1:2]) %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence)

phylo_plant_extant = drop.tip(phylo_big, setdiff(phylo_big$tip.label,
                                                 unique(df.extant.650$species)))

mat_plant_extant = as.matrix(comm_plant_extant[, colnames(comm_plant_extant) %in% phylo_plant_extant$tip.label])
storage.mode(mat_plant_extant) = "numeric"
row.names(mat_plant_extant) = sort(as.numeric(unlist(comm_plant_extant[,'Region_id'])))

#phy_turn_plant_extant = calcu_phy_turn_multiple(tree = phylo_plant_extant,
#                                          x = mat_plant_extant,
#                                      block_size = 40000)

phy_sor_plant_extant = calcu_phy_sor_multiple(tree = phylo_plant_extant,
                                             x = mat_plant_extant,
                                             block_size = 40000)

phy_jac_plant_extant = calcu_phy_jac_multiple(tree = phylo_plant_extant,
                                             x = mat_plant_extant,
                                             block_size = 40000)

phy_rlb_plant_extant = calcu_phy_rlb_multiple(tree = phylo_plant_extant,
                                             x = mat_plant_extant,
                                             block_size = 40000)

#phy_turn_plant_extant[which(phy_turn_plant_extant < 0)]
phy_sor_plant_extant[which(phy_sor_plant_extant < 0)]
phy_jac_plant_extant[which(phy_jac_plant_extant < 0)]
phy_rlb_plant_extant[which(phy_rlb_plant_extant < 0)]
phy_rlb_plant_extant[c(1:4),c(1:4)]

#save(phy_turn_plant_extant,
#    file = 'results/primary_results/distances_beta/phy_turn_plant_extant.rdata')
save(phy_sor_plant_extant,
     file = 'results/primary_results/distances_beta/phy_sor_plant_extant.rdata')
save(phy_jac_plant_extant,
     file = 'results/primary_results/distances_beta/phy_jac_plant_extant.rdata')
save(phy_rlb_plant_extant,
     file = 'results/primary_results/distances_beta/phy_rlb_plant_extant.rdata')


##2.3 natives #####
colnames(df.native.650)
df.native.650$presence = 1
df.native.650 = arrange(df.native.650, df.native.650$Region_id)

df.native.650 = df.native.650 %>% filter(species %in% phylo_big$tip.label)

mat_plant_extant[which(mat_plant_extant > 0)] = 0
mat_plant_native = mat_plant_extant

region_ids = match(df.native.650$Region_id, rownames(mat_plant_native))
sps = match(df.native.650$species, colnames(mat_plant_native))

mat_plant_native[cbind(region_ids, sps)] = df.native.650$presence

#phy_turn_plant_native = calcu_phy_turn_multiple(tree = phylo_plant_extant,
#                                          x = mat_plant_native,
#                                      block_size = 40000)

phy_sor_plant_native = calcu_phy_sor_multiple(tree = phylo_plant_extant,
                                             x = mat_plant_native,
                                             block_size = 40000)

phy_jac_plant_native = calcu_phy_jac_multiple(tree = phylo_plant_extant,
                                             x = mat_plant_native,
                                             block_size = 40000)

phy_rlb_plant_native = calcu_phy_rlb_multiple(tree = phylo_plant_extant,
                                             x = mat_plant_native,
                                             block_size = 40000)

#phy_turn_plant_native[which(phy_turn_plant_native < 0)]
phy_sor_plant_native[which(phy_sor_plant_native < 0)]
phy_jac_plant_native[which(phy_jac_plant_native < 0)]
phy_rlb_plant_native[which(phy_rlb_plant_native < 0)]
phy_rlb_plant_native[c(1:4),c(1:4)]

#save(phy_turn_plant_native,
#    file = 'results/primary_results/distances_beta/phy_turn_plant_native.rdata')
save(phy_sor_plant_native,
     file = 'results/primary_results/distances_beta/phy_sor_plant_native.rdata')
save(phy_jac_plant_native,
     file = 'results/primary_results/distances_beta/phy_jac_plant_native.rdata')
save(phy_rlb_plant_native,
     file = 'results/primary_results/distances_beta/phy_rlb_plant_native.rdata')

print('Finishing calculation of phy_beta for plants')

## 2.6 Patitioning delta ED into 5 possible ways ----
colnames(df.natu.650)
df.natu.650$presence = 1
df.natu.650 = df.natu.650 %>% filter(species %in% phylo_big$tip.label)

colnames(df.native.650)
df.native.650$presence = 1
df.native.650 = df.native.650 %>% filter(species %in% phylo_big$tip.label)

df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1

load('results/primary_results/distances_beta/phy_turn_plant_native.rdata')

turnover_plant_native_mat = phy_turn_plant_native
rm(phy_turn_plant_native)
region_pairs = combn(colnames(turnover_plant_native_mat), 2)

mat_plant_extant[which(mat_plant_extant > 0)] = 0

turnover_plant_path6_mat = turnover_plant_native_mat
turnover_plant_path3_5_mat = turnover_plant_native_mat
turnover_plant_path4_7_mat = turnover_plant_native_mat

#rm(mat_plant_extant)

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
  
  mat_plant_path6 = mat_plant_extant
  mat_plant_path3_5 = mat_plant_extant
  mat_plant_path4_7 = mat_plant_extant
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else {
    
    if (length(natu_1_sps_ori_2) > 0 |
        length(natu_2_sps_ori_1) > 0) {
    
       if(turnover_plant_native_mat[region1, region2] == 0) { ## if regional pairs have no different native species,
      #then whatever two regions spread species, which did not change their phy_turnover!
      turnover_plant_path6_mat[region1, region2] = 0
      turnover_plant_path6_mat[region2, region1] = 0
      
       } else {
      
      path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                          natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
      
      sp_plant_path6 = rbind(nati_region1, 
                             nati_region2
                             ,natu_region1 %>% filter(species %in% path6_sp),
                             natu_region2 %>% filter(species %in% path6_sp)
                             )
      
      
      region_ids = match(sp_plant_path6$Region_id, rownames(mat_plant_path6))
      sps = match(sp_plant_path6$species, colnames(mat_plant_path6))
      
      mat_plant_path6[cbind(region_ids, sps)] = sp_plant_path6$presence
      
      phy_turn_plant_path6 = calcu_phy_turn_pair(
        tree = phylo_plant_extant,
        x = mat_plant_path6[c(region1, region2),])
      
      turnover_plant_path6_mat[region1, region2] = phy_turn_plant_path6
      turnover_plant_path6_mat[region2, region1] = phy_turn_plant_path6
       }
    }
      
      if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_plant_path3_5 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(species %in% path3_5_sp),
                               natu_region2 %>% filter(species %in% path3_5_sp))
      
      
      region_ids = match(sp_plant_path3_5$Region_id, rownames(mat_plant_path3_5))
      sps = match(sp_plant_path3_5$species, colnames(mat_plant_path3_5))
      
      mat_plant_path3_5[cbind(region_ids, sps)] = sp_plant_path3_5$presence
      
      phy_turn_plant_path3_5 = calcu_phy_turn_pair(
        tree = phylo_plant_extant,
        x = mat_plant_path3_5[c(region1, region2),])
      
      turnover_plant_path3_5_mat[region1, region2] = phy_turn_plant_path3_5
      turnover_plant_path3_5_mat[region2, region1] = phy_turn_plant_path3_5
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_plant_path4_7 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(species %in% path4_7_sp),
                               natu_region2 %>% filter(species %in% path4_7_sp))
      
      
      region_ids = match(sp_plant_path4_7$Region_id, rownames(mat_plant_path4_7))
      sps = match(sp_plant_path4_7$species, colnames(mat_plant_path4_7))
      
      mat_plant_path4_7[cbind(region_ids, sps)] = sp_plant_path4_7$presence
      
      phy_turn_plant_path4_7 = calcu_phy_turn_pair(
        tree = phylo_plant_extant,
        x = mat_plant_path4_7[c(region1, region2),])
      
      turnover_plant_path4_7_mat[region1, region2] = phy_turn_plant_path4_7
      turnover_plant_path4_7_mat[region2, region1] = phy_turn_plant_path4_7
      
    }
    
  }
  
  print(paste(i, 'in', ncol(region_pairs)))
}

turnover_plant_delta_path6_mat = log((turnover_plant_path6_mat+0.001) /
                                       (turnover_plant_native_mat+0.001))
which(turnover_plant_delta_path6_mat > 0)


save(turnover_plant_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/plants/turnover_plant_path3_5_mat.rdata')

save(turnover_plant_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/plants/turnover_plant_path6_mat.rdata')

save(turnover_plant_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/plants/turnover_plant_path4_7_mat.rdata')


