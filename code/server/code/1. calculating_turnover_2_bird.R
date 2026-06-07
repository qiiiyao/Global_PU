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
df_trans = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")


#3. Bird turnover: calculation & mapping ####
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

##3.1 exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
exotic_distri_data$presence = 1


##3.4 extant #####
colnames(all_distri_data_c)
all_distri_data_c$presence = 1
all_distri_data_c2 = all_distri_data_c %>%
  group_by(RegionID, ScientificName) %>%
  summarise(presence = max(presence), .groups = "drop")

all_distri_data_c2 = arrange(all_distri_data_c2, all_distri_data_c2$RegionID)
comm_bird_extant = all_distri_data_c2 %>% 
  #filter(RegionID %in% RegionIDs[1:2]) %>% 
  complete(RegionID, ScientificName, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'presence')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = presence)
comm_bird_extant = as_tibble(comm_bird_extant)

phylo_bird_extant = drop.tip(phy_data, setdiff(phy_data$tip.label,
                                                 unique(all_distri_data_c$ScientificName)))

mat_bird_extant = data.matrix(comm_bird_extant[,colnames(comm_bird_extant) %in% phylo_bird_extant$tip.label])
storage.mode(mat_bird_extant) = "numeric"
row.names(mat_bird_extant) = sort(as.numeric(unlist(comm_bird_extant[,'RegionID'])))

#phy_turn_bird_extant = calcu_phy_turn_multiple(tree = phylo_bird_extant,
#                                          x = mat_bird_extant,
#                                      block_size = 40000)

phy_sor_bird_extant = calcu_phy_sor_multiple(tree = phylo_bird_extant,
                                               x = mat_bird_extant,
                                               block_size = 40000)

phy_jac_bird_extant = calcu_phy_jac_multiple(tree = phylo_bird_extant,
                                               x = mat_bird_extant,
                                               block_size = 40000)

phy_rlb_bird_extant = calcu_phy_rlb_multiple(tree = phylo_bird_extant,
                                               x = mat_bird_extant,
                                               block_size = 40000)

#phy_turn_bird_extant[which(phy_turn_bird_extant < 0)]
phy_sor_bird_extant[which(phy_sor_bird_extant < 0)]
phy_jac_bird_extant[which(phy_jac_bird_extant < 0)]
phy_rlb_bird_extant[which(phy_rlb_bird_extant < 0)]
phy_rlb_bird_extant[c(1:4),c(1:4)]

#save(phy_turn_bird_extant,
#    file = 'results/primary_results/distances_beta/phy_turn_bird_extant.rdata')
save(phy_sor_bird_extant,
     file = 'results/primary_results/distances_beta/phy_sor_bird_extant.rdata')
save(phy_jac_bird_extant,
     file = 'results/primary_results/distances_beta/phy_jac_bird_extant.rdata')
save(phy_rlb_bird_extant,
     file = 'results/primary_results/distances_beta/phy_rlb_bird_extant.rdata')


##3.3 natives #####
colnames(native_distri_data)
native_distri_data$presence = 1
native_distri_data = arrange(native_distri_data, native_distri_data$RegionID)

native_distri_data = native_distri_data %>% filter(ScientificName %in% phy_data$tip.label)

mat_bird_extant[which(mat_bird_extant > 0)] = 0
mat_bird_native = mat_bird_extant

RegionIDs = match(native_distri_data$RegionID, rownames(mat_bird_native))
sps = match(native_distri_data$ScientificName, colnames(mat_bird_native))

mat_bird_native[cbind(RegionIDs, sps)] = native_distri_data$presence

#phy_turn_bird_native = calcu_phy_turn_multiple(tree = phylo_bird_extant,
#                                          x = mat_bird_native,
#                                      block_size = 40000)

phy_sor_bird_native = calcu_phy_sor_multiple(tree = phylo_bird_extant,
                                               x = mat_bird_native,
                                               block_size = 40000)

phy_jac_bird_native = calcu_phy_jac_multiple(tree = phylo_bird_extant,
                                               x = mat_bird_native,
                                               block_size = 40000)

phy_rlb_bird_native = calcu_phy_rlb_multiple(tree = phylo_bird_extant,
                                               x = mat_bird_native,
                                               block_size = 40000)

#phy_turn_bird_native[which(phy_turn_bird_native < 0)]
phy_sor_bird_native[which(phy_sor_bird_native < 0)]
phy_jac_bird_native[which(phy_jac_bird_native < 0)]
phy_rlb_bird_native[which(phy_rlb_bird_native < 0)]
phy_rlb_bird_native[c(1:4),c(1:4)]

#save(phy_turn_bird_native,
#    file = 'results/primary_results/distances_beta/phy_turn_bird_native.rdata')
save(phy_sor_bird_native,
     file = 'results/primary_results/distances_beta/phy_sor_bird_native.rdata')
save(phy_jac_bird_native,
     file = 'results/primary_results/distances_beta/phy_jac_bird_native.rdata')
save(phy_rlb_bird_native,
     file = 'results/primary_results/distances_beta/phy_rlb_bird_native.rdata')

print('Finishing calculation of phy_beta for birds')


## 3.6 Patitioning delta ED into 5 possible ways ----
colnames(exotic_distri_data)
exotic_distri_data = exotic_distri_data %>% filter(ScientificName %in% phy_data$tip.label)
exotic_distri_data$presence = 1

colnames(native_distri_data)
native_distri_data = native_distri_data %>% filter(ScientificName %in% phy_data$tip.label)
native_distri_data$presence = 1
load('results/primary_results/distances_beta/phy_turn_bird_native.rdata')

turnover_bird_native_mat = phy_turn_bird_native
rm(phy_turn_bird_native)
region_pairs = combn(colnames(turnover_bird_native_mat), 2)

mat_bird_extant[which(mat_bird_extant > 0)] = 0

turnover_bird_path6_mat = turnover_bird_native_mat
turnover_bird_path3_5_mat = turnover_bird_native_mat
turnover_bird_path4_7_mat = turnover_bird_native_mat

#rm(mat_bird_extant)

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
  
  mat_bird_path6 = mat_bird_extant
  mat_bird_path3_5 = mat_bird_extant
  mat_bird_path4_7 = mat_bird_extant
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else {
    
    if (length(natu_1_sps_ori_2) > 0 |
        length(natu_2_sps_ori_1) > 0) {
      
      if(turnover_bird_native_mat[region1, region2] == 0) { ## if regional pairs have no different native species,
        #then whatever two regions spread species, which did not change their phy_turnover!
        turnover_bird_path6_mat[region1, region2] = 0
        turnover_bird_path6_mat[region2, region1] = 0
        
      } else {
        
        path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
        
        sp_bird_path6 = rbind(nati_region1, 
                                nati_region2
                                ,natu_region1 %>% filter(ScientificName %in% path6_sp),
                                natu_region2 %>% filter(ScientificName %in% path6_sp)
        )
        
        
        RegionIDs = match(sp_bird_path6$RegionID, rownames(mat_bird_path6))
        sps = match(sp_bird_path6$ScientificName, colnames(mat_bird_path6))
        
        mat_bird_path6[cbind(RegionIDs, sps)] = sp_bird_path6$presence
        
        phy_turn_bird_path6 = calcu_phy_turn_pair(
          tree = phylo_bird_extant,
          x = mat_bird_path6[c(region1, region2),])
        
        turnover_bird_path6_mat[region1, region2] = phy_turn_bird_path6
        turnover_bird_path6_mat[region2, region1] = phy_turn_bird_path6
      }
    }
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_bird_path3_5 = rbind(nati_region1, 
                                nati_region2,
                                natu_region1 %>% filter(ScientificName %in% path3_5_sp),
                                natu_region2 %>% filter(ScientificName %in% path3_5_sp))
      
      
      RegionIDs = match(sp_bird_path3_5$RegionID, rownames(mat_bird_path3_5))
      sps = match(sp_bird_path3_5$ScientificName, colnames(mat_bird_path3_5))
      
      mat_bird_path3_5[cbind(RegionIDs, sps)] = sp_bird_path3_5$presence
      
      phy_turn_bird_path3_5 = calcu_phy_turn_pair(
        tree = phylo_bird_extant,
        x = mat_bird_path3_5[c(region1, region2),])
      
      turnover_bird_path3_5_mat[region1, region2] = phy_turn_bird_path3_5
      turnover_bird_path3_5_mat[region2, region1] = phy_turn_bird_path3_5
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_bird_path4_7 = rbind(nati_region1, 
                                nati_region2,
                                natu_region1 %>% filter(ScientificName %in% path4_7_sp),
                                natu_region2 %>% filter(ScientificName %in% path4_7_sp))
      
      
      RegionIDs = match(sp_bird_path4_7$RegionID, rownames(mat_bird_path4_7))
      sps = match(sp_bird_path4_7$ScientificName, colnames(mat_bird_path4_7))
      
      mat_bird_path4_7[cbind(RegionIDs, sps)] = sp_bird_path4_7$presence
      
      phy_turn_bird_path4_7 = calcu_phy_turn_pair(
        tree = phylo_bird_extant,
        x = mat_bird_path4_7[c(region1, region2),])
      
      turnover_bird_path4_7_mat[region1, region2] = phy_turn_bird_path4_7
      turnover_bird_path4_7_mat[region2, region1] = phy_turn_bird_path4_7
      
    }
    
  }
  
  print(paste(i, 'in', ncol(region_pairs)))
}

turnover_bird_delta_path6_mat = log((turnover_bird_path6_mat+0.001) /
                                        (turnover_bird_native_mat+0.001))
which(turnover_bird_delta_path6_mat > 0)


save(turnover_bird_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/birds/turnover_bird_path3_5_mat.rdata')

save(turnover_bird_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/birds/turnover_bird_path6_mat.rdata')

save(turnover_bird_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/birds/turnover_bird_path4_7_mat.rdata')


