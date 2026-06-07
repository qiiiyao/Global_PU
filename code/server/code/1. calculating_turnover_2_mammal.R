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


#1. Mammal turnover: calculation & mapping----
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")
is.rooted(spec_phy.3)


##1.1 exotics-----
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)
sp_overlap_dat_1$presence = 1

##1.4 extant #####
sp_dis_5$presence = 1
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
colnames(sp_mammal_all)

sp_mammal_all2 = sp_mammal_all %>%
  group_by(RegionID, Binomial) %>%
  summarise(presence = max(presence), .groups = "drop")

sp_mammal_all = arrange(sp_mammal_all, sp_mammal_all$RegionID)
comm_mammal_extant = sp_mammal_all2 %>% 
  #filter(RegionID %in% RegionIDs[1:2]) %>% 
  complete(RegionID, Binomial, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(names_from = Binomial,
              values_from = presence)
comm_mammal_extant = as_tibble(comm_mammal_extant)

phylo_mammal_extant = drop.tip(spec_phy.3, setdiff(spec_phy.3$tip.label,
                                                 unique(sp_mammal_all$Binomial)))

mat_mammal_extant = data.matrix(comm_mammal_extant[,
                                colnames(comm_mammal_extant) %in% phylo_mammal_extant$tip.label])
storage.mode(mat_mammal_extant) = "numeric"
row.names(mat_mammal_extant) = sort(as.numeric(unlist(comm_mammal_extant[,'RegionID'])))

#phy_turn_mammal_extant = calcu_phy_turn_multiple(tree = phylo_mammal_extant,
   #                                          x = mat_mammal_extant,
    #                                      block_size = 40000)

phy_sor_mammal_extant = calcu_phy_sor_multiple(tree = phylo_mammal_extant,
                                                 x = mat_mammal_extant,
                                                 block_size = 40000)

phy_jac_mammal_extant = calcu_phy_jac_multiple(tree = phylo_mammal_extant,
                                               x = mat_mammal_extant,
                                               block_size = 40000)

phy_rlb_mammal_extant = calcu_phy_rlb_multiple(tree = phylo_mammal_extant,
                                               x = mat_mammal_extant,
                                               block_size = 40000)

#phy_turn_mammal_extant[which(phy_turn_mammal_extant < 0)]
phy_sor_mammal_extant[which(phy_sor_mammal_extant < 0)]
phy_jac_mammal_extant[which(phy_jac_mammal_extant < 0)]
phy_rlb_mammal_extant[which(phy_rlb_mammal_extant < 0)]
phy_rlb_mammal_extant[c(1:4),c(1:4)]

#save(phy_turn_mammal_extant,
#    file = 'results/primary_results/distances_beta/phy_turn_mammal_extant.rdata')
save(phy_sor_mammal_extant,
     file = 'results/primary_results/distances_beta/phy_sor_mammal_extant.rdata')
save(phy_jac_mammal_extant,
     file = 'results/primary_results/distances_beta/phy_jac_mammal_extant.rdata')
save(phy_rlb_mammal_extant,
     file = 'results/primary_results/distances_beta/phy_rlb_mammal_extant.rdata')


##1.3 natives #####
colnames(sp_dis_5)
sp_dis_5$presence = 1
sp_dis_5 = arrange(sp_dis_5, sp_dis_5$RegionID)

sp_dis_5 = sp_dis_5 %>% filter(Binomial %in% spec_phy.3$tip.label)

mat_mammal_extant[which(mat_mammal_extant > 0)] = 0
mat_mammal_native = mat_mammal_extant

RegionIDs = match(sp_dis_5$RegionID, rownames(mat_mammal_native))
sps = match(sp_dis_5$Binomial, colnames(mat_mammal_native))

mat_mammal_native[cbind(RegionIDs, sps)] = sp_dis_5$presence

#phy_turn_mammal_native = calcu_phy_turn_multiple(tree = phylo_mammal_extant,
#                                          x = mat_mammal_native,
#                                      block_size = 40000)

phy_sor_mammal_native = calcu_phy_sor_multiple(tree = phylo_mammal_extant,
                                               x = mat_mammal_native,
                                               block_size = 40000)

phy_jac_mammal_native = calcu_phy_jac_multiple(tree = phylo_mammal_extant,
                                               x = mat_mammal_native,
                                               block_size = 40000)

phy_rlb_mammal_native = calcu_phy_rlb_multiple(tree = phylo_mammal_extant,
                                               x = mat_mammal_native,
                                               block_size = 40000)

#phy_turn_mammal_native[which(phy_turn_mammal_native < 0)]
phy_sor_mammal_native[which(phy_sor_mammal_native < 0)]
phy_jac_mammal_native[which(phy_jac_mammal_native < 0)]
phy_rlb_mammal_native[which(phy_rlb_mammal_native < 0)]
phy_rlb_mammal_native[c(1:4),c(1:4)]

#save(phy_turn_mammal_native,
#    file = 'results/primary_results/distances_beta/phy_turn_mammal_native.rdata')
save(phy_sor_mammal_native,
     file = 'results/primary_results/distances_beta/phy_sor_mammal_native.rdata')
save(phy_jac_mammal_native,
     file = 'results/primary_results/distances_beta/phy_jac_mammal_native.rdata')
save(phy_rlb_mammal_native,
     file = 'results/primary_results/distances_beta/phy_rlb_mammal_native.rdata')

print('Finishing calculation of phy_beta for mammals')

## 1.6 Patitioning delta ED into 5 possible ways ----
colnames(sp_overlap_dat_2)
sp_overlap_dat_2 = sp_overlap_dat_2 %>% filter(Binomial %in% spec_phy.3$tip.label)

colnames(sp_dis_6)
sp_dis_6 = sp_dis_6 %>% filter(Binomial %in% spec_phy.3$tip.label)

sp_mammal_all = rbind(sp_dis_6, sp_overlap_dat_2)
colnames(sp_mammal_all)
sp_mammal_all$presence = 1

load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')

turnover_mammal_native_mat = phy_turn_mammal_native
rm(phy_turn_mammal_native)
region_pairs = combn(colnames(turnover_mammal_native_mat), 2)

mat_mammal_extant[which(mat_mammal_extant > 0)] = 0

turnover_mammal_path6_mat = turnover_mammal_native_mat
turnover_mammal_path3_5_mat = turnover_mammal_native_mat
turnover_mammal_path4_7_mat = turnover_mammal_native_mat

#rm(mat_mammal_extant)

for (i in 1:ncol(region_pairs)) {
  #i = 1
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
  
  mat_mammal_path6 = mat_mammal_extant
  mat_mammal_path3_5 = mat_mammal_extant
  mat_mammal_path4_7 = mat_mammal_extant
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else {
    
    if (length(natu_1_sps_ori_2) > 0 |
        length(natu_2_sps_ori_1) > 0) {
      
      if(turnover_mammal_native_mat[region1, region2] == 0) { ## if regional pairs have no different native species,
        #then whatever two regions spread species, which did not change their phy_turnover!
        turnover_mammal_path6_mat[region1, region2] = 0
        turnover_mammal_path6_mat[region2, region1] = 0
        
      } else {
        
        path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
        
        sp_mammal_path6 = rbind(nati_region1, 
                               nati_region2
                               ,natu_region1 %>% filter(Binomial %in% path6_sp),
                               natu_region2 %>% filter(Binomial %in% path6_sp)
        )
        
        
        RegionIDs = match(sp_mammal_path6$RegionID, rownames(mat_mammal_path6))
        sps = match(sp_mammal_path6$Binomial, colnames(mat_mammal_path6))
        
        mat_mammal_path6[cbind(RegionIDs, sps)] = sp_mammal_path6$presence
        
        phy_turn_mammal_path6 = calcu_phy_turn_pair(
          tree = phylo_mammal_extant,
          x = mat_mammal_path6[c(region1, region2),])
        
        turnover_mammal_path6_mat[region1, region2] = phy_turn_mammal_path6
        turnover_mammal_path6_mat[region2, region1] = phy_turn_mammal_path6
      }
    }
    
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_mammal_path3_5 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(Binomial %in% path3_5_sp),
                               natu_region2 %>% filter(Binomial %in% path3_5_sp))
      
      
      RegionIDs = match(sp_mammal_path3_5$RegionID, rownames(mat_mammal_path3_5))
      sps = match(sp_mammal_path3_5$Binomial, colnames(mat_mammal_path3_5))
      
      mat_mammal_path3_5[cbind(RegionIDs, sps)] = sp_mammal_path3_5$presence
      
      phy_turn_mammal_path3_5 = calcu_phy_turn_pair(
        tree = phylo_mammal_extant,
        x = mat_mammal_path3_5[c(region1, region2),])
      
      turnover_mammal_path3_5_mat[region1, region2] = phy_turn_mammal_path3_5
      turnover_mammal_path3_5_mat[region2, region1] = phy_turn_mammal_path3_5
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_mammal_path4_7 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(Binomial %in% path4_7_sp),
                               natu_region2 %>% filter(Binomial %in% path4_7_sp))
      
      
      RegionIDs = match(sp_mammal_path4_7$RegionID, rownames(mat_mammal_path4_7))
      sps = match(sp_mammal_path4_7$Binomial, colnames(mat_mammal_path4_7))
      
      mat_mammal_path4_7[cbind(RegionIDs, sps)] = sp_mammal_path4_7$presence
      
      phy_turn_mammal_path4_7 = calcu_phy_turn_pair(
        tree = phylo_mammal_extant,
        x = mat_mammal_path4_7[c(region1, region2),])
      
      turnover_mammal_path4_7_mat[region1, region2] = phy_turn_mammal_path4_7
      turnover_mammal_path4_7_mat[region2, region1] = phy_turn_mammal_path4_7
      
    }
    
  }
  
  print(paste(i, 'in', ncol(region_pairs)))
}

log((phy_turn_mammal_path6+0.001) /
      (turnover_mammal_native_mat[region1, region2]+0.001))

log((phy_turn_mammal_path4_7+0.001) /
      (turnover_mammal_native_mat[region1, region2]+0.001))

turnover_mammal_delta_path6_mat = log((turnover_mammal_path6_mat+0.001) /
                                        (turnover_mammal_native_mat+0.001))


turnover_mammal_delta_path6_mat = log((turnover_mammal_path6_mat+0.001) /
                                       (turnover_mammal_native_mat+0.001))
which(turnover_mammal_delta_path6_mat > 0)


save(turnover_mammal_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path3_5_mat.rdata')

save(turnover_mammal_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path6_mat.rdata')

save(turnover_mammal_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/mammals/turnover_mammal_path4_7_mat.rdata')



## 1.7 Patitioning delta PU into 5 possible ways ----
load("results/primary_results/distances_beta/phy_turn_mammal_native.rdata")
load("results/primary_results/distances_beta/phy_turn_mammal_extant.rdata")

colnames(sp_overlap_dat_2)
sp_overlap_dat_2 = sp_overlap_dat_2 %>% filter(Binomial %in% spec_phy.3$tip.label)

colnames(sp_dis_6)
sp_dis_6 = sp_dis_6 %>% filter(Binomial %in% spec_phy.3$tip.label)

sp_mammal_all = rbind(sp_dis_6, sp_overlap_dat_2)
colnames(sp_mammal_all)
sp_mammal_all$presence = 1

load('results/primary_results/distances_beta/phy_turn_mammal_native.rdata')

region_pairs = combn(colnames(turnover_mammal_native_mat), 2)

mat_mammal_extant[which(mat_mammal_extant > 0)] = 0

regions = colnames(turnover_mammal_native_mat)

phy_turn_mammal_path1 = phy_turn_mammal_native
phy_turn_mammal_path2 = phy_turn_mammal_native
phy_turn_mammal_path3 = phy_turn_mammal_native
phy_turn_mammal_path4 = phy_turn_mammal_native
phy_turn_mammal_path5 = phy_turn_mammal_native

gc()
islands = c('41', '517')

sp_dis_6 = sp_dis_6 %>% filter(Binomial %in% spec_phy.3$tip.label)
sp_dis_6 = sp_dis_6 %>% filter(Binomial %in% spec_phy.3$tip.label)
sp_overlap_dat_2 = sp_overlap_dat_2 %>% filter(Binomial %in% spec_phy.3$tip.label)

delta_pu_mammal_dat = data.frame()

plan(multisession, workers = round(parallel::detectCores()*0.2))

system.time({
  delta_pu_mammal_list = future_map(
    regions,
    function(x) {
      #x = regions[1]
      region1 = x
      
      nati_region1 = sp_dis_6 %>% filter(RegionID == region1)
      natu_region1 = sp_overlap_dat_2 %>% filter(RegionID == region1)
      
      nati_region1_sps = nati_region1$Binomial
      natu_region1_sps = setdiff(natu_region1$Binomial, nati_region1_sps)
      natu_region1_splist = split(natu_region1, natu_region1$Binomial)
      
      natu_from_region1 = sp_overlap_dat_2 %>% filter(Binomial %in% nati_region1_sps)
      
      natu_no_region1 = sp_overlap_dat_2 %>% filter(RegionID != region1 &  
                                                      !(Binomial %in% nati_region1_sps))
      natu_no_region1_splist = split(natu_no_region1, natu_no_region1$Binomial)
      
      mat_mammal_path1 = mat_mammal_extant
      mat_mammal_path2 = mat_mammal_extant
      mat_mammal_path3 = mat_mammal_extant
      mat_mammal_path4 = mat_mammal_extant
      mat_mammal_path5 = mat_mammal_extant
      mat_mammal_path2_4 = mat_mammal_extant
      mat_mammal_path3_5 = mat_mammal_extant
      
      ## Pathway1: Naturalizing species that are native in focal region (R1) to remaining regions
      if (nrow(natu_from_region1) > 0) {
        
        sp_mammal_path1 = rbind(sp_dis_6, 
                                natu_from_region1)
        
        RegionIDs = match(sp_mammal_path1$RegionID, rownames(mat_mammal_path1))
        sps = match(sp_mammal_path1$Binomial, colnames(mat_mammal_path1))
        
        mat_mammal_path1[cbind(RegionIDs, sps)] = sp_mammal_path1$presence
        
        phy_uniq_mammal_path1 = calcu_phy_turn_focal_mean(
          tree = phylo_mammal_extant,
          x = mat_mammal_path1,
          focal_region = region1)
        
        delta_pu_mammal_path1 = phy_uniq_mammal_path1 - mean(phy_turn_mammal_native[,region1])
        
        
      } else {
        delta_pu_mammal_path1 = 0
      }
      
      ## Pathways 2&4
      if (length(natu_region1_splist) > 0) {
        
        phy_uniq_mammal_path2_4_splist = future_map(
          natu_region1_splist,
          function(x) {
            sp_mammal_path = rbind(sp_dis_6, 
                                   x)
            
            RegionIDs = match(sp_mammal_path$RegionID, rownames(mat_mammal_path3_5))
            sps = match(sp_mammal_path$Binomial, colnames(mat_mammal_path3_5))
            
            mat_mammal_path2_4[cbind(RegionIDs, sps)] = sp_mammal_path$presence
            
            phy_uniq_mammal_path = calcu_phy_turn_focal_mean(
              tree = phylo_mammal_extant,
              x = mat_mammal_path2_4,
              focal_region = region1)
            delta_pu = phy_uniq_mammal_path - mean(phy_turn_mammal_native[,region1])
            dat = data.frame(species = unique(x$Binomial),
                             delta_pu = delta_pu)
            return(dat)
          },
          .options = furrr_options(seed = TRUE)
        )
        
        #mean(phy_turn_mammal_extant[,region1])
        
        phy_uniq_mammal_path2_4_dat = data.table::rbindlist(phy_uniq_mammal_path2_4_splist)
        
        ## Pathway 2: Naturalizing species originating from common regions to R1
        mammal_path2_sps = phy_uniq_mammal_path2_4_dat %>% filter(delta_pu < 0) %>% 
          pull(species) %>% unique()
        natu_region1_path2 = natu_region1 %>% filter(Binomial %in% mammal_path2_sps)
        sp_mammal_path2 = rbind(sp_dis_6, 
                                natu_region1_path2)
        
        RegionIDs = match(sp_mammal_path2$RegionID, rownames(mat_mammal_path2))
        sps = match(sp_mammal_path2$Binomial, colnames(mat_mammal_path2))
        
        mat_mammal_path2[cbind(RegionIDs, sps)] = sp_mammal_path2$presence
        
        phy_uniq_mammal_path2 = calcu_phy_turn_focal_mean(
          tree = phylo_mammal_extant,
          x = mat_mammal_path2,
          focal_region = region1)
        delta_pu_mammal_path2 = phy_uniq_mammal_path2 - mean(phy_turn_mammal_native[,region1])
        
        
        ## Pathway 4: Naturalizing species originating from unique regions to R1
        mammal_path4_sps = phy_uniq_mammal_path2_4_dat %>% filter(delta_pu > 0) %>% 
          pull(species) %>% unique()
        natu_region1_path4 = natu_region1 %>% filter(Binomial %in% mammal_path4_sps)
        sp_mammal_path4 = rbind(sp_dis_6, 
                                natu_region1_path4)
        
        RegionIDs = match(sp_mammal_path4$RegionID, rownames(mat_mammal_path4))
        sps = match(sp_mammal_path4$Binomial, colnames(mat_mammal_path4))
        
        mat_mammal_path4[cbind(RegionIDs, sps)] = sp_mammal_path4$presence
        
        phy_uniq_mammal_path4 = calcu_phy_turn_focal_mean(
          tree = phylo_mammal_extant,
          x = mat_mammal_path4,
          focal_region = region1)
        delta_pu_mammal_path4 = phy_uniq_mammal_path4 - mean(phy_turn_mammal_native[,region1])
        
        
      } else {
        delta_pu_mammal_path2 = 0
        delta_pu_mammal_path4 = 0
      }
      
      ## Pathways 3&5
      if (length(natu_no_region1_splist) > 0) {
        
        phy_uniq_mammal_path3_5_splist = future_map(
          natu_no_region1_splist,
          function(x) {
            sp_mammal_path = rbind(sp_dis_6, 
                                   x)
            
            RegionIDs = match(sp_mammal_path$RegionID, rownames(mat_mammal_path3_5))
            sps = match(sp_mammal_path$Binomial, colnames(mat_mammal_path3_5))
            
            mat_mammal_path3_5[cbind(RegionIDs, sps)] = sp_mammal_path$presence
            
            phy_uniq_mammal_path = calcu_phy_turn_focal_mean(
              tree = phylo_mammal_extant,
              x = mat_mammal_path3_5,
              focal_region = region1)
            delta_pu = phy_uniq_mammal_path - mean(phy_turn_mammal_native[,region1])
            dat = data.frame(species = unique(x$Binomial),
                             delta_pu = delta_pu)
            return(dat)
          },
          .options = furrr_options(seed = TRUE)
        )
        
        #mean(phy_turn_mammal_extant[,region1])
        
        phy_uniq_mammal_path3_5_dat = data.table::rbindlist(phy_uniq_mammal_path3_5_splist)
        
        ## Pathway3: Naturalizing close relatives of R1 native to other regions to other remaining ones
        mammal_path3_sps = phy_uniq_mammal_path3_5_dat %>% filter(delta_pu < 0) %>% 
          pull(species) %>% unique()
        natu_region1_path3 = natu_no_region1 %>% filter(Binomial %in% mammal_path3_sps)
        sp_mammal_path3 = rbind(sp_dis_6, 
                                natu_region1_path3)
        
        RegionIDs = match(sp_mammal_path3$RegionID, rownames(mat_mammal_path3))
        sps = match(sp_mammal_path3$Binomial, colnames(mat_mammal_path3))
        
        mat_mammal_path3[cbind(RegionIDs, sps)] = sp_mammal_path3$presence
        
        phy_uniq_mammal_path3 = calcu_phy_turn_focal_mean(
          tree = phylo_mammal_extant,
          x = mat_mammal_path3,
          focal_region = region1)
        delta_pu_mammal_path3 = phy_uniq_mammal_path3 - mean(phy_turn_mammal_native[,region1])
        
        
        ## Pathway5: Naturalizing distant relatives of R1 native to other regions to other remaining ones
        mammal_path5_sps = phy_uniq_mammal_path3_5_dat %>% filter(delta_pu > 0) %>% 
          pull(species) %>% unique()
        natu_region1_path5 = natu_no_region1 %>% filter(Binomial %in% mammal_path5_sps)
        sp_mammal_path5 = rbind(sp_dis_6, 
                                natu_region1_path5)
        
        RegionIDs = match(sp_mammal_path5$RegionID, rownames(mat_mammal_path5))
        sps = match(sp_mammal_path5$Binomial, colnames(mat_mammal_path5))
        
        mat_mammal_path5[cbind(RegionIDs, sps)] = sp_mammal_path5$presence
        
        phy_uniq_mammal_path5 = calcu_phy_turn_focal_mean(
          tree = phylo_mammal_extant,
          x = mat_mammal_path5,
          focal_region = region1)
        delta_pu_mammal_path5 = phy_uniq_mammal_path5 - mean(phy_turn_mammal_native[,region1])
        
        
      } else {
        delta_pu_mammal_path3 = 0
        delta_pu_mammal_path5 = 0
      }
      
      delta_pu_mammal_dat1 = data.frame(region = region1,
                                        delta_pu_path1 = delta_pu_mammal_path1,
                                        delta_pu_path2 = delta_pu_mammal_path2,
                                        delta_pu_path3 = delta_pu_mammal_path3,
                                        delta_pu_path4 = delta_pu_mammal_path4,
                                        delta_pu_path5 = delta_pu_mammal_path5)
    },
    .options = furrr_options(seed = TRUE)
  )}
)

delta_pu_mammal_dat = data.table::rbindlist(delta_pu_mammal_list)

save(delta_pu_mammal_dat,
     file = 'results/primary_results/null_models/delta_pu_mammal_dat.rdata')


