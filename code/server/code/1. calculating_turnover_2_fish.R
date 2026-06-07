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


#5. Fish turnover: calculation & mapping ####
load("data/Fishes/data/my_phy.rdata")
is.rooted(phylo)
load("data/Fishes/data/my_data_used_final.rdata")

colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')


basin_trans = st_read("data/Fishes/data/Basin042017_3119_eck4/Basin042017_3119_eck4.shp")
colnames(basin_trans)[which(colnames(basin_trans) == 'BasinName')] = 'X1.Basin.Name'
basin_mainland = basin_trans %>% filter(!(Island == 1 & 
                                            Area < 5e3))

#plot(basin_mainland$geometry)
#plot(basin_trans$geometry)
basin_mainlands = as.character(sort(unique(basin_mainland$RegionID)))


##5.1 exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
data.used_final_exotics$presence = 1


##5.3 extant #####
data.used_final$presence = 1

comm_fish_extant = data.used_final %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(basin_trans[,c('X1.Basin.Name')], 
            #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
            by = 'X1.Basin.Name') %>% 
  arrange('X1.Basin.Name') 

mat_fish_extant = as.matrix(comm_fish_extant[, colnames(comm_fish_extant) %in% phylo$tip.label])
storage.mode(mat_fish_extant) = "numeric"
row.names(mat_fish_extant) = sort(as.character(unlist(comm_fish_extant[,'X1.Basin.Name'])))

mat_fish_extant = mat_fish_extant[sort(unique(data.used_final_natives$X1.Basin.Name)),]


#phy_turn_fish_extant = calcu_phy_turn_multiple(tree = phylo,
#                                          x = mat_fish_extant,
#                                      block_size = 40000)
phy_sims_fish_core_extant = calcu_phy_sim_multiple_core(tree = phylo,
                                            x = mat_fish_extant,
                                            block_size = 40000)

save(phy_sims_fish_core_extant,
     file = 'results/primary_results/distances_beta/phy_sims_fish_core_extant.rdata')

load('results/primary_results/distances_beta/phy_sims_fish_core_extant.rdata')

pd_region = phy_sims_fish_core_extant$pd_region
pd_tot_mat = phy_sims_fish_core_extant$total
  
# 1. PD_i + PD_j
pd_sum_mat = outer(pd_region, pd_region, "+")

# 2. shared branch length
shared_mat = pd_sum_mat - pd_tot_mat

# 3_1. min(PD_i, PD_j)
pd_min_mat = outer(pd_region, pd_region, pmin)
dimnames(pd_min_mat) = dimnames(pd_tot_mat)

# 3_2. asymmetric matrix(PD_i, PD_j)
pd_focal_mat = matrix(data = rep(pd_region, length(pd_region)), 
                      nrow = length(pd_region), ncol = length(pd_region))
# rows are focal regions, whereas columns are compared region
dimnames(pd_focal_mat) = dimnames(pd_tot_mat)

# 4_1. phylogenetic turnover
phy_turn_fish_extant = 1 - (shared_mat / pd_min_mat)
phy_turn_fish_extant[which(phy_turn_fish_extant < 0 & phy_turn_fish_extant > -1e-12)] = 0
diag(phy_turn_fish_extant) = 0

# 4_2. phylogenetic sorensen
phy_sor_fish_extant = 1 - ((2 * shared_mat) / pd_sum_mat)
phy_sor_fish_extant[which(phy_sor_fish_extant < 0 & phy_sor_fish_extant > -1e-12)] = 0
diag(phy_sor_fish_extant) = 0

# 4_3. phylogenetic jaccard
phy_jac_fish_extant = 1 - (shared_mat / pd_tot_mat)
phy_jac_fish_extant[which(phy_jac_fish_extant < 0 & phy_jac_fish_extant > -1e-12)] = 0
diag(phy_jac_fish_extant) = 0

# 4_4. phylogenetic ruggiero
phy_rlb_fish_extant = 1 - (shared_mat / pd_focal_mat)
phy_rlb_fish_extant[which(phy_rlb_fish_extant < 0 & phy_rlb_fish_extant > -1e-12)] = 0
diag(phy_rlb_fish_extant) = 0

save(phy_turn_fish_extant,
    file = 'results/primary_results/distances_beta/phy_turn_fish_extant.rdata')
save(phy_sor_fish_extant,
     file = 'results/primary_results/distances_beta/phy_sor_fish_extant.rdata')
save(phy_jac_fish_extant,
     file = 'results/primary_results/distances_beta/phy_jac_fish_extant.rdata')
save(phy_rlb_fish_extant,
     file = 'results/primary_results/distances_beta/phy_rlb_fish_extant.rdata')



##5.2 natives #####
colnames(data.used_final_natives)
data.used_final_natives$presence = 1
data.used_final_natives = arrange(data.used_final_natives, data.used_final_natives$X1.Basin.Name)

data.used_final_natives = data.used_final_natives %>% filter(valid_names %in% phylo$tip.label)

mat_fish_extant[which(mat_fish_extant > 0)] = 0
mat_fish_native = mat_fish_extant

X1.Basin.Name = match(data.used_final_natives$X1.Basin.Name, rownames(mat_fish_native))
sps = match(data.used_final_natives$valid_names, colnames(mat_fish_native))

mat_fish_native[cbind(X1.Basin.Name, sps)] = data.used_final_natives$presence

#phy_turn_fish_native = calcu_phy_turn_multiple(tree = phylo,
#                                          x = mat_fish_native,
#                                      block_size = 40000)
phy_sims_fish_core_native = calcu_phy_sim_multiple_core(tree = phylo,
                                                        x = mat_fish_native,
                                                        block_size = 40000)

save(phy_sims_fish_core_native,
     file = 'results/primary_results/distances_beta/phy_sims_fish_core_native.rdata')

load('results/primary_results/distances_beta/phy_sims_fish_core_native.rdata')

pd_region = phy_sims_fish_core_native$pd_region
pd_tot_mat = phy_sims_fish_core_native$total

# 1. PD_i + PD_j
pd_sum_mat = outer(pd_region, pd_region, "+")

# 2. shared branch length
shared_mat = pd_sum_mat - pd_tot_mat

# 3_1. min(PD_i, PD_j)
pd_min_mat = outer(pd_region, pd_region, pmin)
dimnames(pd_min_mat) = dimnames(pd_tot_mat)

# 3_2. asymmetric matrix(PD_i, PD_j)
pd_focal_mat = matrix(data = rep(pd_region, length(pd_region)), 
                      nrow = length(pd_region), ncol = length(pd_region))
# rows are focal regions, whereas columns are compared region
dimnames(pd_focal_mat) = dimnames(pd_tot_mat)

# 4_1. phylogenetic turnover
phy_turn_fish_native = 1 - (shared_mat / pd_min_mat)
phy_turn_fish_native[which(phy_turn_fish_native < 0 & phy_turn_fish_native > -1e-12)] = 0
diag(phy_turn_fish_native) = 0

# 4_2. phylogenetic sorensen
phy_sor_fish_native = 1 - ((2 * shared_mat) / pd_sum_mat)
phy_sor_fish_native[which(phy_sor_fish_native < 0 & phy_sor_fish_native > -1e-12)] = 0
diag(phy_sor_fish_native) = 0

# 4_3. phylogenetic jaccard
phy_jac_fish_native = 1 - (shared_mat / pd_tot_mat)
phy_jac_fish_native[which(phy_jac_fish_native < 0 & phy_jac_fish_native > -1e-12)] = 0
diag(phy_jac_fish_native) = 0

# 4_4. phylogenetic ruggiero
phy_rlb_fish_native = 1 - (shared_mat / pd_focal_mat)
phy_rlb_fish_native[which(phy_rlb_fish_native < 0 & phy_rlb_fish_native > -1e-12)] = 0
diag(phy_rlb_fish_native) = 0

save(phy_turn_fish_native,
     file = 'results/primary_results/distances_beta/phy_turn_fish_native.rdata')
save(phy_sor_fish_native,
     file = 'results/primary_results/distances_beta/phy_sor_fish_native.rdata')
save(phy_jac_fish_native,
     file = 'results/primary_results/distances_beta/phy_jac_fish_native.rdata')
save(phy_rlb_fish_native,
     file = 'results/primary_results/distances_beta/phy_rlb_fish_native.rdata')



print('Finishing calculation of phy_beta for fishes')

## 5.5 Patitioning delta ED into 5 possible ways ----
colnames(data.used_final_exotics)
data.used_final_exotics$presence = 1
data.used_final_exotics = data.used_final_exotics %>% filter(valid_names %in% phylo$tip.label)

colnames(data.used_final_natives)
data.used_final_natives$presence = 1
data.used_final_natives = data.used_final_natives %>% filter(valid_names %in% phylo$tip.label)

data.used_final = rbind(data.used_final_natives, data.used_final_exotics)
colnames(data.used_final)
data.used_final$presence = 1

load('results/primary_results/distances_beta/phy_turn_fish_native.rdata')

turnover_fish_native_mat = phy_turn_fish_native
rm(phy_turn_fish_native)
region_pairs = combn(colnames(turnover_fish_native_mat), 2)

mat_fish_extant[which(mat_fish_extant > 0)] = 0

turnover_fish_path6_mat = turnover_fish_native_mat
turnover_fish_path3_5_mat = turnover_fish_native_mat
turnover_fish_path4_7_mat = turnover_fish_native_mat

#rm(mat_fish_extant)

for (i in 1:ncol(region_pairs)) {
  #i = 252210
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
  
  mat_fish_path6 = mat_fish_extant
  mat_fish_path3_5 = mat_fish_extant
  mat_fish_path4_7 = mat_fish_extant
  
  if (length(natu_region1_sps) + length(natu_region2_sps) == 0) {
    print('no natu sps for both regions')
    
  } else  {
    
    if (length(natu_1_sps_ori_2) > 0 |
        length(natu_2_sps_ori_1) > 0) {
      
       if(turnover_fish_native_mat[region1, region2] == 0) { ## if regional pairs have no different native species,
          #then whatever two regions spread species, which did not change their phy_turnover!
           turnover_fish_path6_mat[region1, region2] = 0
           turnover_fish_path6_mat[region2, region1] = 0
      
       } else {
      
      path6_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_2)],
                          natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_1)]))
      
      sp_fish_path6 = rbind(nati_region1, 
                             nati_region2
                             ,natu_region1 %>% filter(valid_names %in% path6_sp),
                             natu_region2 %>% filter(valid_names %in% path6_sp)
      )
      
      
      X1.Basin.Names = match(sp_fish_path6$X1.Basin.Name, rownames(mat_fish_path6))
      sps = match(sp_fish_path6$valid_names, colnames(mat_fish_path6))
      
      mat_fish_path6[cbind(X1.Basin.Names, sps)] = sp_fish_path6$presence
      
      phy_turn_fish_path6 = calcu_phy_turn_pair(
        tree = phylo,
        x = mat_fish_path6[c(region1, region2),])
      
      turnover_fish_path6_mat[region1, region2] = phy_turn_fish_path6
      turnover_fish_path6_mat[region2, region1] = phy_turn_fish_path6
      }   
   }
  
    if (length(natu_1_sps_ori_3) > 0 &
        length(natu_2_sps_ori_3) > 0) {
      
      path3_5_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_fish_path3_5 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(valid_names %in% path3_5_sp),
                               natu_region2 %>% filter(valid_names %in% path3_5_sp))
      
      
      X1.Basin.Names = match(sp_fish_path3_5$X1.Basin.Name, rownames(mat_fish_path3_5))
      sps = match(sp_fish_path3_5$valid_names, colnames(mat_fish_path3_5))
      
      mat_fish_path3_5[cbind(X1.Basin.Names, sps)] = sp_fish_path3_5$presence
      
      phy_turn_fish_path3_5 = calcu_phy_turn_pair(
        tree = phylo,
        x = mat_fish_path3_5[c(region1, region2),])
      
      turnover_fish_path3_5_mat[region1, region2] = phy_turn_fish_path3_5
      turnover_fish_path3_5_mat[region2, region1] = phy_turn_fish_path3_5
      
      
    } else if (length(natu_1_sps_ori_3) > 0 |
               length(natu_2_sps_ori_3) > 0) {
      
      path4_7_sp = unique(c(natu_region1_sps[which(natu_region1_sps %in% natu_1_sps_ori_3)],
                            natu_region2_sps[which(natu_region2_sps %in% natu_2_sps_ori_3)]))
      
      sp_fish_path4_7 = rbind(nati_region1, 
                               nati_region2,
                               natu_region1 %>% filter(valid_names %in% path4_7_sp),
                               natu_region2 %>% filter(valid_names %in% path4_7_sp))
      
      
      X1.Basin.Names = match(sp_fish_path4_7$X1.Basin.Name, rownames(mat_fish_path4_7))
      sps = match(sp_fish_path4_7$valid_names, colnames(mat_fish_path4_7))
      
      mat_fish_path4_7[cbind(X1.Basin.Names, sps)] = sp_fish_path4_7$presence
      
      phy_turn_fish_path4_7 = calcu_phy_turn_pair(
        tree = phylo,
        x = mat_fish_path4_7[c(region1, region2),])
      
      turnover_fish_path4_7_mat[region1, region2] = phy_turn_fish_path4_7
      turnover_fish_path4_7_mat[region2, region1] = phy_turn_fish_path4_7
      
    }
  }
  print(paste(i, 'in', ncol(region_pairs)))
}

turnover_fish_delta_path6_mat = log((turnover_fish_path6_mat+0.001) /
                                      (turnover_fish_native_mat+0.001))
which(turnover_fish_delta_path6_mat > 0)


save(turnover_fish_path3_5_mat,
     file = 'results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path3_5_mat.rdata')

save(turnover_fish_path6_mat,
     file = 'results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path6_mat.rdata')

save(turnover_fish_path4_7_mat,
     file = 'results/primary_results/partitioning_possible_ways/fishes/turnover_fish_path4_7_mat.rdata')




