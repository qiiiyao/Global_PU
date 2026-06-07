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

library(furrr)
library(future)

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("~/my_pc/Global_ED")
source('code/functions/calculating_phy_turnover_func_2.R')

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

##1.1 exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
exotic_distri_data$presence = 1


##1.4 extant #####
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


## 1.7 Patitioning delta PU into 5 possible ways ----
load("results/primary_results/distances_beta/phy_turn_bird_native.rdata")
load("results/primary_results/distances_beta/phy_turn_bird_extant.rdata")

colnames(exotic_distri_data)
exotic_distri_data = exotic_distri_data %>% filter(ScientificName %in% phy_data$tip.label)
exotic_distri_data$presence = 1

colnames(native_distri_data)
native_distri_data = native_distri_data %>% filter(ScientificName %in% phy_data$tip.label)
native_distri_data$presence = 1

mat_bird_extant[which(mat_bird_extant > 0)] = 0

regions = colnames(phy_turn_bird_native)

gc()
islands = c('41', '517')

plan(multisession, workers = round(parallel::detectCores()*0.2))

system.time({
  delta_pu_bird_list = future_map(
    regions,
    function(x) {
      #x = regions[1]
      region1 = x
      
      nati_region1 = native_distri_data %>% filter(RegionID == region1)
      natu_region1 = exotic_distri_data %>% filter(RegionID == region1)
      
      nati_region1_sps = nati_region1$ScientificName
      natu_region1_sps = setdiff(natu_region1$ScientificName, nati_region1_sps)
      natu_region1_splist = split(natu_region1, natu_region1$ScientificName)
      
      natu_from_region1 = exotic_distri_data %>% filter(ScientificName %in% nati_region1_sps)
      
      natu_no_region1 = exotic_distri_data %>% filter(RegionID != region1 &  
                                                      !(ScientificName %in% nati_region1_sps))
      natu_no_region1_splist = split(natu_no_region1, natu_no_region1$ScientificName)
      
      mat_bird_path1 = mat_bird_extant
      mat_bird_path2 = mat_bird_extant
      mat_bird_path3 = mat_bird_extant
      mat_bird_path4 = mat_bird_extant
      mat_bird_path5 = mat_bird_extant
      mat_bird_path2_4 = mat_bird_extant
      mat_bird_path3_5 = mat_bird_extant
      
      ## Pathway1: Naturalizing species that are native in focal region (R1) to remaining regions
      if (nrow(natu_from_region1) > 0) {
        
        sp_bird_path1 = rbind(native_distri_data, 
                                natu_from_region1)
        
        RegionIDs = match(sp_bird_path1$RegionID, rownames(mat_bird_path1))
        sps = match(sp_bird_path1$ScientificName, colnames(mat_bird_path1))
        
        mat_bird_path1[cbind(RegionIDs, sps)] = sp_bird_path1$presence
        
        phy_uniq_bird_path1 = calcu_phy_turn_focal_mean(
          tree = phylo_bird_extant,
          x = mat_bird_path1,
          focal_region = region1)
        
        delta_pu_bird_path1 = phy_uniq_bird_path1 - mean(phy_turn_bird_native[,region1])
        
        
      } else {
        delta_pu_bird_path1 = 0
      }
      
      ## Pathways 2&4
      if (length(natu_region1_splist) > 0) {
        
        phy_uniq_bird_path2_4_splist = future_map(
          natu_region1_splist,
          function(x) {
            sp_bird_path = rbind(native_distri_data, 
                                   x)
            
            RegionIDs = match(sp_bird_path$RegionID, rownames(mat_bird_path3_5))
            sps = match(sp_bird_path$ScientificName, colnames(mat_bird_path3_5))
            
            mat_bird_path2_4[cbind(RegionIDs, sps)] = sp_bird_path$presence
            
            phy_uniq_bird_path = calcu_phy_turn_focal_mean(
              tree = phylo_bird_extant,
              x = mat_bird_path2_4,
              focal_region = region1)
            delta_pu = phy_uniq_bird_path - mean(phy_turn_bird_native[,region1])
            dat = data.frame(species = unique(x$ScientificName),
                             delta_pu = delta_pu)
            return(dat)
          },
          .options = furrr_options(seed = TRUE)
        )
        
        #mean(phy_turn_bird_extant[,region1])
        
        phy_uniq_bird_path2_4_dat = data.table::rbindlist(phy_uniq_bird_path2_4_splist)
        
        ## Pathway 2: Naturalizing species originating from common regions to R1
        bird_path2_sps = phy_uniq_bird_path2_4_dat %>% filter(delta_pu < 0) %>% 
          pull(species) %>% unique()
        natu_region1_path2 = natu_region1 %>% filter(ScientificName %in% bird_path2_sps)
        sp_bird_path2 = rbind(native_distri_data, 
                                natu_region1_path2)
        
        RegionIDs = match(sp_bird_path2$RegionID, rownames(mat_bird_path2))
        sps = match(sp_bird_path2$ScientificName, colnames(mat_bird_path2))
        
        mat_bird_path2[cbind(RegionIDs, sps)] = sp_bird_path2$presence
        
        phy_uniq_bird_path2 = calcu_phy_turn_focal_mean(
          tree = phylo_bird_extant,
          x = mat_bird_path2,
          focal_region = region1)
        delta_pu_bird_path2 = phy_uniq_bird_path2 - mean(phy_turn_bird_native[,region1])
        
        
        ## Pathway 4: Naturalizing species originating from unique regions to R1
        bird_path4_sps = phy_uniq_bird_path2_4_dat %>% filter(delta_pu > 0) %>% 
          pull(species) %>% unique()
        natu_region1_path4 = natu_region1 %>% filter(ScientificName %in% bird_path4_sps)
        sp_bird_path4 = rbind(native_distri_data, 
                                natu_region1_path4)
        
        RegionIDs = match(sp_bird_path4$RegionID, rownames(mat_bird_path4))
        sps = match(sp_bird_path4$ScientificName, colnames(mat_bird_path4))
        
        mat_bird_path4[cbind(RegionIDs, sps)] = sp_bird_path4$presence
        
        phy_uniq_bird_path4 = calcu_phy_turn_focal_mean(
          tree = phylo_bird_extant,
          x = mat_bird_path4,
          focal_region = region1)
        delta_pu_bird_path4 = phy_uniq_bird_path4 - mean(phy_turn_bird_native[,region1])
        
        
      } else {
        delta_pu_bird_path2 = 0
        delta_pu_bird_path4 = 0
      }
      
      ## Pathways 3&5
      if (length(natu_no_region1_splist) > 0) {
        
        phy_uniq_bird_path3_5_splist = future_map(
          natu_no_region1_splist,
          function(x) {
            sp_bird_path = rbind(native_distri_data, 
                                   x)
            
            RegionIDs = match(sp_bird_path$RegionID, rownames(mat_bird_path3_5))
            sps = match(sp_bird_path$ScientificName, colnames(mat_bird_path3_5))
            
            mat_bird_path3_5[cbind(RegionIDs, sps)] = sp_bird_path$presence
            
            phy_uniq_bird_path = calcu_phy_turn_focal_mean(
              tree = phylo_bird_extant,
              x = mat_bird_path3_5,
              focal_region = region1)
            delta_pu = phy_uniq_bird_path - mean(phy_turn_bird_native[,region1])
            dat = data.frame(species = unique(x$ScientificName),
                             delta_pu = delta_pu)
            return(dat)
          },
          .options = furrr_options(seed = TRUE)
        )
        
        #mean(phy_turn_bird_extant[,region1])
        
        phy_uniq_bird_path3_5_dat = data.table::rbindlist(phy_uniq_bird_path3_5_splist)
        
        ## Pathway3: Naturalizing close relatives of R1 native to other regions to other remaining ones
        bird_path3_sps = phy_uniq_bird_path3_5_dat %>% filter(delta_pu < 0) %>% 
          pull(species) %>% unique()
        natu_region1_path3 = natu_no_region1 %>% filter(ScientificName %in% bird_path3_sps)
        sp_bird_path3 = rbind(native_distri_data, 
                                natu_region1_path3)
        
        RegionIDs = match(sp_bird_path3$RegionID, rownames(mat_bird_path3))
        sps = match(sp_bird_path3$ScientificName, colnames(mat_bird_path3))
        
        mat_bird_path3[cbind(RegionIDs, sps)] = sp_bird_path3$presence
        
        phy_uniq_bird_path3 = calcu_phy_turn_focal_mean(
          tree = phylo_bird_extant,
          x = mat_bird_path3,
          focal_region = region1)
        delta_pu_bird_path3 = phy_uniq_bird_path3 - mean(phy_turn_bird_native[,region1])
        
        
        ## Pathway5: Naturalizing distant relatives of R1 native to other regions to other remaining ones
        bird_path5_sps = phy_uniq_bird_path3_5_dat %>% filter(delta_pu > 0) %>% 
          pull(species) %>% unique()
        natu_region1_path5 = natu_no_region1 %>% filter(ScientificName %in% bird_path5_sps)
        sp_bird_path5 = rbind(native_distri_data, 
                                natu_region1_path5)
        
        RegionIDs = match(sp_bird_path5$RegionID, rownames(mat_bird_path5))
        sps = match(sp_bird_path5$ScientificName, colnames(mat_bird_path5))
        
        mat_bird_path5[cbind(RegionIDs, sps)] = sp_bird_path5$presence
        
        phy_uniq_bird_path5 = calcu_phy_turn_focal_mean(
          tree = phylo_bird_extant,
          x = mat_bird_path5,
          focal_region = region1)
        delta_pu_bird_path5 = phy_uniq_bird_path5 - mean(phy_turn_bird_native[,region1])
        
        
      } else {
        delta_pu_bird_path3 = 0
        delta_pu_bird_path5 = 0
      }
      
      delta_pu_bird_dat1 = data.frame(region = region1,
                                        delta_pu_path1 = delta_pu_bird_path1,
                                        delta_pu_path2 = delta_pu_bird_path2,
                                        delta_pu_path3 = delta_pu_bird_path3,
                                        delta_pu_path4 = delta_pu_bird_path4,
                                        delta_pu_path5 = delta_pu_bird_path5)
    },
    .options = furrr_options(seed = TRUE)
  )}
)

delta_pu_bird_dat = data.table::rbindlist(delta_pu_bird_list)

save(delta_pu_bird_dat,
     file = 'results/primary_results/null_models/delta_pu_bird_dat.rdata')


