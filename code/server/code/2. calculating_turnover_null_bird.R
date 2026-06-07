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
load("code/FYI/Fan_et_al_2023_NC/data.for.shp.plots.4.Rdata")

# load the background data for plotting the world map plot
df_trans = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")


load("results/primary_results/distances_TDWG/geodistances_grid_7.RDATA")

distances_tdwg = as.matrix(geodistances)
rm(geodistances)

#1. Bird turnover: calculation & mapping ####
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


#2. Build the data according to the null model that relocates exotic species from neighboring regions----
# a. find the regions with the lowest 10% number of exotics 
exotic_distri_data_num = exotic_distri_data %>% group_by(RegionID) %>% 
  summarise_at(vars(ScientificName), length) %>% 
  complete(RegionID = df_trans$RegionID, fill = list(ScientificName = 0)) %>%
  arrange(ScientificName) 

top10per_posi = round(length(unique(exotic_distri_data_num$ScientificName)) * 0.1)
top10per = sort(unique(exotic_distri_data_num$ScientificName))[1:top10per_posi]
top10per_regions = exotic_distri_data_num %>% filter(ScientificName %in% top10per) %>%
  pull(RegionID)

## b. Precompute candidate exotic species pools for each focal region
top10per_regions = as.character(sort(top10per_regions))

region_pool_l = lapply(top10per_regions, function(region_i) {
  
  # Identify the 10 nearest neighboring regions
  neighbor_regions_1 = names(sort(distances_tdwg[, region_i])[1:11])
  neighbor_regions_2 = setdiff(neighbor_regions_1, region_i)
  
  # Exotic species occurring in neighboring regions
  neighbor_bird_exotic = all_distri_data_c %>% 
    filter(RegionID %in% neighbor_regions_2,
           SpStatus != "Native") %>% 
    pull(ScientificName) %>% 
    unique()
  
  # Exotic species already present in the focal region
  focal_bird_exotic = all_distri_data_c %>% 
    filter(RegionID == region_i,
           SpStatus != "Native") %>% 
    pull(ScientificName) %>% 
    unique()
  
  # Native species occurring in neighboring regions
  neighbor_bird_native = all_distri_data_c %>% 
    filter(RegionID %in% neighbor_regions_1,
           SpStatus == "Native") %>% 
    pull(ScientificName) %>% 
    unique()
  
  # Restrict the species pool to species that are:
  # (1) non-native in at least one neighboring region
  # (2) native in none of the neighboring regions
  # (3) absent from the focal region
  candidate_species = neighbor_bird_exotic %>%
    setdiff(neighbor_bird_native) %>%
    setdiff(focal_bird_exotic)
  
  # Skip regions with no eligible species
  if (length(candidate_species) == 0) {
    return(NULL)
  }
  
  # Randomly relocate 10% of eligible exotic species
  num = max(1, round(length(candidate_species) * 0.1))
  
  # Store all candidate records for later sampling
  candidate_dat = all_distri_data_c %>% 
    filter(RegionID %in% neighbor_regions_2,
           ScientificName %in% candidate_species)
  
  list(
    region_i = region_i,
    candidate_species = candidate_species,
    candidate_dat = candidate_dat,
    num = num
  )
})

names(region_pool_l) = top10per_regions
region_pool_l = region_pool_l[!sapply(region_pool_l, is.null)]

## c. Construct null communities
exotic_bird_null_l = vector("list", 1000)

for (r in seq_len(1000)) {
  #r = 1
  add_l = lapply(region_pool_l, function(pool) {
    
    # Randomly select species from the candidate pool
    select_bird_exotic = sample(
      pool$candidate_species,
      size = pool$num
    )
    
    # Extract corresponding records
    select_bird_exotic_dat = pool$candidate_dat %>% 
      filter(ScientificName %in% select_bird_exotic) %>% 
      slice_head(n = pool$num)
    
    # Relocate selected species to the focal region
    select_bird_exotic_dat$RegionID = pool$region_i
    
    select_bird_exotic_dat
  })
  
  # Combine relocated species with the original dataset
  add_dat = data.table::rbindlist(add_l)
  sp_overlap_dat_null = rbind(
    exotic_distri_data,
    add_dat
  )
  
  exotic_bird_null_l[[r]] = sp_overlap_dat_null
}


save(exotic_bird_null_l,
     file = 'results/primary_results/null_models/exotic_bird_null_l.rdata')

## d. calculating phylogenentic Simpson similarity for constructed null communities
load("results/primary_results/null_models/exotic_bird_null_l.rdata")

## Precompute fixed region and species sets
# All species appearing in the observed and null datasets
all_species = unique(c(native_distri_data$ScientificName,
                       unlist(lapply(exotic_bird_null_l, function(x) unique(x$ScientificName)))))
all_species = intersect(all_species, phy_data$tip.label)

# All regions appearing in the observed and null datasets
all_regions = unique(c(native_distri_data$RegionID,
                       unlist(lapply(exotic_bird_null_l, function(x) unique(x$RegionID)))))
all_regions = sort(as.numeric(all_regions))

## 2. Construct the fixed phylogeny once
# Drop species that never appear in any observed or null community
phylo_bird_extant_fixed = drop.tip(phy_data,
                                     setdiff(phy_data$tip.label, all_species))

# Use the tree tip order as the matrix column order
species_order = phylo_bird_extant_fixed$tip.label

## 3. Construct the baseline community matrix once
native_distri_data_base = native_distri_data %>%
  filter(ScientificName %in% species_order) %>%
  group_by(RegionID, ScientificName) %>%
  summarise(presence = max(presence), .groups = "drop")

native_mat = matrix(0, nrow = length(all_regions), ncol = length(species_order),
                    dimnames = list(all_regions, species_order))

idx_region = match(as.numeric(native_distri_data_base$RegionID), all_regions)
idx_species = match(native_distri_data_base$ScientificName, species_order)

native_mat[cbind(idx_region, idx_species)] = native_distri_data_base$presence
storage.mode(native_mat) = "numeric"

plan(multisession, workers = round(parallel::detectCores()*0.5))

system.time({
  phy_uni_bird_extant_null = future_map(
    exotic_bird_null_l,
    function(x) {
      
      x2 = x %>%
        filter(ScientificName %in% species_order) %>%
        group_by(RegionID, ScientificName) %>%
        summarise(presence = max(presence), .groups = "drop")
      
      mat_bird_extant = native_mat
      
      idx_region = match(as.numeric(x2$RegionID), all_regions)
      idx_species = match(x2$ScientificName, species_order)
      
      mat_bird_extant[cbind(idx_region, idx_species)] = pmax(
        mat_bird_extant[cbind(idx_region, idx_species)],
        x2$presence
      )
      
      calcu_phy_sim_multiple_all(
        tree = phylo_bird_extant_fixed,
        x = mat_bird_extant,
        block_size = 40000
      )
    },
    .options = furrr_options(seed = TRUE)
  )}
)


save(phy_uni_bird_extant_null,
     file = 'results/primary_results/null_models/phy_uni_bird_extant_null.rdata')


#3. Calculating the delta_PU for the 1000 null models that relocates exotic species from neighboring regions ----
load('results/primary_results/null_models/phy_uni_bird_extant_null.rdata')
load("results/primary_results/distances_beta/phy_turn_bird_native.rdata")
turn_bird_native = data.frame(RegionID = rownames(phy_turn_bird_native),
                                turn = rowMeans(phy_turn_bird_native))
turn_bird_native$RegionID = as.integer(turn_bird_native$RegionID)

system.time({
  phy_turn_bird_delta_null = future_map(
    phy_uni_bird_extant_null,
    function(x) {
      #x = phy_uni_bird_extant_null[[1]]
      phy_turn_bird_all = x$simpson
      turn_bird_all = data.frame(RegionID = rownames(phy_turn_bird_all),
                                   turn = rowMeans(phy_turn_bird_all))
      turn_bird_all$RegionID = as.integer(turn_bird_all$RegionID)
      
      colnames(turn_bird_all)[which(colnames(turn_bird_all) == 'turn')] = 'all_turn'
      PU_bird_delta = turn_bird_native %>% 
        left_join(turn_bird_all, by = 'RegionID')
      PU_bird_delta$delta_PU = log(PU_bird_delta$all_turn/
                                       PU_bird_delta$turn)
      return(PU_bird_delta)
      
    },
    .options = furrr_options(seed = TRUE)
  )}
)

phy_turn_bird_delta_null_dat = data.table::rbindlist(phy_turn_bird_delta_null)

phy_turn_bird_delta_null_dat$rep = rep(c(1:1000),
                                         each = length(unique(phy_turn_bird_delta_null_dat$RegionID)))

phy_turn_bird_delta_null_summary = phy_turn_bird_delta_null_dat %>% group_by(RegionID) %>% 
  summarise(Q1 = quantile(delta_PU, 0.25, na.rm = TRUE),
            median = quantile(delta_PU, 0.50, na.rm = TRUE),  # median
            Q3 = quantile(delta_PU, 0.75, na.rm = TRUE), .groups = "drop")

hist(phy_turn_bird_delta_null_summary$Q1)
hist(phy_turn_bird_delta_null_summary$median)
hist(phy_turn_bird_delta_null_summary$Q3)

save(phy_turn_bird_delta_null_summary,
     file = 'results/primary_results/null_models/phy_turn_bird_delta_null_summary.rdata')


