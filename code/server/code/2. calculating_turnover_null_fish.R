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

load("results/primary_results/distances_basins/geodistances_grid_7.RDATA")
distances_basin = as.matrix(geodistances)
rm(geodistances)


#1. Fish turnover: calculation & mapping ####
load("data/Fishes/data/my_phy.rdata")
is.rooted(phylo)
load("data/Fishes/data/my_data_used_final.rdata")

colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final$presence = 1
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


#2. Build the data according to the null model that relocates exotic valid_names from neighboring regions----
# a. find the regions with the lowest 10% number of exotics 
data.used_final_exotics_num = data.used_final_exotics %>% group_by(X1.Basin.Name) %>% 
  summarise_at(vars(valid_names), length) %>% 
  complete(X1.Basin.Name = basin_trans$X1.Basin.Name, fill = list(valid_names = 0)) %>%
  arrange(valid_names) 

top10per_posi = round(length(unique(data.used_final_exotics_num$valid_names)) * 0.1)
top10per = sort(unique(data.used_final_exotics_num$valid_names))[1:top10per_posi]
top10per_regions = data.used_final_exotics_num %>% filter(valid_names %in% top10per) %>%
  pull(X1.Basin.Name)

## b. Precompute candidate exotic valid_names pools for each focal region
top10per_regions = as.character(sort(top10per_regions))

plan(multicore, workers = 4)

system.time({
  region_pool_l = future_map(
    top10per_regions,
    function(region_i) {
      
      #region_i = top10per_regions[1]
      # Identify the 10 nearest neighboring regions
      neighbor_regions_1 = names(sort(distances_basin[, region_i])[1:11])
      neighbor_regions_2 = setdiff(neighbor_regions_1, region_i)
      
      # Exotic valid_names occurring in neighboring regions
      neighbor_bird_exotic = data.used_final %>% 
        filter(X1.Basin.Name %in% neighbor_regions_2,
               X3.Native.Exotic.Status  != "native") %>% 
        pull(valid_names) %>% 
        unique()
      
      # Exotic valid_names already present in the focal region
      focal_bird_exotic = data.used_final %>% 
        filter(X1.Basin.Name == region_i,
               X3.Native.Exotic.Status != "native") %>% 
        pull(valid_names) %>% 
        unique()
      
      # Native valid_names occurring in neighboring regions
      neighbor_bird_native = data.used_final %>% 
        filter(X1.Basin.Name %in% neighbor_regions_1,
               X3.Native.Exotic.Status == "native") %>% 
        pull(valid_names) %>% 
        unique()
      
      # Restrict the valid_names pool to valid_names that are:
      # (1) non-native in at least one neighboring region
      # (2) native in none of the neighboring regions
      # (3) absent from the focal region
      candidate_valid_names = neighbor_bird_exotic %>%
        setdiff(neighbor_bird_native) %>%
        setdiff(focal_bird_exotic)
      
      # Skip regions with no eligible valid_names
      if (length(candidate_valid_names) == 0) {
        return(NULL)
      }
      
      # Randomly relocate 10% of eligible exotic valid_names
      num = max(1, round(length(candidate_valid_names) * 0.1))
      
      # Store all candidate records for later sampling
      candidate_dat = data.used_final %>% 
        filter(X1.Basin.Name %in% neighbor_regions_2,
               valid_names %in% candidate_valid_names)
      
      list(region_i = region_i,
           candidate_valid_names = candidate_valid_names,
           candidate_dat = candidate_dat,
           num = num)
    },
    .options = furrr_options(seed = TRUE)
  )}
)

names(region_pool_l) = top10per_regions
region_pool_l = region_pool_l[!sapply(region_pool_l, is.null)]


## c. Construct null communities
system.time({
  exotic_fish_null_l = future_map(
    seq_len(100),
    function(repl) {
      add_l = lapply(region_pool_l, function(pool) {
        
        # Randomly select valid_names from the candidate pool
        select_fish_exotic = sample(pool$candidate_valid_names, size = pool$num)
        
        # Extract corresponding records
        select_fish_exotic_dat = pool$candidate_dat %>% 
          filter(valid_names %in% select_fish_exotic) %>% 
          slice_head(n = pool$num)
        
        # Relocate selected valid_names to the focal region
        select_fish_exotic_dat$X1.Basin.Name = pool$region_i
        
        select_fish_exotic_dat
      })
      
      # Combine relocated valid_names with the original dataset
      add_dat = data.table::rbindlist(add_l)
      sp_overlap_dat_null = rbind(data.used_final_exotics, add_dat)
      sp_overlap_dat_null$repl = repl
      
      return(sp_overlap_dat_null)
    },
    .options = furrr_options(seed = TRUE)
  )}
)

rm(region_pool_l)


## d. calculating phylogenentic Simpson similarity for constructed null communities
## Precompute fixed region and valid_names sets
# All valid_names appearing in the observed and null datasets
all_valid_names = unique(c(data.used_final_natives$valid_names,
                       unlist(lapply(exotic_fish_null_l, function(x) unique(x$valid_names)))))
all_valid_names = intersect(all_valid_names, phylo$tip.label)

# All regions appearing in the observed and null datasets
all_regions = unique(c(data.used_final_natives$X1.Basin.Name,
                       unlist(lapply(exotic_fish_null_l, function(x) unique(x$X1.Basin.Name)))))
all_regions = sort(all_regions)

## 2. Construct the fixed phylogeny once
# Drop valid_names that never appear in any observed or null community
phylo_fish_extant_fixed = drop.tip(phylo,
                                    setdiff(phylo$tip.label, all_valid_names))

# Use the tree tip order as the matrix column order
valid_names_order = phylo_fish_extant_fixed$tip.label

## 3. Construct the baseline community matrix once
data.used_final_natives_base = data.used_final_natives %>%
  filter(valid_names %in% valid_names_order) %>%
  group_by(X1.Basin.Name, valid_names) %>%
  summarise(presence = max(presence), .groups = "drop")

native_mat = matrix(0, nrow = length(all_regions), ncol = length(valid_names_order),
                    dimnames = list(all_regions, valid_names_order))

idx_region = match(data.used_final_natives_base$X1.Basin.Name, all_regions)
idx_valid_names = match(data.used_final_natives_base$valid_names, valid_names_order)

native_mat[cbind(idx_region, idx_valid_names)] = data.used_final_natives_base$presence
storage.mode(native_mat) = "numeric"


options(future.globals.maxSize = 2 * 1024^3)  # 2 GB

plan(multicore, workers = 4)

batch_id = split(seq_along(exotic_fish_null_l),
                  ceiling(seq_along(exotic_fish_null_l) / 5))

phy_uni_fish_extant_null = vector("list", length(exotic_fish_null_l))

for (b in seq_along(batch_id)) {
  
  ids = batch_id[[b]]
  
  system.time({phy_uni_fish_extant_null[ids] = future_map(
    exotic_fish_null_l[ids],
    function(x) {
      
      x2 = x %>%
        filter(valid_names %in% valid_names_order) %>%
        group_by(X1.Basin.Name, valid_names) %>%
        summarise(presence = max(presence), .groups = "drop")
      
      mat_fish_extant = native_mat
      
      idx_region = match(x2$X1.Basin.Name, all_regions)
      idx_valid_names = match(x2$valid_names, valid_names_order)
      
      mat_fish_extant[cbind(idx_region, idx_valid_names)] = pmax(
        mat_fish_extant[cbind(idx_region, idx_valid_names)],
        x2$presence)
      
      res = calcu_phy_sim_multiple_all(
        tree = phylo_fish_extant_fixed,
        x = mat_fish_extant,
        block_size = 40000)
      
      rm(mat_fish_extant, x2)
      gc()
      
      return(res)
    },
    .options = furrr_options(seed = TRUE, scheduling = 1)
  )})
  
  gc()
}

save(phy_uni_fish_extant_null,
     file = "results/primary_results/null_models/phy_uni_fish_extant_null.rdata")


#3. Calculating the delta_PU for the 1000 null models that relocates exotic valid_names from neighboring regions ----
load('results/primary_results/null_models/phy_uni_fish_extant_null.rdata')
load("results/primary_results/distances_beta/phy_turn_fish_native.rdata")
turn_fish_native = data.frame(X1.Basin.Name = rownames(phy_turn_fish_native),
                               turn = rowMeans(phy_turn_fish_native))

system.time({
  phy_turn_fish_delta_null = future_map(
    phy_uni_fish_extant_null,
    function(x) {
      #x = phy_uni_fish_extant_null[[1]]
      phy_turn_fish_all = x$simpson
      turn_fish_all = data.frame(X1.Basin.Name = rownames(phy_turn_fish_all),
                                  turn = rowMeans(phy_turn_fish_all))
      
      colnames(turn_fish_all)[which(colnames(turn_fish_all) == 'turn')] = 'all_turn'
      PU_fish_delta = turn_fish_native %>% 
        left_join(turn_fish_all, by = 'X1.Basin.Name')
      PU_fish_delta$delta_PU = log(PU_fish_delta$all_turn/
                                      PU_fish_delta$turn)
      return(PU_fish_delta)
      
    },
    .options = furrr_options(seed = TRUE)
  )}
)

phy_turn_fish_delta_null_dat = data.table::rbindlist(phy_turn_fish_delta_null)

phy_turn_fish_delta_null_dat$rep = rep(seq_len(length(phy_turn_fish_delta_null)),
                                        each = length(unique(phy_turn_fish_delta_null_dat$X1.Basin.Name)))

phy_turn_fish_delta_null_summary = phy_turn_fish_delta_null_dat %>% group_by(X1.Basin.Name) %>% 
  summarise(Q1 = quantile(delta_PU, 0.25, na.rm = TRUE),
            median = quantile(delta_PU, 0.50, na.rm = TRUE),  # median
            Q3 = quantile(delta_PU, 0.75, na.rm = TRUE), .groups = "drop")

hist(phy_turn_fish_delta_null_summary$Q1)
hist(phy_turn_fish_delta_null_summary$median)
hist(phy_turn_fish_delta_null_summary$Q3)

save(phy_turn_fish_delta_null_summary,
     file = 'results/primary_results/null_models/phy_turn_fish_delta_null_summary.rdata')



