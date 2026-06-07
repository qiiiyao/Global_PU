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
shp.glonaf.trans = st_read("data/Plants/shp_glonaf_new_eck4.shp")

load("results/primary_results/distances_GIFT/geodistances_grid_7.RDATA")
distances_gift = as.matrix(geodistances)
rm(geodistances)


#2. Plant turnover: calculation & mapping ----
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

#2. Build the data according to the null model that relocates exotic species from neighboring regions----
# a. find the regions with the lowest 10% number of exotics 
df.natu.650_num = df.natu.650 %>% group_by(Region_id) %>% 
  summarise_at(vars(species), length) %>% 
  complete(Region_id = shp.glonaf.trans$Region_id, fill = list(species = 0)) %>%
  arrange(species) 

top10per_posi = round(length(unique(df.natu.650_num$species)) * 0.1)
top10per = sort(unique(df.natu.650_num$species))[1:top10per_posi]
top10per_regions = df.natu.650_num %>% filter(species %in% top10per) %>%
  pull(Region_id)

## b. Precompute candidate exotic species pools for each focal region
top10per_regions = as.character(sort(top10per_regions))

plan(multicore, workers = 2)

system.time({
  region_pool_l = future_map(
    top10per_regions,
    function(region_i) {
      
      #region_i = top10per_regions[1]
      # Identify the 10 nearest neighboring regions
      neighbor_regions_1 = names(sort(distances_gift[, region_i])[1:11])
      neighbor_regions_2 = setdiff(neighbor_regions_1, region_i)
      
      # Exotic species occurring in neighboring regions
      neighbor_bird_exotic = df.extant.650 %>% 
        filter(Region_id %in% neighbor_regions_2,
               status  != "Native") %>% 
        pull(species) %>% 
        unique()
      
      # Exotic species already present in the focal region
      focal_bird_exotic = df.extant.650 %>% 
        filter(Region_id == region_i,
               status != "Native") %>% 
        pull(species) %>% 
        unique()
      
      # Native species occurring in neighboring regions
      neighbor_bird_native = df.extant.650 %>% 
        filter(Region_id %in% neighbor_regions_1,
               status == "Native") %>% 
        pull(species) %>% 
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
      candidate_dat = df.extant.650 %>% 
        filter(Region_id %in% neighbor_regions_2,
               species %in% candidate_species)
      
      list(region_i = region_i,
        candidate_species = candidate_species,
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
  exotic_plant_null_l = future_map(
    seq_len(100),
    function(repl) {
      add_l = lapply(region_pool_l, function(pool) {
        
        # Randomly select species from the candidate pool
        select_plant_exotic = sample(pool$candidate_species, size = pool$num)
        
        # Extract corresponding records
        select_plant_exotic_dat = pool$candidate_dat %>% 
          filter(species %in% select_plant_exotic) %>% 
          slice_head(n = pool$num)
        
        # Relocate selected species to the focal region
        select_plant_exotic_dat$Region_id = pool$region_i
        
        select_plant_exotic_dat
      })
      
      # Combine relocated species with the original dataset
      add_dat = data.table::rbindlist(add_l)
      sp_overlap_dat_null = rbind(df.natu.650, add_dat)
      sp_overlap_dat_null$repl = repl
      
      return(sp_overlap_dat_null)
    },
    .options = furrr_options(seed = TRUE)
  )}
)

## d. calculating phylogenentic Simpson similarity for constructed null communities
## Precompute fixed region and species sets
# All species appearing in the observed and null datasets
all_species = unique(c(df.native.650$species,
                       unlist(lapply(exotic_plant_null_l, function(x) unique(x$species)))))
all_species = intersect(all_species, phylo_big$tip.label)

# All regions appearing in the observed and null datasets
all_regions = unique(c(df.native.650$Region_id,
                       unlist(lapply(exotic_plant_null_l, function(x) unique(x$Region_id)))))
all_regions = sort(as.numeric(all_regions))

## 2. Construct the fixed phylogeny once
# Drop species that never appear in any observed or null community
phylo_plant_extant_fixed = drop.tip(phylo_big,
                                   setdiff(phylo_big$tip.label, all_species))

# Use the tree tip order as the matrix column order
species_order = phylo_plant_extant_fixed$tip.label

## 3. Construct the baseline community matrix once
df.native.650_base = df.native.650 %>%
  filter(species %in% species_order) %>%
  group_by(Region_id, species) %>%
  summarise(presence = max(presence), .groups = "drop")

native_mat = matrix(0, nrow = length(all_regions), ncol = length(species_order),
                    dimnames = list(all_regions, species_order))

idx_region = match(as.numeric(df.native.650_base$Region_id), all_regions)
idx_species = match(df.native.650_base$species, species_order)

native_mat[cbind(idx_region, idx_species)] = df.native.650_base$presence
storage.mode(native_mat) = "numeric"


options(future.globals.maxSize = 1 * 1024^3)  # 1 GB
plan(multicore, workers = 2)

batch_id = split(seq_along(exotic_plant_null_l),
                 ceiling(seq_along(exotic_plant_null_l) / 4))

phy_uni_plant_extant_null = vector("list", length(exotic_plant_null_l))

for (b in seq_along(batch_id)) {
  
  ids = batch_id[[b]]
  
  system.time({phy_uni_plant_extant_null[ids] = future_map(
    exotic_plant_null_l[ids],
    function(x) {
      
      x2 = x %>%
        filter(species %in% species_order) %>%
        group_by(Region_id, species) %>%
        summarise(presence = max(presence), .groups = "drop")
      
      mat_plant_extant = native_mat
      
      idx_region = match(x2$Region_id, all_regions)
      idx_species = match(x2$species, species_order)
      
      mat_plant_extant[cbind(idx_region, idx_species)] = pmax(
        mat_plant_extant[cbind(idx_region, idx_species)],
        x2$presence)
      
      res = calcu_phy_sim_multiple_all(
        tree = phylo_plant_extant_fixed,
        x = mat_plant_extant,
        block_size = 40000)
      
      rm(mat_plant_extant, x2)
      gc()
      
      return(res)
    },
    .options = furrr_options(seed = TRUE, scheduling = 1)
  )})
  
  gc()
}

## Save intermediate results after each batch
save(phy_uni_plant_extant_null,
     file = "results/primary_results/null_models/phy_uni_plant_extant_null.rdata")


#3. Calculating the delta_PU for the 1000 null models that relocates exotic species from neighboring regions ----
load('results/primary_results/null_models/phy_uni_plant_extant_null.rdata')
load("results/primary_results/distances_beta/phy_turn_plant_native.rdata")
turn_plant_native = data.frame(Region_id = rownames(phy_turn_plant_native),
                              turn = rowMeans(phy_turn_plant_native))
turn_plant_native$Region_id = as.integer(turn_plant_native$Region_id)

system.time({
  phy_turn_plant_delta_null = future_map(
    phy_uni_plant_extant_null,
    function(x) {
      #x = phy_uni_plant_extant_null[[1]]
      phy_turn_plant_all = x$simpson
      turn_plant_all = data.frame(Region_id = rownames(phy_turn_plant_all),
                                 turn = rowMeans(phy_turn_plant_all))
      turn_plant_all$Region_id = as.integer(turn_plant_all$Region_id)
      
      colnames(turn_plant_all)[which(colnames(turn_plant_all) == 'turn')] = 'all_turn'
      PU_plant_delta = turn_plant_native %>% 
        left_join(turn_plant_all, by = 'Region_id')
      PU_plant_delta$delta_PU = log(PU_plant_delta$all_turn/
                                     PU_plant_delta$turn)
      return(PU_plant_delta)
      
    },
    .options = furrr_options(seed = TRUE)
  )}
)

phy_turn_plant_delta_null_dat = data.table::rbindlist(phy_turn_plant_delta_null)

phy_turn_plant_delta_null_dat$rep = rep(seq_len(length(phy_turn_plant_delta_null)),
                                       each = length(unique(phy_turn_plant_delta_null_dat$Region_id)))

phy_turn_plant_delta_null_summary = phy_turn_plant_delta_null_dat %>% group_by(Region_id) %>% 
  summarise(Q1 = quantile(delta_PU, 0.25, na.rm = TRUE),
            median = quantile(delta_PU, 0.50, na.rm = TRUE),  # median
            Q3 = quantile(delta_PU, 0.75, na.rm = TRUE), .groups = "drop")

hist(phy_turn_plant_delta_null_summary$Q1)
hist(phy_turn_plant_delta_null_summary$median)
hist(phy_turn_plant_delta_null_summary$Q3)

save(phy_turn_plant_delta_null_summary,
     file = 'results/primary_results/null_models/phy_turn_plant_delta_null_summary.rdata')



