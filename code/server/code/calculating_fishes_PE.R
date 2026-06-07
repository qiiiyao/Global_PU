### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("~/my_pc/Multi_taxa/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})

setwd("~/my_pc/Multi_taxa")
source('code/calculating_PE_func.R')

#### Fish PE: calculation & mapping ####
load("/home/yaoqi/my_pc/Multi_taxa/data/Fishes/data/my_phy.rdata")
load("/home/yaoqi/my_pc/Multi_taxa/data/Fishes/data/my_data_used_final.rdata")
load("data/Fishes/data/Basin042017_3119.rdata")

df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
colnames(df_trans)[which(colnames(df_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% dplyr::filter(X3.Native.Exotic.Status == 
                                                       'exotic')
data.used_final_natives = data.used_final %>% dplyr::filter(X3.Native.Exotic.Status == 
                                                       'native')

##### natives #####
data.used_final_natives$presence = 1

comm_fish_native = data.used_final_natives %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_trans[,c('X1.Basin.Name', 'Surf_area')], 
                   #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
                   by = 'X1.Basin.Name') %>% 
  relocate(Surf_area, .after = X1.Basin.Name)

pe_fish_native = calcu_PE_RPE_parallel(
  Tree = phylo %>% 
    keep.tip(unique(data.used_final_natives$valid_names)),
  Comm = comm_fish_native[,unique(data.used_final_natives$valid_names)],
  Area = comm_fish_native$Surf_area)

pe_fish_native_1 = cbind(Surf_area = comm_fish_native$Surf_area,
                    pe_fish_native)
pe_fish_native_sf = df_trans %>% left_join(pe_fish_native_1, by = 'Surf_area')
colnames(pe_fish_native_sf)
save(pe_fish_native_sf, file = 'results/primary_results/pe_fish_native_sf.rdata')

##### extant #####
data.used_final$presence = 1

comm_fish_extant = data.used_final %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  left_join(df_trans[,c('X1.Basin.Name', 'Surf_area')], 
            #we calculate PE by multiplying fishes PD by their living basins' drainage surface area
            by = 'X1.Basin.Name') %>% 
  relocate(Surf_area, .after = X1.Basin.Name)

pe_fish_extant = calcu_PE_RPE_parallel(
  Tree = phylo %>% 
    keep.tip(unique(data.used_final$valid_names)),
  Comm = comm_fish_extant[,unique(data.used_final$valid_names)],
  Area = comm_fish_extant$Surf_area)

pe_fish_extant_1 = cbind(Surf_area = comm_fish_extant$Surf_area,
                         pe_fish_extant)
pe_fish_extant_sf = df_trans %>% left_join(pe_fish_extant_1, by = 'Surf_area')
colnames(pe_fish_extant_sf)
save(pe_fish_extant_sf, file = 'results/primary_results/pe_fish_all_sf.rdata')


