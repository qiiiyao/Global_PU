### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("/home/yaoqi/my_pc/Multi_taxa/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'scico')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})

setwd("/home/yaoqi/my_pc/Multi_taxa")
source('code/calculating_PE_func.R')
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)

#### Bird PE: calculation & mapping ####
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

native_distri_data = all_distri_data_c %>% dplyr::filter(SpStatus == 'Native')
exotic_distri_data = all_distri_data_c %>% dplyr::filter(SpStatus == 'alien')
phy_data = ape::read.tree("data/Birds/data/Phylogenetic_Birds.tre")

df.trans = st_transform(df, crs = "+proj=eck4") 
df.trans$area = st_area(df.trans)
sum(as.numeric(df.trans$area) == 0)

##### natives #####
colnames(native_distri_data)
native_distri_data$present = 1
comm_bird_native = native_distri_data %>% 
  complete(RegionID, ScientificName, fill = list(present = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'present')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = present,
              values_fn = mean) %>% 
  left_join(df.trans[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID)

pe_bird_native = calcu_PE_RPE_parallel(Tree = phy_data %>% keep.tip(unique(native_distri_data$ScientificName)),
                                        Comm = as.matrix(comm_bird_native[,unique(native_distri_data$ScientificName)]),
                                        Area = comm_bird_native$area)
pe_bird_native_1 = cbind(RegionID = comm_bird_native$RegionID,
                          pe_bird_native)
pe_bird_native_sf = df.trans %>% left_join(pe_bird_native_1, by = 'RegionID')
colnames(pe_bird_native_sf)
save(pe_bird_native_sf, file = 'results/primary_results/pe_bird_native_sf.rdata')

##### extant #####
colnames(all_distri_data_c)
all_distri_data_c$present = 1
comm_bird_extant = all_distri_data_c %>% 
  complete(RegionID, ScientificName, fill = list(present = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'present')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = present,
              values_fn = mean) %>% 
  left_join(df.trans[,c('RegionID', 'area')],
            by = 'RegionID') %>% 
  relocate(area, .after = RegionID)

pe_bird_extant = calcu_PE_RPE_parallel(Tree = phy_data %>% 
  keep.tip(unique(all_distri_data_c$ScientificName)),
                                       Comm = as.matrix(comm_bird_extant[,unique(all_distri_data_c$ScientificName)]),
                                       Area = comm_bird_extant$area)
pe_bird_extant_1 = cbind(RegionID = comm_bird_extant$RegionID,
                         pe_bird_extant)
pe_bird_extant_sf = df.trans %>% left_join(pe_bird_extant_1, by = 'RegionID')
colnames(pe_bird_extant_sf)
save(pe_bird_extant_sf, file = 'results/primary_results/pe_bird_all_sf.rdata')

