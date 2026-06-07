### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
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

setwd("D:/R projects/Global_ED")



#### Transfer this data into type of TDWG ####
TDWG.glonaf.new = st_read("data/Plants/TDWG_glonaf_new_eck4.shp")
TDWG.glonaf.new_dat = TDWG.glonaf.new %>% st_drop_geometry()
TDWG.glonaf.new_dat = TDWG.glonaf.new_dat[,c('Region_id', 'RegionID')] %>% 
  dplyr::group_by(RegionID) %>% 
  dplyr::distinct(Region_id)
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")

df.native.650_TDWG = df.native.650 %>% left_join(TDWG.glonaf.new_dat,
                                                 by = 'Region_id',
                                                 relationship = 'many-to-many') %>% 
  dplyr::group_by(RegionID) %>% 
  dplyr::distinct(family, genus, status, duplicated, species)

df.natu.650_TDWG = df.natu.650 %>% left_join(TDWG.glonaf.new_dat,
                                                 by = 'Region_id',
                                                 relationship = 'many-to-many') %>% 
  dplyr::group_by(RegionID) %>% 
  dplyr::distinct(family, genus, status, duplicated, species)

colnames(df.natu.650_TDWG) == colnames(df.native.650_TDWG)

df.native.natu.species.650.nonative_TDWG = list(df.native.650_TDWG,
                                                df.natu.650_TDWG)

sort(unique(df.native.650_TDWG$RegionID))
sort(unique(df.natu.650_TDWG$RegionID))

save(df.native.natu.species.650.nonative_TDWG,
file = 'data/Plants/data/df.native.natu.species.650.nonative_TDWG.Rdata')

