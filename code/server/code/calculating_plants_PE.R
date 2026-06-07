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
df = st_read("data/TDWG4/TDWG4_newTibet.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)

#### Plant PE: calculation & mapping ####
load("/home/yaoqi/my_pc/Multi_taxa/data/Plants/data/shp.651.Rdata")
load("/home/yaoqi/my_pc/Multi_taxa/data/Plants/data/phylo.fake.species.653.Rdata")
load("/home/yaoqi/my_pc/Multi_taxa/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)

colnames(df.native.650) == colnames(df.natu.650)

##### extant #####
df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1
comm_plant_extant = df.extant.650 %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id)

pe_plant_extant = calcu_PE_RPE_parallel(Tree = phylo %>% keep.tip(unique(df.extant.650$species)),
                                        Comm = as.matrix(comm_plant_extant[,unique(df.extant.650$species)]),
                                        Area = comm_plant_extant$area)
pe_plant_extant_1 = cbind(Region_id = comm_plant_extant$Region_id,
                          pe_plant_extant)
pe_plant_extant_sf = shp.glonaf.trans %>% left_join(pe_plant_extant_1, by = 'Region_id')
colnames(pe_plant_extant_sf)

save(pe_plant_extant_sf, file = 'results/primary_results/pe_plant_all_sf.rdata')

