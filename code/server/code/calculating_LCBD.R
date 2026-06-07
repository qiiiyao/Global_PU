### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("/home/yaoqi/my_pc/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'scico', 'ggplot2', 'gridExtra')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("/home/yaoqi/my_pc/Global_ED")
source('code/functions/calculating_LCBD_func.R')

# load the background data for plotting the world map plot
df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
load("code/FYI/Codes_and_Data_Fan_et_al_2023/data.for.shp.plots.4.Rdata")


# define color gradients
scico::scico_palette_show()

colors1 = rev(scico::scico(n=8, palette = "lajolla"))
colors2 = scico::scico(n=10, palette = "vik")[1:8] 
colors3 = scico::scico(n=8, begin = 0, end = 0.4, palette = "bam")
colors4 = scico::scico(n=8, palette = "bam")
colors5 = scico::scico(n=8, begin = 0, end = 0.8, palette = "bam")

#plot(1:8, col = colors2, pch = 19, cex = 5)

# define a function to make legend of color gradients
legend.func = function(mycolors, mylabels) {
  group = rep("cc", 8)
  condition = letters[1:8]
  value = rep(1, 8)
  df.legend = data.frame(group, condition, value)
  mycolors.corrected = rev(mycolors)
  ggplot(df.legend, aes(fill = condition, y = value, x = group)) +
    geom_bar(position = "stack", stat = "identity", color = "white") +
    scale_fill_manual(values = mycolors.corrected) +
    theme_classic() +
    theme(
      legend.position = "none", aspect.ratio = 0.03,
      axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(size = 7, color = "black")
    ) +
    scale_y_continuous(breaks = 0:8, labels = mylabels) +
    coord_flip() +
    xlab("")
}


#### Plant ED: calculation & mapping ####
load("/home/yaoqi/my_pc/Global_ED/data/Plants/data/shp.651.Rdata")
load('code/FYI/Codes_and_Data_Cai_et_al_2024/data/checklist_restricted.RDATA')
#load("/home/yaoqi/my_pc/Global_ED/data/Plants/data/phylo.fake.species.653.Rdata")
phylo_big = phytools::read.newick("code/FYI/Codes_and_Data_Cai_et_al_2024/data/v0.1/ALLOTB.tre")
load("/home/yaoqi/my_pc/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)
checklist_restricted$work_species = stringr::str_replace_all(checklist_restricted$work_species,
 " ",
  "_")  # changes spaces to underline

colnames(df.native.650) == colnames(df.natu.650)

##1.3 Remove apomictic species----
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

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
colnames(df.natu.650)
df.natu.650$presence = 1

##### natives #####
colnames(df.native.650)
df.native.650$presence = 1
phylo_plant_native = drop.tip(phylo_big, setdiff(phylo_big$tip.label, unique(df.native.650$species)))
length(unique(df.native.650$species))
comm_plant_native = df.native.650 %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

LCBD_plant_native = calcu_LCBD_parallel(Tree = phylo_plant_native,
                                        Comm = comm_plant_native,
                                        Region_posi = which(colnames(comm_plant_native) == 'Region_id'))

save(LCBD_plant_native,
     file = 'results/primary_results/LCBD_plant_native.rdata')

#### native plants mapping 



##### extant #####
df.extant.650 = rbind(df.native.650, df.natu.650)
colnames(df.extant.650)
df.extant.650$presence = 1
phylo_plant_extant = drop.tip(phylo_big, setdiff(phylo_big$tip.label, df.extant.650$species))

comm_plant_extant = df.extant.650 %>% 
  complete(Region_id, species, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence) %>% 
  left_join(shp.glonaf.trans[,c('Region_id', 'area')],
            by = 'Region_id') %>% 
  relocate(area, .after = Region_id) %>% 
  arrange('Region_id') 

LCBD_plant_extant = calcu_LCBD_parallel(Tree = phylo_plant_extant,
                                        Comm = comm_plant_extant,
                                        Region_posi = which(colnames(comm_plant_extant) == 'Region_id'))
save(LCBD_plant_extant,
     file = 'results/primary_results/LCBD_plant_extant.rdata')

#### extant plants mapping 


#### Bird ED: calculation & mapping ####
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

df.trans = st_transform(df, crs = "+proj=eck4") 
df.trans$area = st_area(df.trans)
sum(as.numeric(df.trans$area) == 0)

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
exotic_distri_data$present = 1

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
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID') 

LCBD_bird_native = calcu_LCBD_parallel(Tree = phy_data,
                                       Comm = comm_bird_native,
                                       Region_posi = which(colnames(comm_bird_native) == 'RegionID'))

save(LCBD_bird_native, file = 'results/primary_results/LCBD_bird_native.rdata')

#### native birds mapping 



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
  relocate(area, .after = RegionID) %>% 
  arrange('RegionID') 

LCBD_bird_extant = calcu_LCBD_parallel(Tree = phy_data,
                                      Comm = comm_bird_extant,
                                      Region_posi = which(colnames(comm_bird_extant) == 'RegionID'))

save(LCBD_bird_extant, file = 'results/primary_results/LCBD_bird_extant.rdata')

#### extant birds mapping 


#### Fish ED: calculation & mapping ####
load("/home/yaoqi/my_pc/Global_ED/data/Fishes/data/my_phy.rdata")
load("/home/yaoqi/my_pc/Global_ED/data/Fishes/data/my_data_used_final.rdata")
load("/home/yaoqi/my_pc/Global_ED/data/Fishes/data/Basin042017_3119.rdata")
df_trans = st_transform(df, crs = "+proj=eck4") 
df_trans$area = st_area(df_trans)
#save(df, file = "data/Fishes/data/Basin042017_3119.rdata")
colnames(df_trans)[which(colnames(df_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
data.used_final_exotics$presence = 1

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
  relocate(Surf_area, .after = X1.Basin.Name)%>% 
  arrange('X1.Basin.Name') 

LCBD_fish_native = calcu_LCBD_parallel(
  Tree = phylo,
  Comm = comm_fish_native,
  Region_posi = which(colnames(comm_fish_native) == 'X1.Basin.Name'))

save(LCBD_fish_native, file = 'results/primary_results/LCBD_fish_native.rdata')


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
  relocate(Surf_area, .after = X1.Basin.Name) %>% 
  arrange('X1.Basin.Name') 

LCBD_fish_extant = calcu_LCBD_parallel(
  Tree = phylo,
  Comm = comm_fish_extant,
  Region_posi = which(colnames(comm_fish_extant) == 'X1.Basin.Name'))

save(LCBD_fish_extant, file = 'results/primary_results/LCBD_fish_extant.rdata')