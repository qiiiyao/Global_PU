rm(list = ls())

library(tidyverse)
library(conflicted)
library(FishPhyloMaker)
library(ape)

#load data
setwd("D:/R projects/Multi_taxa")
occ_drainage = read.csv("data/Fishes/data/all_data_global/Monde_Download_0621.csv",
                        sep = ';',
                        header = T)
occ_drainage_old = read.csv("data/Fishes/data/Occurrence_Table.csv",
                        sep = ';',
                        header = T)
drainage_basins = read.csv("data/Fishes/data/Drainage_Basins_Table.csv",
                           sep = ";")
length(unique(occ_drainage$X6.Fishbase.Valid.Species.Name))
  
#1.The generation of valid species names
occ_drainage$X2.Species.Name.in.Source  =  gsub("[.]", "_", occ_drainage$X2.Species.Name.in.Source)
occ_drainage$X6.Fishbase.Valid.Species.Name  =  gsub("[.]", "_", occ_drainage$X6.Fishbase.Valid.Species.Name)

species_list = unique(occ_drainage$X6.Fishbase.Valid.Species.Name)
all_taxa_names = FishTaxaMaker(data = species_list, allow.manual.insert = TRUE)

#2.Further cleaning of the valid species names file Taxon_data_FishPhyloMaker
all_taxa_names$Taxon_data_FishPhyloMaker$o  =  gsub("/.*", replacement = " ", all_taxa_names$Taxon_data_FishPhyloMaker$o) # removing slash in Perciformes
all_taxa_names$Taxon_data_FishPhyloMaker = all_taxa_names$Taxon_data_FishPhyloMaker[-which(is.na(all_taxa_names$Taxon_data_FishPhyloMaker$s)),
                                                                                    ]#remove the row including NA


#3.Details of the process for determining valid species names
# (1) pull and export the duplicated names of the valid names in the All_info_fishbase 
all_taxa_names$All_info_fishbase[all_taxa_names$All_info_fishbase$valid_names%>%duplicated(),]%>%
  drop_na(valid_names)%>%
  pull(valid_names)%>%
  write.csv(file="data/Fishes/data/duplicated names.csv")

# (2) pull and export the 5 unidentified names of the valid names, but not NA, in the Taxon_data_FishPhyloMaker 
all_taxa_names$All_info_fishbase[!(gsub(" ","_",
                                        all_taxa_names$All_info_fishbase$valid_names)%in%all_taxa_names$Taxon_data_FishPhyloMaker$s),]%>%
  drop_na(valid_names)%>%
  pull(valid_names)%>%
  write.csv(file="data/Fishes/data/not identified names.csv")

# (3) extract the previous invalid species name and replacing the converted species names via FishTaxaMaker(allow.manual.insert=T)
all_taxa_names$Taxon_data_FishPhyloMaker[all_taxa_names$Taxon_data_FishPhyloMaker$s %in%
                                           gsub(" ","_",all_taxa_names$Species_not_in_Fishbase),]%>%
  write.csv(file="data/Fishes/data/na to valid names.csv")


#4.Generate valid species names for the 5 species not identified in 3.(2), and integrate them with the Taxon_data_FishPhyloMaker.
not_identified_names = all_taxa_names$All_info_fishbase[!(gsub(" ",
                                                               "_",
                                                               all_taxa_names$All_info_fishbase$valid_names)%in%
                                                            all_taxa_names$Taxon_data_FishPhyloMaker$s),]%>%
  drop_na(valid_names)%>%
  pull(valid_names)
not_identified_names = gsub(" ", "_", not_identified_names)
not_identified_taxa = FishTaxaMaker(not_identified_names)
not_identified_df = not_identified_taxa$Taxon_data_FishPhyloMaker# the dataframe for 5 unidentified names of the valid names
generated_df = all_taxa_names$Taxon_data_FishPhyloMaker# the automatic dataframe for all the valid names

df_for_phylomaker = rbind(not_identified_df,generated_df)# total of 14868 species

save(df_for_phylomaker, file = 'data/Fishes/data/df_for_phylomaker.rdata')
save(all_taxa_names, file = 'data/Fishes/data/all_taxa_names.rdata')

#5.Construction of a global fish phylogenetic tree
load('data/Fishes/data/df_for_phylomaker.rdata')
source('code/FishPhyloMaker_0.2.0/FishPhyloMaker/R/internal_filter_rank.R')
source('code/FishPhyloMaker_0.2.0/FishPhyloMaker/R/internal_treedata_modif.R')
source('code/FishPhyloMaker_0.2.0/FishPhyloMaker/R/internal_user_opt_printCatFamily.R')
source('code/FishPhyloMaker_0.2.0/FishPhyloMaker/R/internal_user_opt_printCat.R')
source('code/FishPhyloMaker_0.2.0/FishPhyloMaker/R/fishPhyloMaker_yq.R')

#class_fish = whichFishAdd(data = df_for_phylomaker) 
phylo_all_spp = FishPhyloMaker_yq(data = df_for_phylomaker, 
                                 insert.base.node = TRUE, 
                                 return.insertions = TRUE)
# using FishPhyloMaker_yq modifying some bugs in the original function

### the ultimately constructed phy_tree includes 14682 species, as 186 species 
# have not been inserted in the tree.
phylo_all_spp$Insertions_data %>% dplyr::filter(insertions=="Not_inserted")%>%
  pull(s)%>%
  write.csv(file="data/Fishes/data/not inserted to fish tree.csv")

# extract the phy_tree
phylo = phylo_all_spp$Phylogeny
save(phylo, file = 'data/Fishes/data/my_phy.rdata')


#6.Calculation of pairwise phylogenetic distances
library(ape)
#distance_all = cophenetic.phylo(phylo)
#distance_all = distance_all/max(distance_all)
#tibble(distance_all)

#7.Matching occurrence data with phylogenetic data
load("D:/R projects/Multi_taxa/data/Fishes/data/my_phy.rdata")
load("D:/R projects/Multi_taxa/data/Fishes/data/all_taxa_names.rdata")

# update the species valid names in the occurence data to set them consistent with phy_tree
phylo.valid.names.df = all_taxa_names$All_info_fishbase %>% select(1:2) %>% tibble()
phylo.valid.names.df$valid_names = gsub(" ", "_", phylo.valid.names.df$valid_names)
phylo.valid.names.df_na = phylo.valid.names.df %>% 
  dplyr::filter(is.na(valid_names)) %>% 
  mutate(valid_names = user_spp) # replace the name of 5 unidentified valid
# species names to the modified valid species names
phylo.valid.names.df_residue = phylo.valid.names.df %>%
  dplyr::filter(!is.na(valid_names))

phylo.valid.names.df.final = rbind(phylo.valid.names.df_residue, phylo.valid.names.df_na)

data.occurrence = left_join(occ_drainage, phylo.valid.names.df.final,
                           by=c("X6.Fishbase.Valid.Species.Name"="user_spp"))%>%
  tibble() # combine the valid names of species to the occurrence data

data.all = left_join(data.occurrence, drainage_basins) %>% 
  select(X2.Country,everything()) # combine the occurrence data to the basin information

### Discard the countries that no exotics present
country.used = data.all %>% 
  dplyr::filter(X3.Native.Exotic.Status=="exotic") %>% 
  pull(X2.Country) %>% unique() # 123 countries of 143 countries has exotics

data.used.0 = data.all %>% dplyr::filter(X2.Country%in%country.used) %>%
  group_by(X1.Basin.Name) %>%
  distinct(valid_names,.keep_all = T) %>%
  #select(c(1,2,4,9)) %>%
  arrange(X2.Country,desc(X3.Native.Exotic.Status),X1.Basin.Name) %>%
  ungroup() # after discard, there are *** species totally

# exclude the 184 species that not in the phy_tree, which include 6 exotics that 
# distributed in the 4 countries.
species.exclude = (data.used.0$valid_names %>% 
                     unique)[!((data.used.0$valid_names%>%unique) %in% phylo$tip.label)]

data.exclude = data.used.0 %>% 
  dplyr::filter(valid_names %in% species.exclude)

data.exclude %>% 
  summarise(n.country=n_distinct(X2.Country),
            n.basin=n_distinct(X1.Basin.Name),
            n.status=n_distinct(X3.Native.Exotic.Status),
            n.species=n_distinct(valid_names)) # list of excluded species 

data.exclude %>% dplyr::filter(X3.Native.Exotic.Status == "exotic")# list of excluded 4 exotics species 

data.used.1 = data.used.0 %>% dplyr::filter(!(valid_names %in% species.exclude))
data.used.1$valid_names %>% unique() %>% length()# length of retained 14618 species
data.used.1 %>% print(n=100)

# correct the exotic status of Cyprinus_carpio in China to native
data.used = data.used.1 %>% 
  mutate(X3.Native.Exotic.Status = if_else(X2.Country=="China"&valid_names=="Cyprinus_carpio",
                                           "native",X3.Native.Exotic.Status))


#--------------------Analyze the introduced alien species

#8.# find the exotics originally defined by the different basin, 
# which means even in the same country, but in the different basin
data.distinct = data.used %>% group_by(X2.Country,X3.Native.Exotic.Status) %>%
  distinct(valid_names,.keep_all = T)
data.distinct = data.distinct %>% ungroup()
data.distinct %>% print(n=100)

data.dup = data.distinct %>% group_by(X2.Country) %>% 
  mutate(dup=duplicated(valid_names)) %>% 
  dplyr::filter(dup) # translocated species are both natives and exotics in a country
data.dup = data.dup %>% unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)
data.dup%>%print(n=100)

# add a new column: "X3.Native.Exotic.Status_country" to define the exotics 
# by the country, not the basin
data.used$X3.Native.Exotic.Status_country = data.used$X3.Native.Exotic.Status

data.used_final = data.used %>%
  unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)

data.used_final[data.used_final$new.col %in% data.dup$new.col,]$X3.Native.Exotic.Status_country = "native"

data.used_final%>%group_by(X2.Country,X3.Native.Exotic.Status_country)%>%
  distinct(valid_names,.keep_all = T)%>%ungroup()%>%
  group_by(X2.Country) %>% mutate(dup=duplicated(valid_names))%>%dplyr::filter(dup)
# check that there is no duplicated species that are both natives and exotics in a country

country.used.final = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status_country=="exotic") %>% 
  pull(X2.Country)%>%unique()# only 120 countries have exotics according to X3.Native.Exotic.Status_country

data.used_final# final data frame
data.used_final %>% print(n=200)

# summarise 
data.used_final %>% 
  summarise(n.country = n_distinct(X2.Country),
            n.basin = n_distinct(X1.Basin.Name),
            n.speces = n_distinct(valid_names))

data.used_final %>% 
  group_by(X3.Native.Exotic.Status) %>% 
  summarise(n.country=n_distinct(X2.Country),
            n.basin=n_distinct(X1.Basin.Name),
            n.speces=n_distinct(valid_names))

save(data.used_final, file = 'data/Fishes/data/my_data_used_final.rdata')
