### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
library(reshape2)
#library(plyr)
library(ape)
library(MASS)
library(vegan)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
#library(betareg)
library(lme4)
#library(glmmTMB)
library(INLA)
library(sf)

setwd("~/my_pc/Global_ED")

# define color gradients
colors_4d = c('#91989f', '#efbb24', '#1b813e', '#398fb7')
colors_10d = c('#91989f',
               colorRampPalette(c('#FFF9B0', '#FFE066', '#FFC300', '#FF9E00', '#FF7A00', '#E85C00'))(6),
               colorRampPalette(c('#9BD4E4', '#398FB7', '#004A78'))(3))

df = st_read("data/TDWG4/TDWG4_newTibet_eck4.shp")


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
phy_bird = read.tree("data/Birds/data/Phylogenetic_Birds.tre")


df.trans = st_transform(df, crs = "+proj=eck4") 
df.trans$area = st_area(df.trans)
sum(as.numeric(df.trans$area) == 0)


load('results/primary_results/distances_beta/ED_bird_native.rdata')

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species

bird_sps_ED = phyloregion::evol_distinct(tree = phy_bird, type = "fair.proportion")
bird_sps_ED = data.frame(species = names(bird_sps_ED), 
                         ED = bird_sps_ED)

exotic_distri_data$present = 1
comm_bird_exotic = exotic_distri_data %>% 
  complete(ScientificName, RegionID, fill = list(present = 0)) %>%
  dplyr::select(c('RegionID', 'ScientificName', 'present')) %>% 
  pivot_wider(names_from = ScientificName,
              values_from = present,
              values_fn = mean) %>% 
  complete(RegionID = unique(df$RegionID)) %>% ## assume that 
  ##  the absence regions have no naturalized aliens
  mutate(
    across(everything(), \(x) {
      if (is.factor(x)) {
        fct_expand(x, "0") %>%  tidyr::replace_na("0")
      } else if (is.character(x)) {
        tidyr::replace_na(x, "0")
      } else if (inherits(x, "units")) {
        unit = units(x)
        tidyr::replace_na(x, set_units(0, unit, mode = "standard"))
      } else {
        tidyr::replace_na(x, 0)
      }
    })
  ) %>% 
  dplyr::arrange('RegionID') %>% 
  pivot_longer(cols = Acanthis_flammea:Zosterops_natalis, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(bird_sps_ED, by = 'species') %>% 
  rename(c('exotic_ED' = 'ED'))


mean_ED_natives_bird = cbind(RegionID = df.trans$RegionID,
                             ED_bird_native$df)

sp_bird_exo_nat = comm_bird_exotic %>% 
  left_join(mean_ED_natives_bird, by = 'RegionID') %>% 
  rename(c('mean_native_ED' = 'mean_ED')) %>% 
  rename(c('mean_native_EDR' = 'mean_EDR'))


sp_bird_exo_nat$delta_exo_nat = log(sp_bird_exo_nat$exotic_ED/
                                      sp_bird_exo_nat$mean_native_ED)

sp_bird_exo_nat = sp_bird_exo_nat %>% filter(!is.na(delta_exo_nat))

### get predict data for mod_bird_sp_inv_ed_inla
lincombs.bird.data.estab.delta_e_n_single = data.frame(delta_exo_nat=seq(min(sp_bird_exo_nat$delta_exo_nat),
                                                                         max(sp_bird_exo_nat$delta_exo_nat),
                                                                         length=100))

lincombs.bird.matrix.estab.delta_e_n_single=model.matrix(~delta_exo_nat,
                                                         data=lincombs.bird.data.estab.delta_e_n_single)
lincombs.bird.matrix.estab.delta_e_n_single=as.data.frame(lincombs.bird.matrix.estab.delta_e_n_single)
lincombs.bird.estab.delta_e_n_single=inla.make.lincombs(lincombs.bird.matrix.estab.delta_e_n_single)

inla.model_lincombs.bird.estab.delta_e_n_single = inla(presence ~ delta_exo_nat
                                                       + f(RegionID, model="iid")  
                                                       + f(species, delta_exo_nat, model="iid")
                                                       ,
                                                       data = sp_bird_exo_nat,
                                                       family = "zeroinflatedbinomial1",
                                                       control.compute = list(dic = TRUE,
                                                                              waic = TRUE, 
                                                                              cpo = TRUE,
                                                                              config = TRUE),
                                                       control.predictor = list(compute = TRUE),
                                                       quantiles = c(0.025, 0.5, 0.975),
                                                       lincomb = lincombs.bird.estab.delta_e_n_single)

lincombs.bird.posterior.estab.delta_e_n_single = inla.model_lincombs.bird.estab.delta_e_n_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.bird.estab.delta_e_n_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.bird.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.bird.data.estab.delta_e_n_single$predicted.value=unlist(lapply(lincombs.bird.posterior.estab.delta_e_n_single,
                                                                        function(x)inla.emarginal(fun=plogis,x)))
lincombs.bird.data.estab.delta_e_n_single$lower=unlist(lapply(lincombs.bird.posterior.estab.delta_e_n_single,
                                                              function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.bird.data.estab.delta_e_n_single$upper=unlist(lapply(lincombs.bird.posterior.estab.delta_e_n_single,
                                                              function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

lincombs.bird.data.estab.delta_e_n_single_ran2 = list(summary = inla.model_lincombs.bird.estab.delta_e_n_single$summary.fixed,
                                                      prediction = lincombs.bird.data.estab.delta_e_n_single)

#lincombs.bird.data.estab.delta_e_n_single
save(lincombs.bird.data.estab.delta_e_n_single_ran2, 
     file = 'results/primary_results/lincombs.bird.data.estab.delta_e_n_single_ran2.rdata')

rm(sp_bird_exo_nat)

#### Fish ED: calculation & mapping ####
load("data/Fishes/data/my_phy.rdata")
phy_fish = phylo

load("data/Fishes/data/my_data_used_final.rdata")
df_trans = st_read("data/Fishes/data/Basin042017_3119_eck4/Basin042017_3119_eck4.shp")

#save(df, file = "data/Fishes/data/Basin042017_3119.rdata")
colnames(df_trans)[which(colnames(df_trans) == 'BasinName')] = 'X1.Basin.Name'
colnames(data.used_final)
str(data.used_final)
unique(data.used_final$X3.Native.Exotic.Status) # define exotics according to different basin, not countries
data.used_final_exotics = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'exotic')
data.used_final_natives = data.used_final %>% 
  dplyr::filter(X3.Native.Exotic.Status == 'native')

load('results/primary_results/distances_beta/ED_fish_native.rdata')

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
fish_sps_ED = phyloregion::evol_distinct(tree = phy_fish, type = "fair.proportion")
fish_sps_ED = data.frame(species = names(fish_sps_ED), 
                         ED = fish_sps_ED)

data.used_final_exotics$presence = 1

comm_fish_exotic = data.used_final_exotics %>% 
  complete(valid_names, X1.Basin.Name, fill = list(presence = 0)) %>%
  dplyr::select(c('X1.Basin.Name', 'valid_names', 'presence')) %>% 
  pivot_wider(names_from = valid_names,
              values_from = presence,
              values_fn = mean) %>% 
  complete(X1.Basin.Name = unique(df$BasinName)) %>% ## assume that 
  ##  the absence regions have no naturalized aliens
  mutate(
    across(everything(), \(x) {
      if (is.factor(x)) {
        fct_expand(x, "0") %>%  tidyr::replace_na("0")
      } else if (is.character(x)) {
        tidyr::replace_na(x, "0")
      } else if (inherits(x, "units")) {
        unit = units(x)
        tidyr::replace_na(x, set_units(0, unit, mode = "standard"))
      } else {
        tidyr::replace_na(x, 0)
      }
    })
  )  %>% 
  dplyr::arrange('X1.Basin.Name') %>% 
  pivot_longer(cols = Abbottina_rivularis:Zacco_platypus, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(fish_sps_ED, by = 'species') %>% 
  rename(c('exotic_ED' = 'ED'))

mean_ED_natives_fish = cbind(Basin.Name = sort(unique(data.used_final_natives$X1.Basin.Name)),
                             ED_fish_native$df)

sp_fish_exo_nat = comm_fish_exotic %>% 
  left_join(mean_ED_natives_fish, by = join_by('X1.Basin.Name' == 'Basin.Name')) %>% 
  rename(c('mean_native_ED' = 'mean_ED')) %>% 
  rename(c('mean_native_EDR' = 'mean_EDR'))

sp_fish_exo_nat$delta_exo_nat = log(sp_fish_exo_nat$exotic_ED/
                                      sp_fish_exo_nat$mean_native_ED)

sp_fish_exo_nat = sp_fish_exo_nat %>% filter(!is.na(delta_exo_nat))


### get predict data for mod_fish_sp_inv_ed_inla
lincombs.fish.data.estab.delta_e_n_single = data.frame(delta_exo_nat=seq(min(sp_fish_exo_nat$delta_exo_nat),
                                                                         max(sp_fish_exo_nat$delta_exo_nat),
                                                                         length=100))

lincombs.fish.matrix.estab.delta_e_n_single=model.matrix(~delta_exo_nat,
                                                         data=lincombs.fish.data.estab.delta_e_n_single)
lincombs.fish.matrix.estab.delta_e_n_single=as.data.frame(lincombs.fish.matrix.estab.delta_e_n_single)
lincombs.fish.estab.delta_e_n_single=inla.make.lincombs(lincombs.fish.matrix.estab.delta_e_n_single)

inla.model_lincombs.fish.estab.delta_e_n_single = inla(presence ~ delta_exo_nat
                                                       + f(X1.Basin.Name, model="iid")  
                                                       + f(species, delta_exo_nat, model="iid")
                                                       ,
                                                       data = sp_fish_exo_nat,
                                                       family = "zeroinflatedbinomial1",
                                                       control.compute = list(dic = TRUE,
                                                                              waic = TRUE, 
                                                                              cpo = TRUE,
                                                                              config = TRUE),
                                                       control.predictor = list(compute = TRUE),
                                                       quantiles = c(0.025, 0.5, 0.975),
                                                       lincomb = lincombs.fish.estab.delta_e_n_single)

lincombs.fish.posterior.estab.delta_e_n_single = inla.model_lincombs.fish.estab.delta_e_n_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.fish.estab.delta_e_n_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.fish.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.fish.data.estab.delta_e_n_single$predicted.value=unlist(lapply(lincombs.fish.posterior.estab.delta_e_n_single,
                                                                        function(x)inla.emarginal(fun=plogis,x)))
lincombs.fish.data.estab.delta_e_n_single$lower=unlist(lapply(lincombs.fish.posterior.estab.delta_e_n_single,
                                                              function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.fish.data.estab.delta_e_n_single$upper=unlist(lapply(lincombs.fish.posterior.estab.delta_e_n_single,
                                                              function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

lincombs.fish.data.estab.delta_e_n_single_ran2 = list(summary = inla.model_lincombs.fish.estab.delta_e_n_single$summary.fixed,
                                                      prediction = lincombs.fish.data.estab.delta_e_n_single)

#lincombs.fish.data.estab.delta_e_n_single
save(lincombs.fish.data.estab.delta_e_n_single_ran2, 
     file = 'results/primary_results/lincombs.fish.data.estab.delta_e_n_single_ran2.rdata')


