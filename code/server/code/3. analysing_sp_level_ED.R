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

#### Mammal ED: calculation & mapping ####
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")

load('results/primary_results/distances_beta/ED_mammal_native.rdata')

phy_mammal = spec_phy.3

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
sp_overlap_dat_1$Binomial = gsub(' ', '_', sp_overlap_dat_1$Binomial)
sp_overlap_dat_1$presence = 1

mammal_sps_ED = phyloregion::evol_distinct(tree = phy_mammal, type = "fair.proportion")
mammal_sps_ED = data.frame(species = names(mammal_sps_ED), 
                           ED = mammal_sps_ED)

comm_mammal_exotic = sp_overlap_dat_1 %>% 
  sf::st_drop_geometry() %>% 
  complete(Binomial, RegionID, fill = list(presence = 0)) %>%
  dplyr::select(c('RegionID', 'Binomial', 'presence')) %>% 
  pivot_wider(names_from = Binomial,
              values_from = presence,
              values_fn = mean) %>% 
  complete(RegionID = unique(df$RegionID)) %>% ## assume that 
  ## the absence regions have no naturalized aliens
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
  pivot_longer(cols = Alces_alces:Wallabia_bicolor, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(mammal_sps_ED, by = 'species') %>% 
  plyr::rename(c('ED' = 'exotic_ED')) %>% 
  dplyr::arrange('RegionID')




mean_ED_natives_mammal = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                                                    ED_mammal_native$df)

sp_mammal_exo_nat = comm_mammal_exotic %>% 
  left_join(mean_ED_natives_mammal, by = 'RegionID') %>% 
  rename(c('mean_native_ED' = 'mean_ED')) %>% 
  rename(c('mean_native_EDR' = 'mean_EDR'))
  

sp_mammal_exo_nat$delta_exo_nat = log(sp_mammal_exo_nat$exotic_ED/
                                        sp_mammal_exo_nat$mean_native_ED)
nrow(sp_mammal_exo_nat %>% 
  filter(presence == 1)) / nrow(sp_mammal_exo_nat)

sp_mammal_exo_nat = sp_mammal_exo_nat %>% filter(!is.na(delta_exo_nat))

plot(sp_mammal_exo_nat$delta_exo_nat, sp_mammal_exo_nat$presence)

nrow(sp_mammal_exo_nat %>% filter(presence == 0))/nrow(sp_mammal_exo_nat)

mod_exotic_ED = lme4::glmer(presence ~ exotic_ED 
                   + (1|RegionID) 
                   #+ (1|species)
                   , 
                   family = binomial, data = sp_mammal_exo_nat)
summary(mod_exotic_ED)



sp_mammal_exo_nat$Region_species = paste(sp_mammal_exo_nat$RegionID,
                                         sp_mammal_exo_nat$species,
                                         sep = '_')

lme_delta_ED = glmer(presence ~ delta_exo_nat + (1|RegionID) + (delta_exo_nat | species), 
                   data = sp_mammal_exo_nat,
                   family = binomial(link = "cloglog"),
                   nAGQ=0, 
                   control=glmerControl(optimizer = "nloptwrap"))
summary(lme_delta_ED)

inla_delta_ED = inla(presence ~ delta_exo_nat + f(RegionID, model="iid") + 
                  f(species, delta_exo_nat, model="iid"),
     family = "zeroinflatedbinomial1",
     data = sp_mammal_exo_nat 
     #Ntrials=Ntrials
     )
summary(inla_delta_ED)


### get predict data for mod_mammal_sp_inv_ed_inla
lincombs.mammal.data.estab.delta_e_n_single = data.frame(delta_exo_nat=seq(min(sp_mammal_exo_nat$delta_exo_nat),
                                                                      max(sp_mammal_exo_nat$delta_exo_nat),
                                                                      length=100))

lincombs.mammal.matrix.estab.delta_e_n_single=model.matrix(~delta_exo_nat,
                                                          data=lincombs.mammal.data.estab.delta_e_n_single)
lincombs.mammal.matrix.estab.delta_e_n_single=as.data.frame(lincombs.mammal.matrix.estab.delta_e_n_single)
lincombs.mammal.estab.delta_e_n_single=inla.make.lincombs(lincombs.mammal.matrix.estab.delta_e_n_single)

inla.model_lincombs.mammal.estab.delta_e_n_single = inla(presence ~ delta_exo_nat
                                                        + f(RegionID, model="iid")  
                                                        #+ f(species, delta_exo_nat, model="iid")
                                                        ,
                                                        data = sp_mammal_exo_nat,
                                                        family = "zeroinflatedbinomial1",
                                                        control.compute = list(dic = TRUE,
                                                                               waic = TRUE, 
                                                                               cpo = TRUE,
                                                                               config = TRUE),
                                                        control.predictor = list(compute = TRUE),
                                                        quantiles = c(0.025, 0.5, 0.975),
                                                        lincomb = lincombs.mammal.estab.delta_e_n_single)

lincombs.mammal.posterior.estab.delta_e_n_single = inla.model_lincombs.mammal.estab.delta_e_n_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.mammal.estab.delta_e_n_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.mammal.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.mammal.data.estab.delta_e_n_single$predicted.value=unlist(lapply(lincombs.mammal.posterior.estab.delta_e_n_single,
                                                                         function(x)inla.emarginal(fun=plogis,x)))
lincombs.mammal.data.estab.delta_e_n_single$lower=unlist(lapply(lincombs.mammal.posterior.estab.delta_e_n_single,
                                                               function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.mammal.data.estab.delta_e_n_single$upper=unlist(lapply(lincombs.mammal.posterior.estab.delta_e_n_single,
                                                               function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))


#lincombs.mammal.data.estab.delta_e_n_single
save(lincombs.mammal.data.estab.delta_e_n_single, 
     file = 'results/primary_results/fitted_mods/lincombs.mammal.data.estab.delta_e_n_single.rdata')


##### Draw logistic curve for mammal establishment ~ delta_exo_nat #####
load("results/primary_results/fitted_mods/lincombs.mammal.data.estab.delta_e_n_single.rdata")

(estab.mammal.delta_ED_single.logistic=ggplot(data=lincombs.mammal.data.estab.delta_e_n_single,
                                    aes(x = delta_exo_nat,
                                        y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_2d[1],alpha=0.2)+
    geom_line(color=colors_2d[1],size=1.2,
              linetype = 1)+
    #scale_color_gradient(low = turbo(4)[2],
                 #        high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    #geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
            #   color = colors_4d[1],
            #   shape=1,
            #   alpha = 0.2,
            #   position=position_jitter(height=0.02))+
    labs(x=' ', y='Nutralization probability')+
    annotate(geom="text",x=min(lincombs.mammal.data.estab.delta_e_n_single$delta_exo_nat)+
               (max(lincombs.mammal.data.estab.delta_e_n_single$delta_exo_nat)-
                  min(lincombs.mammal.data.estab.delta_e_n_single$delta_exo_nat))*0.5,
             y=c(max(lincombs.mammal.data.estab.delta_e_n_single$upper),
                 max(lincombs.mammal.data.estab.delta_e_n_single$upper)*0.90),
             label = c("italic(beta)[Delta*'ED'] == '0.10'",
                       "'95% CI' == '[0.04, 0.17]'"),
             parse=T,size=4)+
    ggtitle('Mammals')+
    theme_test()+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
      plot.margin=unit(c(0.4,0.4,0.4,0.4),units="lines")
    ) 
)

rm(inla.model_lincombs.mammal.estab.delta_e_n_single)


#### Plant ED: calculation & mapping ####
gc()
load("data/Plants/data/df.native.natu.species.650.nonative.Rdata")
load('results/primary_results/distances_beta/ED_plant_native.rdata')
phylo_big = phytools::read.newick("code/FYI/Cai_et_al_2024_NEE/data/v0.1/ALLOTB.tre")
shp.glonaf.new = st_read("data/Plants/shp_glonaf_new_eck4.shp")

phy_plant = phylo_big
shp.glonaf.trans = shp.glonaf.new 
rm(shp.glonaf.new)

colnames(df.native.650) == colnames(df.natu.650)

##### exotics #####
#1. create the species presence/absence matrix with rows as sampling regions and columns as species
plant_sps_ED = phyloregion::evol_distinct(tree = phy_plant, type = "fair.proportion")
plant_sps_ED = data.frame(species = names(plant_sps_ED), 
                           ED = plant_sps_ED)


colnames(df.natu.650)
df.natu.650$presence = 1
comm_plant_exotic = df.natu.650 %>% 
  complete(species, Region_id, fill = list(presence = 0)) %>%
  dplyr::select(c('Region_id', 'species', 'presence')) %>% 
  pivot_wider(names_from = species,
              values_from = presence,
              values_fn = mean) %>% 
  complete(Region_id = unique(shp.glonaf.trans$Region_id)) %>% ## assume that 
  ## the absence regions have no naturalized aliens
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
  dplyr::arrange('Region_id') %>% 
  pivot_longer(cols = Abelia_chinensis:Zygophyllum_fabago, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(plant_sps_ED, by = 'species') %>% 
  rename(c('ED' = 'exotic_ED'))


mean_ED_natives_plant = cbind(Region_id = sort(unique(df.native.650$Region_id)),
                              ED_plant_native$df)

sp_plant_exo_nat = comm_plant_exotic %>% 
  left_join(mean_ED_natives_plant, by = 'Region_id') %>% 
  rename(c('mean_ED' = 'mean_native_ED')) %>% 
  rename(c('mean_EDR' = 'mean_native_EDR'))

gc()

sp_plant_exo_nat$delta_exo_nat = log(sp_plant_exo_nat$exotic_ED/
                                        sp_plant_exo_nat$mean_native_ED)

nrow(sp_plant_exo_nat %>% 
       filter(presence == 1)) / nrow(sp_plant_exo_nat)

plot(sp_plant_exo_nat$delta_exo_nat, sp_plant_exo_nat$presence)


mod_plant = lme4::glmer(presence ~ delta_exo_nat + (1|Region_id), 
                   family = binomial, data = sp_plant_exo_nat)
summary(mod_plant)

mod_plant_sp_inv_delta_ed_inla = inla(presence ~ delta_exo_nat + f(Region_id, model="iid")+ 
                                        f(species, delta_exo_nat, model="iid"),
                 family = "zeroinflatedbinomial1",
                 data = sp_plant_exo_nat,
                 verbose = T
                 #Ntrials=Ntrials
)

summary(mod_plant_sp_inv_delta_ed_inla)

save(mod_plant_sp_inv_delta_ed_inla,
     file = 'results/primary_results/mod_plant_sp_inv_delta_ed_inla.rdata')
summary(mod_plant_inla)



##### Draw logistic curve for plant establishment ~ delta_exo_nat #####
load("results/primary_results/lincombs.plant.data.estab.delta_e_n_single_ran2.rdata")
lincombs.plant.data.estab.delta_e_n_single_ran2$summary[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.plant.data.estab.delta_e_n_single = lincombs.plant.data.estab.delta_e_n_single_all$prediction
  
load("~/my_pc/Global_ED/results/primary_results/fitted_mods/lincombs.plant.data.estab.delta_e_n_single_all.rdata")
lincombs.plant.data.estab.delta_e_n_single_all$summary[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.plant.data.estab.delta_e_n_single = lincombs.plant.data.estab.delta_e_n_single_all$prediction

(estab.plant.delta_ED_single.logistic=ggplot(data=lincombs.plant.data.estab.delta_e_n_single,
                                              aes(x = delta_exo_nat,
                                                  y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_2d[2],alpha=0.2)+
    geom_line(color=colors_2d[2],size=1.2,
              linetype = 1)+
    #scale_color_gradient(low = turbo(4)[2],
    #        high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    #geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
    #   color = colors_4d[1],
    #   shape=1,
    #   alpha = 0.2,
    #   position=position_jitter(height=0.02))+
    labs(x=' ', y='Nutralization probability')+
    annotate(geom="text",x=min(lincombs.plant.data.estab.delta_e_n_single$delta_exo_nat)+
               (max(lincombs.plant.data.estab.delta_e_n_single$delta_exo_nat)-
                  min(lincombs.plant.data.estab.delta_e_n_single$delta_exo_nat))*0.5,
             y=c(max(lincombs.plant.data.estab.delta_e_n_single$upper),
                 max(lincombs.plant.data.estab.delta_e_n_single$upper)*0.90),
             label = c("italic(beta)[Delta*'ED'] == '-0.20'",
                       "'95% CI' == '[-0.21, -0.19]'"),
             parse=T,size=4)+
    ggtitle('Plants')+
    theme_test()+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
      plot.margin=unit(c(0.4,0.4,0.4,0.4),units="lines")
    ) 
)



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


load('results/primary_results/ED_bird_native.rdata')

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
  arrange('RegionID') %>% 
  pivot_longer(cols = Acanthis_flammea:Zosterops_natalis, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(bird_sps_ED, by = 'species') %>% 
  rename('exotic_ED' = 'ED')


mean_ED_natives_bird = cbind(RegionID = df.trans$RegionID,
                              ED_bird_native$df)

sp_bird_exo_nat = comm_bird_exotic %>% 
  left_join(mean_ED_natives_bird, by = 'RegionID') %>% 
  rename('mean_native_ED' = 'mean_ED') %>% 
  rename('mean_native_EDR' = 'mean_EDR')


sp_bird_exo_nat$delta_exo_nat = log(sp_bird_exo_nat$exotic_ED/
                                       sp_bird_exo_nat$mean_native_ED)
nrow(sp_bird_exo_nat %>% 
       filter(presence == 1)) / nrow(sp_bird_exo_nat)

plot(sp_bird_exo_nat$delta_exo_nat, sp_bird_exo_nat$presence)


mod_bird = glmer(presence ~ delta_exo_nat + (1|RegionID) + (delta_exo_nat | species), 
                 data = sp_bird_exo_nat,
                 family = binomial(link = "cloglog"),
                 nAGQ=0, 
                 control=glmerControl(optimizer = "nloptwrap"))
summary(mod_bird)

library(INLA)

mod_bird_sp_inv_ed_inla = inla(presence ~ delta_exo_nat + f(RegionID, model="iid") + 
                                 f(species, delta_exo_nat, model="iid"),
                               family = "zeroinflatedbinomial1",
                                data = sp_bird_exo_nat,
                                verbose = F,
                               Ntrials = rep(1, nrow(sp_bird_exo_nat)),
                               control.family = list(
                                 hyper = list(theta = list(prior = "logitbeta",
                                                                   param = c(1, 1)))
                               )
                                #Ntrials=Ntrials
)

save(mod_bird_sp_inv_ed_inla,
     file = 'results/primary_results/mod_bird_sp_inv_ed_inla.rdata')
summary(mod_bird_sp_inv_ed_inla)



##### Draw logistic curve for bird establishment ~ delta_exo_nat #####
load("results/primary_results/lincombs.bird.data.estab.delta_e_n_single_all.rdata")
lincombs.bird.data.estab.delta_e_n_single_all$summary[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.bird.data.estab.delta_e_n_single = lincombs.bird.data.estab.delta_e_n_single_all$prediction

(estab.bird.delta_ED_single.logistic=ggplot(data=lincombs.bird.data.estab.delta_e_n_single,
                                             aes(x = delta_exo_nat,
                                                 y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_2d[2],alpha=0.2)+
    geom_line(color=colors_2d[2],size=1.2,
              linetype = 1)+
    #scale_color_gradient(low = turbo(4)[2],
    #        high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    #geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
    #   color = colors_4d[1],
    #   shape=1,
    #   alpha = 0.2,
    #   position=position_jitter(height=0.02))+
    labs(x='log (ED_exotic / MED_native)', y='Nutralization probability')+
    annotate(geom="text",x=min(lincombs.bird.data.estab.delta_e_n_single$delta_exo_nat)+
               (max(lincombs.bird.data.estab.delta_e_n_single$delta_exo_nat)-
                  min(lincombs.bird.data.estab.delta_e_n_single$delta_exo_nat))*0.5,
             y=c(max(lincombs.bird.data.estab.delta_e_n_single$upper),
                 max(lincombs.bird.data.estab.delta_e_n_single$upper)*0.90),
             label = c("italic(beta)[Delta*'ED'] == '-0.22'",
                       "'95% CI' == '[-0.28, -0.16]'"),
             parse=T,size=4)+
    ggtitle('Birds')+
    theme_test()+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
      plot.margin=unit(c(0.4,0.4,0.4,0.4),units="lines")
    ) 
)




#### Fish ED: calculation & mapping ####
load("D:/R projects/Global_ED/data/Fishes/data/my_phy.rdata")
phy_fish = phylo

load("D:/R projects/Global_ED/data/Fishes/data/my_data_used_final.rdata")
df = st_read("data/Fishes/data/Basin042017_3119.shp")
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

load('results/primary_results/ED_fish_native.rdata')

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
  arrange('X1.Basin.Name') %>% 
  pivot_longer(cols = Abbottina_rivularis:Zacco_platypus, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(fish_sps_ED, by = 'species') %>% 
  rename('exotic_ED' = 'ED')

mean_ED_natives_fish = cbind(Basin.Name = sort(unique(data.used_final_natives$X1.Basin.Name)),
                             ED_fish_native$df)

sp_fish_exo_nat = comm_fish_exotic %>% 
  left_join(mean_ED_natives_fish, by = join_by('X1.Basin.Name' == 'Basin.Name')) %>% 
  rename('mean_native_ED' = 'mean_ED') %>% 
  rename('mean_native_EDR' = 'mean_EDR')

sp_fish_exo_nat$delta_exo_nat = log(sp_fish_exo_nat$exotic_ED/
                                    sp_fish_exo_nat$mean_native_ED)
nrow(sp_fish_exo_nat %>% 
       filter(presence == 1)) / nrow(sp_fish_exo_nat)

plot(sp_fish_exo_nat$delta_exo_nat, sp_fish_exo_nat$presence)


mod_fish = lme4::glmer(presence ~ scale(delta_exo_nat) + (1|X1.Basin.Name)
                       + (delta_exo_nat | species), 
                       data = sp_fish_exo_nat,
                       family = binomial(link = "cloglog"),
                       nAGQ=0, 
                       control=glmerControl(optimizer = "nloptwrap"))
summary(mod_fish)

library(INLA)

mod_fish_sp_inv_ed_inla = inla(presence ~ delta_exo_nat + f(X1.Basin.Name, model="iid"),
                               family = "zeroinflatedbinomial1",
                               data = sp_fish_exo_nat,
                               verbose = F
                               #Ntrials=Ntrials
)

save(mod_fish_sp_inv_ed_inla,
    file = 'results/primary_results/mod_fish_sp_inv_ed_inla.rdata')
summary(mod_fish_sp_inv_ed_inla)



##### Draw logistic curve for fish establishment ~ delta_exo_nat #####
load("results/primary_results/lincombs.fish.data.estab.delta_e_n_single_all.rdata")
lincombs.fish.data.estab.delta_e_n_single_all$summary[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.fish.data.estab.delta_e_n_single = lincombs.fish.data.estab.delta_e_n_single_all$prediction

(estab.fish.delta_ED_single.logistic=ggplot(data=lincombs.fish.data.estab.delta_e_n_single,
                                            aes(x = delta_exo_nat,
                                                y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_2d[2],alpha=0.2)+
    geom_line(color=colors_2d[2],size=1.2,
              linetype = 1)+
    #scale_color_gradient(low = turbo(4)[2],
    #        high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
    #   color = colors_4d[1],
    #   shape=1,
    #   alpha = 0.2,
    #   position=position_jitter(height=0.02))+
    labs(x='log (ED_exotic / MED_native)', y='Nutralization probability')+
    annotate(geom="text",x=min(lincombs.fish.data.estab.delta_e_n_single$delta_exo_nat)+
               (max(lincombs.fish.data.estab.delta_e_n_single$delta_exo_nat)-
                  min(lincombs.fish.data.estab.delta_e_n_single$delta_exo_nat))*0.5,
             y=c(max(lincombs.fish.data.estab.delta_e_n_single$upper),
                 max(lincombs.fish.data.estab.delta_e_n_single$upper)*0.90),
             label = c("italic(beta)[Delta*'ED'] == '-0.48'",
                       "'95% CI' == '[-0.51, -0.45]'"),
             parse=T,size=4)+
    ggtitle('Fishes')+
    theme_test()+
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, 
                                hjust = 0.5),
      plot.margin=unit(c(0.4,0.4,0.4,0.4),units="lines")
    ) 
)


#### Merge and export the plots ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

##### Mean ED #####
plot_width = 0.48
plot_height = 0.48  # Square-ish proportions

figs_sp_delta_ED_estab = ggdraw() +
  # Top row
  draw_plot(estab.mammal.delta_ED_single.logistic,   x = 0.00, y = 0.50, width = plot_width, height = plot_height) +
  draw_plot(estab.plant.delta_ED_single.logistic,   x = 0.50, y = 0.50, width = plot_width, height = plot_height) +
  # Bottom row
  draw_plot(estab.bird.delta_ED_single.logistic, x = 0.00, y = 0.00, width = plot_width, height = plot_height) +
  draw_plot(estab.fish.delta_ED_single.logistic,  x = 0.50, y = 0.00, width = plot_width, height = plot_height) +
  
  # Labels in top-left corner of each plot
  draw_plot_label(
    label = c("a", "b", "c", "d"),
    size = 14,
    x = c(0.01, 0.51, 0.01, 0.51),
    y = c(0.97, 0.97, 0.47, 0.47)
  )


emf('figures/figs_sp_delta_ED_estab.emf',
 height=25 * 0.6, width=25 * 0.6, coordDPI = 600, 
 emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
 #ensure text looks correct on the viewing system
 units = 'cm')
figs_sp_delta_ED_estab
dev.off() #turn off device and finalize file

emf('figures/figs_sp_delta_ED_estab_PPT.emf',
    height=25 * 0.6, width=40 * 0.6, coordDPI = 600, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    #ensure text looks correct on the viewing system
    units = 'cm')
figs_sp_delta_ED_estab
dev.off() #turn off device and finalize file



#emf('figures/figs_ED_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_ED_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_ED_nat_extant = ggdraw() +
  draw_plot(plants_native_ED_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_ED_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_ED_map_all, x = 0, y = 0.53, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_ED_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_ED_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_ED_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_ED_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_ED_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_ED_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_ED_nat_extant
dev.off() #turn off device and finalize file


# Compare the PE patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.015         # small horizontal gap between columns

figs_ED_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_ED_map_all,  x = 0,                 y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_ED_map_all,  x = plot_width + gap,  y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_delta_ED_map_all,   x = 2*(plot_width + gap), y = 0.75, width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_ED_map_all,   x = 0,                 y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_extant_ED_map_all,   x = plot_width + gap,  y = 0.5,  width = plot_width, height = plot_height) +
  draw_plot(birds_delta_ED_map_all,    x = 2*(plot_width + gap), y = 0.5,  width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_ED_map_all, x = 0,                 y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_ED_map_all, x = plot_width + gap,  y = 0.25, width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_ED_map_all,  x = 2*(plot_width + gap), y = 0.25, width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_ED_map_all,   x = 0,                 y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_ED_map_all,   x = plot_width + gap,  y = 0,    width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_ED_map_all,   x = 2*(plot_width + gap), y = 0,    width = plot_width, height = plot_height) +
  
  # Labels a–l
  draw_plot_label(
    label = letters[1:12],
    size = 13,
    x = rep(c(0.01, plot_width + gap + 0.01, 2*(plot_width + gap) + 0.01), 4),
    y = rep(c(0.75 + plot_height - 0.015,
              0.5 + plot_height - 0.015,
              0.25 + plot_height - 0.015,
              0 + plot_height - 0.015), each = 3)
  )

png(filename = 'figures/figs_ED_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_ED_nat_extant_delta
dev.off() #turn off device and finalize file



##### Mean EDR #####
plot_width = 0.5
plot_height = 0.20

figs_EDR_nat_exo = ggdraw() +
  draw_plot(plants_native_EDR_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_exotic_EDR_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_EDR_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_exotic_EDR_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_EDR_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_exotic_EDR_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_EDR_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_exotic_EDR_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_EDR_nat_exo.png',
    #plot = figs_EDR_nat_exo,
    height=40, width=25, # setting emfPlusFontToPath=TRUE to 
    res = 300,
    # ensure text looks correct on the viewing system
    units = 'cm')
figs_EDR_nat_exo
dev.off() #turn off device and finalize file


#emf('figures/figs_EDR_nat_exo.emf',
# height=40, width=25, coordDPI = 1, 
# emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
# ensure text looks correct on the viewing system
# units = 'cm')
#figs_EDR_nat_exo
#dev.off() #turn off device and finalize file

# Compare the PE patterns of only natives to including native + exotics
figs_EDR_nat_extant = ggdraw() +
  draw_plot(plants_native_EDR_map_all, x = 0, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(plants_extant_EDR_map_all, x = 0.5, y = 0.75, width = plot_width, height = plot_height) +
  draw_plot(birds_native_EDR_map_all, x = 0, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(birds_extant_EDR_map_all, x = 0.5, y = 0.54, width = plot_width, height = plot_height) +
  draw_plot(mammals_native_EDR_map_all, x = 0, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_EDR_map_all, x = 0.5, y = 0.33, width = plot_width, height = plot_height) +
  draw_plot(fishes_native_EDR_map_all, x = 0, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_EDR_map_all, x = 0.5, y = 0.12, width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f", "g", "h"),
    size = 14,
    x = c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52, 0.02, 0.52),
    y = c(0.75+plot_height-0.02, 0.75+plot_height-0.02,
          0.54+plot_height-0.02, 0.54+plot_height-0.02,
          0.33+plot_height-0.02, 0.33+plot_height-0.02,
          0.12+plot_height-0.02, 0.12+plot_height-0.02)
  )


png(filename = 'figures/figs_EDR_nat_extant.png',
    height=40, width=25, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)
figs_EDR_nat_extant
dev.off() #turn off device and finalize file


# Compare the PE patterns of natives, native + exotics, and their differences

# Adjusted dimensions: small horizontal gap
plot_width = 0.31   # < 1/3 to allow spacing
plot_height = 0.24  # keep a compact vertical layout
gap = 0.015         # small horizontal gap between columns

figs_EDR_nat_extant_delta = ggdraw() +
  # Row 1 – Plants
  draw_plot(plants_native_EDR_map_all,  x = 0,                 y = 0.75,
            width = plot_width, height = plot_height) +
  draw_plot(plants_extant_EDR_map_all,  x = plot_width + gap,  y = 0.75,
            width = plot_width, height = plot_height) +
  draw_plot(plants_delta_EDR_map_all,   x = 2*(plot_width + gap), y = 0.75,
            width = plot_width, height = plot_height) +
  
  # Row 2 – Birds
  draw_plot(birds_native_EDR_map_all,   x = 0,                 y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(birds_extant_EDR_map_all,   x = plot_width + gap,  y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(birds_delta_EDR_map_all,    x = 2*(plot_width + gap), y = 0.5,
            width = plot_width, height = plot_height) +
  
  # Row 3 – Mammals
  draw_plot(mammals_native_EDR_map_all, x = 0,                 y = 0.25,
            width = plot_width, height = plot_height) +
  draw_plot(mammals_extant_EDR_map_all, x = plot_width + gap,  y = 0.25,
            width = plot_width, height = plot_height) +
  draw_plot(mammals_delta_EDR_map_all,  x = 2*(plot_width + gap), y = 0.25,
            width = plot_width, height = plot_height) +
  
  # Row 4 – Fishes
  draw_plot(fishes_native_EDR_map_all,   x = 0,                 y = 0,
            width = plot_width, height = plot_height) +
  draw_plot(fishes_extant_EDR_map_all,   x = plot_width + gap,  y = 0,
            width = plot_width, height = plot_height) +
  draw_plot(fishes_delta_EDR_map_all,   x = 2*(plot_width + gap), y = 0,
            width = plot_width, height = plot_height) +
  
  # Labels a–l
  draw_plot_label(
    label = letters[1:12],
    size = 13,
    x = rep(c(0.01, plot_width + gap + 0.01, 2*(plot_width + gap) + 0.01), 4),
    y = rep(c(0.75 + plot_height - 0.015,
              0.5 + plot_height - 0.015,
              0.25 + plot_height - 0.015,
              0 + plot_height - 0.015), each = 3)
  )

png(filename = 'figures/figs_EDR_nat_extant_delta.png',
    height=40, width=25*1.5, units = 'cm', # setting emfPlusFontToPath=TRUE to 
    res = 600)

figs_EDR_nat_extant_delta
dev.off() #turn off device and finalize file

