### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
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

setwd("D:/R projects/Global_ED")
source('code/functions/calculating_ED_func.R')

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
colors5 = scico::scico(n=7, begin = 0, end = 0.4, palette = "bam")
colors5[length(colors5)+1] = colors4[length(colors4)-2]

colors_2d = scico::scico(n=2, begin = 0.2, end = 0.7, palette = "bam")

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


#### Mammal ED: calculation & mapping ####
load("data/Mammals/Results_data/native/sp_nati.rdata")
load("data/Mammals/Results_data/alien/sp_alien.rdata")
load("data/Mammals/Results_data/all_phy.rdata")

load('results/primary_results/ED_mammal_native.rdata')

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
  arrange('RegionID') %>% 
  pivot_longer(cols = Alces_alces:Wallabia_bicolor, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(mammal_sps_ED, by = 'species') %>% 
  rename('exotic_ED' = 'ED')



mean_ED_natives_mammal = cbind(RegionID = sort(unique(sp_dis_5$Region.ID)),
                               ED_mammal_native$df)

sp_mammal_exo_nat = comm_mammal_exotic %>% 
  left_join(mean_ED_natives_mammal, by = 'RegionID') %>% 
  rename('mean_native_ED' = 'mean_ED') %>% 
  rename('mean_native_EDR' = 'mean_EDR')


sp_mammal_exo_nat$delta_exo_nat = log(sp_mammal_exo_nat$exotic_ED/
                                        sp_mammal_exo_nat$mean_native_ED)
nrow(sp_mammal_exo_nat %>% 
       filter(presence == 1)) / nrow(sp_mammal_exo_nat)

sp_mammal_exo_nat = sp_mammal_exo_nat %>% filter(!is.na(delta_exo_nat))

library(INLA)

#mod_mammal_sp_exotic_ed_inla = inla(presence ~ exotic_ED
              #            + f(RegionID, model="iid")
               #           ,
                #          family = "zeroinflatedbinomial1",
                #          data = sp_mammal_exo_nat 
#)

#summary(mod_mammal_sp_exotic_ed_inla)


### get predict data for mod_mammal_sp_exotic_ed_inla
lincombs.mammal.data.estab.exotic_ED_single = data.frame(exotic_ED=seq(min(sp_mammal_exo_nat$exotic_ED),
                                                                       max(sp_mammal_exo_nat$exotic_ED),
                                                                           length=100))

lincombs.mammal.matrix.estab.exotic_ED_single=model.matrix( ~ exotic_ED,
                                                           data=lincombs.mammal.data.estab.exotic_ED_single)
lincombs.mammal.matrix.estab.exotic_ED_single=as.data.frame(lincombs.mammal.matrix.estab.exotic_ED_single)
lincombs.mammal.estab.exotic_ED_single=inla.make.lincombs(lincombs.mammal.matrix.estab.exotic_ED_single)

inla.model_lincombs.mammal.estab.exotic_ED_single = inla(presence ~ exotic_ED
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
                                                         lincomb = lincombs.mammal.estab.exotic_ED_single)

lincombs.mammal.posterior.estab.exotic_ED_single = inla.model_lincombs.mammal.estab.exotic_ED_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.mammal.estab.exotic_ED_single$summary.fixed[c(1,3,5)] %>% round(3)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.mammal.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
f_posi = which(!is.infinite(sapply(lincombs.mammal.posterior.estab.exotic_ED_single,
       function(x)sum(x))))
lincombs.mammal.data.estab.exotic_ED_single = data.frame(exotic_ED = lincombs.mammal.data.estab.exotic_ED_single[f_posi,])
lincombs.mammal.data.estab.exotic_ED_single$predicted.value=unlist(lapply(lincombs.mammal.posterior.estab.exotic_ED_single[f_posi],
                                                                          function(x)inla.emarginal(fun=plogis,x)))
lincombs.mammal.data.estab.exotic_ED_single$lower=unlist(lapply(lincombs.mammal.posterior.estab.exotic_ED_single[f_posi],
                                                                function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.mammal.data.estab.exotic_ED_single$upper=unlist(lapply(lincombs.mammal.posterior.estab.exotic_ED_single[f_posi],
                                                                function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))


lincombs.mammal.data.estab.exotic_ED_single_all = list(summary = inla.model_lincombs.mammal.estab.exotic_ED_single$summary.fixed,
                                                       prediction = lincombs.mammal.data.estab.exotic_ED_single)
#lincombs.mammal.data.estab.exotic_ED_single
save(lincombs.mammal.data.estab.exotic_ED_single_all, 
     file = 'results/primary_results/lincombs.mammal.data.estab.exotic_ED_single_all.rdata')

##### Draw logistic curve for mammal establishment ~ exotic_ED #####
load("results/primary_results/lincombs.mammal.data.estab.exotic_ED_single_all.rdata")
lincombs.mammal.data.estab.exotic_ED_single_all$summary[c(1,3,5)]%>%round(3)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.mammal.data.estab.exotic_ED_single = lincombs.mammal.data.estab.exotic_ED_single_all$prediction


(estab.mammal.exotic_ED_single.logistic=ggplot(data=lincombs.mammal.data.estab.exotic_ED_single,
                                              aes(x = exotic_ED,
                                                  y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_2d[2],alpha=0.2)+
    geom_line(color=colors_2d[2],size=1.2,
              linetype = 2)+
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
    annotate(geom="text",x=min(lincombs.mammal.data.estab.exotic_ED_single$exotic_ED)+
               (max(lincombs.mammal.data.estab.exotic_ED_single$exotic_ED)-
                  min(lincombs.mammal.data.estab.exotic_ED_single$exotic_ED))*0.5,
             y=c(max(lincombs.mammal.data.estab.exotic_ED_single$upper),
                 max(lincombs.mammal.data.estab.exotic_ED_single$upper)*0.90),
             label = c("italic(beta)[Exotic_ED] == '-0.006'",
                       "'95% CI' == '[-0.014, 0.001]'"),
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
# exotic_ED   -0.006     -0.014      0.001

#### Plant ED: calculation & mapping ####
load("D:/R projects/Global_ED/data/Plants/data/shp.651.Rdata")
load("D:/R projects/Global_ED/data/Plants/data/phylo.fake.species.653.Rdata")
phy_plant = phylo
load("D:/R projects/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
load('results/primary_results/ED_plant_native.rdata')

shp.glonaf.trans = st_transform(shp.glonaf.new, crs = "+proj=eck4") 
shp.glonaf.trans$area = st_area(shp.glonaf.trans)
sum(as.numeric(shp.glonaf.trans$area) == 0)

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
  arrange('Region_id') %>% 
  pivot_longer(cols = Abelia_chinensis:Zygophyllum_fabago, 
               names_to = 'species', values_to = 'presence') %>% 
  left_join(plant_sps_ED, by = 'species') %>% 
  rename('exotic_ED' = 'ED')


mean_ED_natives_plant = cbind(Region_id = sort(unique(df.native.650$Region_id)),
                              ED_plant_native$df)

sp_plant_exo_nat = comm_plant_exotic %>% 
  left_join(mean_ED_natives_plant, by = 'Region_id') %>% 
  rename('mean_native_ED' = 'mean_ED') %>% 
  rename('mean_native_EDR' = 'mean_EDR')


sp_plant_exo_nat$delta_exo_nat = log(sp_plant_exo_nat$exotic_ED/
                                       sp_plant_exo_nat$mean_native_ED)
nrow(sp_plant_exo_nat %>% 
       filter(presence == 1)) / nrow(sp_plant_exo_nat)

#plot(sp_plant_exo_nat$exotic_ED, sp_plant_exo_nat$presence)

mod_plant_sp_exotic_ed_inla = inla(presence ~ exotic_ED + f(Region_id, model="iid")  
                                   #+ f(species, exotic_ED, model="iid")
                                   ,
                                   family = "zeroinflatedbinomial1",
                                   data = sp_plant_exo_nat,
                                   verbose = F
                                   #Ntrials=Ntrials
)



##### Draw logistic curve for plant establishment ~ delta_exo_nat #####
load("results/primary_results/lincombs.plant.data.estab.exotic_ED_single_all.rdata")
lincombs.plant.data.estab.exotic_ED_single_all$summary[c(1,3,5)]%>%round(3)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.plant.data.estab.exotic_ED_single = lincombs.plant.data.estab.exotic_ED_single$prediction

(estab.plant.exotic_ED_single.logistic=ggplot(data=lincombs.plant.data.estab.exotic_ED_single,
                                             aes(x = exotic_ED,
                                                 y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_2d[2], alpha=0.2)+
    geom_line(color=colors_2d[2], size=1.2,
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
    annotate(geom="text",x=min(lincombs.plant.data.estab.exotic_ED_single$exotic_ED)+
               (max(lincombs.plant.data.estab.exotic_ED_single$exotic_ED)-
                  min(lincombs.plant.data.estab.exotic_ED_single$exotic_ED))*0.5,
             y=c(max(lincombs.plant.data.estab.exotic_ED_single$upper),
                 max(lincombs.plant.data.estab.exotic_ED_single$upper)*0.90),
             label = c("italic(beta)[Exotic_ED] == '-0.021'",
                       "'95% CI' == '[-0.022, -0.021]'"),
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


mean_ED_natives_bird = cbind(RegionID = df_trans$RegionID,
                             ED_bird_native$df)

sp_bird_exo_nat = comm_bird_exotic %>% 
  left_join(mean_ED_natives_bird, by = 'RegionID') %>% 
  rename('mean_native_ED' = 'mean_ED') %>% 
  rename('mean_native_EDR' = 'mean_EDR')


sp_bird_exo_nat$delta_exo_nat = log(sp_bird_exo_nat$exotic_ED/
                                      sp_bird_exo_nat$mean_native_ED)
nrow(sp_bird_exo_nat %>% 
       filter(presence == 1)) / nrow(sp_bird_exo_nat)

mod_bird_sp_exotic_ed_inla = inla(presence ~ exotic_ED + f(RegionID, model="iid")
                                  #+ f(species, exotic_ED, model="iid")
                                  ,
                                  family = "zeroinflatedbinomial1",
                                  data = sp_bird_exo_nat,
                                  verbose = F
                                  #Ntrials=Ntrials
)



##### Draw logistic curve for bird establishment ~ delta_exo_nat #####
load("results/primary_results/lincombs.bird.data.estab.exotic_ED_single_all.rdata")
lincombs.bird.data.estab.exotic_ED_single_all$summary[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.bird.data.estab.exotic_ED_single = lincombs.bird.data.estab.exotic_ED_single_all$prediction

(estab.bird.exotic_ED_single.logistic=ggplot(data=lincombs.bird.data.estab.exotic_ED_single,
                                            aes(x = exotic_ED,
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
    labs(x='Exotic_ED', y='Nutralization probability')+
    annotate(geom="text",x=min(lincombs.bird.data.estab.exotic_ED_single$exotic_ED)+
               (max(lincombs.bird.data.estab.exotic_ED_single$exotic_ED)-
                  min(lincombs.bird.data.estab.exotic_ED_single$exotic_ED))*0.5,
             y=c(max(lincombs.bird.data.estab.exotic_ED_single$upper),
                 max(lincombs.bird.data.estab.exotic_ED_single$upper)*0.90),
             label = c("italic(beta)[Exotic_ED] == '-0.04'",
                       "'95% CI' == '[-0.05, -0.03]'"),
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

#plot(sp_fish_exo_nat$exotic_ED, sp_fish_exo_nat$presence)

mod_fish_sp_exotic_ed_inla = inla(presence ~ exotic_ED + f(X1.Basin.Name, model="iid") 
                                  # + f(species, exotic_ED, model="iid")
                                  ,
                                  family = "zeroinflatedbinomial1",
                                  data = sp_fish_exo_nat,
                                  verbose = F
                                  #Ntrials=Ntrials
)


##### Draw logistic curve for fish establishment ~ delta_exo_nat #####
load("results/primary_results/lincombs.fish.data.estab.exotic_ED_single_all.rdata")
lincombs.fish.data.estab.exotic_ED_single_all$summary[c(1,3,5)]%>%round(3)##Extracting effects and confidence intervals for prediction curves from raw data
lincombs.fish.data.estab.exotic_ED_single = lincombs.fish.data.estab.exotic_ED_single_all$prediction

(estab.fish.exotic_ED_single.logistic=ggplot(data=lincombs.fish.data.estab.exotic_ED_single,
                                            aes(x = exotic_ED,
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
    labs(x='Exotic_ED', y='Nutralization probability')+
    annotate(geom="text",x=min(lincombs.fish.data.estab.exotic_ED_single$exotic_ED)+
               (max(lincombs.fish.data.estab.exotic_ED_single$exotic_ED)-
                  min(lincombs.fish.data.estab.exotic_ED_single$exotic_ED))*0.5,
             y=c(max(lincombs.fish.data.estab.exotic_ED_single$upper),
                 max(lincombs.fish.data.estab.exotic_ED_single$upper)*0.90),
             label = c("italic(beta)[Exotic_ED] == '-0.020'",
                       "'95% CI' == '[-0.022, -0.018]'"),
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

##### Exotic_ED #####
plot_width = 0.48
plot_height = 0.48  # Square-ish proportions

figs_sp_exotic_ED_estab = ggdraw() +
  # Top row
  draw_plot(estab.mammal.exotic_ED_single.logistic,
            x = 0.00, y = 0.50, width = plot_width, height = plot_height) +
  draw_plot(estab.plant.exotic_ED_single.logistic,
            x = 0.50, y = 0.50, width = plot_width, height = plot_height) +
  # Bottom row
  draw_plot(estab.bird.exotic_ED_single.logistic,
            x = 0.00, y = 0.00, width = plot_width, height = plot_height) +
  draw_plot(estab.fish.exotic_ED_single.logistic,
            x = 0.50, y = 0.00, width = plot_width, height = plot_height) +
  
  # Labels in top-left corner of each plot
  draw_plot_label(
    label = c("a", "b", "c", "d"),
    size = 14,
    x = c(0.01, 0.51, 0.01, 0.51),
    y = c(0.97, 0.97, 0.47, 0.47)
  )


emf('figures/figs_sp_exotic_ED_estab.emf',
    height=25 * 0.6, width=25 * 0.6, coordDPI = 600, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    #ensure text looks correct on the viewing system
    units = 'cm')
figs_sp_exotic_ED_estab
dev.off() #turn off device and finalize file

emf('figures/figs_sp_exotic_ED_estab_PPT.emf',
    height=25 * 0.6, width=40 * 0.6, coordDPI = 600, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    #ensure text looks correct on the viewing system
    units = 'cm')
figs_sp_exotic_ED_estab
dev.off() #turn off device and finalize file


