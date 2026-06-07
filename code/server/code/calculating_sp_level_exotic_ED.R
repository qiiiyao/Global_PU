### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#### Package management ####
# list of packages
# Manually install the PhyloMeasures package
#install.packages("~/my_pc/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
#install.packages("INLA",repos=c(getOption("repos"),
# INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
rm(list = ls())
requirements = c("PhyloMeasures", "phangorn", "foreach", "dplyr", "ape", 
                 'tidyr', 'sf', 'raster', 'INLA')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})
library(INLA)
setwd("~/my_pc/Global_ED")
source('code/calculating_ED_func.R')

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
load("~/my_pc/Global_ED/data/Plants/data/shp.651.Rdata")
load("~/my_pc/Global_ED/data/Plants/data/phylo.fake.species.653.Rdata")
phy_plant = phylo
load("~/my_pc/Global_ED/data/Plants/data/df.native.natu.species.650.nonative.Rdata")
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

save(mod_plant_sp_exotic_ed_inla,
     file = 'results/primary_results/mod_plant_sp_exotic_ed_inla.rdata')
#load('results/primary_results/mod_plant_sp_exotic_ed_inla.rdata')
summary(mod_plant_sp_exotic_ed_inla)

### get predict data for mod_plant_sp_exotic_ed_inla
lincombs.plant.data.estab.exotic_ED_single = data.frame(exotic_ED=seq(min(sp_plant_exo_nat$exotic_ED),
                                                                      max(sp_plant_exo_nat$exotic_ED),
                                                                      length=100))

lincombs.plant.matrix.estab.exotic_ED_single=model.matrix( ~ exotic_ED,
                                                          data=lincombs.plant.data.estab.exotic_ED_single)
lincombs.plant.matrix.estab.exotic_ED_single=as.data.frame(lincombs.plant.matrix.estab.exotic_ED_single)
lincombs.plant.estab.exotic_ED_single=inla.make.lincombs(lincombs.plant.matrix.estab.exotic_ED_single)

inla.model_lincombs.plant.estab.exotic_ED_single = inla(presence ~ exotic_ED
                                                        + f(Region_id, model="iid")  
                                                        #+ f(species, exotic_ED, model="iid")
                                                        ,
                                                        data = sp_plant_exo_nat,
                                                        family = "zeroinflatedbinomial1",
                                                        control.compute = list(dic = TRUE,
                                                                               waic = TRUE, 
                                                                               cpo = TRUE,
                                                                               config = TRUE),
                                                        control.predictor = list(compute = TRUE),
                                                        quantiles = c(0.025, 0.5, 0.975),
                                                        lincomb = lincombs.plant.estab.exotic_ED_single)

lincombs.plant.posterior.estab.exotic_ED_single = inla.model_lincombs.plant.estab.exotic_ED_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.plant.estab.exotic_ED_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.plant.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.plant.data.estab.exotic_ED_single$predicted.value=unlist(lapply(lincombs.plant.posterior.estab.exotic_ED_single,
                                                                         function(x)inla.emarginal(fun=plogis,x)))
lincombs.plant.data.estab.exotic_ED_single$lower=unlist(lapply(lincombs.plant.posterior.estab.exotic_ED_single,
                                                               function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.plant.data.estab.exotic_ED_single$upper=unlist(lapply(lincombs.plant.posterior.estab.exotic_ED_single,
                                                               function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

lincombs.plant.data.estab.exotic_ED_single_all = list(summary = inla.model_lincombs.plant.estab.exotic_ED_single$summary.fixed,
                                                  prediction = lincombs.plant.data.estab.exotic_ED_single)
#lincombs.plant.data.estab.exotic_ED_single
save(lincombs.plant.data.estab.exotic_ED_single_all, 
     file = 'results/primary_results/lincombs.plant.data.estab.exotic_ED_single_all.rdata')





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
save(mod_bird_sp_exotic_ed_inla,
     file = 'results/primary_results/mod_bird_sp_exotic_ed_inla.rdata')
#load('results/primary_results/mod_bird_sp_exotic_ed_inla.rdata')
#summary(mod_bird_sp_exotic_ed_inla)





### get predict data for mod_bird_sp_exotic_ed_inla
lincombs.bird.data.estab.exotic_ED_single = data.frame(exotic_ED=seq(min(sp_bird_exo_nat$exotic_ED),
                                                                      max(sp_bird_exo_nat$exotic_ED),
                                                                      length=100))

lincombs.bird.matrix.estab.exotic_ED_single=model.matrix(~exotic_ED,
                                                          data=lincombs.bird.data.estab.exotic_ED_single)
lincombs.bird.matrix.estab.exotic_ED_single=as.data.frame(lincombs.bird.matrix.estab.exotic_ED_single)
lincombs.bird.estab.exotic_ED_single=inla.make.lincombs(lincombs.bird.matrix.estab.exotic_ED_single)

inla.model_lincombs.bird.estab.exotic_ED_single = inla(presence ~ exotic_ED
                                                        + f(RegionID, model="iid")  
                                                        #+ f(species, exotic_ED, model="iid")
                                                        ,
                                                        data = sp_bird_exo_nat,
                                                        family = "zeroinflatedbinomial1",
                                                        control.compute = list(dic = TRUE,
                                                                               waic = TRUE, 
                                                                               cpo = TRUE,
                                                                               config = TRUE),
                                                        control.predictor = list(compute = TRUE),
                                                        quantiles = c(0.025, 0.5, 0.975),
                                                        lincomb = lincombs.bird.estab.exotic_ED_single)

lincombs.bird.posterior.estab.exotic_ED_single = inla.model_lincombs.bird.estab.exotic_ED_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.bird.estab.exotic_ED_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.bird.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.bird.data.estab.exotic_ED_single$predicted.value=unlist(lapply(lincombs.bird.posterior.estab.exotic_ED_single,
                                                                         function(x)inla.emarginal(fun=plogis,x)))
lincombs.bird.data.estab.exotic_ED_single$lower=unlist(lapply(lincombs.bird.posterior.estab.exotic_ED_single,
                                                               function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.bird.data.estab.exotic_ED_single$upper=unlist(lapply(lincombs.bird.posterior.estab.exotic_ED_single,
                                                               function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

lincombs.bird.data.estab.exotic_ED_single_all = list(summary = inla.model_lincombs.bird.estab.exotic_ED_single$summary.fixed,
                                                  prediction = lincombs.bird.data.estab.exotic_ED_single)
#lincombs.bird.data.estab.exotic_ED_single
save(lincombs.bird.data.estab.exotic_ED_single_all, 
     file = 'results/primary_results/lincombs.bird.data.estab.exotic_ED_single_all.rdata')





#### Fish ED: calculation & mapping ####
load("~/my_pc/Global_ED/data/Fishes/data/my_phy.rdata")
phy_fish = phylo

load("~/my_pc/Global_ED/data/Fishes/data/my_data_used_final.rdata")
load("~/my_pc/Global_ED/data/Fishes/data/Basin042017_3119.rdata")
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
save(mod_fish_sp_exotic_ed_inla,
    file = 'results/primary_results/mod_fish_sp_exotic_ed_inla.rdata')
load('results/primary_results/mod_fish_sp_exotic_ed_inla.rdata')
summary(mod_fish_sp_exotic_ed_inla)


### get predict data for mod_fish_sp_exotic_ed_inla
lincombs.fish.data.estab.exotic_ED_single = data.frame(exotic_ED=seq(min(sp_fish_exo_nat$exotic_ED),
                                                                      max(sp_fish_exo_nat$exotic_ED),
                                                                      length=100))

lincombs.fish.matrix.estab.exotic_ED_single=model.matrix( ~ exotic_ED,
                                                          data=lincombs.fish.data.estab.exotic_ED_single)
lincombs.fish.matrix.estab.exotic_ED_single=as.data.frame(lincombs.fish.matrix.estab.exotic_ED_single)
lincombs.fish.estab.exotic_ED_single=inla.make.lincombs(lincombs.fish.matrix.estab.exotic_ED_single)

inla.model_lincombs.fish.estab.exotic_ED_single = inla(presence ~ exotic_ED
                                                        + f(X1.Basin.Name, model="iid")  
                                                        #+ f(species, exotic_ED, model="iid")
                                                        ,
                                                        data = sp_fish_exo_nat,
                                                        family = "zeroinflatedbinomial1",
                                                        control.compute = list(dic = TRUE,
                                                                               waic = TRUE, 
                                                                               cpo = TRUE,
                                                                               config = TRUE),
                                                        control.predictor = list(compute = TRUE),
                                                        quantiles = c(0.025, 0.5, 0.975),
                                                        lincomb = lincombs.fish.estab.exotic_ED_single)

lincombs.fish.posterior.estab.exotic_ED_single = inla.model_lincombs.fish.estab.exotic_ED_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.fish.estab.exotic_ED_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.fish.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.fish.data.estab.exotic_ED_single$predicted.value=unlist(lapply(lincombs.fish.posterior.estab.exotic_ED_single,
                                                                         function(x)inla.emarginal(fun=plogis,x)))
lincombs.fish.data.estab.exotic_ED_single$lower=unlist(lapply(lincombs.fish.posterior.estab.exotic_ED_single,
                                                               function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.fish.data.estab.exotic_ED_single$upper=unlist(lapply(lincombs.fish.posterior.estab.exotic_ED_single,
                                                               function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

lincombs.fish.data.estab.exotic_ED_single_all = list(summary = inla.model_lincombs.fish.estab.exotic_ED_single$summary.fixed,
                                                  prediction = lincombs.fish.data.estab.exotic_ED_single)
#lincombs.fish.data.estab.exotic_ED_single
save(lincombs.fish.data.estab.exotic_ED_single_all, 
     file = 'results/primary_results/lincombs.fish.data.estab.exotic_ED_single_all.rdata')

