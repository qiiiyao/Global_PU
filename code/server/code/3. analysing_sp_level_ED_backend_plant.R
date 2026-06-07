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
  rename(c('exotic_ED' = 'ED'))


mean_ED_natives_plant = cbind(Region_id = sort(unique(df.native.650$Region_id)),
                              ED_plant_native$df)

sp_plant_exo_nat = comm_plant_exotic %>% 
  left_join(mean_ED_natives_plant, by = 'Region_id') %>% 
  rename(c('mean_native_ED' = 'mean_ED')) %>% 
  rename(c('mean_native_EDR' = 'mean_EDR'))

gc()

sp_plant_exo_nat$delta_exo_nat = log(sp_plant_exo_nat$exotic_ED/
                                       sp_plant_exo_nat$mean_native_ED)

sp_plant_exo_nat = sp_plant_exo_nat %>% filter(!is.na(delta_exo_nat))


### get predict data for mod_plant_sp_inv_ed_inla
lincombs.plant.data.estab.delta_e_n_single = data.frame(delta_exo_nat=seq(min(sp_plant_exo_nat$delta_exo_nat),
                                                                          max(sp_plant_exo_nat$delta_exo_nat),
                                                                          length=100))

lincombs.plant.matrix.estab.delta_e_n_single=model.matrix(~delta_exo_nat,
                                                          data=lincombs.plant.data.estab.delta_e_n_single)
lincombs.plant.matrix.estab.delta_e_n_single=as.data.frame(lincombs.plant.matrix.estab.delta_e_n_single)
lincombs.plant.estab.delta_e_n_single=inla.make.lincombs(lincombs.plant.matrix.estab.delta_e_n_single)

inla.model_lincombs.plant.estab.delta_e_n_single = inla(presence ~ delta_exo_nat
                                                        + f(Region_id, model="iid")  
                                                        + f(species, delta_exo_nat, model="iid")
                                                        ,
                                                        data = sp_plant_exo_nat,
                                                        family = "zeroinflatedbinomial1",
                                                        control.compute = list(dic = TRUE,
                                                                               waic = TRUE, 
                                                                               cpo = TRUE,
                                                                               config = TRUE),
                                                        control.predictor = list(compute = TRUE),
                                                        quantiles = c(0.025, 0.5, 0.975),
                                                        lincomb = lincombs.plant.estab.delta_e_n_single)

lincombs.plant.posterior.estab.delta_e_n_single = inla.model_lincombs.plant.estab.delta_e_n_single$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.plant.estab.delta_e_n_single$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.plant.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.plant.data.estab.delta_e_n_single$predicted.value=unlist(lapply(lincombs.plant.posterior.estab.delta_e_n_single,
                                                                         function(x)inla.emarginal(fun=plogis,x)))
lincombs.plant.data.estab.delta_e_n_single$lower=unlist(lapply(lincombs.plant.posterior.estab.delta_e_n_single,
                                                               function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.plant.data.estab.delta_e_n_single$upper=unlist(lapply(lincombs.plant.posterior.estab.delta_e_n_single,
                                                               function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

lincombs.plant.data.estab.delta_e_n_single_ran2 = list(summary = inla.model_lincombs.plant.estab.delta_e_n_single$summary.fixed,
                                                       prediction = lincombs.plant.data.estab.delta_e_n_single)
#lincombs.plant.data.estab.delta_e_n_single
save(lincombs.plant.data.estab.delta_e_n_single_ran2, 
     file = 'results/primary_results/lincombs.plant.data.estab.delta_e_n_single_ran2.rdata')


rm(sp_plant_exo_nat)