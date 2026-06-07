#Lirong Cai
#Email:lirong.cai18@gmail.com

#In this script, we modeled geographic patterns of phylogenetic and taxonomic turnover

#0 loading----------------------------
rm(list = ls())
library(reshape2)
library(plyr)
library(betareg)
library(MASS)
#library(gdm)
library(vegan)
library(rgeos)
library(sjPlot)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

setwd("D:/R projects/Global_ED")

colors_4d = c('#91989f', '#efbb24', '#1b813e', '#398fb7')
colors_10d = c('#91989f',
  colorRampPalette(c('#FFF9B0', '#FFE066', '#FFC300', '#FF9E00', '#FF7A00', '#E85C00'))(6),
  colorRampPalette(c('#9BD4E4', '#398FB7', '#004A78'))(3))

mat_to_long = function(mat, value_name) {
  stopifnot(is.matrix(mat))
  
  df = as.data.frame(as.table(mat))
  colnames(df) = c("region_i", "region_j", value_name)
  
  df
}


#1. modelling for Mammals and Birds: TDWG4----------------------------
# load trade flow value matrix
load('results/primary_results/distances_TDWG/tradeflow_mat.RDATA')
load("results/primary_results/predictors_TDWG/RegionID_direct_trade.rdata")

#load individual climate distances
load('results/primary_results/distances_TDWG/clim_dist.rdata')

#load beta diversity
load('results/primary_results/distances_beta/LCBD_mammal_all.rdata')
load('results/primary_results/distances_beta/LCBD_mammal_native.rdata')
load('results/primary_results/distances_beta/LCBD_bird_extant.rdata')
load('results/primary_results/distances_beta/LCBD_bird_native.rdata')
turnover_mammal_extant_mat = LCBD_mammal_all$beta_mat
turnover_mammal_native_mat = LCBD_mammal_native$beta_mat
turnover_mammal_delta_mat = log((turnover_mammal_extant_mat+0.001) /
                                  (turnover_mammal_native_mat+0.001))

turnover_bird_extant_mat = LCBD_bird_extant$beta_mat
turnover_bird_native_mat = LCBD_bird_native$beta_mat
turnover_bird_delta_mat = log((turnover_bird_extant_mat+0.001) /
                              (turnover_bird_native_mat+0.001))

dim(turnover_mammal_delta_mat)
dim(turnover_bird_delta_mat)

clim_mat = as.matrix(clim_dist)
df_trade = mat_to_long(tradeflow_mat, "trade_flow")
df_clim  = mat_to_long(clim_mat, "clim_dist")

df_mammal = mat_to_long(turnover_mammal_delta_mat, "delta_turnover_mammal")
df_bird   = mat_to_long(turnover_bird_delta_mat,   "delta_turnover_bird")

df_trade_all = df_trade %>%
  left_join(df_clim,  by = c("region_i", "region_j")) %>%
  left_join(df_mammal, by = c("region_i", "region_j")) %>%
  left_join(df_bird,   by = c("region_i", "region_j"))

df_trade_all = df_trade_all %>%
  filter(region_i != region_j)

df_trade_all = df_trade_all %>%
  mutate(
    region_i = as.numeric(region_i),
    region_j = as.numeric(region_j),
    pair_id = paste(
      pmin(region_i, region_j),
      pmax(region_i, region_j),
      sep = "_"
    )
  ) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  dplyr::select(-pair_id)

hist(df_trade_all$delta_turnover_mammal)
hist(df_trade_all$trade_flow)
hist(df_trade_all$clim_dist)
str(df_trade_all)

##1.1 All regions with direct and relative population * trade data----
df_trade_all = df_trade_all %>%
  dplyr::filter(
    !is.na(delta_turnover_mammal),
    !is.na(trade_flow),
    !is.na(clim_dist),
    is.finite(delta_turnover_mammal),
    is.finite(trade_flow),
    is.finite(clim_dist)
  )
df_trade_direct = df_trade_all %>% filter(region_i %in% RegionID_direct_trade & 
                                            region_j %in% RegionID_direct_trade)

df_trade_all$modi_trade_flow = scale(log(df_trade_all$trade_flow+0.001))
df_trade_all$modi_clim_dist = scale(log(df_trade_all$clim_dist+0.001))
hist(df_trade_all$modi_trade_flow)
hist(df_trade_all$modi_clim_dist)

lm_trade_clim_mammal = lm(delta_turnover_mammal ~ modi_trade_flow + modi_clim_dist,
                  data = df_trade_all)
summary(lm_trade_clim_mammal)

lm_trade_clim_bird = lm(delta_turnover_bird ~ modi_trade_flow + modi_clim_dist,
                          data = df_trade_all)
summary(lm_trade_clim_bird)

##1.2 if just retain regions with direct trade data----
df_trade_direct$modi_trade_flow = as.numeric(scale(log10(df_trade_direct$trade_flow+0.001)))
df_trade_direct$modi_clim_dist = as.numeric(scale(log10(df_trade_direct$clim_dist+0.001)))

table(df_trade_direct$trade_flow)
range(df_trade_direct$clim_dist)
hist(df_trade_direct$modi_trade_flow)
hist(df_trade_direct$modi_clim_dist)
hist(df_trade_direct$modi_trade_flow)
hist(df_trade_direct$modi_clim_dist)
str(df_trade_direct)

lm_trade_clim_mammal_direct = lm(delta_turnover_mammal ~ modi_trade_flow + modi_clim_dist,
                          data = df_trade_direct)
summary(lm_trade_clim_mammal_direct)
sjPlot::plot_model(lm_trade_clim_mammal_direct, 
           type = 'pred')

lm_trade_clim_bird_direct = lm(delta_turnover_bird ~ modi_trade_flow + modi_clim_dist,
                        data = df_trade_direct)
summary(lm_trade_clim_bird_direct)
plot_model(lm_trade_clim_bird_direct, 
           type = 'pred')


plot_trade_clim_mammal_direct = plot_model(lm_trade_clim_mammal_direct, type = 'est',
                                          terms = c('modi_trade_flow',
                                                    'modi_clim_dist'),
                                          show.p = T) + 
  geom_hline(yintercept = 0, linetype = 2, color = 'black')+
  scale_fill_manual(values=c(colors_4d[4], colors_4d[4]))+
  scale_color_manual(values=c(colors_4d[4], colors_4d[4]))+
  scale_y_continuous(limits = c(-0.05, 0.003)) + 
  scale_x_discrete(labels = c('Climate distance',
                              'Trade flow'))+
  labs(x = '', y = 'Coefficient estimate')+
  ggtitle('Mammals')+
  theme(#plot.title = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.25, 0.8),
    plot.title = element_text(hjust = 0.5))

pred_trade_clim_mammal_direct = get_model_data(lm_trade_clim_mammal_direct, type = 'est',
                                           terms = c('modi_trade_flow',
                                                     'modi_clim_dist')) 


plot_trade_clim_bird_direct = plot_model(lm_trade_clim_bird_direct, type = 'est',
                                           terms = c('modi_trade_flow',
                                                     'modi_clim_dist'),
                                           show.p = T) + 
  geom_hline(yintercept = 0, linetype = 2, color = 'black')+
  scale_fill_manual(values=c(colors_4d[4], colors_4d[4]))+
  scale_color_manual(values=c(colors_4d[4], colors_4d[4]))+
  scale_y_continuous(limits = c(-0.05, 0.003)) + 
  scale_x_discrete(labels = c('Climate distance',
                              'Trade flow'))+
  labs(x = '', y = 'Coefficient estimate')+
  ggtitle('Birds')+
  theme(#plot.title = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.25, 0.8),
    plot.title = element_text(hjust = 0.5))

pred_trade_clim_bird_direct = get_model_data(lm_trade_clim_bird_direct, type = 'est',
                                               terms = c('modi_trade_flow',
                                                         'modi_clim_dist')) 

#2. modelling for Plants: GIFT----------------------------
# load trade flow value matrix
load('results/primary_results/distances_GIFT/tradeflow_mat.RDATA')
load("results/primary_results/predictors_GIFT/Region_id_direct_trade.rdata")

#load individual climate distances
load('results/primary_results/distances_GIFT/clim_dist.rdata')

#load beta diversity
load('results/primary_results/distances_beta/LCBD_plant_extant.rdata')
load('results/primary_results/distances_beta/LCBD_plant_native.rdata')
turnover_plant_extant_mat = LCBD_plant_extant$beta_mat
turnover_plant_native_mat = LCBD_plant_native$beta_mat
turnover_plant_delta_mat = log((turnover_plant_extant_mat+0.001) /
                                  (turnover_plant_native_mat+0.001))

dim(turnover_plant_delta_mat)


clim_mat = as.matrix(clim_dist)
df_trade = mat_to_long(tradeflow_mat, "trade_flow")
df_clim  = mat_to_long(clim_mat, "clim_dist")

df_plant = mat_to_long(turnover_plant_delta_mat, "delta_turnover_plant")

df_trade_all = df_trade %>%
  left_join(df_clim,  by = c("region_i", "region_j")) %>%
  left_join(df_plant, by = c("region_i", "region_j")) 

df_trade_all = df_trade_all %>%
  filter(region_i != region_j)

df_trade_all = df_trade_all %>%
  mutate(
    region_i = as.numeric(region_i),
    region_j = as.numeric(region_j),
    pair_id = paste(
      pmin(region_i, region_j),
      pmax(region_i, region_j),
      sep = "_"
    )
  ) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  dplyr::select(-pair_id)

hist(df_trade_all$delta_turnover_plant)
hist(df_trade_all$trade_flow)
hist(df_trade_all$clim_dist)
str(df_trade_all)

##2.1 All regions with direct and relative population * trade data----
df_trade_all = df_trade_all %>%
  dplyr::filter(
    !is.na(delta_turnover_plant),
    !is.na(trade_flow),
    !is.na(clim_dist),
    is.finite(delta_turnover_plant),
    is.finite(trade_flow),
    is.finite(clim_dist)
  )
df_trade_direct = df_trade_all %>% filter(region_i %in% Region_id_direct_trade & 
                                            region_j %in% Region_id_direct_trade)

df_trade_all$modi_trade_flow = scale(log(df_trade_all$trade_flow+0.001))
df_trade_all$modi_clim_dist = scale(log(df_trade_all$clim_dist+0.001))
hist(df_trade_all$modi_trade_flow)
hist(df_trade_all$modi_clim_dist)

lm_trade_clim_plant = lm(delta_turnover_plant ~ modi_trade_flow + modi_clim_dist,
                          data = df_trade_all)
summary(lm_trade_clim_plant)


##2.2 if just retain regions with direct trade data----
df_trade_direct$modi_trade_flow = as.numeric(scale(log10(df_trade_direct$trade_flow+0.001)))
df_trade_direct$modi_clim_dist = as.numeric(scale(log10(df_trade_direct$clim_dist+0.001)))

table(df_trade_direct$trade_flow)
range(df_trade_direct$clim_dist)
hist(df_trade_direct$modi_trade_flow)
hist(df_trade_direct$modi_clim_dist)
hist(df_trade_direct$modi_trade_flow)
hist(df_trade_direct$modi_clim_dist)
str(df_trade_direct)

lm_trade_clim_plant_direct = lm(delta_turnover_plant ~ modi_trade_flow + modi_clim_dist ,
                                 data = df_trade_direct)
summary(lm_trade_clim_plant_direct)
sjPlot::plot_model(lm_trade_clim_plant_direct, 
                   type = 'pred')

plot_trade_clim_plant_direct = plot_model(lm_trade_clim_plant_direct, type = 'est',
           terms = c('modi_trade_flow',
                     'modi_clim_dist'),
           show.p = T) + 
  geom_hline(yintercept = 0, linetype = 2, color = 'black')+
  scale_fill_manual(values=c(colors_4d[4], colors_4d[4]))+
  scale_color_manual(values=c(colors_4d[4], colors_4d[4]))+
  scale_y_continuous(limits = c(-0.05, 0.003)) + 
  scale_x_discrete(labels = c('Climate distance',
                              'Trade flow'))+
  labs(x = '', y = 'Coefficient estimate')+
  ggtitle('Plants')+
  theme(#plot.title = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.25, 0.8),
    plot.title = element_text(hjust = 0.5))

pred_trade_clim_plant_direct = get_model_data(lm_trade_clim_plant_direct, type = 'est',
                                             terms = c('modi_trade_flow',
                                                       'modi_clim_dist')) 

#3. Exporting all taxon column figures----------------------------
pred_trade_clim_plant_direct$taxa = 'Plants'
pred_trade_clim_mammal_direct$taxa = 'Mammals'
pred_trade_clim_bird_direct$taxa = 'Birds'

pred_trade_clim_all_direct = rbind(pred_trade_clim_plant_direct, 
                                   pred_trade_clim_mammal_direct,
                                   pred_trade_clim_bird_direct) %>%
                           as_tibble()

pred_trade_clim_all_direct$taxa = factor(pred_trade_clim_all_direct$taxa,
                                      levels = c("Plants",
                                                 "Mammals",
                                                 "Birds"))

pred_trade_clim_all_direct$term = factor(pred_trade_clim_all_direct$term,
                                         levels = c("modi_clim_dist",
                                                    "modi_trade_flow"))
pred_trade_clim_all_direct$sig = ifelse(pred_trade_clim_all_direct$p.value < 0.05, 1, 0)
pred_trade_clim_all_direct$sig = factor(pred_trade_clim_all_direct$sig,
                          levels = c(1,
                                     0))

ggthemr::ggthemr(palette = "fresh", layout = "clean")
plot_trade_clim_all_direct = ggplot()+
  facet_wrap(vars(taxa), scales = "free")+
  geom_point(mapping = aes(y = term, x = estimate, shape = sig),
             data = pred_trade_clim_all_direct,
             size = 2.5,
             color = 'black')+
  geom_linerange(mapping = aes(y = term, xmin = conf.low, xmax = conf.high),
                 data = pred_trade_clim_all_direct) + 
  geom_vline(xintercept = 0, linetype = 2, color = colors_4d[1])+
  scale_shape_manual(values = c(16, 1))+
  scale_x_continuous(limits = c(-0.04, 0.01), 
                     breaks = c(-0.04, -0.02, 0)) + 
  scale_y_discrete(labels = c('Climate distance',
                              'Trade flow'))+
  labs(x = 'Coefficient estimate', y = '')+
  theme(#plot.title = element_blank(),
    legend.position = 'None',
    legend.position.inside = c(0.25, 0.8),
    plot.title = element_text(hjust = 0.5))



#4. Export figures for paper----
library(devEMF)
library(cowplot)

emf('figures/plot_trade_clim_all_direct.emf',
    width = 25*0.8, height = 10*0.8, coordDPI = 600, 
    units = 'cm',
    emfPlusFontToPath=TRUE # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
)
plot_trade_clim_all_direct
dev.off() #turn off device and finalize file

