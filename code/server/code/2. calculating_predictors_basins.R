### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#0 loading required packages-------------------
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
rm(list = ls())
requirements = c("foreach", "dplyr", 'tidyr', 'sf', 'terra', 'raster', 'scico')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  } 
})

setwd("~/my_pc/Global_ED")

# load the background maps for quantifying the predictors
df = st_read("data/Fishes/data/Basin042017_3119/Basin042017_3119.shp")
df = sf::st_make_valid(df)  


#3. trade flow----
gc()
library(countrycode)
library(stringr)
library(stringi)
library(tibble)
library(data.table)
library(dplyr)

df = st_read("data/Fishes/data/Basin042017_3119/Basin042017_3119.shp")
df = sf::st_make_valid(df)  
trade_0 = read.csv("data/trade_flow/BACI_HS22_Y2022_V202501.csv", header = T,
                   sep = ',')
country_code = read.csv("data/trade_flow/country_codes_V202501.csv", header = T,
                        sep = ',')
#product_code = read.csv("data/trade_flow/product_codes_HS22_V202501.csv", header = T,
#                        sep = ',')
population = terra::rast("data/population/gpw_v4_population_count_rev11_2020_30_min.tif")

str(trade_0)
str(country_code)
str(df)

df = df %>% st_drop_geometry()

df = df %>%
  mutate(iso_code = countrycode(Country, origin = "country.name", destination = "iso3c"))

intersect(sort(unique(df$iso_code)), sort(unique(country_code$country_iso3)))
no_trade_id = setdiff(sort(unique(df$iso_code)), sort(unique(country_code$country_iso3)))


df_country_code = df %>% 
  left_join(country_code, by = join_by('iso_code' == 'country_iso3'))

df_no_country_code = df_country_code %>% filter(is.na(country_code))
unique(df_no_country_code$iso_code)
unique(df_no_country_code$Country)

## 3.1 Join country/region names of basins and trade data----

# modify some ISO_code to make some region in the TYDG match country_code 
df_2 = arrange(df, df$iso_code)
df_2$my_iso_code = df_2$iso_code
df_2$my_iso_code[df_2$Country == "Puerto Rico"] = "PUS"
df_2$my_iso_code[df_2$Country == "Taiwan"] = "CHN"

df_country_code = df_2 %>% left_join(country_code,
                                     by = join_by('my_iso_code' == 'country_iso3')) %>% 
  filter(!is.na(country_code))
table(sort(df_country_code$Region_id))
length(unique(df_country_code$BasinName))/length(df$BasinName)
# 0.991664 of basins has trade data
df_no_country_code = df_2 %>% left_join(country_code,
                                        by = join_by('my_iso_code' == 'country_iso3')) %>% 
  filter(is.na(country_code))
unique(df_no_country_code$Country)

basins_no_trade = sort(unique(df_no_country_code$BasinName))

#st_write(df_no_country_code, 
#         'results/primary_results/predictors_basins/df_no_country_code.shp')
save(basins_no_trade,
     file = 'results/primary_results/predictors_basins/basins_no_trade.rdata')


## 3.2 Join country/region names of basins and population data----
# make projection coincidence
df = st_read("data/Fishes/data/Basin042017_3119/Basin042017_3119.shp")
df_for_population = sf::st_transform(df, crs = sf::st_crs(population))

# convert to SpatVector and do zonal
v_for_population = terra::vect(df_for_population)
# could be clip to high speed
population_crop = terra::crop(population, v_for_population, snap="out")

# zonal: by=Region_id, fun="mean"
z_for_population = cbind(BasinName = v_for_population$BasinName,
                         terra::extract(population_crop,
                                        v_for_population, fun = "sum", na.rm = TRUE))
colnames(z_for_population)[which(colnames(z_for_population) == 'gpw_v4_population_count_rev11_2020_30_min')] = 'population_count_2020'

df_country_code_population = df_country_code %>% left_join(z_for_population[,c('BasinName',
                                                                               'population_count_2020')],
                                                           by = 'BasinName')
rm(population_crop, v_for_population)

## 3.3 Construct trade flow matrix----
#distance for tradeflow among regions
df_country_code_population = df_country_code_population %>% st_drop_geometry()
df_country_code_population$rel_popu = NA

if(sum(is.na(df_country_code_population$population_count_2020)) > 0){
  df_country_code_population[is.na(df_country_code_population$population_count_2020),]$population_count_2020 = 0
}

list_country_code_population = split(df_country_code_population,
                                     df_country_code_population$country_code)
list_country_code_population = lapply(list_country_code_population,
                                      function(x){
                                        #x = list_country_code_population[['76']]
                                        if (nrow(x) == 1){x$rel_popu = 1
                                        } else {
                                          x$rel_popu = x$population_count_2020/sum(x$population_count_2020,
                                                                                   na.rm = T)}
                                        
                                        
                                        return(x)
                                      })
df_country_code_population = data.table::rbindlist(list_country_code_population)

df_country_code_population_direct = df_country_code_population %>% filter(rel_popu == 1)

basin_direct_trade = sort(df_country_code_population_direct$BasinName)
save(basin_direct_trade,
     file = 'results/primary_results/predictors_basins/basin_direct_trade.rdata')

### 3.3.1 Symmetrical matrix for regional ED----
tradeflow_mat_s = matrix(NA, nrow = length(unique(df_country_code_population$BasinName)),
                         ncol = length(unique(df_country_code_population$BasinName)))

rownames(tradeflow_mat_s) = sort(unique(df_country_code_population$BasinName))
colnames(tradeflow_mat_s) = sort(unique(df_country_code_population$BasinName))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_s)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_s)){
    #j = 1
    region_i = rownames(tradeflow_mat_s)[i]
    region_j = rownames(tradeflow_mat_s)[j]
    
    country_code_i = df_country_code_population %>% filter(BasinName == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(BasinName == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
      pull(rel_popu) %>% unique()
    
    flow = (sum(tr[J(country_code_i, country_code_j), v], 
                na.rm = T) + sum(tr[J(country_code_j, country_code_i), v], 
                                 na.rm = T)) * rel_popu_i * rel_popu_j
    tradeflow_mat_s[j,i] = flow
  }
}


tradeflow_mat_s[upper.tri(tradeflow_mat_s,
                          diag=FALSE)] = t(tradeflow_mat_s)[upper.tri(tradeflow_mat_s,
                                                                      diag=FALSE)]
tradeflow_mat_s[1:5,1:5]

tradeflow_mat_s = tradeflow_mat_s[sort(as.character(row.names(tradeflow_mat_s))),
                                  sort(as.character(colnames(tradeflow_mat_s)))]
tradeflow_mat_s[1:10,1:10]
colSums(tradeflow_mat_s, na.rm = T)


save(tradeflow_mat_s,
     file="results/primary_results/distances_basins/tradeflow_mat_s.RDATA")



### 3.3.1 Asymmetrical matrix for sp ED----
tradeflow_mat_a = matrix(NA, nrow = length(unique(df_country_code_population$BasinName)),
                         ncol = length(unique(df_country_code_population$BasinName)))

rownames(tradeflow_mat_a) = sort(unique(df_country_code_population$BasinName))
colnames(tradeflow_mat_a) = sort(unique(df_country_code_population$BasinName))

tr = as.data.table(trade_0)[, .(v = sum(v, na.rm=TRUE)), by=.(i,j)] ##v: Values in thousand USD
setkey(tr, i, j)

for (i in 1:nrow(tradeflow_mat_a)){
  print(i)
  # i = 42
  for (j in i:nrow(tradeflow_mat_a)){
    #j = 1
    region_i = rownames(tradeflow_mat_a)[i]
    region_j = rownames(tradeflow_mat_a)[j]
    
    country_code_i = df_country_code_population %>% filter(BasinName == region_i) %>% 
      pull(country_code) %>% unique()
    country_code_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
      pull(country_code) %>% unique()
    
    rel_popu_i = df_country_code_population %>% filter(BasinName == region_i)%>% 
      pull(rel_popu) %>% unique()
    rel_popu_j = df_country_code_population %>% filter(BasinName == region_j)%>% 
      pull(rel_popu) %>% unique()
    
    flow_ji = sum(tr[J(country_code_i, country_code_j), v], 
                  na.rm = T)  * rel_popu_j
    tradeflow_mat_a[j,i] = flow_ji
    
    flow_ij = sum(tr[J(country_code_j, country_code_i), v], 
                  na.rm = T)  * rel_popu_j
    tradeflow_mat_a[i,j] = flow_ij
  }
}


tradeflow_mat_a[1:5,1:5]

tradeflow_mat_a = tradeflow_mat_a[sort(as.character(row.names(tradeflow_mat_a))),
                                  sort(as.character(colnames(tradeflow_mat_a)))]
tradeflow_mat_a[1:10,1:10]
colSums(tradeflow_mat_a, na.rm = T)

save(tradeflow_mat_a,
     file="results/primary_results/distances_basins/tradeflow_mat_a.RDATA")





