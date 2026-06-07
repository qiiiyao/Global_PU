import fiona
import geopandas as gpd
import pandas as pd
import numpy as np
import os

df = gpd.read_file('D:/Globe taxa distribution/TDWG4/TDWG4_newTibet.shp')

# Get all the layers from the .gdb file 
layers = fiona.listlayers('D:/Globe taxa distribution/Birds/DistributionData/BOTW.gdb')

gdf_sp_name = gpd.read_file('D:/Globe taxa distribution/Birds/DistributionData/BOTW.gdb',
    layer=layers[1])

## Check crs
sp.crs
df.crs ## Same coordinate reference systems!

## clip
#sp_no_overlap_1 = pd.DataFrame()
#sp_overlap_1 = pd.DataFrame()
for x in range(0,17414):
   #x = 17413
   sp = gpd.read_file(r'D:/Globe taxa distribution/Birds/DistributionData/BOTW.gdb',
    layer=layers[0], rows=slice(x,(x+1)))
   clipped_shp = sp.overlay(df, how='intersection',
                            keep_geom_type=False)
   if (len(clipped_shp) < 1):
     print(x)
     #sp_no_overlap_1 = pd.concat([sp, sp_no_overlap_1])
     sp.to_file('D:/Globe taxa distribution/Birds/tranformed data/sp_no_overlap_'+ str(x) +'.shp',  
                   driver='ESRI Shapefile')
   else:
     #sp_overlap_1 = pd.concat([clipped_shp, sp_overlap_1])
     clipped_shp_dat = clipped_shp.iloc[0:len(clipped_shp),0:(len(clipped_shp.columns)-1)]
     clipped_shp_dat.to_csv('D:/Globe taxa distribution/Birds/tranformed data/sp_overlap_' + str(x) + '.csv')


## Merge non-overlap shp files
file_path = 'D:/Globe taxa distribution/Birds/tranformed data/'  # file path
lst = os.listdir(file_path)
print(lst)

lst_1 = [x for x in lst if x.endswith('.shp')]
print(lst_1)
len(lst_1)
sp_no_overlap = pd.DataFrame()
#sp_overlap = pd.DataFrame()
for x in lst_1:
   #x = lst_1[1]
   sp = gpd.read_file("D:/Globe taxa distribution/Birds/tranformed data/"+x)
   sp_no_overlap = pd.concat([sp, sp_no_overlap])

len(np.unique(sp_no_overlap.loc[:,'sci_name']))
sp_no_overlap.to_file('D:/Globe taxa distribution/Birds/results/sp_no_overlap.shp',
                   driver='ESRI Shapefile')

### overlap species
lst_2 = [x for x in lst if x.endswith('.csv')]
print(lst_2)
len(lst_2)

sp_overlap = pd.DataFrame()
for x in lst_2:
   #x = lst_2[1]
   sp = gpd.read_file("D:/Globe taxa distribution/Birds/tranformed data/"+x)
   sp_overlap = pd.concat([sp, sp_overlap])
len(sp_overlap)
len(np.unique(sp_overlap.loc[:,'sci_name']))
sp_overlap.to_csv('D:/Globe taxa distribution/Birds/results/sp_overlap.csv')

### All species name in overlap and non-overlap
len(np.unique(np.concatenate((sp_overlap.loc[:,'sci_name'].to_numpy(),
                         sp_no_overlap.loc[:,'sci_name'].to_numpy()),
                          axis=None)))



