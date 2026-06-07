import fiona
import geopandas as gpd
import pandas as pd
import numpy as np
import os

df = gpd.read_file('D:/Globe taxa distribution/GADM 401/gadm401shp/gadm401.shp')
sp_nonoverlap = gpd.read_file(r'D:/Globe taxa distribution/Birds/results/sp_no_overlap.shp')
## Check crs
sp_nonoverlap.crs
df.crs ## Same coordinate reference systems!

## Clip
sp_no_overlap_1 = pd.DataFrame()
sp_overlap_1 = pd.DataFrame()
for x in range(0,len(sp_nonoverlap)):
   #x = 17413
   sp = gpd.read_file(r'D:/Globe taxa distribution/Birds/results/sp_no_overlap.shp',
    rows=slice(x,(x+1)))
   clipped_shp = sp.overlay(df, how='intersection',
                            keep_geom_type=False)
   if (len(clipped_shp) < 1):
     print(x)
     sp_no_overlap_1 = pd.concat([sp, sp_no_overlap_1])
     #sp.to_file('D:/Globe taxa distribution/Birds/tranformed data/sp_no_overlap_'+ str(x) +'.shp',  
     #              driver='ESRI Shapefile')
   else:
     sp_overlap_1 = pd.concat([clipped_shp, sp_overlap_1])
     #clipped_shp_dat = clipped_shp.iloc[0:len(clipped_shp),0:(len(clipped_shp.columns)-1)]
     #clipped_shp_dat.to_csv('D:/Globe taxa distribution/Birds/tranformed data/sp_overlap_' + str(x) + '.csv')

len(sp_no_overlap_1) ### Still 19 species cannot overlap!!!
