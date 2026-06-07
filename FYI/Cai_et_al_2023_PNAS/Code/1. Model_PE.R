#Lirong Cai
#Email:lirong.cai18@gmail.com

# Model phylogenetic endemism (PE)
# Here we take PE calculated based on range size of species as the area of regions where a species occur (i.e. PE.area.merged), based on the dataset with missing species added to the phylogeny and excluding apomictic taxa as an example

# 0 loading------------------------------------------------
rm(list=ls())
library(dplyr)
#load data
setwd('D:/R projects/Multi_taxa/code/FYI/Codes_and_Data_Cai_et_al_2023')
load("Data/Data.RData")

#load function
source("Code/Function_autocov_sphere.R")


#1 transform data for models-----
summary(Data)
Data<-Data%>% dplyr::select(entity_ID,PE.area.merged,longitude,latitude, Area,
                            Elev , Soildiv , MAT , MAP, TS , PS , 
                            LengthGrow , SLMP, VT_LGM, TempStability_LGM, 
                            TempAnomaly_midPliocene)
#remove regions without complete data
Data<-Data[complete.cases(Data),]#818

#transform response variable
Data$PE.area.merged<-log10(Data$PE.area.merged)
hist(Data$PE.area.merged,xlab = "log10(PE.area)")

#log-transform predictors having skewed distributions
Data$Area <- log10(Data$Area)
Data$Elev <- log10(Data$Elev)
Data$Soildiv<-log10(Data$Soildiv+30)
Data$MAP<-log10(Data$MAP+700)
Data$TS<-log10(Data$TS+10)
Data$PS<-log10(Data$PS)
Data$VT_LGM <- log10(Data$VT_LGM)
Data$TempStability_LGM <- log10(Data$TempStability_LGM+0.1)
Data$TempAnomaly_midPliocene<- log10(Data$TempAnomaly_midPliocene+1)
#center and scale continuous predictors
Data[,c(5:(ncol(Data)))] <- scale(Data[,c(5:(ncol(Data)))])

# coordinates for RAC
coords <- cbind(Data$longitude, Data$latitude)


#2 Non-spatial model----
#define formula
formulaString <- PE.area.merged ~ Area + Elev + Soildiv + MAT + MAP + TS + PS + 
  LengthGrow + SLMP + VT_LGM + TempStability_LGM + 
  TempAnomaly_midPliocene+Area:SLMP + Elev:SLMP + Soildiv:SLMP + 
  MAT:SLMP + MAP:SLMP + TS:SLMP + PS:SLMP + 
  LengthGrow:SLMP +VT_LGM:SLMP + TempStability_LGM:SLMP + 
  TempAnomaly_midPliocene:SLMP

LM.m<- lm(formula=formulaString,data=Data)
# plot correlogram
correlogram <- ncf::correlog(Data$longitude, Data$latitude,resid(LM.m),increment=200, resamp = 300, latlon=TRUE)
plot(correlogram)

#3 Spatial model (RAC)----
# define formula of RAC
formula_RAC <- PE.area.merged ~ Area + Elev + Soildiv + MAT + MAP + TS + PS + 
  LengthGrow + SLMP + VT_LGM + TempStability_LGM + 
  TempAnomaly_midPliocene+Area:SLMP + Elev:SLMP + Soildiv:SLMP + 
  MAT:SLMP + MAP:SLMP + TS:SLMP + PS:SLMP + 
  LengthGrow:SLMP +VT_LGM:SLMP + TempStability_LGM:SLMP + 
  TempAnomaly_midPliocene:SLMP+ rac

rac<- autocov_sphere(z=resid(LM.m), xy=coords, type = "one", zero.policy = TRUE, style = "W", longlat=T) 
RAC.PEA.rmapom <- lm(formula_RAC , Data)
summary(RAC.PEA.rmapom)

# plotting a correlogram
correlogram_RAC <- ncf::correlog(Data$longitude, Data$latitude,resid(RAC.PEA.rmapom),increment=200, resamp = 300, latlon=TRUE)
plot(correlogram_RAC)

#save(RAC.PEA.rmapom,file="data/results/Models/RAC_PEA_rmapom.Rdata")
