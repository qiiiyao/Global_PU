#Lirong Cai
#Email:lirong.cai18@gmail.com

# Model standardized effect size of relative phylogenetic endemism (RPE.ses)
# Here we take RPE.ses calculated based on range size of species as the area of regions where a species occur (RPE.ses.area.merged), based on the dataset with missing species added to the phylogeny and excluding apomictic taxa as an example

#0 loading------------------------------------------------------------------------------
rm(list = ls())
library(dotwhisker)
library(dplyr)
library(ggpubr)
library(spdep)
library(dbscan)

#load function
source(".../Code/Function_autocov_sphere.R")
#load data
load(".../Data/Data.RData")

Data<- Data %>% dplyr::select(entity_ID,geo_entity,RPE.ses.area.merged,type.area.merged,Geo_type,
                                             latitude,longitude,Elev,
                                             VT_LGM,TempStability_LGM,TempAnomaly_midPliocene)

#1. multiple mean comparisons----------------------------------
Data$type.area.merged <- ordered(Data$type.area.merged,  levels = c("super-endemic","mixed","paleo","neo","notsignificant"))

##1.1 Kruskal-Wallis test----
kruskal.test(Elev~ type.area.merged, data = Data) 
kruskal.test(TempStability_LGM~ type.area.merged, data = Data) 
kruskal.test(VT_LGM~ type.area.merged, data = Data)
kruskal.test(TempAnomaly_midPliocene~ type.area.merged, data = Data)
##1.2 Multiple pairwise-comparison between groups----
pairwise.wilcox.test(Data$Elev, Data$type.area.merged,alternative="two.sided",p.adjust.method = "holm")
pairwise.wilcox.test(Data$TempStability_LGM, Data$type.area.merged,alternative="two.sided",p.adjust.method = "holm")
pairwise.wilcox.test(Data$VT_LGM, Data$type.area.merged,alternative="two.sided",p.adjust.method = "holm")
pairwise.wilcox.test(Data$TempAnomaly_midPliocene, Data$type.area.merged,alternative="two.sided",p.adjust.method = "holm")


#2 Non-spatial model------------------
##2.1 data for modelling----
Data.lm<-Data%>%select(entity_ID,RPE.ses.area.merged,type.area.merged,Geo_type,
                        latitude,longitude,Elev,TempStability_LGM)%>%filter(type.area.merged!="notsignificant")#141
#remove Heterogenous islands
table(Data.lm$Geo_type)
Data.lm<-Data.lm[which(Data.lm$Geo_type!="Heterogenous origins"),]
Data.lm$Geo_type<-as.factor(Data.lm$Geo_type)
Data.lm$Geo_type <- relevel(Data.lm$Geo_type,ref = "Mainland")

hist(Data.lm$RPE.ses,xlab = "RPE.ses",cex.lab=2, cex.axis=2)#138

#log-transform predictors having skewed distributions
Data.lm$Elev <- log(Data.lm$Elev)
Data.lm$TempStability_LGM <- log(Data.lm$TempStability_LGM)
summary(Data.lm)

#center and scale continuous predictors
Data.lm[,c(7:ncol(Data.lm))] <- scale(Data.lm[,c(7:ncol(Data.lm))])

##2.2 fit models----
#get formula
formulaString <- RPE.ses.area.merged~Geo_type+Elev+TempStability_LGM

LM.RPEses.area<- lm(formula=formulaString,data=Data.lm)
summary(LM.RPEses.area)

##2.3  plot correlogram----
cor.LM<- ncf::correlog(x=Data.lm$longitude,y=Data.lm$latitude,  
                                 z=resid(LM.RPEses.area),
                                 latlon=TRUE, resamp=300, increment=200)
plot(cor.LM)


#3 Spatial models (RAC)-----------------------------

##3.1 fit RAC model----
#get formula
formula_RAC <-RPE.ses.area.merged~Geo_type+Elev+TempStability_LGM+rac
#create autovariate
rac<- autocov_sphere(z=resid(LM.RPEses.area), xy=cbind(Data.lm$longitude, Data.lm$latitude), 
                     type = "one", zero.policy = TRUE, style = "B", longlat=T) 
RAC.RPEses.area <- lm(formula_RAC , Data.lm)
summary(RAC.RPEses.area)
#save(RAC.RPEses.area,file="data/results/Models/RAC_RPEarea_rmapom.Rdata")

##3.2 plot correlogram----
cor.RAC<- ncf::correlog(x=Data.lm$longitude,y=Data.lm$latitude, 
                                 z=resid(RAC.RPEses.area),
                                 latlon=TRUE, resamp=300, increment=200)
plot(cor.RAC)