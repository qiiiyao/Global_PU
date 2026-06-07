#Lirong Cai
#Email:lirong.cai18@gmail.com


#We run 1000 null model randomizations using the Scientific Compute Cluster at GWDG, the joint data center of Max Planck Society for the Advancement of Science (MPG) and University of Göttingen.
#The results obtained were used for CANAPE and calculating standardized effect size of relative phylogenetic endemism.

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#set private library path 
.libPaths( c( .libPaths(), "...") )
.libPaths()

#
library(ape)
library(PhyloMeasures)
library(picante)
library(phangorn)
library(dplyr)
library(vegan)
gc()
#load data
load(".../species_comm.RData")#Comm
load(".../Area.RData")#Area
load(".../Tree.RData")#Tree

#trees
Tree.orig<-Tree
Tree.scale<-Tree
Tree.comp<-Tree
Tree.orig2<-Tree
Tree.scale2<-Tree
Tree.comp2<-Tree

# randomize comm
RandComm<-simulate(nullmodel(Comm, "curveball"), nsim = 1, burnin = 100000,thin = 0)[, , 1]

RandComm.df<-as.data.frame(RandComm)
if(!(all(colnames(RandComm.df)==colnames(RandComm))&all(rownames(RandComm.df)==rownames(RandComm)))){
  stop("error in converting matrix to data.frame")
}

#
totallength_edge<- sum(Tree$edge.length)
count_edge<-length(Tree$edge.length)
#
for(i in 1:length(Tree$edge.length))
{

  clade=Tree$tip.label[Descendants(Tree,Tree$edge[i,2],"tip")[[1]]]
  
  if(length(clade) > 1)
  {
    sp_occur<- RandComm.df %>%dplyr::select(all_of(clade)) %>% rowSums()
    cladeArea =sum(Area[which(sp_occur>0)])  
    cladecounts = length(which(sp_occur>0)) 
  }
  if(length(clade) == 1)
  {
    cladeArea = sum(Area[which(RandComm[,which(colnames(RandComm) %in% clade)]>0)])
    cladecounts = sum(RandComm[,which(colnames(RandComm) %in% clade)])
  }
  
  #1.weighted by area
  #a. original tree
  Tree.orig$edge.length[i] = Tree.orig$edge.length[i]/cladeArea
  #b. scaled tree: scale edge lengths as fraction of the total tree length
  Tree.scale$edge.length[i] = Tree.scale$edge.length[i]/totallength_edge/cladeArea #clade(cladeArea) should be non-zero; so taxa which not occur in the communities shoud be dropped from phylogeny 
  #c. comparison tree: retains the actual tree topology but makes all branches of equal length(here, use 1).
  Tree.comp$edge.length[i] = 1/count_edge/cladeArea
  
  #2.weighted by region counts
  #a. original tree
  Tree.orig2$edge.length[i] = Tree.orig2$edge.length[i]/cladecounts
  #b. scaled tree: scale edge lengths as fraction of the total tree length
  Tree.scale2$edge.length[i] = Tree.scale2$edge.length[i]/totallength_edge/cladecounts #clade should be non-zero; so taxa which not occur in the communities shoud be dropped from phylogeny
  #c. comparison tree (also scaled): retains the actual tree topology but makes all branches of equal length.
  Tree.comp2$edge.length[i] = 1/count_edge/cladecounts
}
#calculate metrics
#1.weighted by area
PE.area = pd.query(Tree.orig, RandComm)
PE.numerator.area = pd.query(Tree.scale, RandComm)
PE.denominator.area = pd.query(Tree.comp,RandComm)
RPE.area=PE.numerator.area/PE.denominator.area

#2.weighted by region counts
PE.count = pd.query(Tree.orig2, RandComm)
PE.numerator.count= pd.query(Tree.scale2, RandComm)
PE.denominator.count = pd.query(Tree.comp2,RandComm)
RPE.count=PE.numerator.count/PE.denominator.count

out<- data.frame(entity_ID=row.names(Comm), PE.area=PE.area,PE.numerator.area=PE.numerator.area, 
                 PE.denominator.area=PE.denominator.area, RPE.area=RPE.area,
                 PE.count=PE.count,PE.numerator.count=PE.numerator.count, 
                 PE.denominator.count=PE.denominator.count, RPE.count=RPE.count)

save(out, file=file.path(paste0("..../out",  args[1], ".Rdata")))
write.csv(out,file=file.path(paste0("..../out",  args[1], ".csv")))
rm(list = ls())
gc()

