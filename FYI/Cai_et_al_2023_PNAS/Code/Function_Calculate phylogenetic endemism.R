#Lirong Cai
#Email:lirong.cai18@gmail.com

library(PhyloMeasures)
library(phangorn)
library(foreach)
library(dplyr)

#Tree: phylogenetic tree
#Comm: species presence/absence matrix with rows as sampling regions and columns as species
#Area: a vector representing area size of each sampling region  
 
#calculate PE and RPE weighted by region area-------------------------------------------------

calcu_PE_RPE_parallel<-function(Tree, Comm, Area){
  
  parallel::detectCores()
  n.cores <- parallel::detectCores()-1
  my.cluster <- parallel::makeCluster(n.cores,  type = "PSOCK",outfile="") 
  doParallel::registerDoParallel(cl = my.cluster)
  Tree.orig<-Tree
  Tree.scale<-Tree
  Tree.comp<-Tree
  Comm.df<-as.data.frame(Comm)
  if(!(all(colnames(Comm.df)==colnames(Comm))&all(rownames(Comm.df)==rownames(Comm)))){
    stop("error in converting matrix to data.frame")
  }
  
  totallength_edge<- sum(Tree$edge.length)
  count_edge<-length(Tree$edge.length)
  
  
    cladeArea_table <- foreach (i=1:length(Tree$edge.length),.combine='c',.packages = c("ape","phangorn","dplyr"),.inorder=T) %dopar%{
      
      clade=Tree$tip.label[Descendants(Tree,Tree$edge[i,2],"tip")[[1]]]
      if(length(clade) > 1)
      {
        #cladeArea = sum(Area[which(apply(Comm[,which(colnames(Comm) %in% clade)],1,sum)>0)])
        sp_occur<- Comm.df %>%select(all_of(clade)) %>% rowSums()
        cladeArea =sum(Area[which(sp_occur>0)])  
        
      }
      if(length(clade) == 1)
      {
        cladeArea = sum(Area[which(Comm[,which(colnames(Comm) %in% clade)]>0)])
      }
      return(cladeArea)
      
    }

  
  # Stop cluster
  parallel::stopCluster(cl = my.cluster)
  closeAllConnections()
  gc()
  
  #a. original tree
  Tree.orig$edge.length = Tree.orig$edge.length/cladeArea_table
  #b. scaled tree:scale edge lengths as fraction of the total tree length
  Tree.scale$edge.length = Tree.scale$edge.length/totallength_edge/cladeArea_table 
  #c. The comparison tree retains the actual tree topology but makes all branches of equal length.
  Tree.comp$edge.length= rep(1, count_edge)/count_edge/cladeArea_table

  
  PE = pd.query(Tree.orig, Comm)
  PE.numerator = pd.query(Tree.scale, Comm)
  PE.denominator = pd.query(Tree.comp,Comm)
  RPE=PE.numerator/PE.denominator
  out<- data.frame(PE=PE,PE.numerator=PE.numerator, PE.denominator=PE.denominator, RPE=RPE)
  return(out)
  
}



#calculate PE and RPE weighted by region counts -------------------------------------------------
calcu_PEcount_parallel<-function(Tree, Comm){
  
  parallel::detectCores()
  n.cores <- parallel::detectCores()-1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK",outfile="") 
  doParallel::registerDoParallel(cl = my.cluster)

  Tree.orig<-Tree
  Tree.scale<-Tree
  Tree.comp<-Tree
  Comm.df<-as.data.frame(Comm)
  if(!(all(colnames(Comm.df)==colnames(Comm))&all(rownames(Comm.df)==rownames(Comm)))){
    stop("error in converting matrix to data.frame")
  }
  
  totallength_edge<- sum(Tree$edge.length)
  count_edge<-length(Tree$edge.length)
  
  
  cladecounts_table <- foreach (i=1:length(Tree$edge.length),.combine='c',.packages = c("ape","phangorn","dplyr"),.inorder=T) %dopar%{
    
    clade=Tree$tip.label[Descendants(Tree,Tree$edge[i,2],"tip")[[1]]]
    if(length(clade) > 1)
    {
      sp_occur<- Comm.df %>%select(all_of(clade)) %>% rowSums()
      cladecounts = length(which(sp_occur>0)) #count regions which species occurred
    }
    if(length(clade) == 1)
    {
      cladecounts = sum(Comm[,which(colnames(Comm) %in% clade)])
    }
    
    return(cladecounts)
    
  }


  
  # Stop cluster
  parallel::stopCluster(cl = my.cluster)
  closeAllConnections()
  gc()
  
  #a. original tree
  Tree.orig$edge.length = Tree.orig$edge.length/cladecounts_table
  #b. scaled tree:scale edge lengths as fraction of the total tree length
  Tree.scale$edge.length = Tree.scale$edge.length/totallength_edge/cladecounts_table 
  #c. The comparison tree retains the actual tree topology but makes all branches of equal length.
  Tree.comp$edge.length= rep(1, count_edge)/count_edge/cladecounts_table
  
  
  PE = pd.query(Tree.orig, Comm)
  PE.numerator = pd.query(Tree.scale, Comm)
  PE.denominator = pd.query(Tree.comp,Comm)
  RPE=PE.numerator/PE.denominator
  out<- data.frame(PE=PE,PE.numerator=PE.numerator, PE.denominator=PE.denominator, RPE=RPE)
  return(out)
  
}

