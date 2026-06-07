### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
#rm(list = ls())
# Package management
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
requirements = c("PhyloMeasures", "phangorn", "parallel", "dplyr",
                 "ape", 'phyloregion', 'adespatial')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})

setwd("D:/R projects/Global_ED")

#Load phylogenetic Tree of bird families from package "ape"
#data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
#comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
#for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
#colnames(comm) = bird.families$tip.label

#Use query function to calculate pd values for each community
#pd.query(bird.families,comm)

#Tree: phylogenetic Tree
#Comm: species presence/absence matrix with rows as sampling regions and columns as species
#Area: a vector representing area size of each sampling region  

#calculate PE and RPE weighted by region area-------------------------------------------------
#Tree = spec_phy.3
#Comm = comm_mammal_exotic
#Region_posi = which(colnames(comm_mammal_exotic) == 'RegionID')

#phyloregion::evol_distinct(x, type="fair.proportion")

calcu_turnover_parallel = function(Tree, Comm, Region_posi) {
  
  parallel::detectCores()
  n.cores = parallel::detectCores()-1
  my.cluster = parallel::makeCluster(n.cores, type = "PSOCK", outfile="") 
  doParallel::registerDoParallel(cl = my.cluster)
  Tree.orig = Tree
  Tree.scale = Tree
  Tree.comp = Tree
  Comm.df = as.data.frame(Comm)
  if(!(all(colnames(Comm.df)==colnames(Comm))&all(rownames(Comm.df)==rownames(Comm)))){
    stop("error in converting matrix to data.frame")
  }
  
  
  Comm_1 = as.matrix(Comm[,which(colnames(Comm) %in% Tree$tip.label)])
  row.names(Comm_1) = unlist(Comm[,Region_posi])
  
  all(Tree$tip.label %in% colnames(Comm))
  
  #common branch length
  cbl = cbl.query(Tree, Comm_1, standardize = FALSE)
  
  if (sum(cbl[upper.tri(cbl)] == 0) > 0){
    replaces = cbind(which(cbl == 0, arr.ind = T), NA)
    replaces[,3] = replaces[,1] - replaces[,2]
    replaces = as_tibble(replaces)
    replaces = replaces %>% filter(V3 != 0)
    
    #minimum of total branch length between a pair of communities
    PD = pd.query(Tree, Comm_1)
    PD_min = sapply(PD, function(x) sapply(PD, function(y) min(x,y)))
    #phylogenetic turnover (Simpson dissimilarity index)
    phylo_sim_all_0 = 1-(cbl/PD_min)
    # PD_min = common branches (a) + minimum unique branch of one region (min(b,c))
    
    phylo_sim_all_0[as.matrix(replaces[,c(1,2)])] = 1-0
    
  } else {
    #minimum of total branch length between a pair of communities
    PD = pd.query(Tree, Comm_1)
    PD_min = sapply(PD, function(x) sapply(PD, function(y) min(x,y)))
    #phylogenetic turnover (Simpson dissimilarity index)
    phylo_sim_all_0 = 1-(cbl/PD_min)
    # PD_min = common branches (a) + minimum unique branch of one region (min(b,c))
    
  }
  
  colnames(phylo_sim_all_0) = rownames(Comm_1)
  rownames(phylo_sim_all_0) = rownames(Comm_1)
  diag(phylo_sim_all_0) = 0
  
  
  out = phylo_sim_all_0
  return(out)
  
}





calcu_turnover_simple = function(Tree, Comm) {
  
  Comm_1 = as.matrix(Comm[,which(colnames(Comm) %in% Tree$tip.label)])
  #common branch length
  cbl = cbl.query(Tree, Comm_1, standardize = FALSE)
  
  if (sum(cbl[upper.tri(cbl)] == 0) > 0){
    replaces = cbind(which(cbl == 0, arr.ind = T), NA)
    replaces[,3] = replaces[,1] - replaces[,2]
    replaces = as_tibble(replaces)
    replaces = replaces %>% filter(V3 != 0)
    
    #minimum of total branch length between a pair of communities
    PD = pd.query(Tree, Comm_1)
    PD_min = sapply(PD, function(x) sapply(PD, function(y) min(x,y)))
    #phylogenetic turnover (Simpson dissimilarity index)
    phylo_sim_all_0 = 1-(cbl/PD_min)
    # PD_min = common branches (a) + minimum unique branch of one region (min(b,c))
    
    phylo_sim_all_0[as.matrix(replaces[,c(1,2)])] = 1-0
    
  } else {
    #minimum of total branch length between a pair of communities
    PD = pd.query(Tree, Comm_1)
    PD_min = sapply(PD, function(x) sapply(PD, function(y) min(x,y)))
    #phylogenetic turnover (Simpson dissimilarity index)
    phylo_sim_all_0 = 1-(cbl/PD_min)
    # PD_min = common branches (a) + minimum unique branch of one region (min(b,c))
    
  }
  
  colnames(phylo_sim_all_0) = rownames(Comm_1)
  rownames(phylo_sim_all_0) = rownames(Comm_1)
  diag(phylo_sim_all_0) = 0
  
  out = phylo_sim_all_0
  return(out)
  
}




