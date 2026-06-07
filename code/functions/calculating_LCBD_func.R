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

calcu_LCBD_parallel = function(Tree, Comm, Region_posi) {
  
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
  system.time({cbl = cbl.query(Tree, Comm_1, standardize = FALSE)})
  #minimum of total branch length between a pair of communities
  PD = pd.query(Tree, Comm_1)
  PD_min = sapply(PD, function(x) sapply(PD, function(y) min(x,y)))
  #phylogenetic turnover (Simpson dissimilarity index)
  phylo_sim_all_0 = 1-(cbl/PD_min)
  colnames(phylo_sim_all_0) = rownames(Comm_1)
  rownames(phylo_sim_all_0) = rownames(Comm_1)
  
  # Discard rows and columns that are all NA
  non_na_rows = rowSums(!is.na(phylo_sim_all_0)) > 0
  non_na_cols = colSums(!is.na(phylo_sim_all_0)) > 0
  phylo_sim_all = phylo_sim_all_0[non_na_rows, non_na_cols]
  
  phylo_sim_all = as.dist(phylo_sim_all)
  
  #3 Local contribution to phylogenetic beta diversity----------
  dissimilarities=phylo_sim_all
  
  ade4::is.euclid(dissimilarities)
  dissimilarities[dissimilarities <0] = 0
  LCBD_simp=LCBD.comp((dissimilarities), sqrt.D = TRUE, save.D = T)
  summary(LCBD_simp$LCBD)
  
  #if LCBD have negative values, use Lingoes transformation
  if (min(LCBD_simp$LCBD)<0){
    dissimilarities=ade4::lingoes(dissimilarities,print =T,cor.zero =F)
    ade4::is.euclid(dissimilarities)
    LCBD_simp=LCBD.comp(dissimilarities, sqrt.D = F, save.D = T)
  }
  
  summary(LCBD_simp$LCBD)
  #save data
  LCBD_simp_geo = data.frame(RegionID = as.numeric(attr(dissimilarities, "Labels")),
                             LCBD = LCBD_simp$LCBD, 
                             SSI = LCBD_simp$LCBD * LCBD_simp$beta['SStotal'])
  
  out = list(LCBD_simp_geo = LCBD_simp_geo,
             beta_mat = phylo_sim_all_0)
  return(out)
  
}

