### Set up R environments 
### Part of the code adapted from Cai_et_al_2023_PNAS
rm(list = ls())
# Package management
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
requirements = c("PhyloMeasures", "phangorn", "parallel", "dplyr",
                 "ape", 'phyloregion')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})

setwd("D:/R projects/Global_ED")

#Load phylogenetic tree of bird families from package "ape"
#data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
#comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
#for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
#colnames(comm) = bird.families$tip.label

#Use query function to calculate pd values for each community
#pd.query(bird.families,comm)

#Tree: phylogenetic tree
#Comm: species presence/absence matrix with rows as sampling regions and columns as species
#Area: a vector representing area size of each sampling region  

#calculate PE and RPE weighted by region area-------------------------------------------------
#Tree = spec_phy.3
#Comm = comm_mammal_exotic[,unique(sp_overlap_dat_1$Binomial)]
#Area = comm_mammal_exotic$area

#phyloregion::evol_distinct(x, type="fair.proportion")

calcu_mean_ED_parallel = function(Tree, Comm, Area){
  
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
  
  all_ED_vec = phyloregion::evol_distinct(tree = Tree, type = "fair.proportion")
  all_ED_names = names(all_ED_vec)
  all_ED_df = data.frame(species = all_ED_names, ED = all_ED_vec)
  sps = colnames(Comm.df)
  
  clusterExport(my.cluster, varlist = c("all_ED_df", 'sps', 'Area'), envir = environment())

  mean_ED = parApply(cl = my.cluster, Comm.df, 1, 
                     function(x){
                       #x = Comm.df[2,]
                       selected_sps = sps[which(x == 1)]
                       sps_ED = dplyr::filter(all_ED_df, species %in% selected_sps)
                       med = mean(sps_ED$ED)
                     })
  
  ED_list = parApply(cl = my.cluster, Comm.df, 1, 
                         function(x){
                           #x = Comm.df[2,]
                           selected_sps = sps[which(x == 1)]
                           sps_ED = dplyr::filter(all_ED_df, species %in% selected_sps)
                           return(sps_ED)
                         })
  
  
  sps_range_sizes = parApply(cl = my.cluster, Comm.df, 2, 
                     function(x){
                       #x = Comm.df[2,]
                       selected_areas = Area[which(x == 1)]
                       sps_range_size = sum(selected_areas)
                     })
  
  all_EDR = (all_ED_df %>% dplyr::filter(species %in% names(sps_range_sizes)) %>% 
                                    arrange(species) %>% 
                pull(ED))/sps_range_sizes[sort(names(sps_range_sizes))]
  
  clusterExport(my.cluster, varlist = c("all_ED_df", 'sps', 'Area', 'all_EDR'), envir = environment())
  mean_EDR = parApply(cl = my.cluster, Comm.df, 1, 
                     function(x){
                       #x = Comm.df[2,]
                       selected_sps = sps[which(x == 1)]
                       sps_EDR = all_EDR[selected_sps]
                       medr = mean(sps_EDR)
                     })
  
  EDR_list = parApply(cl = my.cluster, Comm.df, 1, 
                      function(x){
                        #x = Comm.df[2,]
                        selected_sps = sps[which(x == 1)]
                        sps_EDR = all_EDR[selected_sps]
                      })
  
  # Stop cluster
  parallel::stopCluster(cl = my.cluster)
  closeAllConnections()
  gc()
  
  out_df = data.frame(mean_ED = mean_ED,
                      mean_EDR = mean_EDR)
  out = list(df = out_df, 
             ED_list = ED_list,
             EDR_list = EDR_list)
  return(out)
  
}

