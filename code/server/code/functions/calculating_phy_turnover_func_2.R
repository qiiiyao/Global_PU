library(phyloregion)

#1. Phylogenetic Simpson index----

#tree = phylo_mammal_extant
#x = mat_mammal_native[c('41', '84'),]

calcu_phy_turn_pair = function(tree, x){
  
  com.tot.pair = colSums(x)>0
  com.all = rbind(x, com.tot.pair)
  
  com.sparse = dense2sparse(com.all)
  pd.all = phyloregion::PD(x = com.sparse, tree) # PD for each community of the community matrix
  
  shared = sum(pd.all[1:2]) - pd.all[3]  # a
  
  min.pd = min(pd.all) # min(b,c)
  
  #phylogenetic turnover (Simpson dissimilarity index)
  phylo_sim_all_0 = 1 - (shared/min.pd)
  
  return(phylo_sim_all_0)
  
} # end of function


#tree = phylo_mammal_extant
#x = mat_mammal_path2
#focal_region = region1

calcu_phy_turn_focal_mean = function(tree, x, focal_region) {
  
  ## Ensure focal_region can be either row name or row index
  if (is.character(focal_region)) {
    focal_id = match(focal_region, rownames(x))
  } else {
    focal_id = focal_region
  }
  
  ## Identify other regions
  other_ids = setdiff(seq_len(nrow(x)), focal_id)
  
  ## PD of all single regions
  pd_single = phyloregion::PD(
    x = dense2sparse(x),
    tree
  )
  
  pd_focal = pd_single[focal_id]
  pd_other = pd_single[other_ids]
  
  ## Construct focal + other union communities only
  x_other = x[other_ids, , drop = FALSE]
  x_focal = x[focal_id, , drop = FALSE]
  
  x_union = x_other
  focal_sps = which(x_focal[1, ] > 0)
  x_union[, focal_sps] = 1
  
  ## PD of focal + each other region
  pd_union = phyloregion::PD(
    x = dense2sparse(x_union),
    tree
  )
  
  ## Shared phylogenetic diversity
  shared_pd = pd_focal + pd_other - pd_union
  
  ## Simpson phylogenetic turnover
  phy_turn = 1 - (shared_pd / pmin(pd_focal, pd_other))
  
  phy_turn[which(phy_turn < 0 & phy_turn > -1e-12)] = 0
  
  ## Return mean turnover of focal region against all other regions
  mean_phy_turn = mean(c(phy_turn,0), na.rm = TRUE) # plus focal region vs. itself
  
  return(mean_phy_turn)
}


#tree = phylo_plant_native
#x = mat_plant_native[1:100,]
#block_size = 400

calcu_phy_turn_multiple = function(tree, x, block_size){
  
  pd_tot_mat = matrix(nrow = nrow(x),
                      ncol = nrow(x))
  rownames(pd_tot_mat) = rownames(x)
  colnames(pd_tot_mat) = rownames(x)
  
  region_pairs = combn(rownames(x), 2)
  
  sep_pairs = c(seq(0, ncol(region_pairs), block_size), ncol(region_pairs))
  sep_region_pairs = data.frame(start = (sep_pairs[1:(length(sep_pairs)-1)]+1),
                                end = sep_pairs[2:length(sep_pairs)])
  
  for (b in 1:nrow(sep_region_pairs)) {
    #b = 1
    pair_sub = region_pairs[,sep_region_pairs[b,1]:sep_region_pairs[b,2]]
    #pair_sub = region_pairs[,1:100]
    
    ## union communities
    com_union = (x[pair_sub[1, ], , drop = FALSE] +
                 x[pair_sub[2, ], , drop = FALSE]) > 0
    
    storage.mode(com_union) = "numeric"
    
    colnames(com_union) = colnames(x)
    
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  
  pd_region = pd_region[rownames(pd_tot_mat)]
  
  # 1. PD_i + PD_j
  pd_sum_mat = outer(pd_region, pd_region, "+")
  
  # 2. shared branch length
  shared_mat = pd_sum_mat - pd_tot_mat
  
  # 3. min(PD_i, PD_j)
  pd_min_mat = outer(pd_region, pd_region, pmin)
  dimnames(pd_min_mat) = dimnames(pd_tot_mat)
  
  # 4. phylogenetic turnover
  phy_turn_mat = 1 - (shared_mat / pd_min_mat)
  
  phy_turn_mat[which(phy_turn_mat < 0 & phy_turn_mat > -1e-12)] = 0
  
  diag(phy_turn_mat) = 0
  
  out = phy_turn_mat
  return(out)
  
} # end of function



#2. Phylogenetic Sorensen index----
calcu_phy_sor_pair = function(tree, x){
  
  com.tot.pair = colSums(x)>0
  com.all = rbind(x, com.tot.pair)
  
  com.sparse = dense2sparse(com.all)
  pd.all = phyloregion::PD(x = com.sparse, tree) # PD for each community of the community matrix
  
  sum.all.com = sum(pd.all[1:2])
  shared = sum.all.com - pd.all[3]  # a
  
  #phylogenetic sorensen (Simpson dissimilarity index)
  phylo_sim_all_0 = 1 - ((2 * shared)/sum.all.com)
  
  return(phylo_sim_all_0)
  
} # end of function

#tree = phylo_plant_native
#x = mat_plant_native[1:100,]
#block_size = 400

calcu_phy_sor_multiple = function(tree, x, block_size){
  
  pd_tot_mat = matrix(nrow = nrow(x),
                      ncol = nrow(x))
  rownames(pd_tot_mat) = rownames(x)
  colnames(pd_tot_mat) = rownames(x)
  
  region_pairs = combn(rownames(x), 2)
  
  sep_pairs = c(seq(0, ncol(region_pairs), block_size), ncol(region_pairs))
  sep_region_pairs = data.frame(start = (sep_pairs[1:(length(sep_pairs)-1)]+1),
                                end = sep_pairs[2:length(sep_pairs)])
  
  for (b in 1:nrow(sep_region_pairs)) {
    #b = 1
    pair_sub = region_pairs[,sep_region_pairs[b,1]:sep_region_pairs[b,2]]
    #pair_sub = region_pairs[,1:100]
    
    ## union communities
    com_union = (x[pair_sub[1, ], , drop = FALSE] +
                   x[pair_sub[2, ], , drop = FALSE]) > 0
    
    storage.mode(com_union) = "numeric"
    
    colnames(com_union) = colnames(x)
    
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  
  pd_region = pd_region[rownames(pd_tot_mat)]
  
  # 1. PD_i + PD_j
  pd_sum_mat = outer(pd_region, pd_region, "+")
  
  # 2. shared branch length
  shared_mat = pd_sum_mat - pd_tot_mat
  
  # 4. phylogenetic sorover
  phy_sor_mat = 1 - ((2 * shared_mat) / pd_sum_mat)
  
  phy_sor_mat[which(phy_sor_mat < 0 & phy_sor_mat > -1e-12)] = 0
  
  diag(phy_sor_mat) = 0
  
  out = phy_sor_mat
  return(out)
  
} # end of function



#3. Phylogenetic Jaccard index----
calcu_phy_jac_pair = function(tree, x){
  
  com.tot.pair = colSums(x)>0
  com.all = rbind(x, com.tot.pair)
  
  com.sparse = dense2sparse(com.all)
  pd.all = phyloregion::PD(x = com.sparse, tree) # PD for each community of the community matrix
  
  shared = sum(pd.all[1:2]) - pd.all[3]  # a
  
  #phylogenetic jaccard
  phylo_sim_all_0 = 1 - (shared/pd.all[3])
  
  return(phylo_sim_all_0)
  
} # end of function

#tree = phylo_plant_native
#x = mat_plant_native[1:100,]
#block_size = 400

calcu_phy_jac_multiple = function(tree, x, block_size){
  
  pd_tot_mat = matrix(nrow = nrow(x),
                      ncol = nrow(x))
  rownames(pd_tot_mat) = rownames(x)
  colnames(pd_tot_mat) = rownames(x)
  
  region_pairs = combn(rownames(x), 2)
  
  sep_pairs = c(seq(0, ncol(region_pairs), block_size), ncol(region_pairs))
  sep_region_pairs = data.frame(start = (sep_pairs[1:(length(sep_pairs)-1)]+1),
                                end = sep_pairs[2:length(sep_pairs)])
  
  for (b in 1:nrow(sep_region_pairs)) {
    #b = 1
    pair_sub = region_pairs[,sep_region_pairs[b,1]:sep_region_pairs[b,2]]
    #pair_sub = region_pairs[,1:100]
    
    ## union communities
    com_union = (x[pair_sub[1, ], , drop = FALSE] +
                   x[pair_sub[2, ], , drop = FALSE]) > 0
    
    storage.mode(com_union) = "numeric"
    
    colnames(com_union) = colnames(x)
    
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  
  pd_region = pd_region[rownames(pd_tot_mat)]
  
  # 1. PD_i + PD_j
  pd_sum_mat = outer(pd_region, pd_region, "+")
  
  # 2. shared branch length
  shared_mat = pd_sum_mat - pd_tot_mat
  
  # 4. phylogenetic jaccard
  phy_jac_mat = 1 - (shared_mat / pd_tot_mat)
  
  phy_jac_mat[which(phy_jac_mat < 0 & phy_jac_mat > -1e-12)] = 0
  
  diag(phy_jac_mat) = 0
  
  out = phy_jac_mat
  return(out)
  
} # end of function





#4. Phylogenetic Ruggiero index----
calcu_phy_rlb_pair = function(tree, x){
  
  com.tot.pair = colSums(x)>0
  com.all = rbind(x, com.tot.pair)
  
  com.sparse = dense2sparse(com.all)
  pd.all = phyloregion::PD(x = com.sparse, tree) # PD for each community of the community matrix
  
  shared = sum(pd.all[1:2]) - pd.all[3]  # a
  
  #phylogenetic ruggiero (Ruggiero dissimilarity index)
  phylo_sim_all_0 = 1 - (shared/pd.all[1:2])
  
  names(phylo_sim_all_0) = rownames(x)
  return(phylo_sim_all_0)
  
} # end of function

#tree = phylo_plant_native
#x = mat_plant_native[1:100,]
#block_size = 400

calcu_phy_rlb_multiple = function(tree, x, block_size){
  
  pd_tot_mat = matrix(nrow = nrow(x),
                      ncol = nrow(x))
  rownames(pd_tot_mat) = rownames(x)
  colnames(pd_tot_mat) = rownames(x)
  
  region_pairs = combn(rownames(x), 2)
  
  sep_pairs = c(seq(0, ncol(region_pairs), block_size), ncol(region_pairs))
  sep_region_pairs = data.frame(start = (sep_pairs[1:(length(sep_pairs)-1)]+1),
                                end = sep_pairs[2:length(sep_pairs)])
  
  for (b in 1:nrow(sep_region_pairs)) {
    #b = 1
    pair_sub = region_pairs[,sep_region_pairs[b,1]:sep_region_pairs[b,2]]
    #pair_sub = region_pairs[,1:100]
    
    ## union communities
    com_union = (x[pair_sub[1, ], , drop = FALSE] +
                   x[pair_sub[2, ], , drop = FALSE]) > 0
    
    storage.mode(com_union) = "numeric"
    
    colnames(com_union) = colnames(x)
    
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  
  pd_region = pd_region[rownames(pd_tot_mat)]
  
  # 1. PD_i + PD_j
  pd_sum_mat = outer(pd_region, pd_region, "+")
  
  # 2. shared branch length
  shared_mat = pd_sum_mat - pd_tot_mat
  
  # 3. asymmetric matrix(PD_i, PD_j)
  pd_focal_mat = matrix(data = rep(pd_region, length(pd_region)), 
                        nrow = length(pd_region), ncol = length(pd_region))
  # rows are focal regions, whereas columns are compared region
  dimnames(pd_focal_mat) = dimnames(pd_tot_mat)
  
  # 4. phylogenetic ruggiero
  phy_rlb_mat = 1 - (shared_mat / pd_focal_mat)
  
  phy_rlb_mat[which(phy_rlb_mat < 0 & phy_rlb_mat > -1e-12)] = 0
  
  diag(phy_rlb_mat) = 0
  
  out = phy_rlb_mat
  return(out)
  
} # end of function


#5. Phylogenentic beta all index ----


calcu_phy_sim_multiple_all = function(tree, x, block_size){
  
  pd_tot_mat = matrix(nrow = nrow(x),
                      ncol = nrow(x))
  rownames(pd_tot_mat) = rownames(x)
  colnames(pd_tot_mat) = rownames(x)
  
  region_pairs = combn(rownames(x), 2)
  
  sep_pairs = c(seq(0, ncol(region_pairs), block_size), ncol(region_pairs))
  sep_region_pairs = data.frame(start = (sep_pairs[1:(length(sep_pairs)-1)]+1),
                                end = sep_pairs[2:length(sep_pairs)])
  
  for (b in 1:nrow(sep_region_pairs)) {
    #b = 1
    pair_sub = region_pairs[,sep_region_pairs[b,1]:sep_region_pairs[b,2]]
    #pair_sub = region_pairs[,1:100]
    
    ## union communities
    com_union = (x[pair_sub[1, ], , drop = FALSE] +
                   x[pair_sub[2, ], , drop = FALSE]) > 0
    
    storage.mode(com_union) = "numeric"
    
    colnames(com_union) = colnames(x)
    
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  
  pd_region = pd_region[rownames(pd_tot_mat)]
  
  # 1. PD_i + PD_j
  pd_sum_mat = outer(pd_region, pd_region, "+")
  
  # 2. shared branch length
  shared_mat = pd_sum_mat - pd_tot_mat
  
  # 3_1. min(PD_i, PD_j)
  pd_min_mat = outer(pd_region, pd_region, pmin)
  dimnames(pd_min_mat) = dimnames(pd_tot_mat)
  
  # 3_2. asymmetric matrix(PD_i, PD_j)
  pd_focal_mat = matrix(data = rep(pd_region, length(pd_region)), 
                        nrow = length(pd_region), ncol = length(pd_region))
  # rows are focal regions, whereas columns are compared region
  dimnames(pd_focal_mat) = dimnames(pd_tot_mat)
  
  # 4_1. phylogenetic turnover
  phy_turn_mat = 1 - (shared_mat / pd_min_mat)
  phy_turn_mat[which(phy_turn_mat < 0 & phy_turn_mat > -1e-12)] = 0
  diag(phy_turn_mat) = 0
  
  # 4_2. phylogenetic sorensen
  phy_sor_mat = 1 - ((2 * shared_mat) / pd_sum_mat)
  phy_sor_mat[which(phy_sor_mat < 0 & phy_sor_mat > -1e-12)] = 0
  diag(phy_sor_mat) = 0
  
  # 4_3. phylogenetic jaccard
  phy_jac_mat = 1 - (shared_mat / pd_tot_mat)
  phy_jac_mat[which(phy_jac_mat < 0 & phy_jac_mat > -1e-12)] = 0
  diag(phy_jac_mat) = 0
  
  # 4_4. phylogenetic ruggiero
  phy_rlb_mat = 1 - (shared_mat / pd_focal_mat)
  phy_rlb_mat[which(phy_rlb_mat < 0 & phy_rlb_mat > -1e-12)] = 0
  diag(phy_rlb_mat) = 0
  
  out = list(simpson = phy_turn_mat, sor = phy_sor_mat,
             jac = phy_jac_mat, rlb = phy_rlb_mat)
  return(out)
  
} # end of function



#6. phylogenetic beta sim core ----

calcu_phy_sim_multiple_core = function(tree, x, block_size){
  
  pd_tot_mat = matrix(nrow = nrow(x),
                      ncol = nrow(x))
  rownames(pd_tot_mat) = rownames(x)
  colnames(pd_tot_mat) = rownames(x)
  
  region_pairs = combn(rownames(x), 2)
  
  sep_pairs = c(seq(0, ncol(region_pairs), block_size), ncol(region_pairs))
  sep_region_pairs = data.frame(start = (sep_pairs[1:(length(sep_pairs)-1)]+1),
                                end = sep_pairs[2:length(sep_pairs)])
  
  for (b in 1:nrow(sep_region_pairs)) {
    #b = 1
    pair_sub = region_pairs[,sep_region_pairs[b,1]:sep_region_pairs[b,2]]
    #pair_sub = region_pairs[,1:100]
    
    ## union communities
    com_union = (x[pair_sub[1, ], , drop = FALSE] +
                   x[pair_sub[2, ], , drop = FALSE]) > 0
    
    storage.mode(com_union) = "numeric"
    
    colnames(com_union) = colnames(x)
    
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  
  pd_region = pd_region[rownames(pd_tot_mat)]
  
  out = list(pd_region = pd_region, total = pd_tot_mat)
  return(out)
  
} # end of function




# Paired matrix to distance matrix conversion (utility function) ----

dist.mat = function(com,pair) {
  
  ncom = nrow(com)
  distmat = matrix(nrow=ncom,ncol=ncom,0,dimnames=list(rownames(com),rownames(com)))
  st = c(0,cumsum(seq(ncom-1,2)))+1
  end = cumsum(seq(ncom-1,1))
  for (i in 1:(ncom-1)) distmat[i,(ncom:(seq(1,ncom)[i]))]=c(pair[end[i]:st[i]],0)
  distmat = as.dist(t(distmat))
  return(distmat)
  
} # end of function dist.mat

