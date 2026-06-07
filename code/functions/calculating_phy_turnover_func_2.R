library(phyloregion)

#tree = phylo
#x = mat_fish_path3_5[c(region1, region2),]

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



calcu_phy_turn_simple = function(tree, x, Region_posi){
  
  com = as.matrix(x[, colnames(x) %in% tree$tip.label])
  storage.mode(com) = "numeric"
  row.names(com) = unlist(x[,Region_posi])
  
  com = x
  # pariwise comparisons
  combin = combn(rownames(com),2) # table with all pairs
  
  pd = phyloregion::PD(x = com, tree) # PD for each community of the community matrix
  names(pd) = rownames(com)
  
  com.tot.pair = t(apply(combin,2,function(x) colSums(com[x,])>0))
  com.tot.pair[which(com.tot.pair == T)] = 1
  
  if(nrow(com.tot.pair) == 1) {
    com.tot.pair = rbind(com.tot.pair, com.tot.pair)
    pd.tot.pair = unique(phyloregion::PD(x = com.tot.pair, tree))
    
  } else{
    #pd.tot.pair = pdnew(com.tot.pair[c(1:3),c(1:3)],tree)[,"PD"]  # PD of the two communities combined
    pd.tot.pair = phyloregion::PD(x = com.tot.pair, tree)
  }
  
  sum.pd.pair = apply(combin,2,function(x) sum(pd[x])) # Sum of PD for each community, separetely
  
  shared = sum.pd.pair - pd.tot.pair  # a
  shared = dist.mat(com, shared)
  
  min.not.shared = apply(pd.tot.pair - t(combn(pd, 2)), 1, min) # min(b,c)
  min.not.shared = dist.mat(com, min.not.shared)
  
  #phylogenetic turnover (Simpson dissimilarity index)
  phylo_sim_all_0 = min.not.shared/(min.not.shared + shared)
  phylo_sim_all_0 = as.matrix(phylo_sim_all_0)
  
  #which(phylo_sim_all_1 < 0)
  colnames(phylo_sim_all_0) = rownames(com)
  rownames(phylo_sim_all_0) = rownames(com)
  
  out = phylo_sim_all_0
  return(out)
  
} # end of function



#tree = spec_phy.3
#x = comm_mammal_path2[,2:ncol(comm_mammal_path2)]
#block_size = 50000

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
    
    com_union = dense2sparse(com_union)
    pd_union = phyloregion::PD(com_union, tree)
    pd_union = as.numeric(pd_union)
    
    i1 = pair_sub[1, ]
    i2 = pair_sub[2, ]
    
    pd_tot_mat[cbind(i1, i2)] = pd_union
    pd_tot_mat[cbind(i2, i1)] = pd_union
    
  }
  x = dense2sparse(x)
  pd_region = phyloregion::PD(x, tree) # PD for each community of the community matrix
  names(pd_region) = rownames(pd_tot_mat)
  
  # 1. PD_i + PD_j
  pd_sum_mat = outer(pd_region, pd_region, "+")
  
  # 2. shared branch length
  shared_mat = pd_sum_mat - pd_tot_mat
  
  # 3. min(PD_i, PD_j)
  pd_min_mat = outer(pd_region, pd_region, pmin)
  dimnames(pd_min_mat) = dimnames(pd_tot_mat)
  
  # 4. phylogenetic turnover
  phy_turn_mat = 1 - shared_mat / pd_min_mat
  
  phy_turn_mat[which(phy_turn_mat < 0 & phy_turn_mat > -1e-12)] = 0
  
  diag(phy_turn_mat) = 0
  
  out = phy_turn_mat
  return(out)
  
} # end of function



############ Paired matrix to distance matrix conversion (utility function) #######################

dist.mat = function(com,pair) {
  
  ncom = nrow(com)
  distmat = matrix(nrow=ncom,ncol=ncom,0,dimnames=list(rownames(com),rownames(com)))
  st = c(0,cumsum(seq(ncom-1,2)))+1
  end = cumsum(seq(ncom-1,1))
  for (i in 1:(ncom-1)) distmat[i,(ncom:(seq(1,ncom)[i]))]=c(pair[end[i]:st[i]],0)
  distmat = as.dist(t(distmat))
  return(distmat)
  
} # end of function dist.mat

