library(picante)
library(ape)
library(dplyr)
library(sf)
library(ggplot2)
library(ggpmisc)
library(cowplot)

setwd("E:/Doctoral/World_Database/World_Birds/Final_data/")

phylo_tre <- read.tree("Phylogenetic_Birds.tre")
species_distribution <- read.csv("Distribution_data_note.csv")
species_distribution$Species <- gsub(" ", "_", species_distribution$ScientificName)

# No consideration of presence, origin and Seasonality
alien_distribution <- subset(species_distribution, SpStatus == "alien")
native_distribution <- subset(species_distribution, SpStatus == "Native")

alien_distribution_unique <- alien_distribution %>% distinct(RegionID, ScientificName, .keep_all = TRUE)
native_distribution_unique <- native_distribution %>% distinct(RegionID, ScientificName, .keep_all = TRUE)

distribution_unique <- rbind(alien_distribution_unique, native_distribution_unique)
Num_alien <- alien_distribution_unique %>% count(RegionID)
Num_native <- native_distribution_unique %>% count(RegionID)
Num_all <- distribution_unique %>% count(RegionID)

Num_species <- merge(Num_all, Num_native, by = "RegionID", all = TRUE)
Num_species <- merge(Num_species, Num_alien, by = "RegionID", all = TRUE)
colnames(Num_species) <- c("RegionID", "Num_all", "Num_native", "Num_alien")

# Consideration of presence, origin and Seasonality
species_distribution_breed <- subset(species_distribution, !(seasonal %in% c("3", "4", "5")))
species_distribution_breed <- subset(species_distribution_breed, !(presence %in% c("3", "4", "5", "6", "7")))

alien_distribution_breed <- subset(species_distribution_breed, SpStatus == "alien")
native_distribution_breed <- subset(species_distribution_breed, SpStatus == "Native")

alien_distribution_breed_unique <- alien_distribution_breed %>% distinct(RegionID, ScientificName, .keep_all = TRUE)
native_distribution_breed_unique <- native_distribution_breed %>% distinct(RegionID, ScientificName, .keep_all = TRUE)

distribution_breed_unique <- rbind(alien_distribution_breed_unique, native_distribution_breed_unique)
Num_alien <- alien_distribution_breed_unique %>% count(RegionID)
Num_native <- native_distribution_breed_unique %>% count(RegionID)
Num_all <- distribution_breed_unique %>% count(RegionID)

Num_species_breed <- merge(Num_all, Num_native, by = "RegionID", all = TRUE)
Num_species_breed <- merge(Num_species_breed, Num_alien, by = "RegionID", all = TRUE)
colnames(Num_species_breed) <- c("RegionID", "Num_all", "Num_native", "Num_alien")

Region.ID1 <- unique(alien_distribution_unique$RegionID)
Region.ID2 <- unique(native_distribution_unique$RegionID)
Region.ID <- intersect(Region.ID1, Region.ID2)
# Region.ID <- setdiff(Region.ID, c(596, 442, 569, 453, 189, 517, 41)) #only have alien species

Index_results <- data.frame(RegionID=Region.ID,mpd.beta=NA,mntd.beta=NA, mpd.alpha.natu=NA, 
                            mpd.alpha.native=NA,mntd.alpha.natu=NA, mntd.alpha.native=NA, pd.natu=NA,
                            pd.native=NA, SR.natu=NA, SR.native=NA)

# Function
mpdnull <- function(data.natu,data.na,results){
  #calcute beta.mpd & beta.mntd between data.naturalized and data.native
  #and alpha.mpd & alpha.mptd of data.naturalized and data.native
  #and faith's PD
  data.pd=rbind(data.na,data.natu)
  # convert data format
  data.pd$naturalized=0
  data.pd$naturalized[which(data.pd$SpStatus %in% c('alien'))]=1
  data.pd$native=0
  data.pd$native[which(data.pd$SpStatus=='Native')]=1
  data.pd=data.pd[,c("Species","naturalized","native")] %>% unique()
  cnames=data.pd$Species
  data.pd2=data.pd[,-1]
  data.pd2=t(data.pd2)
  colnames(data.pd2)=cnames
  re.temp=match.comm.dist(comm = data.pd2,dis = results$distancematrix)
  data.pd2=re.temp$comm
  distm=re.temp$dist
  # calculate mpd
  mpd.beta=as.matrix(comdist(data.pd2,distm,abundance.weighted = F))[1,2]
  mntd.beta=as.matrix(comdistnt(data.pd2,distm,abundance.weighted = F))[1,2]
  mpd.alpha=mpd(data.pd2,distm,abundance.weighted = F)
  mntd.alpha=mntd(data.pd2,distm,abundance.weighted = F)
  re.temp=match.phylo.comm(phy = results$tree,comm = data.pd2)
  data.pd2=re.temp$comm
  mytree=re.temp$phy
  # calculate pd
  # Attempt to calculate PD using mytree
  pd.alpha <- tryCatch({
    pd(samp = data.pd2, tree = mytree, include.root = TRUE)
  }, error = function(e) {
    message("Error with mytree, switching to phylo_tre: ", e$message)
    
    # Try the second option with phylo_tre
    tryCatch({
      pd(samp = data.pd2, tree = phylo_tre, include.root = TRUE)
    }, error = function(e) {
      message("Error with phylo_tre, switching to include.root = FALSE: ", e$message)
      
      # Modify include.root to FALSE if both trees fail
      pd(samp = data.pd2, tree = phylo_tre, include.root = FALSE)
    })
  })
  
  pdlist=data.frame(mpd.beta=mpd.beta,mntd.beta=mntd.beta,
                    mpd.alpha.natu=mpd.alpha[1],mpd.alpha.native=mpd.alpha[2],
                    mntd.alpha.natu=mntd.alpha[1],mntd.alpha.native=mntd.alpha[2],
                    pd.natu=pd.alpha["naturalized","PD"],pd.native=pd.alpha["native","PD"],
                    SR.natu=pd.alpha["naturalized","SR"],SR.native=pd.alpha["native","SR"])
  return(pdlist)
}

for (k in 1:length(Region.ID)){
  re_id=Region.ID[k]
  data_na=subset(native_distribution_breed_unique,native_distribution_breed_unique$RegionID == re_id)
  data_natu=subset(alien_distribution_breed_unique, (alien_distribution_breed_unique$RegionID == re_id))
  
  if (length(data_natu$Species)>0){
    species_to_keep <- c(data_na$Species, data_natu$Species)
    species_to_remove <- setdiff(phylo_tre$tip.label, species_to_keep)
    tree=drop.tip(phylo_tre,species_to_remove)
    
    distm = cophenetic.phylo(tree)
    results=list(tree=tree, distancematrix=distm)
    
    Index_results[k,2:11]=unlist(mpdnull(data.natu = data_natu,data.na = data_na,results = results))
  }else{
    Index_results[k,2:11]=NA
  }
  
  print(k)
}

write.csv(Index_results, file = "Birds_index.csv")

#plot x = lat, y = mpd
library(sf)
tdwg <- st_read("E:/Doctoral/Task/TDWG4/TDWG4_newTibet.shp")
tdwg$Lat_abs <- abs(tdwg$Lat)
tdwg <- tdwg[,c("RegionID", "Lon", "Lat", "Lat_abs")]

Birds_index <- left_join(Index_results, tdwg, by = "RegionID")

plot(Birds_index$Lat_abs, Birds_index$mpd.beta)

# Create the plot with linear regression line and confidence interval
# Define the function to create the plot
myplot <- function(data, x, y, xlab, ylab) {
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(color = "gray", size = 1.5) +  # Scatter plot with gray points and adjusted size
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +  # Linear regression line with 95% confidence interval
    stat_poly_eq(aes(label = paste(..rr.label..)),  # Add R-squared value to the plot
                 formula = y ~ x, parse = TRUE, size = 3) +  # Use the formula y ~ x for linear regression
    labs(x = xlab, y = ylab) +  # Axis labels
    theme_classic() +  # Change to a classic theme without gridlines
    theme(panel.grid = element_blank())  # Remove gridlines
}

p1 <- myplot(Birds_index, "Lat_abs", "mpd.beta", "Absolute Latitude (°)", "MPD")
p2 <- myplot(Birds_index, "Lat_abs", "mntd.beta", "Absolute Latitude (°)", "MNTD")

pdf(paste0('figures',"/Fig.MPD.LAT.pdf",sep=""),width = 7, height = 3)
ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(p2, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("a","b"), size = 10,
                  x = c(0,0.5), y = c(1,1))
dev.off()

# null model
# function
#### Function
samplempd.full.natu <- function(numbers){
  speciesnumber=length(data_natu$Species)
  speciespool$num=c(1:length(speciespool$Species))
  data.sam.id=sample(speciespool$num,size = speciesnumber,replace = F)
  data.sam=speciespool[data.sam.id,-length(speciespool)]
  mpdmntd=mpdnull(data.natu = data.sam,data.na = data_na,results = results)
}

match.comm.dist.fsy=function (comm, dis){
  res <- list()
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to match distance matrix and community data")
  }
  disclass <- dis
  dis <- as.matrix(dis)
  distaxa <- rownames(dis)
  if (is.null(distaxa)) {
    warning("Distance matrix lacks taxa names, these are required to match community and distance matrix. Data are returned unsorted. Assuming that distance matrix and community data taxa columns are in the same order!")
    if (inherits(disclass, "dist")) {
      return(list(comm = comm, dist = as.dist(dis)))
    }
    else {
      return(list(comm = comm, dist = dis))
    }
  }
  if (!all(distaxa %in% commtaxa)) {
    print("Dropping taxa from the distance matrix because they are not present in the community data:")
    print(setdiff(distaxa, commtaxa))
    dis <- dis[intersect(distaxa, commtaxa), intersect(distaxa, 
                                                       commtaxa)]
    distaxa <- rownames(dis)
  }
  if (any(!(commtaxa %in% distaxa))) {
    print("Dropping taxa from the community because they are not present in the distance matrix:")
    print(setdiff(commtaxa, distaxa))
    res$comm <- comm[, intersect(commtaxa, distaxa)]
  }
  else {
    res$comm <- comm
  }
  if (inherits(disclass, "dist")) {
    res$dist <- as.dist(dis[colnames(res$comm), colnames(res$comm)])
  }
  else {
    res$dist <- dis[colnames(res$comm), colnames(res$comm)]
  }
  return(res)
}

# data process
#species.pool - all naturalized species except for the natives in the focal region
alien_sp_pool_global <- alien_distribution_unique[!duplicated(alien_distribution_unique$Species), ]

# data.na.all.new=data.na.APnew.618.intree.nonative.gbif
# data.natu.all.new=data.natu.APnew.618.intree.nonative.gbif

all_results <- data.frame()

for (k in 1:length(Region.ID)){
  re_id=Region.ID[k]
  
  results.null <- data.frame(RegionID=rep(re_id,1000),mpd.beta=NA,mntd.beta=NA, mpd.alpha.natu=NA, 
                              mpd.alpha.native=NA,mntd.alpha.natu=NA, mntd.alpha.native=NA, pd.natu=NA,
                              pd.native=NA, SR.natu=NA, SR.native=NA)
  
  data_na=subset(native_distribution_unique,native_distribution_unique$RegionID == re_id)
  data_natu=subset(alien_distribution_unique,alien_distribution_unique$RegionID == re_id)
  
  #SES.MPD
  #delete native species from species pool-native
  sp.pool.sub=alien_sp_pool_global
  # sp.pool.sub$status='naturalized'
  # data_na=data_na[,c("species","genus","family","group","status")]
  data_both=rbind(data_na,sp.pool.sub)
  id=duplicated(data_both[,c("Species")])
  data_both=data_both[!id,]
  speciespool=subset(data_both,SpStatus == 'alien')
  
  #
  species_to_keep <- c(data_both$Species)
  species_to_remove <- setdiff(phylo_tre$tip.label, species_to_keep)
  tree=drop.tip(phylo_tre,species_to_remove)
  distm = cophenetic.phylo(tree)
  results=list(tree=tree, distancematrix=distm)
  
  #
  for (i in 1:1000){
    results.null[i,2:11]=samplempd.full.natu(i)
  }
  
  # After each iteration, append the results for this RegionID to the all_results data frame
  all_results = rbind(all_results, results.null)
  
  print(k)
}

save(all_results, Index_results, file = "Birds_index_results.rda")

#############
# Continent #
#############
continent <- read.csv("E:/Doctoral/World_Database/World_Birds/TDWG_RegionID.csv")
continent$ContinentID <- as.numeric(continent$ContinentID)
alien_sp_pool_continent <- left_join(alien_distribution_unique, continent, by = "RegionID")

# data.na.all.new=data.na.APnew.618.intree.nonative.gbif
# data.natu.all.new=data.natu.APnew.618.intree.nonative.gbif


# parallel computing
process_region <- function(k) {
  re_id <- Region.ID[k]
  
  results.null <- data.frame(RegionID = rep(re_id, 1000), mpd.beta = NA, mntd.beta = NA, 
                             mpd.alpha.natu = NA, mpd.alpha.native = NA, mntd.alpha.natu = NA, 
                             mntd.alpha.native = NA, pd.natu = NA, pd.native = NA, 
                             SR.natu = NA, SR.native = NA)
  
  data_na <- subset(native_distribution_unique, native_distribution_unique$RegionID == re_id)
  data_natu <- subset(alien_distribution_unique, alien_distribution_unique$RegionID == re_id)
  
  # SES.MPD
  Contin.ID <- continent %>% filter(RegionID == re_id) %>% select(ContinentID)
  Contin.ID <- as.numeric(Contin.ID)
  sp.pool.sub <- alien_sp_pool_continent %>% filter(ContinentID == Contin.ID)
  sp.pool.sub <- subset(sp.pool.sub, select = -ContinentID)
  sp.pool.sub <- sp.pool.sub[!duplicated(sp.pool.sub$Species), ]
  
  data_both <- rbind(data_na, sp.pool.sub)
  id <- duplicated(data_both[, c("Species")])
  data_both <- data_both[!id, ]
  speciespool <- subset(data_both, SpStatus == 'alien')
  
  species_to_keep <- c(data_both$Species)
  species_to_remove <- setdiff(phylo_tre$tip.label, species_to_keep)
  tree <- drop.tip(phylo_tre, species_to_remove)
  distm <- cophenetic.phylo(tree)
  results <- list(tree = tree, distancematrix = distm)
  
  for (i in 1:1000) {
    results.null[i, 2:11] <- samplempd.full.natu(i)
  }
  
  return(results.null)
}

library(parallel)

# 使用并行计算
cl <- makeCluster(detectCores() - 1)  # 使用所有可用的核心数，减去1
# Load necessary libraries on each worker node
clusterEvalQ(cl, {
  library(dplyr)
  library(ape)  # Load ape for the drop.tip function
  library(picante)
})

# Export the necessary function and data objects to the workers
clusterExport(cl, c("Region.ID", "native_distribution_unique", "alien_distribution_unique", 
                    "continent", "alien_sp_pool_continent", "phylo_tre", "samplempd.full.natu",
                    "data_natu", "data_na", "speciespool", "mpdnull", "match.comm.dist", "results"))  # Export the mpdnull function

# Use parLapply for parallel computation
all_results_list <- parLapply(cl, 1:4, process_region)

# Close the cluster after use
stopCluster(cl)

# Combine the results
all_results <- do.call(rbind, all_results_list)

# Save the results
save(all_results, file = "Birds_index_results_global.rda")