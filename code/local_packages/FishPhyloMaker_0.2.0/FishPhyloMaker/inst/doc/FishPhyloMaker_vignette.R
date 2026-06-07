## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install_pkg, echo=TRUE, eval=FALSE---------------------------------------
#  devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main")

## ----read_data, eval=FALSE, echo=TRUE-----------------------------------------
#  library(FishPhyloMaker)
#  data(neotropical_comm)
#  data_comm <- neotropical_comm[, -c(1, 2)] # removing Latitude and Longitude

## ----taxon_data, echo=TRUE, eval=FALSE----------------------------------------
#  taxon_data <- FishTaxaMaker(data_comm, allow.manual.insert = TRUE)
#  Characidae
#  Characiformes
#  Characidae
#  Characiformes
#  Characidae
#  Characiformes
#  Loricariidae
#  Siluriformes
#  Characidae
#  Characiformes
#  Cichlidae
#  Cichliformes
#  Crenuchidae
#  Characiformes
#  Gymnotidae
#  Gymnotiformes
#  Loricariidae
#  Siluriformes
#  Loricariidae
#  Siluriformes
#  Loricariidae
#  Siluriformes
#  Loricariidae
#  Siluriformes
#  Heptapteridae
#  Siluriformes
#  Characidae
#  Characiformes
#  Loricariidae
#  Siluriformes
#  Characidae
#  Characiformes

## ----phylo_make, eval=FALSE, echo=TRUE----------------------------------------
#  phylo_fish_streams <- FishPhyloMaker(data = taxon_data,
#                                       return.insertions = TRUE,
#                                       insert.base.node = TRUE,
#                                       progress.bar = TRUE)

## ----plot_phylo, eval=FALSE, echo=TRUE----------------------------------------
#  library(ggtree)
#  tree.PR<- phylo_fish_streams$Phylogeny
#  
#  tree.PR <- ape::makeNodeLabel(tree.PR)
#  phylo <- tree.PR
#  
#  rm.famNames <- which(table(taxon_dataPR$f) == 1) # monotipic families
#  names.fam <- setdiff(unique(taxon_dataPR$f), names(rm.famNames)) # removing monotipic families from the names
#  
#  for (i in 1:length(names.fam)) {
#    set <- subset(taxon_dataPR, f == names.fam[i])
#    phylo <- ape::makeNodeLabel(phylo, "u", nodeList = list(Fam_name = set$s))
#  
#    phylo$node.label[which(phylo$node.label ==
#                             "Fam_name") ] <- paste(set$f[1])
#  }
#  
#  pos.node <- unlist(lapply(names.fam, function(x){
#    which(phylo$node.label == x) + length(phylo$tip.label)
#  }))
#  
#  df.phylo <- data.frame(Fam.names = names.fam,
#                         node.number = pos.node)
#  
#  plot.base <- ggtree(phylo) + theme_tree2()
#  plot1 <- revts(plot.base) + scale_x_continuous(labels=abs)
#  
#  
#  PR.PG <- plot1 + geom_hilight(data = df.phylo, aes(node = node.number, fill = Fam.names),
#                        alpha = .6) +
#    scale_fill_viridis(discrete = T, name = "Family names")

