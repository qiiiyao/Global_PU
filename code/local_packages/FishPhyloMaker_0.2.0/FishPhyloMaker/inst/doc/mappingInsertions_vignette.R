## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FishPhyloMaker)
data("spp_afrotropic")

## ----formatData, echo=TRUE, eval=FALSE----------------------------------------
#  taxon_data <- FishTaxaMaker(data = spp_afrotropic, allow.manual.insert = TRUE)

## ----makingPhylo, echo=TRUE, eval=FALSE---------------------------------------
#  phylo_fish_Afrotropics <- FishPhyloMaker(data = taxon_data$Taxon_data_FishPhyloMaker,
#                                           return.insertions = TRUE,
#                                           insert.base.node = TRUE,
#                                           progress.bar = TRUE)

## ----familyNames, echo=TRUE, eval=FALSE---------------------------------------
#  library(phytools)
#  library(ggtree)
#  library(ggplot2)
#  insertions_org <- phylo_fish_Afrotropics$Insertions_data[match(tree$tip.label, phylo_fish_Afrotropics$Insertions_data$s), ]
#  p.base <- ggtree(tree, layout = "circular", size = .3)  %<+% insertions_org +
#    geom_treescale(x = 0, width = 20, linesize = .5, color = "blue",
#                   fontsize = 0) + #  plot the scale bar
#    annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar
#  
#  p.full <- p.base +
#    geom_tippoint(aes(color = insertions),
#                  size = .5, alpha = .8) +
#    theme(legend.position = "bottom") +
#    guides(color = guide_legend(override.aes = list(size = 2)))  +
#    scale_color_viridis_d(name = NULL, na.translate = F,
#                          labels = c("Congeneric F", "Congeneric", "Family", "Order", "Present"))

