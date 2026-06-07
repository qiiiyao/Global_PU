### Set up R environments 
### This function adapted from Baselga_et_al_2025_Ecol._Lett.
#rm(list = ls())
# Package management
# list of packages
# Manually install the PhyloMeasures package
#install.packages("D:/R projects/Global_ED/code/PhyloMeasures_2.1.tar.gz",
#     repos = NULL, type = "source")
requirements = c("PhyloMeasures", "phangorn", "doParallel", "dplyr",
                 "ape", 'phyloregion', 'adespatial', 'betapart', 'geodist')

# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})

setwd("D:/R projects/Global_ED")

############################################################################
### Function to compute the phylogenetic version of Ruggiero's dissimilarity 
### (index beta.rbl in Koleff et al 2023)
### Based on code from the R package phyloregion (Daru et al. 2020) 

### The function takes:
### a presence/absence table (x) 
### and a phylogenetic tree (phy)

### The function returns a matrix of phylogenetic Ruggiero's dissimilarity
### values, in which focal cells are rows, and all other cells are columns.

### The result is a square matrix because the index is not symmetrical!
### Therefore, dissimilarity between Site 1 and Site 2 can be different 
### from the dissimilarity between Site 2 and Site 1.

phylobeta.rug = function (x, phy) 
{
  if (!is(x, "sparseMatrix")) 
    stop("x needs to be a sparse matrix!")
  if (!identical(sort(colnames(x)), sort(phy$tip.label))) 
    stop("Labels of community matrix and tree differ!")
  
  phylo_community <- function(x, phy) {
    el <- numeric(max(phy$edge))
    el[phy$edge[, 2]] <- phy$edge.length
    x <- x[, phy$tip.label]
    anc <- phangorn::Ancestors(phy, seq_along(phy$tip.label))
    anc <-  mapply(c, seq_along(phy$tip.label), anc, SIMPLIFY=FALSE)
    M <- Matrix::sparseMatrix(as.integer(rep(seq_along(anc), 
                                             lengths(anc))), as.integer(unlist(anc)), x = 1L)
    commphylo <- x %*% M
    commphylo@x[commphylo@x > 1e-8] <- 1
    list(Matrix = commphylo, edge.length = el)
  }    
  x <- phylo_community(x, phy)
  pd_tmp <- (x$Matrix %*% x$edge.length)[, 1]
  l <- length(pd_tmp)
  m <- l - 1L
  SHARED <- Matrix::tcrossprod(x$Matrix, x$Matrix %*% Matrix::Diagonal(x = x$edge.length))
  SHARED <- as.dist(as.matrix(SHARED))
  B <- pd_tmp[rep(1:m, m:1)] - SHARED
  C <- pd_tmp[sequence(m:1) + rep(1:m, m:1)] - SHARED
  not.shared <- as.matrix(B)
  not.shared[lower.tri(not.shared)] <- C
  res <- not.shared / (as.matrix(SHARED) + not.shared)
  diag(res) <- 0
  res
}
