#Lirong Cai
#Email:lirong.cai18@gmail.com

autocov_sphere <- function(z, xy, type="inverse", zero.policy=NULL,
                         style="B", longlat=NULL) {
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  stopifnot(is.vector(z))
  if (type=="one") expo <- 0
  if (type=="inverse") expo <- 1
  if (type=="inverse.squared") expo <- 2
  
  
  if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
  stopifnot(ncol(xy) == 2)
  
  coords <- as.matrix(xy)
  nb_tri <- spdep::tri2nb(coords)
  nb <- graph2nb(soi.graph(nb_tri, coords))

  
  if (any(card(nb) == 0)) warning(paste("With value", nbs,
                                        "some points have no neighbours"))
  
  nbd <- nbdists(nb, xy, longlat=longlat)
  
  if (expo == 0) lw <- nb2listw(nb, style=style, zero.policy=zero.policy)
  else {
    gl <- lapply(nbd, function(x) 1/(x^expo))
    lw <- nb2listw(nb, glist=gl, style=style, zero.policy=zero.policy)
  }
  lag.listw(lw, z, zero.policy=zero.policy) 
  
}
