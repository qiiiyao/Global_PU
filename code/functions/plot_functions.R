
# define a function to make legend of color gradients
legend.func = function(mycolors, mylabels) {
  group = rep("cc", 8)
  condition = letters[1:8]
  value = rep(1, 8)
  df.legend = data.frame(group, condition, value)
  mycolors.corrected = rev(mycolors)
  ggplot(df.legend, aes(fill = condition, y = value, x = group)) +
    geom_bar(position = "stack", stat = "identity", color = "white") +
    scale_fill_manual(values = mycolors.corrected) +
    theme_classic() +
    theme(
      legend.position = "none", aspect.ratio = 0.03,
      axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(size = 7, color = "black")
    ) +
    scale_y_continuous(breaks = 0:8, labels = mylabels) +
    coord_flip() +
    xlab("")
}


## Based on code by Hijmans 
#https://stackoverflow.com/questions/54144269/bivariate-choropleth-map-in-r

# function to make colour matrix
makeCM = function(breaks=10, upperleft, upperright, lowerleft, lowerright) { 
  m = matrix(ncol=breaks, nrow=breaks)
  b = breaks-1
  b = (0:b)/b
  col1 = rgb(colorRamp(c(upperleft, lowerleft))(b), max=255)
  col2 = rgb(colorRamp(c(upperright, lowerright))(b), max=255)
  cm = apply(cbind(col1, col2), 1, function(i) rgb(colorRamp(i)(b), max=255))
  cm[, ncol(cm):1 ]
  
}


makeCM_bilinear = function(breaks = 9,
                           upperleft,   # (-1,+1)
                           upperright,  # (+1,+1)
                           lowerleft,   # (-1,-1)
                           lowerright,  # (+1,-1)
                           lighten = 0.7,   # how much to brighten near center (0–1)
                           power   = 1.5,   # >1 keeps the center lighter for longer
                           neutral_radius = 0 # 0–1 flat light core size (0 = none)
) {
  n = breaks
  x_seq = seq(-1, 1, length.out = n)
  y_seq = seq( 1,-1, length.out = n)
  
  cols_corner = grDevices::col2rgb(c(upperleft, upperright, lowerleft, lowerright))
  M = matrix(NA_character_, nrow = n, ncol = n)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      # [0,1] normalized coordinates
      fx = (x_seq[j] + 1) / 2
      fy = (y_seq[i] + 1) / 2
      
      # bilinear weights
      w_ul = (1 - fx) * fy
      w_ur = fx * fy
      w_ll = (1 - fx) * (1 - fy)
      w_lr = fx * (1 - fy)
      
      # bilinear RGB from corners  (this preserves edge gradients)
      rgb_bi = w_ul*cols_corner[,1] + w_ur*cols_corner[,2] +
        w_ll*cols_corner[,3] + w_lr*cols_corner[,4]
      
      # distance from center (Chebyshev) -> 0 at center, 1 at edges
      r  = max(abs(x_seq[j]), abs(y_seq[i]))
      # central lightening factor: 1 at center -> 0 at edges
      s0 = 1 - r
      s  = pmax(0, (s0 - neutral_radius) / (1 - neutral_radius))^power
      w  = lighten * s
      
      # mix toward white only near the center; edges remain rgb_bi
      rgb_val = (1 - w) * rgb_bi + w * 255
      rgb_val = pmax(pmin(rgb_val, 255), 0)
      
      M[i, j] = grDevices::rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255)
    }
  }
  M
}



# function to plot colour matrix
plotCM_centred = function(cm, xlab="", ylab="", main="", cex.lab=1) {
  #cm = cmat.pe.pu
  #n = cm
  #n = matrix(1:length(cm), nrow=nrow(cm), byrow=TRUE)
  #r = raster(n)
  #cm = cm[, ncol(cm):1 ]
  #image(r, col=cm, axes=FALSE, xlab=xlab, ylab=ylab, main=main, cex.lab=cex.lab)
  
  #cm = cmat.pe.pu
  #xlab="" 
  #ylab="" 
  #main="" 
  #cex.lab=1
  
  #cm = cmat.pe.pu             
  nmat = matrix(1:length(cm), nrow = nrow(cm), byrow = TRUE)
  r = raster::raster(nmat)
  
  cm2 = matrix(as.vector(cm), nrow = nrow(cm), byrow = T)
  
  cols = as.vector(cm2)       

  k = length(cols)          
  brks = seq(0.5, k + 0.5, by = 1)
  
  raster::image(r,
                col   = cols,
                breaks= brks,           
                axes=FALSE, xlab=xlab, ylab=ylab, main=main, cex.lab=cex.lab
  )
  
}

plotCM = function(cm, xlab="", ylab="", main="", cex.lab=1) {
  #cm = cmat.pe.pu
  #n = cm
  n = matrix(1:length(cm), nrow=nrow(cm), byrow=TRUE)
  r = raster(n)
  cm = cm[, ncol(cm):1 ]
  image(r, col=cm, axes=FALSE, xlab=xlab, ylab=ylab, main=main, cex.lab=cex.lab)
  
}



# Function to make bivariate raster, cutting by quantiles
rasterCM = function(x, y, n) {
  q1 = quantile(x, seq(0,1,1/(n)))
  q2 = quantile(y, seq(0,1,1/(n)))
  r1 = cut(x, q1, include.lowest=TRUE)
  r2 = cut(y, q2, include.lowest=TRUE)
  overlay(r1, r2, fun=function(i, j) {
    (j-1) * n + i
  })
}   

rasterCM_centered = function(x, y, n) {
  #x = ED_mammal_delta_raster_moll_mask
  #y = LCBD_mammal_delta_raster_moll_mask
  #n=4
  q1 = quantile(x, seq(0,1,1/(n)), digits = 12)
  q2 = quantile(y, seq(0,1,1/(n)), digits = 12)

  mid = (n/2) + 1  

  q1[mid] = 0
  q2[mid] = 0
  
  q1[1:mid] = seq(min(q1), 0, abs(min(q1)-0)/(n/2))
  q1[mid:length(q1)] = seq(0, max(q1), abs(max(q1)-0)/(n/2))
  
  q2[1:mid] = seq(min(q2), 0, abs(min(q2)-0)/(n/2))
  q2[mid:length(q2)] = seq(0, max(q2), abs(max(q2)-0)/(n/2))
  
  r1 = cut(x, q1, include.lowest=TRUE)
  r2 = cut(y, q2, include.lowest=TRUE)
  
  overlay(r1, r2, fun=function(i, j) {
    (j-1) * n + i
  })
  
}  
