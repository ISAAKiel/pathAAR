################################################################################
## Collection of Functions needed to reconstruct empirical paths
## =============================================================================
## Author of Functions: Oliver Nakoinz & Franziska Faupel
## Version: 04
## Date of last changes: 18.09.2018
## Licence Script: 
################################################################################

#===============================================================================
# distance functions
#===============================================================================
# euclidian distance
edist  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))} 
# euclidian distance weighted with Gauß distribution of static kernel
gau1   <- function(x, sd){stats::dnorm(edist(x), mean=0, sd=sd)}         

#===============================================================================
# Calculating static Kernel
#===============================================================================
makestatkde <- function(sd1, sgdf, df, x = 1, y = 2){
    stest <- 0
   if(stest != sd1) {
      base_kde <- sgdf
      for (i in seq(along=base_kde@data$v)){
        pdist <- cbind(sp::coordinates(base_kde)[i,1],
                       sp::coordinates(base_kde)[i,2],df[,x],df[,y])
        base_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sd1))
      }
    }
    stest <- sd1
    return(base_kde)
  }

#===============================================================================
# Dynamic KDE
#===============================================================================
# 1. KDE Function 
kernel.par <- function(xp,s,int){
  x1 <- asin(-s)
  y1 <- cos(x1)
  x2 <- (-1/s)^0.5
  y2 <- 1/x2
  a <- y2-y1
  b <- x2-x1
  c <- x1/xp
  return(c(xp,a,b,c,int))
}

# 3. KDE Function 
kernel1d <- function(x,kp){
  d <- edist(x)
  xp <- kp[1]*kp[4]
  a  <- kp[2]
  b  <- kp[3]
  c  <- kp[4]
  int  <- kp[5]
  x <- d*c
  if (x <= xp) {y <- cos(x) + a + (de-de*x/xp)}
  if (d > xp)  {y <- 1/(x + b)}
  y <- y*int
  return(y)
}

#===============================================================================
# Path Reconstruction
#===============================================================================

repath <- function(base_kde, nn, sgdf){
  
  #iter <- 2
  #rw     <- 500   # width of raster defined in m
  #tresh  <- 0.05  # treshold, to delete phaths in areas with really low density values (KDE), meaning calculation artefacts 
  #f_sd1  <- 4     # factor defining size if the first kernel, which generate the stucture of dynamic kernel
  #f1     <- 0.2      # factor defining the minimum border of dynamic kernel (raster width) f1*mean(nn)  ## 0.2
  #f2     <- 0.4      # factor defining the maximum border of dynamic kernel f2*mean(nn)
  #f3     <- 0.5      # MinimalinentitÃÂ¤t des Kernels
  #f4     <- 1        # MaximalinentitÃÂ¤t des Kernels
  #s      <- -0.3     # Kernelparameter: incline starting from ponit 1
  #de     <- 0.7      # hight of additional kernell peak
  #sw     <- 12       # width of picture, cm
  #mwin   <- 9        # Mowing-window-size for ridge detection (4,9,16)
  #xp <- 750     # Kernelparameter: x-wert   von Punkt 1
  
  for(i  in 1:iter) {
    
    d2 <- min(base_kde@data$v)
    d1 <- max(base_kde@data$v)
    a  <- -(d2*f1-d1*f2)/(d1-d2)
    b  <- (f1-f2)/(d1-d2)
    factor <- function(x){a+b*x}
    fsd2 <- factor(base_kde@data$v)
    sd2 <- fsd2*nn
    a  <- -(d1*f3-d2*f4)/(d2-d1)
    b  <- (f3-f4)/(d2-d1)
    factor_i <- function(x){a+b*x}
    fint <- factor_i(base_kde@data$v)
    
    dyn_kde <- sgdf
    
    for (i in seq(along=dyn_kde@data$v)){
      sdi <- sd2[i]
      finti <- fint[i]
      kerpar <- kernel.par(sdi,s,finti)
      pdist <- cbind(sp::coordinates(dyn_kde)[i,1],sp::coordinates(dyn_kde)[i,2],ppp_g$x[],ppp_g$y[])
      dyn_kde@data$v[i]<- sum(apply(pdist,1, kernel1d ,kp=kerpar))
    }
    
    ras_ridges <- sgdf
    ras        <- dyn_kde
    ras_ridges@data$v <- 1
    ras@data$v[is.na(ras@data$v)] <- 10000000000
    
    if (mwin==4){
      for (i in 1:(length(ras@data$v)-colums))  {
        ind <- c(i,i+1,i+colums,i+colums+1)
        ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
        ind_min2 <-ind[ind_min1]
        ras_ridges@data$v[ind_min2] <- 0
      }
    }
    
    if (mwin==9){
      for (i in 1:(length(ras@data$v)-(2*colums)-2))  { 
        ind <- c(i,i+1,i+2,i+colums,i+colums+1,i+colums+2,i+(2*colums),i+(2*colums)+1,i+(2*colums)+2)
        ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
        ind_min2 <-ind[ind_min1]
        ras_ridges@data$v[ind_min2] <- 0
      }
    }
    
    if (mwin==16){
      for (i in 1:(length(ras@data$v)-(3*colums)-3))  { 
        ind <- c(i,i+1,i+2,i+3,i+colums,i+colums+1,i+colums+2,i+colums+3,i+(2*colums),i+(2*colums)+1,i+(2*colums)+2,i+(2*colums)+3,i+(3*colums),i+(3*colums)+1,i+(3*colums)+2,i+(3*colums)+3)
        ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
        ind_min2 <-ind[ind_min1]
        ras_ridges@data$v[ind_min2] <- 0
      }
    }
    ras_ridges_corr <- ras_ridges
    tresh_dens <- tresh*max(ras@data$v)
    ras_ridges_corr@data$v[which (ras@data$v < tresh_dens)] <- 0
    
  } 
  
  list_repath <- list(KDE=dyn_kde, ras_dens_ridges=ras_ridges_corr)
  
  return(list_repath)
}
