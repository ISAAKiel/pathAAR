#' Reconstructing Empirical Pathes
#' 
#' description of function
#' 
#' @title repath
#' 
#' @param nodes_x a vector containing metric x coordinates of nodes
#' @param nodes_y a vector containing metric y coordinates of nodes
#' @param nodes_id a vector containing ID for nodes
#' @param features a data.frame containing metric x and y coordinates of features, and feature type. Coordinates are expected to be the first two columns.
#' @param type_col a character string naming the columname containing feature types.
#' @param pre_size numeric, amount of letters, e.g. characters before typenumbers. Defaults to 1.
#' @param method character string, the distance measure to be 
#'   used ("euclidean", "maximum", "manhattan", "canberra", 
#'   "binary" or "minkowski"). Defaults to euclidean distance.
#' @param plotted a Boolean operator defining whether a plot should be created. Defaults to FALSE. Edge widths are scaled by maximum distance values and enlarged by factor 2. 
#' 
#' @return a list containing a graph object of classes tidygraph ("tbl_graph") resp. igraph ("igraph") and the cultural distance matrix.
#'         NOTE: The output igraph object contains zeros for display purpose. The cultural distance matrix (cult_dist_matr) is correct, containing NA for missing values.  
#'
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{}>
#' 
#' @examples
#' 
#'
#' @export 

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
