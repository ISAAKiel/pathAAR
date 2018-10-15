#' Reconstructing Empirical Pathes
#' 
#' description of function
#' 
#' @title repath
#' 
#' @param df
#' @param sgdf
#' 
#' @param x numeric, indicating column number of x coordinates
#' @param y numeric, indicating column number of y coordinates
#'  
#' 
#' @return a list containing a graph object of classes tidygraph ("tbl_graph") resp. igraph ("igraph") and the cultural distance matrix.
#'         NOTE: The output igraph object contains zeros for display purpose. The cultural distance matrix (cult_dist_matr) is correct, containing NA for missing values.  
#'
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#' 
#'
#' @export 

repath <- function(df, sgdf, x = 1, y = 2, rw){
  
  pppm <- spatstat::ppp(df[,1], df[,2], 
                         window = spatstat::owin(
                           xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
                           yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
                           unitname="m"))
  f_sd1 <- sd1gen(pppm)
  num <- length(sgdf@data$v)
  base_kde <- makestatkde(pppm, f_sd1[[1]], sgdf, df, x = 1, y = 2, num=num)
  iter  <- 2
  tresh  <- 0.05  # treshold, to delete phaths in areas with really low density values (KDE), meaning calculation artefacts 
  f1     <- 0.2      # factor defining the minimum border of dynamic kernel (raster width) f1*mean(nn)  ## 0.2
  f2     <- 0.4      # factor defining the maximum border of dynamic kernel f2*mean(nn)
  f3     <- 0.5      # minimal intensity of Kernel
  f4     <- 1        # maximal intensity of Kernel
  s      <- -0.3     # Kernelparameter: incline starting from ponit 1
  mwin   <- 9        # Mowing-window-size for ridge detection (4,9,16)
  nn <- f_sd1[[2]]
  for(i  in 1:iter) {
    
    d2 <- min(base_kde@data$v)
    d1 <- max(base_kde@data$v)
    a  <- -(d2*f1-d1*f2)/(d1-d2)
    b  <- (f1-f2)/(d1-d2)
    factor <- function(x){a+b*x} # linear function to scale density values 
    fsd2 <- factor(base_kde@data$v) 
    sd2 <- fsd2*nn  # scaled density values multiplied by nearest neighbourhood distance
    a  <- -(d1*f3-d2*f4)/(d2-d1)
    b  <- (f3-f4)/(d2-d1)
    factor_i <- function(x){a+b*x} # linear function to scale density values 
    fint <- factor_i(base_kde@data$v)
    
    dyn_kde <- sgdf
    
    for (i in seq(along=dyn_kde@data$v)){
      sdi <- sd2[i]
      finti <- fint[i]
      kerpar <- kernel.par(sdi,s,finti)
      pdist <- cbind(sp::coordinates(dyn_kde)[i,1],sp::coordinates(dyn_kde)[i,2],pppm$x[],pppm$y[])
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

#' Calculating Static KDE
#' 
#' description of function
#' 
#' @title makestatkde
#' 
#' @param pppm PointPattern object of Spatial Data
#' @param f_sd1 factor defining size of the first kernel, which generate the stucture of dynamic kernel
#' @param sgdf grid object to store results
#' @param df data frame, containing coordinates of monuments
#' @param x numeric, indicating column number of x coordinates
#' @param y numeric, indicating column number of y coordinates
#'  
#' 
#' @return SpatialGrid Object  
#'
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#' # Crating Test Data Randomly
#' testmatrix <- data.frame(x = abs(rnorm(100)*50), y = abs(rnorm(100)*50))
#' # Setting geographical frame
#' xmin    <- 0
#' xmax    <- max(testmatrix$x)
#' ymin    <- 0
#' ymax    <- max(testmatrix$y)
#' ext_ai <- extent(xmin, xmax, ymin, ymax)
#' 
#' sv <- (xmax-xmin)/(ymax-ymin)
#' rw     <- 10   # width of raster defined in m
#' 
#' rows  <- round((ymax-ymin)/rw, 0) + 1                                    
#' colums <- round((xmax-xmin)/rw, 0) + 1                                      
#' v <- cbind(1:(colums*rows))                                              
#' df <- data.frame(v)                                                         
#' gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
#' sgdf    <- sp::SpatialGridDataFrame(gt, df)
#' 
#' pppm <- spatstat::ppp(testmatrix[,1], testmatrix[,2], 
#'                        window = spatstat::owin(
#'                        xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
#'                        yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
#'                        unitname="m"))
#'                        
#' base_kde <- makestatkde(pppm, sgdf=sgdf, df=testmatrix, x=1, y=2, num=length(df[,1]))
#'
#' @export

makestatkde <- function(pppm, f_sd1=4, sgdf, df, x = 1, y = 2, num){
  sd1 <- sd1gen(pppm)
  sd1 <- sd1[[1]]
  stest <- 0
  if(stest != sd1) {
    base_kde <- sgdf
    for (i in 1:num){
      pdist <- cbind(sp::coordinates(base_kde)[i,1],
                     sp::coordinates(base_kde)[i,2],df[,x],df[,y])
      base_kde@data$v[i]<- sum(apply(pdist,1,gau1,sd=sd1))
    }
  }
  stest <- sd1
  return(base_kde)
}

#' Factor defining size of first Kernel
#' 
#' factor defining size of the first kernel, which generate the stucture of dynamic kernel
#' 
#' @title sd1gen
#' 
#' @param pppm PointPattern Object of Data
#' @param f_sd1 factor scaling the mean of nearest neighbourhoud 
#'  
#' 
#' @return list containg two numeric  
#'
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#' # Crating Test Data Randomly
#' testmatrix <- data.frame(x = abs(rnorm(100)*50), y = abs(rnorm(100)*50))
#' # Setting geographical frame
#' xmin    <- 0
#' xmax    <- max(testmatrix$x)
#' ymin    <- 0
#' ymax    <- max(testmatrix$y)
#' ext_ai <- extent(xmin, xmax, ymin, ymax)
#' 
#' sv <- (xmax-xmin)/(ymax-ymin)
#' rw     <- 10   # width of raster defined in m
#' 
#' rows  <- round((ymax-ymin)/rw, 0) + 1                                    
#' colums <- round((xmax-xmin)/rw, 0) + 1                                      
#' v <- cbind(1:(colums*rows))                                              
#' df <- data.frame(v)                                                         
#' gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
#' sgdf    <- sp::SpatialGridDataFrame(gt, df)
#' 
#' pppm <- spatstat::ppp(testmatrix[,1], testmatrix[,2], 
#'                        window = spatstat::owin(
#'                        xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
#'                        yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
#'                        unitname="m"))
#'                        
#' f_sd1 <- sd1gen(pppm)
#'
#' @export

sd1gen <- function(pppm, f_sd1 = 4) {
  sd1  <- f_sd1*mean(nndist(pppm))
  nn <- mean(nndist(pppm))
  sd1 <- list(sd1, nn)
  return(sd1)
}


#' Factor defining size of first Kernel
#' 
#' factor defining size of the first kernel, which generate the stucture of dynamic kernel
#' 
#' @title gau1
#' 
#' @param x numeric, vector 
#' @param sd numeric, standard deviation
#'  
#' 
#' @return density value of normal distribution
#'
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#'  
#'  x <- c(2,4,1,5,7,8)
#'  sd <- mean(x)
#'  gau1(x, sd)
#'
#' @export

gau1   <- function(x, sd){
  stats::dnorm(edist(x), mean=0, sd=sd)
  } 

#' Euclidian Distance in multidimensional space
#' 
#' Calculates the euclidian distance in a multidimensional space
#' 
#' @title edist
#' 
#' @param x numeric
#'
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#'  
#'  x <- c(2,4,1,5,7,8)
#'  edist(x)
#' 

edist  <- function(x){
  sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))
} 


#' kernel.par
#' 
#' Description
#' 
#' @title kernel.par
#' 
#' @param xp numeric, scaled density values multiplied by nearest neighbourhood distance
#' @param s numeric, Kernelparameter: incline starting from point 1
#' @param int numeric, scaled density values
#'
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#'  
#'  xp <- 6
#'  s <- -0.3
#'  int <- 7
#'  kernel.par(xp,s,int)
#'

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



#' KDE Function 
#' 
#' Description
#' 
#' @title kernel1d
#' 
#' @param x 
#' @param kp 
#'
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @examples
#'  
#'  x <- c(2,4,1,5,7,8)
#'  edist(x)
#'


kernel1d <- function(x,kp){
  d <- edist(x)
  xp <- kp[1]*kp[4]
  a  <- kp[2]
  b  <- kp[3]
  c  <- kp[4]
  int  <- kp[5]
  x <- d*c
  if (x <= xp) {y <- cos(x) + a + (0.7-0.7*x/xp)} # 0,7 -> hight of additional kernell peak
  if (d > xp)  {y <- 1/(x + b)}
  y <- y*int
  return(y)
}