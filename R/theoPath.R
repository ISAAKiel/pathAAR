#' Calculating a Delauny triangulation and preparing coordinates for rWalk loops
#' 
#' In absence of known contemporary nodes, like settlements, denser monument
#' locations can be used to reconstruct nodes of infrastructre. This function
#' will use the Kernel Density Approach combined with a moving window, in 
#' order to detect nodes in areas with less monuments (from a global perspective).
#' 
#' @title theo_del
#' 
#' @param maxima
#' @param win 
#'  
#' @return 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 

theo_del <- function(maxima, win){
  library("sp")
  # Delaunay Triangulation to detect possible Connections
  ppp_max <- spatstat::ppp(maxima@coords[,1],maxima@coords[,2], window=win) #ppp are needed for delauny function
  max_del <- spatstat::delaunay(ppp_max)
  sp_con <- as(max_del, "SpatialPolygons") #converting tess object to a polygon
  sl_con <- as(sp_con, "SpatialLines")
  con <- ggplot2::fortify(sp_con) # converting polygon to a data frame
  coord_from <- con[which(con$order<4),]
  coord_from <- coord_from[,c(1,2)]
  coord_to <- con[which(con$order>1),]
  coord_to <- coord_to[,c(1,2)]
  con <- cbind(coord_from, coord_to) # creating new dataframe with coordinates of connections
  colnames(con) <- c("xFrom", "yFrom", "xTo", "yTo")
  # loop to detect doubled connections
  for(i in 1:length(con[,1])){ 
    delete <- (id=1:length(con[,1])) #creating a dummy column
    con <- cbind(con, delete) #adding a dummy column
    a <- which(con$xFrom[i] == con$xTo[c(1:length(con[,1]))] & con$yFrom[i] == con$yTo[c(1:length(con[,1]))] ) #detecting point(i) FROM in TO Columns
    b <- which(con$xTo[i] == con$xFrom[c(1:length(con[,1]))] & con$yTo[i] == con$yFrom[c(1:length(con[,1]))] ) #detecting point(i) TO in FROM Columns
    c <- intersect(a,b) # intersecting, to detect doubeled connections
    con  <- con[!con$delete %in% c ,]
    con <- con[,c(1:4)]
  }
  
  # Preparing data to use in loops (is needed to trim the raster object inside the loop properly)
  con$xmin <- ifelse(con$xFrom < con$xTo, con$xFrom, con$xTo) 
  con$xmax <- ifelse(con$xFrom > con$xTo, con$xFrom, con$xTo) 
  con$ymin <- ifelse(con$yFrom < con$yTo, con$yFrom, con$yTo)
  con$ymax <- ifelse(con$yFrom > con$yTo, con$yFrom, con$yTo) 
  # Deleting long conection, which is represented by others but to long for a rWalk calculation
  for(i in 1:length(con[,1])){
    con$dis[i] <- round(sqrt(((con$xmax[i]-con$xmin[i])^2)+((con$ymax[i]-con$ymin[i])^2)))
  }
  
  delete <- (id=1:length(con[,1])) #creating a dummy column
  con <- cbind(con, delete)
  a <- which(con$dis >=100000)
  con <- con[!con$delete %in% a ,]
  
  CONN <- list(con, sl_con)
  plot(sl_con)
  return(CONN)
  
}



#' Local Centers of Infrastructure from Monument Location
#' 
#' In absence of known contemporary nodes, like settlements, denser monument
#' locations can be used to reconstruct nodes of infrastructre. This function
#' will use the Kernel Density Approach combined with a moving window, in 
#' order to detect nodes in areas with less monuments (from a global perspective).
#' 
#' @title localMax
#' 
#' @param df data frame, containing coordinates of path associated features
#' @param x numeric, indicating column number of x coordinates
#' @param y numeric, indicating column number of y coordinates
#' @param r numeric, size of moving window in meter, default = 5000
#'  
#' @return SpatialPoints 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 

localMax <- function(df, x=1, y=2, r){
  sp_g <- sp::SpatialPoints(cbind(df[,1], df[,2]))
  # calculating static KDE
  bb <- sp::bbox(sp_g) # setting bondingbox for KDE
  win <- spatstat::owin(xrange=c(bb[1,1],bb[1,2]), yrange= c(bb[2,1],bb[2,2]), unitname="m")
  ppp_g <- spatstat::ppp(sp_g@coords[,1], sp_g@coords[,2], window=win)
  dens <- density(ppp_g, kernel="gaussian", sigma=r/2, dimyx=c(36,56), 
                  w=win, edge=TRUE, at="pixels") # calculating static KDE
  
  sgdf_dens <- maptools::as.SpatialGridDataFrame.im(dens)
  ras_dens <- sgdf_dens
  
  ras_dens@data$v[which(is.na(ras_dens@data$v))] <- 0
  m <- max(ras_dens@data$v)
  s <- m / 10
  indmax <- c()
  indplan <- c()
  # moving window approach 
  for (i in seq(along=ras_dens@data$v)) {
    x <- sp::coordinates(ras_dens)[i,1]
    y <- sp::coordinates(ras_dens)[i,2]
    z <- ras_dens@data$v[i]
    indx <- which((sp::coordinates(ras_dens)[,1] > x - r) &
                    + (sp::coordinates(ras_dens)[,1] < x + r))
    indy <- which((sp::coordinates(ras_dens)[,2] > y - r) &
                    + (sp::coordinates(ras_dens)[,2] < y + r))
    indxy <- intersect(indx,indy)
    if (max(ras_dens@data[indxy,1]) == z & z > s)
    {indmax[length(indmax)+1] <- i}
    if (stats::sd(ras_dens@data[indxy,1]) == 0)
    {indplan[length(indplan)+1] <- i}
    rm(indx)
    rm(indy)
    rm(indxy)
  }
  # exracting coordinates of local density maxmia
  mn      <- length(indmax)
  mx      <- sp::coordinates(ras_dens)[indmax,1]
  my      <- sp::coordinates(ras_dens)[indmax,2]
  mx2     <- sp::coordinates(ras_dens)[indplan,1]
  my2     <- sp::coordinates(ras_dens)[indplan,2]
  mz      <- ras_dens@data[indmax,1]
  maxima  <- data.frame(cbind(mx,my,mz))
  sp::coordinates(maxima)=~mx+my

  rm(sgdf_dens)
  rm(ras_dens)
  return(maxima)
}


#' Calculates Slope Values out of Elevation
#' 
#' The `hdiff` function will calculate a cost surfcae, where costs refer to the use of
#' wheeled vehicules with a critical slope value of 8%. Here the 
#' adapted function published by Herzog 2012 is used.
#' 
#' @title hdiff
#' 
#' @param s numeric vector
#'  
#' @return numeric 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 

hdiff <- function(x){
  y <- abs(x[2]-x[1])
  return(y)}

#' Costs for Wheeled Vehicules 
#' 
#' The `drive_i` function will calculate a cost surfcae, where costs refer to the use of
#' wheeled vehicules with a critical slope value of 8%. Here the 
#' adapted function published by Herzog 2012 is used.
#' 
#' Note: As gdistance uses resitor instead of cost, the function creates inverse costs
#' 
#' @title drive_i
#' 
#' @param s numeric 
#'  
#' @return numeric 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 

drive_i <- function(s){
  x <- (1/(1 + (s / 8)^2))
  return(x)
  }

#' Costs for Energy Expenditure for Pedestrians
#' 
#' The `name` function will calculate a cost surfcae, where costs refer to energy
#' expenditure for a pedestrian. Initially published by Llobera & ... Here the 
#' adapted function published by Herzog 2012 is used.
#' 
#' Note: As gdistance uses resitor instead of cost, the function creates inverse costs
#' 
#' @title walk_i
#' 
#' @param s numeric 
#'  
#' @return numeric 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 

walk_i <- function(s){
   x <- 1/(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 
        78.199 * s^3 + 93.419 * s^2 + 19.825 * s + 1.64)
   return(x)
   }
