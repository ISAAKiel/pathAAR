################################################################################
## Collection of Functions needed to reconstruct empirical paths
## =============================================================================
## Author of Functions: Oliver Nakoinz & Franziska Faupel
## Version: 04
## Date of last changes: 18.09.2018
## Licence Script: 
################################################################################

#===============================================================================
# cost functions
#===============================================================================
# walk <- function(s){(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3 + 93.419 * s^2 + 19.825 * s + 1.64)} # Herzog 2012)
walk_i <- function(s){1/(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 
                           78.199 * s^3 + 93.419 * s^2 + 19.825 * s + 1.64)}
# drive  <- function(s){(1 + (s / 8)^2)}# Herzog 2012)
drive_i <- function(s){(1/(1 + (s / 8)^2))}
# Cost Surface 
hdiff <- function(x){abs(x[2]-x[1])}
#===============================================================================
# Detecting local density centers using a moving window approach 
#===============================================================================
# KDE to detect local density
localMax <- function(sp_g){
  r <- 15000 # sd of moving window in m
  # calculating static KDE
  bb <- sp::bbox(sp_g) # setting bondingbox for KDE
  win <- spatstat::owin(xrange=c(bb[1,1],bb[1,2]), yrange= c(bb[2,1],bb[2,2]), unitname="m")
  ppp_g <- spatstat::ppp(sp_g@coords[,1], sp_g@coords[,2], window=win)
  dens <- density(ppp_g, kernel="gaussian", sigma=1000, dimyx=c(36,56), 
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
    if (sd(ras_dens@data[indxy,1]) == 0)
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
  
  #plot(maxima) # checking results
  rm(sgdf_dens)
  rm(ras_dens)
  return(maxima)
}

#===============================================================================
# Calculating a Delauny triangulation and preparing coordinates for rWalk loops
#===============================================================================

con_rWalk <- function(maxima, win){
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

#===============================================================================
# Theoretical Paths using cost functions definded by Herzog 2012 
#===============================================================================

theoPath_herzog <- function(emp_ai, ras_ai, con, method, theta){
  
  ras_M1 <- emp_ai
  # Loop to calculate rWalk for all connections
  #-------------------------------------------------------------------------------
  for(i in 1:length(con[,1])){
    xmin <- con[i,5] - p 
    xmax <- con[i,6] + p
    ymin <- con[i,7] - p
    ymax <- con[i,8] + p
    ras <- raster::crop(ras_ai, raster::extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
    tran_hdiff<- gdistance::transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
    slope <- gdistance::geoCorrection(tran_hdiff, type="r", scl=TRUE)  
    adj <- raster::adjacent(ras, cells=1:raster::ncell(ras), direction=8)                     # Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
    cost_walki <- slope                                                   # storing Geo-Information in later TransitionLayer Object
    
    # Different Cost Functions
    if(method == "walk_i"){
    cost_walki[adj] <- walk_i(slope[adj])                          # Calculating cost surface using different functions/equation
    }
    
    if(method == "drive_i"){
      cost_walki[adj] <- drive_i(slope[adj])                          # Calculating cost surface using different functions/equation
    }
    
    cost_walki <- gdistance::geoCorrection(cost_walki, type="r", scl=TRUE)                       # conductivity=cost/dist; time=1/conductivity
    from <- data.frame(con[i, c(1,2)])
    to <- data.frame(con[i, c(3,4)])
    sp::coordinates(from)=~xFrom+yFrom; sp::coordinates(to)=~xTo+yTo; 
    random <- gdistance::passage(cost_walki,from, to, theta = theta, totalNet="total") # rWalk: expectaion for a passage of a cell
    ras_M1 <- raster::mosaic(ras_M1, random, fun=max) # summing up expactation of single rWalk connections
  }
  
  plot(ras_M1)
  
  rm(random); rm(to); rm(from); rm(tran_hdiff); rm(slope); rm(cost_walki)
  
  return(ras_M1)
  
}

#===============================================================================
# Theoretical Paths with additional parameter and 2nd version with factor *1000
#===============================================================================

theoPath_param <- function(emp_ai, ras_ai, ras_view, con, method, theta){
  
  ras_M3 <- emp_ai
  ras_M31000 <- emp_ai
  # Loop to calculate rWalk for all connections
  for(i in 1:length(con[,1])){
    ## Setting area for randomwalk analysis of each connection
    xmin <- con[i,5] - p 
    xmax <- con[i,6] + p
    ymin <- con[i,7] - p
    ymax <- con[i,8] + p
    from <- data.frame(con[i, c(1,2)])
    to <- data.frame(con[i, c(3,4)])
    sp::coordinates(from)=~xFrom+yFrom; sp::coordinates(to)=~xTo+yTo; 
    ## Calculating a cost surface using viewwhed parameters
    ras_v <- raster::crop(ras_view, raster::extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
    adj_v <- raster::adjacent(ras_v, cells=1:raster::ncell(ras_v), directions=8)
    cost_v <- gdistance::transition(ras_v, transitionFunction=max, directions=8, symm=TRUE)
    cost_v <- gdistance::geoCorrection(cost_v, type="r", scl=TRUE)
    cost_v[adj_v] <- max(cost_v[adj_v])
    ## Calculatting a cost surface using Elevation parameters
    ras <- raster::crop(ras_ai, raster::extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
    tran_hdiff<- gdistance::transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
    slope <- gdistance::geoCorrection(tran_hdiff, type="r", scl=TRUE)
    adj <- raster::adjacent(ras, cells=1:raster::ncell(ras), direction=8)# Setting up a adjacent matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
    cost_walki <- slope                                                   # storing Geo-Information in later TransitionLayer Object
    
    # Different Cost Functions
    if(method == "walk_i"){
      cost_walki[adj] <- walk_i(slope[adj])                          # Calculating cost surface using different functions/equation
    }
    
    if(method == "drive_i"){
      cost_walki[adj] <- drive_i(slope[adj])                          # Calculating cost surface using different functions/equation
    }
    
    cost_walki <- gdistance::geoCorrection(cost_walki, type="r", scl=TRUE)                       # conductivity=cost/dist; time=1/conductivity
    ## Adding view-costs to elevation based cost surface
    cost_walkia <- cost_walki + (cost_v*0.1)
    cost_walkib <- cost_walki + (cost_v *100)
    random <- gdistance::passage(cost_walkia,from, to, theta = theta, totalNet="total") # rWalk: expactaion for a passage of a cell
    ras_M3 <- raster::mosaic(ras_M3, random, fun=max) # summing up expactation of single rWalk connections
    random <- gdistance::passage(cost_walkib,from, to, theta = theta*10, totalNet="total") # rWalk: expactaion for a passage of a cell
    ras_M31000 <- raster::mosaic(ras_M31000, random, fun=max) # summing up expactation of single rWalk connections
  }
  
  plot(ras_M3)
  plot(ras_M31000)
  
  rm(random); rm(to); rm(from); rm(tran_hdiff); rm(slope); rm(cost_walkia); rm(cost_walkib); rm(ras_v); rm(adj_v); rm(cost_v)
  
  FAC <- list("param"=ras_M3, "param_1000"=ras_M31000)
  
  return(FAC)
  
}
