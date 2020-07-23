#' Theoretical Paths using cost functions defined by Herzog 2012
#' 
#' If there are no actual parts of a path network known, the `theoPath_herzog` function can be used to reconstruct pathways based on randomised shortest paths connecting known regions with higher densities of sites, e.g. monuments. An underlying cost surface is created by using the cost functions defined by I. Herzog (2012) either for walking or driving. This function is a useful step in the evaluation of reconstructed paths. 
#'
#' 
#' @title theoPath_herzog
#' 
#' @param ras_ai RasterLayer, raster with elevation values
#' @param con data.frame, connections of a Delaunay triangulation as a result of function theo_del or your own method
#' @param method chr, either "walk_i" for pedestrians or "drive_i" for vehicles. For further informations look up the respective functions
#' @param theta numeric, parameter controls randomisation of walk. Lower values equal more exploration of the walker around the shortest path, while if theta approaches zero the walk becomes totally random 
#' @param p numeric, buffer zone around rWalk rasters, used in loop
#' @param type chr, either "c" (default) for least-cost distances or "r" for random walks. As stated by J. van Etten, there is no analytical way as of now to decide for intermediate values of theta which type should be choosen. For further informations see ?gdistance::geoCorrection
#'  
#' @return raster object, with values of summed up expectations of single rWalk connections
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' @author Hendrik Raese <\email{h.raese@@roots.uni-kiel.de}>
#' 
#' @examples 
#' 
#' set.seed(123)
#' 
#' # Creating random test data 
#' testmatrix <- data.frame(
#'  x = abs(runif(100, 1000, 1500)), 
#'  y = abs(runif(100, 1000, 1500)))
#'
#' # Calculate local centres of infrastructure from monument location
#' pd=25
#' sw=10
#' r=100
#' maxima <- localMax(df=testmatrix, r=r, sw=sw, pd=pd)
#' 
#' # Setting geographical frame
#' xmin    <- 0
#' xmax    <- max(testmatrix$x)
#' ymin    <- 0
#' ymax    <- max(testmatrix$y)
#' ext_ai <- raster::extent(xmin, xmax, ymin, ymax)
#' 
#' # Coordinates used to set frame corner for definition of the aspect ratio
#' sv <- (xmax-xmin)/(ymax-ymin)
#' rw     <- 5   # width of raster defined in m
#' 
#' # Definition of frame expansion and defining frame                              
#' rows  <- round((ymax-ymin)/rw, 0) + 1                                  
#' colums <- round((xmax-xmin)/rw, 0) + 1                                     
#' v <- cbind(1:(colums*rows))                                              
#' df <- data.frame(v)                                                         
#' gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
#' sgdf    <- sp::SpatialGridDataFrame(gt, df)
#' 
#' # Initialising observation window for theoretical connections
#' win <- spatstat::owin(c(xmin, xmax),c(ymin, ymax))
#' 
#' # calculating theoretical connections via delaunay triangulation
#' theo_con <- theo_del(maxima,win)
#' 
#' # Setting up an artificial elevation map with random values
#' emap <- sgdf
#' emap@data$v <- sample((50:56), length(emap@data$v), replace=T)
#' ras_emap <- raster(emap)
#' 
#' # Run the function with chosen parameters for method, theta and p
#' theo_run <- theoPath_herzog(ras_ai=ras_emap,
#'                             method="drive_i",
#'                             theo_con[[1]], 
#'                             theta=0.001, 
#'                             p=5, 
#'                             type="r")
#'
#'@export



theoPath_herzog <- function(ras_ai, con, method, theta, p, type="c"){
  
  # Initialise empty raster with same extent and resolution as the elevation raster
  ras_M1 <- ras_ai
  ras_M1@data@values <- NA
  
  # Loop to calculate rWalk for all connections
  #-------------------------------------------------------------------------------
  for(i in 1:length(con[,1])){
    xmin <- con[i,5] - p 
    xmax <- con[i,6] + p
    ymin <- con[i,7] - p
    ymax <- con[i,8] + p
    ras <- raster::crop(ras_ai, raster::extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
    tran_hdiff<- gdistance::transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE) # Attention: If directions > 4 is selected, there will be a need for a geographical correction because the orthogonal and diagonal distances between all connected cells are not the same (1 vs. sqr(2))! 
    slope <- gdistance::geoCorrection(tran_hdiff, type=type, scl=TRUE) # Geographic correction is needed because the TransitionLayer was created with >4 directions. For random walks there should be an North-South and an East-West correction which calls for type="r" correction, while least-cost distances on a LonLat basis can be corrected with type="c". 
    adj <- raster::adjacent(ras, cells=1:raster::ncell(ras), direction=8) # Setting up an adjacency matrix using 8 directions to construct a transition layer to store cost values calculated in the next step. 
    cost_method <- slope # storing Geo-Information in later TransitionLayer Object
    
    # Different Cost Functions
    if(method == "walk_i"){
      cost_method[adj] <- walk_i(slope[adj]) # Calculating cost surface using different functions/equations
    }
    
    if(method == "drive_i"){
      cost_method[adj] <- drive_i(slope[adj]) # Calculating cost surface using different functions/equations
    }
    
    cost_method <- gdistance::geoCorrection(cost_method, type=type, scl=TRUE) # conductivity=cost/dist; time=1/conductivity
    from <- data.frame(con[i, c(1,2)])
    to <- data.frame(con[i, c(3,4)])
    sp::coordinates(from)=~xFrom+yFrom; sp::coordinates(to)=~xTo+yTo; 
    random <- gdistance::passage(cost_method,from, to, theta = theta, totalNet="total") # rWalk: expectation for a passage of a cell
    ras_M1 <- raster::mosaic(ras_M1, random, fun=max) # summing up expectation of single rWalk connections
  }
  
  raster::plot(ras_M1)
  
  rm(random); rm(to); rm(from); rm(tran_hdiff); rm(slope); rm(cost_method)
  
  return(ras_M1)
  
}



#' Theoretical Paths using cost functions defined by Herzog 2012 with additional cost parameter
#' 
#' If there are no actual parts of a path network known, the `theoPath_herzog_param` function can be used to reconstruct pathways based on randomised shortest paths connecting known regions with higher densities of sites, e.g. monuments. An underlying cost surface is created by using the cost functions defined by I. Herzog (2012) either for walking or driving. By supplying an additional parameter (visibility, friction costs, etc.), hypotheses on the influence of different variables can be tested. This function is a useful step in the evaluation of reconstructed paths. 
#'
#' 
#' @title theoPath_param
#' 
#' @param ras_ai RasterLayer, raster with elevation values
#' @param ras_para RasterLayer, raster with additional values (visibility, friction costs, etc.)
#' @param con data.frame, connections of a Delaunay triangulation as a result of function theo_del or your own method
#' @param method chr, either "walk_i" for pedestrians or "drive_i" for vehicles. For further informations look up the respective functions
#' @param theta numeric, parameter controls randomisation of walk. Lower values equal more exploration of the walker around the shortest path, while if theta approaches zero the walk becomes totally random 
#' @param p numeric, buffer zone around rWalk rasters, used in loop
#' @param type chr, either "c" (default) for least-cost distances or "r" for random walks. As stated by J. van Etten, there is no analytical way as of now to decide for intermediate values of theta which type should be choosed. For further informations see ?gdistance::geoCorrection
#'  
#' @return List, two RasterLayer with values of summed up expectations of single rWalk connections. The item param with 0.1x the influence of the supplied parameter and param_1000 with 100x the influence of the supplied parameter and 10x the value of theta.
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' @author Hendrik Raese <\email{h.raese@@roots.uni-kiel.de}>
#' 
#' @examples
#' 
#' set.seed(123)
#' 
#' # Creating random test data 
#' testmatrix <- data.frame(
#'  x = abs(runif(100, 1000, 1500)), 
#'  y = abs(runif(100, 1000, 1500)))
#'
#' # Calculate local centres of infrastructure from monument location
#' pd=25
#' sw=10
#' r=100
#' maxima <- localMax(df=testmatrix, r=r, sw=sw, pd=pd)
#' 
#' # Setting geographical frame
#' xmin    <- 0
#' xmax    <- max(testmatrix$x)
#' ymin    <- 0
#' ymax    <- max(testmatrix$y)
#' ext_ai <- raster::extent(xmin, xmax, ymin, ymax)
#' 
#' # Coordinates used to set frame corner for definition of the aspect ratio
#' sv <- (xmax-xmin)/(ymax-ymin)
#' rw     <- 5   # width of raster defined in m
#' 
#' # Definition of frame expansion and defining frame                              
# 'rows  <- round((ymax-ymin)/rw, 0) + 1                                    
#' colums <- round((xmax-xmin)/rw, 0) + 1                                     
#' v <- cbind(1:(colums*rows))                                              
#' df <- data.frame(v)                                                         
#' gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
#' sgdf    <- sp::SpatialGridDataFrame(gt, df)
#' 
#' # Initialising observation window for theoretical connections
#' win <- owin(c(xmin, xmax),c(ymin, ymax))
#' 
#' # calculating theoretical connections via delaunay triangulation
#' theo_con <- theo_del(maxima,win)
#' 
#' # Setting up an artificial elevation map with random values
#' emap <- sgdf
#' emap@data$v <- sample((50:56), length(emap@data$v), replace=T)
#' ras_emap <- raster(emap)
#' 
#' # Friction Raster for preferring lowlands:
#' para <- ras_sgdf
#' para@data@values[which(para@data@values <0 )] <- 0
#' para@data@values <- para@data@values/ max(para@data@values)
#' raster::plot(para)
#' 
#' # Run the function with chosen parameters for method, theta and p
#' theo_run <- theoPath_param(ras_ai=ras_emap, 
#'                            ras_para=para, 
#'                            method="drive_i",
#'                            theo_con[[1]], 
#'                            theta=0.001, 
#'                            p=5, 
#'                            type="r")
#'
#'@export

theoPath_param <- function(ras_ai, ras_para, con, method, theta, p, type="c"){
  
  # Initialise empty raster with same extent and resolution as the elevation raster
  ras_M3 <- ras_ai
  ras_M31000 <- ras_M3
  ras_M3@data@values <- NA
  ras_M31000@data@values <- NA
  
  # Loop to calculate rWalk for all connections
  for(i in 1:length(con[,1])){
    ## Setting area for random walk analysis of each connection
    xmin <- con[i,5] - p 
    xmax <- con[i,6] + p
    ymin <- con[i,7] - p
    ymax <- con[i,8] + p
    from <- data.frame(con[i, c(1,2)])
    to <- data.frame(con[i, c(3,4)])
    sp::coordinates(from)=~xFrom+yFrom; sp::coordinates(to)=~xTo+yTo; 
    ## Calculating a cost surface using provided additional parameter
    ras_par <- raster::crop(ras_para, raster::extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
    adj_par <- raster::adjacent(ras_par, cells=1:raster::ncell(ras_par), directions=8)
    cost_par <- gdistance::transition(ras_par, transitionFunction=max, directions=8, symm=TRUE)
    cost_par <- gdistance::geoCorrection(cost_par, type=type, scl=TRUE)
    cost_par[adj_par] <- max(cost_par[adj_par])
    ## Calculating a cost surface using Elevation parameters
    ras <- raster::crop(ras_ai, raster::extent(xmin, xmax, ymin,ymax), filename= "corssings" , snap='near', overwrite=TRUE)
    tran_hdiff<- gdistance::transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)
    slope <- gdistance::geoCorrection(tran_hdiff, type=type, scl=TRUE)
    adj <- raster::adjacent(ras, cells=1:raster::ncell(ras), direction=8) # Setting up an adjacency matrix using 8 directions to construct a transition Layer to store cost values calculated in the next step
    cost_method <- slope # storing Geo-Information in later TransitionLayer Object
    
    # Different Cost Functions
    if(method == "walk_i"){
      cost_method[adj] <- walk_i(slope[adj]) # Calculating cost surface using different functions/equations
    }
    
    if(method == "drive_i"){
      cost_method[adj] <- drive_i(slope[adj]) # Calculating cost surface using different functions/equations
    }
    
    cost_method <- gdistance::geoCorrection(cost_method, type=type, scl=TRUE) # conductivity=cost/dist; time=1/conductivity
    ## Adding costs of the additional parameter to the elevation based cost surface
    cost_methoda <- cost_method + (cost_par*0.1)
    cost_methodb <- cost_method + (cost_par *100)
    random <- gdistance::passage(cost_methoda,from, to, theta = theta, totalNet="total") # rWalk: expectation for a passage of a cell
    ras_M3 <- raster::mosaic(ras_M3, random, fun=max) # summing up expectation of single rWalk connections
    random <- gdistance::passage(cost_methodb,from, to, theta = theta*10, totalNet="total") # rWalk: expectation for a passage of a cell
    ras_M31000 <- raster::mosaic(ras_M31000, random, fun=max) # summing up expectation of single rWalk connections
  }
  
  raster::plot(ras_M3)
  raster::plot(ras_M31000)
  
  rm(random); rm(to); rm(from); rm(tran_hdiff); rm(slope); rm(cost_methoda); rm(cost_methodb); rm(ras_par); rm(adj_par); rm(cost_par)
  
  FAC <- list("param"=ras_M3, "param_1000"=ras_M31000)
  
  return(FAC)
  
}



#' Local Centres of Infrastructure from Monument Location
#' 
#' In absence of known contemporary nodes, like settlements, the density of monument
#' locations can be used to reconstruct nodes of infrastructure. This function
#' will use the Kernel Density Approach combined with a moving window, in 
#' order to detect nodes in areas with less monuments (from a global perspective).
#' 
#' @title localMax
#' 
#' @param df data frame, containing coordinates of path associated features
#' @param sw numeric, value used for defining threshold for defining a local maximum. A higher value means more local maxima are stored, default = 10. It's suggested not to change this parameter and instead work with the size of the moving window!
#' @param pd numeric, dimension of pixels for intermediate raster image of the used density function in meter, default = 500
#' @param r numeric, sets the size of moving window and sigma in the used density function in meter, default = 5000
#'  
#' @return SpatialPointsDataFrame with coordinates of local maxima
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' @author Hendrik Raese <\email{h.raese@@roots.uni-kiel.de}>
#' 
#' @examples 
#' set.seed(123)
#'
#' # Creating random test data 
#' testmatrix <- data.frame(
#'  x = abs(runif(100, 1000, 1500)), 
#'  y = abs(runif(100, 1000, 1500)))
#'
#' pd=25
#' sw=10
#' r=100
#' maxima <- localMax(df=testmatrix, r=r, sw=sw, pd=pd)
#' 
#' # Plot the result
#' maxima <- data.frame(maxima)
#' 
#' ggplot2::ggplot() +
#'  ggplot2::geom_point(data= testmatrix, ggplot2::aes(x=x, y=y)) +
#'  ggplot2::geom_point(data=maxima, 
#'                      ggplot2::aes(x=mx, y=my, color="red"), 
#'                      shape=15, 
#'                      size=3, 
#'                      show.legend = FALSE)
#' 
#' @export

localMax <- function(df, r=5000, sw=10, pd=500){
  sp_g <- sp::SpatialPoints(cbind(df[,1], df[,2]))
  # calculating static KDE
  bb <- sp::bbox(sp_g) # setting spatial bounding box for KDE
  win <- spatstat::owin(xrange=c(bb[1,1],bb[1,2]), yrange= c(bb[2,1],bb[2,2]), unitname="m") # calculating observation window
  ppp_g <- spatstat::ppp(sp_g@coords[,1], sp_g@coords[,2], window=win)
  dens <- stats::density(ppp_g, kernel="gaussian", sigma=r/2, dimyx=c((round((win$yrange[2]-win$yrange[1])/pd)),(round((win$xrange[2]-win$xrange[1])/pd))), 
                  w=win, edge=TRUE, at="pixels") # calculating static KDE
  
  sgdf_dens <- maptools::as.SpatialGridDataFrame.im(dens)
  ras_dens <- sgdf_dens
  
  ras_dens@data$v[which(is.na(ras_dens@data$v))] <- 0
  m <- max(ras_dens@data$v)
  s <- m / sw 
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
    if (stats::sd(ras_dens@data[indxy,1]) == 0) #An dieser Stelle teilweise rausgehauen - stattdessen ras_dens@data$v ? Funktioniert aber so, wie es ist 
      if (stats::sd(ras_dens@data$v) == 0)
      {indplan[length(indplan)+1] <- i}
    rm(indx)
    rm(indy)
    rm(indxy)
  }
  # extracting coordinates of local density maxmia
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


#' Calculating a Delaunay triangulation and preparing coordinates for rWalk loops
#' 
#' In absence of known contemporary nodes, like settlements, denser monument
#' locations can be used to reconstruct nodes of infrastructure. This function
#' will use the Kernel Density Approach combined with a moving window, in 
#' order to detect nodes in areas with less monuments (from a global perspective).
#' 
#' @title theo_del
#' 
#' @param maxima SpatialPoints object
#' @param win owin object
#'  
#' @return list with nodes and edges of centers of infrastructure
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 
#' @export

theo_del <- function(maxima, win){
  # Delaunay Triangulation to detect possible Connections
  ppp_max <- spatstat::ppp(maxima@coords[,1],maxima@coords[,2], window=win) #ppp are needed for Delauny function
  max_del <- spatstat::delaunay(ppp_max)
  sp_con <- methods::as(max_del, "SpatialPolygons") #converting tess object to a polygon
  sl_con <- methods::as(sp_con, "SpatialLines")
  
  # according to the documentation the function fortify may be deprecated in the future, maybe it should be replaced with the appropriate function from the 'broom' package as described by ggplot2 documentation. BUT: In the broom package documentation it is stated that development of sp tidier is halted and may be deprecated in the future. They propose changing to sf instead of sp. That would effect the whole package here...  
  con <- ggplot2::fortify(sp_con) # converting polygon to a data frame 
  #con <- broom::tidy(sp_con) # works with warnings and 'group' is chr instead of factor but as only order is used, that is not a problem.
  #con <- sf::st_as_sf(sp_con)  
  
  coord_from <- con[which(con$order<4),]
  coord_from <- coord_from[,c(1,2)]
  coord_to <- con[which(con$order>1),]
  coord_to <- coord_to[,c(1,2)]
  con <- cbind(coord_from, coord_to) # creating new data frame with coordinates of connections
  colnames(con) <- c("xFrom", "yFrom", "xTo", "yTo")
  # loop to detect doubled connections
  for(i in 1:length(con[,1])){ 
    delete <- (id=1:length(con[,1])) #creating a dummy column
    con <- cbind(con, delete) #adding a dummy column
    a <- which(con$xFrom[i] == con$xTo[c(1:length(con[,1]))] & con$yFrom[i] == con$yTo[c(1:length(con[,1]))] ) #detecting point(i) FROM in TO Columns
    b <- which(con$xTo[i] == con$xFrom[c(1:length(con[,1]))] & con$yTo[i] == con$yFrom[c(1:length(con[,1]))] ) #detecting point(i) TO in FROM Columns
    c <- intersect(a,b) # intersecting, to detect doubled connections
    con  <- con[!con$delete %in% c ,]
    con <- con[,c(1:4)]
  }
  
  # Preparing data to use in loops (is needed to trim the raster object inside the loop properly)
  con$xmin <- ifelse(con$xFrom < con$xTo, con$xFrom, con$xTo) 
  con$xmax <- ifelse(con$xFrom > con$xTo, con$xFrom, con$xTo) 
  con$ymin <- ifelse(con$yFrom < con$yTo, con$yFrom, con$yTo)
  con$ymax <- ifelse(con$yFrom > con$yTo, con$yFrom, con$yTo) 
  # Deleting long connection, which is represented by others but to long for a rWalk calculation
  for(i in 1:length(con[,1])){
    con$dis[i] <- round(sqrt(((con$xmax[i]-con$xmin[i])^2)+((con$ymax[i]-con$ymin[i])^2)))
  }
  
  delete <- (id=1:length(con[,1])) #creating a dummy column
  con <- cbind(con, delete)
  a <- which(con$dis >=100000)
  con <- con[!con$delete %in% a ,]
  
  CONN <- list(con, sl_con)
  graphics::plot(sl_con)
  return(CONN)
  
}



#' Calculates Slope Values out of Elevation
#' 
#' The `hdiff` function will calculate slope based on the numeric value of the elevation
#' in adjacent cells
#' 
#' @title hdiff
#' 
#' @param x numeric vector
#'  
#' @return numeric 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' 

hdiff <- function(x){
  y <- abs(x[2]-x[1])
  return(y)}

#' Costs for Wheeled Vehicles 
#' 
#' The `drive_i` function will calculate a cost surface, where costs refer to the use of
#' wheeled vehicles with a critical slope value of 8%. Here the 
#' adapted function published by Herzog 2012 based on Llobera 2007 is used.
#' 
#' Note: As gdistance uses resistor instead of cost, the function creates inverse costs
#' 
#' @title drive_i
#' 
#' @param s numeric 
#'  
#' @return numeric 
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' @author Hendrik Raese <\email{h.raese@@ufg.uni-kiel.de}>
#' 

drive_i <- function(s){
  x <- (1/(1 + (s / (8/100))^2))
  return(x)
}

#' Costs for Energy Expenditure for Pedestrians
#' 
#' The `walk_i` function will calculate a cost surface, where costs refer to energy
#' expenditure for a pedestrian. Initially published by Llobera & ... Here the 
#' adapted function published by Herzog 2012 is used.
#' 
#' Note: As gdistance uses resistor instead of cost, the function creates inverse costs
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
