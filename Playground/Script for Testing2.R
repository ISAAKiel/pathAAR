library(devtools)

install_github("ISAAKiel/pathAAR")
library(pathAAR)

 set.seed(123)
 
 # Creating random test data 
 testmatrix <- data.frame(
  x = abs(runif(100, 1000, 1500)), 
  y = abs(runif(100, 1000, 1500)))

 # Calculate local centres of infrastructure from monument location
 pd=25
 sw=10
 r=100
 maxima <- localMax(df=testmatrix, r=r, sw=sw, pd=pd)
 
 # Setting geographical frame
 xmin    <- min(testmatrix$x)
 xmax    <- max(testmatrix$x)
 ymin    <- min(testmatrix$y)
 ymax    <- max(testmatrix$y)
 ext_ai <- raster::extent(xmin, xmax, ymin, ymax)
 
 # Coordinates used to set frame corner for definition of the aspect ratio
 sv <- (xmax-xmin)/(ymax-ymin)
 rw     <- 5   # width of raster defined in m
 
 # Definition of frame expansion and defining frame                              
 rows  <- round((ymax-ymin)/rw, 0) + 1                                  
 colums <- round((xmax-xmin)/rw, 0) + 1                                     
 v <- cbind(1:(colums*rows))                                              
 df <- data.frame(v)                                                         
 gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
 sgdf    <- sp::SpatialGridDataFrame(gt, df)
 
 # Initialising observation window for theoretical connections
 win <- spatstat::owin(c(xmin, xmax),c(ymin, ymax))
 
 # calculating theoretical connections via delaunay triangulation
 theo_con <- theo_del(maxima,win)
 
 
 # Setting up an artificial elevation map with random values
 emap <- sgdf
 emap@data$v <- sample((50:56), length(emap@data$v), replace=TRUE)
 ras_emap <- raster::raster(emap)
 
 theo_run <- theoPath_herzog(ras_ai=ras_emap,method="drive_i",theo_con[[1]], theta=0.001, p=5, type="r")
 
  # Friction Raster for preferring lowlands:
  para <- ras_emap
  para@data@values[which(para@data@values <0 )] <- 0
  para@data@values <- para@data@values/ max(para@data@values)
  raster::plot(para)
  
  # Run the function with chosen parameters for method, theta and p
  theo_run <- theoPath_param(ras_ai=ras_emap, 
                             ras_para=para, 
                             method="drive_i",
                             theo_con[[1]], 
                             theta=0.001, 
                             p=5, 
                             type="r")
  