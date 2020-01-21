library(ggplot2)
library(ggspatial)
library(sp)
library(spdep)
library(raster)
library(spatstat)
library(sf)

source("1script/repath.R") #Rscript with function repath and asscociated custom functions
source("1script/theoPath.R") #Rscript with functions theo_path and asscociated custom functions

# Loading point data of barrow graves and create sf-PointObject for later plotting with ggplot2
testmatrix <- read.csv("2data/test_points.csv")
sp_map <- st_read("2data/test_points.csv", options=c("X_POSSIBLE_NAMES=x","Y_POSSIBLE_NAMES=y"))
st_crs(sp_map) <- 25833
#sp_test <- sp::SpatialPointsDataFrame(cbind(testmatrix$x, testmatrix$y), proj4string = sp::CRS("+init=epsg:25833"))

# Loading raster data of DEM and create df for later plotting with ggplo2
testrast <- raster::raster("3geodata/Raster/EUDEM_DGM25_UTM33_test_mask.tif")
testmap <- raster::rasterToPoints(testrast)
testmap <- data.frame(testmap)
colnames(testmap) <- c("x","y","value")

for (i in seq(from=250, to=10000,by=250)){
  tryCatch({
    # Find Local Centers of Infrastructure from Monument Location
    print(i)
    maxima <- localMax(testmatrix, r=i)
    proj4string(maxima) <- sp::CRS("+init=epsg:25833")
    mxm <- sf::st_as_sf(maxima,coords = 1:2)
    
    # Plot locations of local maxima
    # mmap <- ggplot()+
    # ggtitle(i)+
    # theme(plot.title = element_text(hjust = 0.5))+
    # geom_point(data=mxm, aes(x=maxima@coords[,1], y=maxima@coords[,2]))
    # plot(maxima@coords[,1] , maxima@coords[,2] )
    # ggsave(filename = paste("5figures/Maxima_Test_Radius",i,".png",sep=""), plot=mmap)
    
    # Setting geographical frame
    xmin    <- min(sp_map$x)
    xmax    <- max(sp_map$x)
    ymin    <- min(sp_map$y)
    ymax    <- max(sp_map$y)
    ext_ai <- raster::extent(xmin, xmax, ymin, ymax)
    
    # Coordinates to set frame corner to define aspect ratio                      
    sv <- (xmax-xmin)/(ymax-ymin)
    rw     <- 25   # width of raster defined in m 
    
    # Defintion of frame expansion and defining frame                              
    rows  <- round((ymax-ymin)/rw, 0) + 1                                    
    colums <- round((xmax-xmin)/rw, 0) + 1                                     
    v <- cbind(1:(colums*rows))                                              
    df <- data.frame(v)                                                         
    gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
    sgdf    <- sp::SpatialGridDataFrame(gt, df) 
    
    # Initialising observation window that in this case equals the reserach area
    win <- spatstat::owin(c(xmin, xmax),c(ymin, ymax))
    
    # calculating theoretical connections via delaunay triangulation
    theo_con <- theo_del(maxima,win)
    
    # Plot theoretical connections
    theo_con_lines <- sf::st_as_sf(theo_con[[2]])
    st_crs(theo_con_lines) <- 25833
    
    #theo_con_lines_map <- ggplot()+
    #  ggtitle(i)+
    #  theme(plot.title = element_text(hjust = 0.5))+
    #  geom_sf(data=sp_map)+
    #  geom_sf(data = theo_con_lines)+
    #  coord_sf(datum = sf::st_crs(25833))
    #ggsave(filename = paste("5figures/Theo_lines_Test",i,".png",sep=""), plot=theo_con_lines_map)
    
    # Initialising empty raster to store summed up expectations of single rWalk connections
    ras_empty <- ras_emap
    ras_empty@data@values <- NA
    
    # Run the function with chosen parameters for method, theta and p
    theo_run <- theoPath_herzog(emp_ai=ras_empty,ras_ai=ras_emap,method="drive_i",theo_con[[1]], theta=0.001, p=1000, type="c")
    
    # Plot the 
    proj4string(theo_run) <- sp::CRS("+init=epsg:25833")
    theo_run_map <- raster::as.data.frame(theo_run, xy=TRUE)
    colnames(theo_run_map) <- c("x", "y", "value")
    
    theo_paths_map <- ggplot() + 
      ggtitle(i)+
      geom_raster(data=theo_run_map, aes(x=x, y=y, fill=value), alpha=0.8) +
      viridis::scale_fill_viridis() +
      geom_point(data = sp_map, aes(x=x, y=y))+
      geom_sf(data = theo_con_lines)+
      coord_sf(datum = sf::st_crs(25833))
    ggsave(filename = paste("5figures/Theo_paths_Test",i,".png",sep=""), plot=theo_paths_map)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
