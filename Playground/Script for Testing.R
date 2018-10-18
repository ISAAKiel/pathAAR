library(sp)
library(raster) 
library(rgdal) 
library(spatstat) 
library(spdep) 
library(maptools)

source("R/repath.R")
source("R/theoPath.R")

# Crating Test Data Randomly
testmatrix <- data.frame(
  x = abs(rnorm(100)*50), 
  y = abs(rnorm(100)*50)
)

#plot(testmatrix$x, testmatrix$y)

maxima <- localMax(testmatrix, r=0.5)

plot(maxima@coords[,1] , maxima@coords[,2] )

# Setting geographical frame
xmin    <- 0
xmax    <- max(testmatrix$x)
ymin    <- 0
ymax    <- max(testmatrix$y)
ext_ai <- extent(xmin, xmax, ymin, ymax)

# Coordinates to set frame corner to define aspect ratio                      
sv <- (xmax-xmin)/(ymax-ymin)
rw     <- 10   # width of raster defined in m

# Defintion of frame expantion and defining frame                              
rows  <- round((ymax-ymin)/rw, 0) + 1 ; rows                                    
colums <- round((xmax-xmin)/rw, 0) + 1 ; colums                                     
v <- cbind(1:(colums*rows))                                              
df <- data.frame(v)                                                         
gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
sgdf    <- sp::SpatialGridDataFrame(gt, df) 


PATH <- repath(df=testmatrix, sgdf, rw=rw)

