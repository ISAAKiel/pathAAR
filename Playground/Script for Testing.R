library(sp)
library(raster) 
library(rgdal) 
library(spatstat) 
library(spdep) 
library(maptools)
#library(broom)
library(sf)
#library(foreach)

source("R/repath.R")
source("R/theoPath.R")

# Crating Test Data Randomly
testmatrix <- data.frame(
  x = abs(rnorm(100)*50), 
  y = abs(rnorm(100)*50)
)

#plot(testmatrix$x, testmatrix$y)

maxima <- localMax(testmatrix, r=10)

plot(maxima@coords[,1] , maxima@coords[,2] )

# Setting geographical frame
xmin    <- 0
xmax    <- max(testmatrix$x)
ymin    <- 0
ymax    <- max(testmatrix$y)
ext_ai <- extent(xmin, xmax, ymin, ymax)

# Coordinates to set frame corner to define aspect ratio                      
sv <- (xmax-xmin)/(ymax-ymin)
rw     <- 5   # width of raster defined in m 
# Rasterweite angepasst -> darf nicht zu grob sein, sonst haut es einen in der Funktion theoPath_herzog bei "tran_hdiff<- gdistance::transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)" raus -> Zellen dürfen nicht direkt nebeneinander liegen


# Defintion of frame expantion and defining frame                              
rows  <- round((ymax-ymin)/rw, 0) + 1 ; rows                                    
colums <- round((xmax-xmin)/rw, 0) + 1 ; colums                                     
v <- cbind(1:(colums*rows))                                              
df <- data.frame(v)                                                         
gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
sgdf    <- sp::SpatialGridDataFrame(gt, df) 

win <- owin(c(xmin, xmax),c(ymin, ymax))

PATH <- repath(df=testmatrix, sgdf, rw=rw)

theo_con <- theo_del(maxima,win)
theo_con2 <- theo_del(maxima,win)

sgdf@data$v <- sample((50:56), length(sgdf@data$v), replace=T)

ras_sgdf <- raster(sgdf)
ras_empty <- ras_sgdf
ras_empty@data@values <- NA
moep <- theoPath_herzog(emp_ai=ras_empty,ras_ai=ras_sgdf,method="drive_i",theo_con[[1]], theta=0.001, p=5)
#plot(sl_con)

emp_ai <- ras_empty
ras_ai <- ras_sgdf 
p <- 5
con <- theo_con[[1]] 
theta <- 0.001
method <- "drive_i"

# Friction Raster for Prefering lowlands:
para <- ras_sgdf
para@data@values[which(para@data@values <0 )] <- 0
para@data@values <- para@data@values/ max(para@data@values)
plot(para)

theo_run <- theoPath_param(emp_ai=ras_empty,ras_ai=ras_ai, ras_para=para, method="drive_i",theo_con[[1]], theta=0.001, p=5)



