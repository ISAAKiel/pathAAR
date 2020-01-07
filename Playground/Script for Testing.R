library(sp)
library(raster) 
library(rgdal) 
library(spatstat) 
library(spdep) 
library(maptools)
#library(broom)
library(sf)
#library(foreach)
library(rgrass7) # requires an installation of GRASS GIS 7.x

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
rw     <- 5   # width of raster defined in m # Rasterweite angepasst -> darf nicht zu grob sein, sonst haut es einen in der Funktion theoPath_herzog bei "tran_hdiff<- gdistance::transition(ras, transitionFunction=hdiff, directions=8, symm=TRUE)" raus -> Zellen dürfen nicht direkt nebeneinander liegen

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

## Create cumulative viewshed (based on script by D. Knitter)

tmp <- tempdir()

# use random elevation map for viewshed calculation
ras_sgdf2 <- ras_sgdf

# Prepare coordinate vectors/object for viewshed calculation

coords_list <- rasterToPoints(x = ras_sgdf2)

str(coords_list)


# transform data from RasterFile to GridTopology
rf <- raster::writeRaster(ras_sgdf2, filename=file.path(tmp, "test.tif"), format="GTiff", overwrite=TRUE)

dem <- readGDAL(file.path(tmp, "test.tif"))

# set arbitrary projected crs for raster
proj4string(dem) <- CRS("+init=epsg:3587")

# initialise empty raster
view <- dem
view@data$band1 <- 0


## parseGRASS("r.viewshed")

# enforce use of sp for working with GRASS 
use_sp()


# Set up GRASS and load data

loc <- initGRASS("C:/Program Files/GRASS GIS 7.8", home=tmp, mapset = "PERMANENT", override = TRUE)

execGRASS("g.proj", flags = c("c"), parameters = list(proj4=dem@proj4string@projargs))              

writeRAST(x = dem,
          vname="dem",
          flags = c("overwrite")
)            

writeRAST(x = view,
          vname=paste0("view_all"),
          flags = c("overwrite")
)               

execGRASS("g.region",
          parameters = list(raster = "dem",
                            res = as.character(dem@grid@cellsize[1])),
          flags = c("p")
)

# Start of the cumulative viewshed analysis 

cumview_dir <- for (i in 1:length(coords_list[,1])){   
  
  execGRASS("r.viewshed",
            flags = c("overwrite","b", "quiet"),
            parameters = list(input = "dem",
                              output = "view_tmp",
                              coordinates = cbind(coords_list[i,1],coords_list[i,2]),
                              max_distance = 15,
                              memory = 5000)
  )
  
  execGRASS("r.mapcalc",
            parameters = list(expression = paste("'view_all' = 'view_all' + view_tmp")),
            flags = c("overwrite", "quiet")
  )
}

cumview_dir <-   readRAST("view_all")

plot(cumview_dir)

