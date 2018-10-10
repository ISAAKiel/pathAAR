##################################################################
##R-Script: Theoretical Paths
## ===============================================================
## Project: MOSAICnet Research School 2018, Bibracte
## Author: Franziska Faupel (O.Nakoinz)
## Version: 01
## Date of last changes: 18.09.2018
## Data (Graves): Gravemounds in Daenisch Wholdd, Bronze Age
## Author of data (Graves): GSDHL Kiel
## Data (SRTM): www.naturalearthdata.com (with contributions from 
##                F. Faupel)
## Author of data (SRTM): Tom Patterson,Nathaniel Vaughn Kelso and 
##                Franziska Faupel
## Purpose: Reconstructing Path System using empirical Data
## Licence Script: GPL (http://www.gnu.org/licenses/gpl-3.0.html)
##################################################################
# Loading Packages which function are used frequently
library(ggplot2)
# Load functions to you namespace! 
source("9praes/Bloc 3/FUN_theoPath.R")
#================================================================
## 1. Load Data 
#================================================================
g <- read.table("2data/BA_SH_gravemounds.csv", header=TRUE,
                sep=',')
g <- g[-c(which(g$X<=4250000),which(g$X>=4350000)),]
sp_g <- sp::SpatialPoints(cbind(g$X, g$Y))
raster::projection(sp_g)  <- "+init=epsg:3035"
# Digital Elevation Modell of Schleswig-Holstein, GER
dem_SH <- sp::read.asciigrid("3geodata/SH_srtm90.asc")
ras_SH <- raster::raster(dem_SH)
ras_SH <- raster::crop(ras_SH, raster::extent(sp_g@bbox))
rm(dem_SH)
# Cumulative viewshed analysis 
ras_view <- rgdal::readGDAL("3geodata/SH_viewshed2500.asc")
ras_view <- raster::raster(ras_view)
# Inversed cumulative viewshed analysis 
ras_viewi <- rgdal::readGDAL("3geodata/SH_viewshed_invers2500.asc")
ras_viewi <- raster::raster(ras_viewi)
# Spread of bog
ras_bog <- rgdal::readGDAL("3geodata/SH_moore.tif")

#================================================================
## 2. Preparing Spatial Objects
#================================================================
# Provide empty rasterobjekt to store random walk values
grid_ai <- as(ras_SH, 'SpatialGridDataFrame')
grid_ai@data$X3geodata.SH_srtm90.asc <- NA
emp_ai <- raster::raster(grid_ai)
rm(grid_ai)
# Friction Raster for Prefering lowlands:
ras_h <- ras_SH
ras_h@data@values[which(ras_h@data@values <0 )] <- 0
ras_h@data@values <- ras_h@data@values/ max(ras_h@data@values)
plot(ras_h)
# Friction Raster for Prefering Highlands
ras_hi <- ras_h
ras_hi@data@values <- (ras_h@data@values *-1)+1
plot(ras_hi)
# Friction Raster to Avoid Bogs
ras_bogi <- raster::raster(ras_bog)
ras_bogi <- raster::crop(ras_bogi , 
                         raster::extent(emp_ai@extent))
ras_bogi <- raster::mosaic(emp_ai, ras_bogi, fun=max, tolerance=0.3)
ras_bogi@data@values[which(ras_bogi@data@values > 0)] <- 1
ras_bogi@data@values <- (ras_bogi@data@values *-1)+1
# Croppig other raster files
ras_view <- raster::crop(ras_view , 
                         raster::extent(emp_ai@extent))
ras_viewi <- raster::crop(ras_viewi , 
                          raster::extent(emp_ai@extent))

# Predefinded Value for rWalk
p <- 10000 # buffer zone around rWalk rasters, used in loop

#================================================================
## 3. Local Maxima from Monument Location
#================================================================
# Or use site locations, please generate a ppp of the coordinates
maxima <- localMax(sp_g)
# Preparing coordinates for rWalk loops
bb <- sp::bbox(sp_g) # setting bondingbox for KDE
win <- spatstat::owin(xrange=c(bb[1,1],bb[1,2]), yrange= 
                        c(bb[2,1],bb[2,2]), unitname="m")
CONN <- con_rWalk(maxima, win)
con <- CONN[[1]]
#================================================================
## 4. Basic Modell
#================================================================
# Asigning values to raster cells which are crossed by direct 
# lines from nbl_jun
## No need to interpolate before, rasterize does it....
ras_b <- raster::rasterize(CONN[[2]], emp_ai, field=10000, 
                           func='length', update=TRUE)
plot(ras_b)

#================================================================
## 5. HERZOG 2012 WALKING --> M1
#================================================================

ras_M1 <- theoPath_herzog(emp_ai, ras_ai=ras_SH,con, 
                          method="walk_i", theta = 0.001)

#================================================================
## 6. HERZOG 2012 DRIVING  --> M2
#================================================================

ras_M2 <- theoPath_herzog(emp_ai, ras_ai=ras_SH,con, 
                          method="drive_i", theta = 0.001)

#================================================================
## 7. Additional Parameter: Viewshed
#================================================================

fac_M3 <- theoPath_param(emp_ai, ras_ai=ras_SH,
                         ras_view=ras_viewi, con, method="walk_i", 
                         theta = 0.001)

#================================================================
## 8. Additional Parameter: Inverse Viewshed
#================================================================

fac_M4 <- theoPath_param(emp_ai, ras_ai=ras_SH, 
                      ras_view=ras_view, con, method = "walk_i", 
                      theta = 0.001)

#================================================================
## 9. Additional Parameter: Prefering Highlands
#================================================================

fac_M5 <- theoPath_param(emp_ai, ras_ai=ras_SH, 
                         ras_view=ras_h, con, method = "walk_i", 
                         theta = 0.001)

#================================================================
## 10. Additional Parameter: Avoiding Bogs
#================================================================

fac_M6 <- theoPath_param(emp_ai, ras_ai=ras_SH, ras_view=ras_bogi, 
                         con, method = "walk_i", theta = 0.001)

#================================================================
## 11. Export
#================================================================

raster::writeRaster(ras_b,"3geodata/raster_basic", format="ascii",
                    overwrite=TRUE)

raster::writeRaster(ras_M1,"3geodata/raster_walk", format="ascii",
                    overwrite=TRUE)

raster::writeRaster(ras_M2,"3geodata/raster_drive", format="ascii",
                    overwrite=TRUE)

raster::writeRaster(fac_M3[[1]],"3geodata/raster_view", 
                    format="ascii",
                    overwrite=TRUE)

raster::writeRaster(fac_M3[[2]],"3geodata/raster_view1000", 
                    format="ascii",
                    overwrite=TRUE)

raster::writeRaster(fac_M4[[1]],"3geodata/raster_viewi", 
                    format="ascii",
                    overwrite=TRUE)
raster::writeRaster(fac_M4[[2]],"3geodata/raster_viewi1000", 
                    format="ascii",
                    overwrite=TRUE)

raster::writeRaster(fac_M5[[1]],"3geodata/raster_hight", 
                    format="ascii",
                    overwrite=TRUE)
raster::writeRaster(fac_M5[[2]],"3geodata/raster_hight1000", 
                    format="ascii",
                    overwrite=TRUE)

raster::writeRaster(fac_M6[[1]],"3geodata/raster_bog", 
                    format="ascii",
                    overwrite=TRUE)
raster::writeRaster(fac_M6[[2]],"3geodata/raster_bog1000", 
                    format="ascii",
                    overwrite=TRUE)

#================================================================
## Report
#================================================================
savehistory(file="6report/random_pathe01.Rhistory")