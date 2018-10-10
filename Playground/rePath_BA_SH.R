##################################################################
## R-Script: Reconstructing Bronze Age Paths in Schleswig-Holstein
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
source("9praes/Bloc 3/FUN_rePath.R")

#================================================================
## 1. Load Data and set arear of Interesst
#================================================================

g <- read.table("../../2data/BA_SH_gravemounds.csv", header=TRUE,
                sep=',')
g <- g[-c(which(g$X<=4250000),which(g$X>=4350000)),]

sp_g <- sp::SpatialPoints(cbind(g$X, g$Y))
raster::projection(sp_g)  <- "+init=epsg:3035"

#==================================================================
# Plotting not needed
library(ggplot2)
ggplot()+
  ggspatial::geom_spatial(data=sp_g)+ 
  theme_void()
#==================================================================

## Setting a Grid to store results
xmin    <- sp_g@bbox[1,1]
xmax    <- sp_g@bbox[1,2]
ymin    <- sp_g@bbox[2,1]
ymax    <- sp_g@bbox[2,2]

#================================================================
## 2. Create Spatial Grid for KDE
#================================================================
# Coordinates to set frame corner to define aspect ratio
sv      <- (xmax-xmin)/(ymax-ymin)
# Defintion of frame expantion and defining frame
rw      <- 500
rows    <- round((ymax-ymin)/rw, 0) + 1 ; rows      # no of rows
colums  <- round((xmax-xmin)/rw, 0) + 1 ; colums   # no of coulmns
df      <- data.frame(v=cbind(1:(colums*rows)))        # vector for all cells
gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
sgdf    <- sp::SpatialGridDataFrame(gt, df, 
          proj4string = sp::CRS(as.character("+init=epsg:3035" ))) # Grid wird projeziert

#==================================================================
# Plotting not needed
df_sgdf <- as.data.frame(sgdf)
colnames(df_sgdf) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data=df_sgdf, aes(x=x, y=y, fill=value), alpha=0.8) +
  viridis::scale_fill_viridis() +
  geom_point(data = sp_g, aes(x=sp_g@coords[,1], y=sp_g@coords[,2]))
#==================================================================

#================================================================
## 3. Static KDE
#================================================================

ppp_g <- spatstat::ppp(sp_g@coords[,1], sp_g@coords[,2], 
                       window = spatstat::owin(
                         xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
                         yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
                         unitname="m"))

# variables depending on nearest neighbourhood distnce
nn   <- mean(spatstat::nndist(ppp_g))
sd1  <- 4*nn  # factor defining size of the first kernel,
              # which generate the stucture of dynamic kernel

# calculate static-KDE

base_kde <- makestatkde(sd1 = sd1, sgdf = sgdf, df = g)

#==================================================================
# Plotting not needed
df_kde <- as.data.frame(base_kde)
colnames(df_kde) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data=df_kde, aes(x=x, y=y, fill=value), alpha=0.8) +
  viridis::scale_fill_viridis() +
  geom_point(data = sp_g, aes(x=sp_g@coords[,1], y=sp_g@coords[,2]))
#==================================================================

maptools::writeAsciiGrid(base_kde,"base_kde")

#================================================================
## 4. Dynamic KDE
#================================================================
# define parameters
iter   <- 2
rw     <- 500      # width of raster defined in m
tresh  <- 0.05     # treshold, to delete paths in areas with low 
#    density values (KDE), meaning calculation 
#    artefacts 
f_sd1  <- 4        # factor defining size of the first kernel,
#    which generates the stucture of dynamic 
#    kernel
f1     <- 0.2      # factor defining the minimum border of dynamic 
#    kernel (raster width) f1*mean(nn)  ## 0.2
f2     <- 0.4      # factor defining the maximum border of dynamic 
#    kernel f2*mean(nn)
f3     <- 0.5      # minimal intensity of Kernel
f4     <- 1        # maximal intensity of Kernel
s      <- -0.3     # Kernelparameter: incline starting from point 1
de     <- 0.7      # hight of additional kernel peak
sw     <- 12       # width of picture, cm
mwin   <- 9        # Mowing-window-size for ridge detection (4,9,16)
xp     <- 750      # Kernelparameter: x-value of point 1

PATH <- repath(base_kde, nn, sgdf)

#maptools::writeAsciiGrid(PATH[[dyn_kde]],"dyn_kde")
dyn_kde <- maptools::readAsciiGrid("../../4ws/dyn_kde")

#==================================================================
# Plotting not needed
df_kde <- as.data.frame(dyn_kde)
colnames(df_kde) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data=df_kde, aes(x=x, y=y, fill=value), alpha=0.8) +
  viridis::scale_fill_viridis() 

#maptools::writeAsciiGrid(PATH[[ras_dens_ridges]],
# "ras_ridges_corr")r
ras_dens_ridges <- maptools::readAsciiGrid("
                      ../../4ws/ras_ridges_corr")

#==================================================================
# Plotting not needed
df_kde <- as.data.frame(ras_dens_ridges)
colnames(df_kde) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data=df_kde, aes(x=x, y=y, fill=value), alpha=0.8) 
#==================================================================

#==================================================================
## 5. Results
#==================================================================
# extracting the coordinates
coord_rid <- rbind(sp::coordinates(ras_dens_ridges)[which
             (ras_dens_ridges@data$ras_ridges_corr==1),])

df <- data.frame(coord_rid, id=1:length(coord_rid[,1]))
sp::coordinates(df)=~s1+s2
df <- sp::remove.duplicates(df, zero=1000, remove.second=TRUE) 
# width of kde 500, so the radius needs to be 1000

rid_cle <- df@coords 

nb_cle <- spdep::graph2nb(spdep::relativeneigh(rid_cle)) 
wts <- seq(1:length(nb_cle)); wts[] <- 1
nbl_cle <- spdep::nb2lines(nb_cle, wts, rid_cle, 
           proj4string=sp::CRS(as.character("+init=epsg:3035" )))

ggplot()+
  ggspatial::geom_spatial(data=nbl_cle, colour = "darkgrey")+ 
  theme_void()

coord_rid <- rbind(sp::coordinates(ras_dens_ridges)[which
             (ras_dens_ridges@data$ras_ridges_corr==1),])

df <- data.frame(coord_rid, id=1:length(coord_rid[,1]))
sp::coordinates(df)=~s1+s2
df <- sp::remove.duplicates(df, zero=2500, remove.second=TRUE)

rid_cle <- df@coords 

nb_cle <- spdep::graph2nb(spdep::relativeneigh(rid_cle)) 
wts <- seq(1:length(nb_cle)); wts[] <- 1
nbl_cle <- spdep::nb2lines(nb_cle, wts, rid_cle, 
            proj4string=sp::CRS(as.character("+init=epsg:3035" )))

ggplot()+
  ggspatial::geom_spatial(data=nbl_cle, colour = "darkgrey")+ 
  theme_void()

#================================================================
## 6. Export
#================================================================

## Points
ddd <- data.frame(1:length(coord_rid[,1]))
coord_rid <-  sp::SpatialPointsDataFrame(coordsnew, ddd)
rgdal::writeOGR(coord_rid, "../../6results", 
                "coord_density_ridges", layer="coord_den_ridges", 
                driver="ESRI Shapefile", overwrite=TRUE)

ddd <- data.frame(1:length(rid_cle[,1]))
rid_cle <-  sp::SpatialPointsDataFrame(rid_cle, ddd)
rgdal::writeOGR(rid_cle, "../../6results", "coord_path", 
                layer="coord_path", driver="ESRI Shapefile", 
                overwrite=TRUE)

## Lines
maptools::writeLinesShape(nbl_cle, "../../6results/lines_path")
#================================================================
## 7. Report
#================================================================
savehistory(file="7report/rePath_SH_BA.Rhistory")