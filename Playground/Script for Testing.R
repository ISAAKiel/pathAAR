library(sp)
library(raster) 
library(rgdal) 
library(spatstat) 
library(spdep) 
library(maptools)

source("R/repath.R")

# Crating Test Data Randomly
testmatrix <- data.frame(
  x = abs(rnorm(100)*50), 
  y = abs(rnorm(100)*50)
)

#plot(testmatrix$x, testmatrix$y)

# Setting geographical frame
xmin    <- 0
xmax    <- max(testmatrix$x)
ymin    <- 0
ymax    <- max(testmatrix$y)
ext_ai <- extent(xmin, xmax, ymin, ymax)

# Coordinates to set frame corner to define aspect ratio                      # Rahmeneckpunkte werden genutzt, das Seitenverhältnis des Ergebnisses zu definieren
sv <- (xmax-xmin)/(ymax-ymin)
rw     <- 10   # width of raster defined in m

# Defintion of frame expantion and defining frame                              # Rahmenausdehnung wird definiert und Rahmen wird erstellt
rows  <- round((ymax-ymin)/rw, 0) + 1 ; rows                                    #Anzahl der Zeilen
colums <- round((xmax-xmin)/rw, 0) + 1 ; colums                                     #Anzahl der Spalten
v <- cbind(1:(colums*rows))                                              #Vektor f alle Gridzellen
df <- data.frame(v)                                                         #Konvertiert den Vektor zu einem Dataframe
gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
sgdf    <- sp::SpatialGridDataFrame(gt, df) # Grid wird projeziert


# This should be done by the repath Function later!
pppm <- spatstat::ppp(testmatrix[,1], testmatrix[,2], 
                      window = spatstat::owin(
                        xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
                        yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
                        unitname="m"))

# sd1 <- sd1gen(pppm)
base_kde <- makestatkde(pppm, sgdf=sgdf, df=testmatrix, x=1, y=2, num=length(df[,1]))


iter   <- 2
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

PATH <- repath(base_kde, sgdf)


dyn_KDE <- repath()

