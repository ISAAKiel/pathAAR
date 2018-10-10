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
v <- cbind(1:(colums*rows))                                              #Vektor fÃÂ¼r alle Gridzellen
df <- data.frame(v)                                                         #Konvertiert den Vektor zu einem Dataframe
gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
sgdf    <- sp::SpatialGridDataFrame(gt, df, 
                                    proj4string = sp::CRS(as.character("+init=epsg:3035" ))) # Grid wird projeziert


# This should be done by the repath Function later!
pppm <- spatstat::ppp(testmatrix[,1], testmatrix[,2], 
                      window = spatstat::owin(
                        xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
                        yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
                        unitname="m"))

f_sd1 <- sd1gen(pppm)
base_kde <- makestatkde(pppm, sgdf, testmatrix)
