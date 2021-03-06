% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theoPath.R
\name{theoPath_herzog}
\alias{theoPath_herzog}
\title{theoPath_herzog}
\usage{
theoPath_herzog(ras_ai, con, method, theta, p, type = "c")
}
\arguments{
\item{ras_ai}{RasterLayer, raster with elevation values}

\item{con}{data.frame, connections of a Delaunay triangulation as a result of function theo_del or your own method}

\item{method}{chr, either "walk_i" for pedestrians or "drive_i" for vehicles. For further informations look up the respective functions}

\item{theta}{numeric, parameter controls randomisation of walk. Lower values equal more exploration of the walker around the shortest path, while if theta approaches zero the walk becomes totally random}

\item{p}{numeric, buffer zone around rWalk rasters, used in loop}

\item{type}{chr, either "c" (default) for least-cost distances or "r" for random walks. As stated by J. van Etten, there is no analytical way as of now to decide for intermediate values of theta which type should be choosen. For further informations see ?gdistance::geoCorrection}
}
\value{
raster object, with values of summed up expectations of single rWalk connections
}
\description{
Theoretical Paths using cost functions defined by Herzog 2012
}
\details{
If there are no actual parts of a path network known, the `theoPath_herzog` function can be used to reconstruct pathways based on randomised shortest paths connecting known regions with higher densities of sites, e.g. monuments. An underlying cost surface is created by using the cost functions defined by I. Herzog (2012) either for walking or driving. This function is a useful step in the evaluation of reconstructed paths.
}
\examples{

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

# calculating theoretical connections via Delaunay triangulation
theo_con <- theo_del(maxima,win)

# Setting up an artificial elevation map with random values
emap <- sgdf
emap@data$v <- sample((50:56), length(emap@data$v), replace=TRUE)
ras_emap <- raster::raster(emap)

# Run the function with chosen parameters for method, theta and p
theo_run <- theoPath_herzog(ras_ai=ras_emap,
                            method="drive_i",
                            theo_con[[1]], 
                            theta=0.001, 
                            p=5, 
                            type="r")

}
\author{
Franziska Faupel <\email{ffaupel@ufg.uni-kiel.de}>

Oliver Nakoinz <\email{oliver.nakoinz.i@gmail.com}>

Hendrik Raese <\email{h.raese@roots.uni-kiel.de}>
}
