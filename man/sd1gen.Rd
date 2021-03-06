% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repath.R
\name{sd1gen}
\alias{sd1gen}
\title{sd1gen}
\usage{
sd1gen(pppm, f_sd1 = 4)
}
\arguments{
\item{pppm}{PointPattern Object of Data}

\item{f_sd1}{factor scaling the mean of nearest neighbourhoud}
}
\value{
list containg two numeric
}
\description{
Factor defining size of first Kernel
}
\details{
Factor defining size of the first Kernel, which generates the stucture of dynamic Kernel
}
\examples{
# Crating Test Data Randomly
testmatrix <- data.frame(x = abs(rnorm(100)*50), y = abs(rnorm(100)*50))
# Setting geographical frame
xmin    <- 0
xmax    <- max(testmatrix$x)
ymin    <- 0
ymax    <- max(testmatrix$y)
ext_ai <- raster::extent(xmin, xmax, ymin, ymax)

sv <- (xmax-xmin)/(ymax-ymin)
rw     <- 10   # width of raster defined in m

rows  <- round((ymax-ymin)/rw, 0) + 1                                    
colums <- round((xmax-xmin)/rw, 0) + 1                                      
v <- cbind(1:(colums*rows))                                              
df <- data.frame(v)                                                         
gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
sgdf    <- sp::SpatialGridDataFrame(gt, df)

pppm <- spatstat::ppp(testmatrix[,1], testmatrix[,2], 
                       window = spatstat::owin(
                       xrange=c(sgdf@bbox[1,1],sgdf@bbox[1,2]),
                       yrange=c(sgdf@bbox[2,1],sgdf@bbox[2,2]), 
                       unitname="m"))
                       
f_sd1 <- sd1gen(pppm)

}
\author{
Franziska Faupel <\email{ffaupel@ufg.uni-kiel.de}>

Oliver Nakoinz <\email{oliver.nakoinz.i@gmail.com}>

Hendrik Raese <\email{h.raese@roots.uni-kiel.de}>
}
