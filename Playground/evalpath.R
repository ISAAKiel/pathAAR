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

#' Theoretical Paths using cost functions defined by Herzog 2012
#' 
#' If there are no actual parts of a path network known, the `theoPath_herzog` function can be used to reconstruct pathways based on randomised shortest paths connecting known regions with higher densities of sites, e.g. monuments. An underlying cost surface is created by using the cost functions defined by I. Herzog (2012) either for walking or driving. This function is a useful step in the evaluation of reconstructed paths. 
#'
#' 
#' @title theoPath_herzog
#' 
#' @param ras_ai RasterLayer, raster with elevation values
#' @param con data.frame, connections of a Delaunay triangulation as a result of function theo_del or your own method
#' @param method chr, either "walk_i" for pedestrians or "drive_i" for vehicles. For further informations look up the respective functions
#' @param theta numeric, parameter controls randomisation of walk. Lower values equal more exploration of the walker around the shortest path, while if theta approaches zero the walk becomes totally random 
#' @param p numeric, buffer zone around rWalk rasters, used in loop
#' @param type chr, either "c" (default) for least-cost distances or "r" for random walks. As stated by J. van Etten, there is no analytical way as of now to decide for intermediate values of theta which type should be choosen. For further informations see ?gdistance::geoCorrection
#'  
#' @return raster object, with values of summed up expectations of single rWalk connections
#'           
#' @author Franziska Faupel <\email{ffaupel@@ufg.uni-kiel.de}>
#' @author Oliver Nakoinz <\email{oliver.nakoinz.i@@gmail.com}>
#' @author Hendrik Raese <\email{h.raese@@roots.uni-kiel.de}>
#' 
#' @examples 
#' 
#' set.seed(123)
#' 
#' # Creating random test data 
#' testmatrix <- data.frame(
#'  x = abs(runif(100, 1000, 1500)), 
#'  y = abs(runif(100, 1000, 1500)))
#'
#' # Calculate local centres of infrastructure from monument location
#' pd=25
#' sw=10
#' r=100
#' maxima <- localMax(df=testmatrix, r=r, sw=sw, pd=pd)
#' 
#' # Setting geographical frame
#' xmin    <- min(testmatrix$x)
#' xmax    <- max(testmatrix$x)
#' ymin    <- min(testmatrix$y)
#' ymax    <- max(testmatrix$y)
#' ext_ai <- raster::extent(xmin, xmax, ymin, ymax)
#' 
#' # Coordinates used to set frame corner for definition of the aspect ratio
#' sv <- (xmax-xmin)/(ymax-ymin)
#' rw     <- 5   # width of raster defined in m
#' 
#' # Definition of frame expansion and defining frame                              
#' rows  <- round((ymax-ymin)/rw, 0) + 1                                  
#' colums <- round((xmax-xmin)/rw, 0) + 1                                     
#' v <- cbind(1:(colums*rows))                                              
#' df <- data.frame(v)                                                         
#' gt      <- sp::GridTopology(c(xmin, ymin), c(rw, rw), c(colums, rows))
#' sgdf    <- sp::SpatialGridDataFrame(gt, df)
#' 
#' # Initialising observation window for theoretical connections
#' win <- spatstat::owin(c(xmin, xmax),c(ymin, ymax))
#' 
#' # calculating theoretical connections via Delaunay triangulation
#' theo_con <- theo_del(maxima,win)
#' 
#' # Setting up an artificial elevation map with random values
#' emap <- sgdf
#' emap@data$v <- sample((50:56), length(emap@data$v), replace=TRUE)
#' ras_emap <- raster::raster(emap)
#' 
#' # Run the function with chosen parameters for method, theta and p
#' theo_run <- theoPath_herzog(ras_ai=ras_emap,
#'                             method="drive_i",
#'                             theo_con[[1]], 
#'                             theta=0.001, 
#'                             p=5, 
#'                             type="r")
#'
#'@export




# Loading Packages which function are used frequently
library(ggplot2)
# Load functions to you namespace! 
source("9praes/Bloc 3/FUN_evalPath.R")
#================================================================
## 1. Load Data 
#================================================================
# Monuments
g <- read.table("2data/BA_SH_gravemounds.csv", header=TRUE,
                sep=',')
g <- g[-c(which(g$X<=4250000),which(g$X>=4350000)),]
sp_g <- sp::SpatialPoints(cbind(g$X, g$Y))
raster::projection(sp_g)  <- "+init=epsg:3035"
# Basic Model
ras_b <- sp::read.asciigrid("3geodata/raster_basic.asc")
ras_b <- raster::raster(ras_b)
# Theo. Model 1
ras_M1 <- sp::read.asciigrid("3geodata/raster_walk.asc")
ras_M1 <- raster::raster(ras_M1)
# Theo. Model 2
ras_M2 <- sp::read.asciigrid("3geodata/raster_drive.asc")
ras_M2 <- raster::raster(ras_M2)
# Theo. Model 3
ras_M3 <- sp::read.asciigrid("3geodata/raster_view.asc")
ras_M3 <- raster::raster(ras_M3)
# Theo. Model 3*1000
ras_M31000 <- sp::read.asciigrid("3geodata/raster_view1000.asc")
ras_M31000 <- raster::raster(ras_M31000)
# Theo. Model 4
ras_M4 <- sp::read.asciigrid("3geodata/raster_viewi.asc")
ras_M4 <- raster::raster(ras_M4)
# Theo. Model 4*1000
ras_M41000 <- sp::read.asciigrid("3geodata/raster_viewi1000.asc")
ras_M41000 <- raster::raster(ras_M41000)
# Theo. Model 5
ras_M5 <- sp::read.asciigrid("3geodata/raster_hight.asc")
ras_M5 <- raster::raster(ras_M5)
# Theo. Model 5*1000
ras_M51000 <- sp::read.asciigrid("3geodata/raster_hight1000.asc")
ras_M51000 <- raster::raster(ras_M51000)
# Theo. Model 6
ras_M6 <- sp::read.asciigrid("3geodata/raster_bog.asc")
ras_M6 <- raster::raster(ras_M6)
# Theo. Model 6*1000
ras_M61000 <- sp::read.asciigrid("3geodata/raster_bog1000.asc")
ras_M61000 <- raster::raster(ras_M61000)
# define evalues
pd <- position_dodge(0.1)
#================================================================
## 2. Deviding Passage Values by Quantile
#================================================================

grid_b <- ras_b
b_q100 <- rbind(sp::coordinates(grid_b)[which(grid_b@data@values == 10000 ),]) 
poi_b_q100 <- data.frame(b_q100)
coordinates(poi_b_q100)=~x+y 
dis_b_q100 <- FNN::get.knnx(sp::coordinates(poi_b_q100), sp::coordinates(sp_g),k=1) 

M1<- quan_theo(ras_M1, sp_g)

M2<- quan_theo(ras_M2, sp_g)

M3<- quan_theo(ras_M3, sp_g)
M31000<- quan_theo(ras_M31000, sp_g)

M4<- quan_theo(ras_M4, sp_g)
M41000<- quan_theo(ras_M41000, sp_g)

M5<- quan_theo(ras_M5, sp_g)
M51000<- quan_theo(ras_M51000, sp_g)

M6<- quan_theo(ras_M6, sp_g)
M61000<- quan_theo(ras_M61000, sp_g)

#================================================================
## 3. Combining distances of all theoretical modells into one DATA FRAME eval
#================================================================

id <- (id=1:length(sp_g@coords[,1]))
# T1
#======================================================
eval <- data.frame(sp_g, id, M1[[4]])
names(eval)[names(eval)== "nn.index"] <- "t1_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t1_q100"
eval <- data.frame(eval, M1[[3]])
names(eval)[names(eval)== "nn.index"] <- "t1_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t1_q50"
eval <- data.frame(eval, M1[[2]])
names(eval)[names(eval)== "nn.index"] <- "t1_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t1_q75"
eval <- data.frame(eval, M1[[1]])
names(eval)[names(eval)== "nn.index"] <- "t1_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t1_q25"
# T2
#=======================================================
eval <- data.frame(eval, M2[[4]])
names(eval)[names(eval)== "nn.index"] <- "t2_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t2_q100"
eval <- data.frame(eval, M2[[3]])
names(eval)[names(eval)== "nn.index"] <- "t2_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t2_q50"
eval <- data.frame(eval, M2[[2]])
names(eval)[names(eval)== "nn.index"] <- "t2_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t2_q75"
eval <- data.frame(eval, M2[[1]])
names(eval)[names(eval)== "nn.index"] <- "t2_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t2_q25"
# T3a
# ============================================================
eval <- data.frame(eval, M3[[4]])
names(eval)[names(eval)== "nn.index"] <- "t3a_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3a_q100"
eval <- data.frame(eval, M3[[3]])
names(eval)[names(eval)== "nn.index"] <- "t3a_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3a_q50"
eval <- data.frame(eval, M3[[2]])
names(eval)[names(eval)== "nn.index"] <- "t3a_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3a_q75"
eval <- data.frame(eval, M3[[1]])
names(eval)[names(eval)== "nn.index"] <- "t3a_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3a_q25"
# T3b
# ============================================================
eval <- data.frame(eval, M31000[[4]])
names(eval)[names(eval)== "nn.index"] <- "t3b_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3b_q100"
eval <- data.frame(eval, M31000[[3]])
names(eval)[names(eval)== "nn.index"] <- "t3b_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3b_q50"
eval <- data.frame(eval, M31000[[2]])
names(eval)[names(eval)== "nn.index"] <- "t3b_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3b_q75"
eval <- data.frame(eval, M31000[[1]])
names(eval)[names(eval)== "nn.index"] <- "t3b_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t3b_q25"
# T4a
# ============================================================
eval <- data.frame(eval, M4[[4]])
names(eval)[names(eval)== "nn.index"] <- "t4a_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4a_q100"
eval <- data.frame(eval, M4[[3]])
names(eval)[names(eval)== "nn.index"] <- "t4a_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4a_q50"
eval <- data.frame(eval, M4[[2]])
names(eval)[names(eval)== "nn.index"] <- "t4a_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4a_q75"
eval <- data.frame(eval, M4[[1]])
names(eval)[names(eval)== "nn.index"] <- "t4a_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4a_q25"
# T4b
# ============================================================
eval <- data.frame(eval, M41000[[4]])
names(eval)[names(eval)== "nn.index"] <- "t4b_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4b_q100"
eval <- data.frame(eval, M41000[[3]])
names(eval)[names(eval)== "nn.index"] <- "t4b_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4b_q50"
eval <- data.frame(eval, M41000[[2]])
names(eval)[names(eval)== "nn.index"] <- "t4b_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4b_q75"
eval <- data.frame(eval, M41000[[1]])
names(eval)[names(eval)== "nn.index"] <- "t4b_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t4b_q25"
# T5a
# ============================================================
eval <- data.frame(eval, M5[[4]])
names(eval)[names(eval)== "nn.index"] <- "t5a_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5a_q100"
eval <- data.frame(eval, M5[[3]])
names(eval)[names(eval)== "nn.index"] <- "t5a_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5a_q50"
eval <- data.frame(eval, M5[[2]])
names(eval)[names(eval)== "nn.index"] <- "t5a_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5a_q75"
eval <- data.frame(eval, M5[[1]])
names(eval)[names(eval)== "nn.index"] <- "t5a_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5a_q25"
# T5b
# ============================================================
eval <- data.frame(eval, M51000[[4]])
names(eval)[names(eval)== "nn.index"] <- "t5b_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5b_q100"
eval <- data.frame(eval, M51000[[3]])
names(eval)[names(eval)== "nn.index"] <- "t5b_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5b_q50"
eval <- data.frame(eval, M51000[[2]])
names(eval)[names(eval)== "nn.index"] <- "t5b_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5b_q75"
eval <- data.frame(eval, M51000[[1]])
names(eval)[names(eval)== "nn.index"] <- "t5b_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t5b_q25"
# T6a
# ============================================================
eval <- data.frame(eval, M6[[4]])
names(eval)[names(eval)== "nn.index"] <- "t6a_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6a_q100"
eval <- data.frame(eval, M6[[3]])
names(eval)[names(eval)== "nn.index"] <- "t6a_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6a_q50"
eval <- data.frame(eval, M6[[2]])
names(eval)[names(eval)== "nn.index"] <- "t6a_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6a_q75"
eval <- data.frame(eval, M6[[1]])
names(eval)[names(eval)== "nn.index"] <- "t6a_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6a_q25"
# T6b
# ============================================================
eval <- data.frame(eval, M61000[[4]])
names(eval)[names(eval)== "nn.index"] <- "t6b_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6b_q100"
eval <- data.frame(eval, M61000[[3]])
names(eval)[names(eval)== "nn.index"] <- "t6b_cell_q50"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6b_q50"
eval <- data.frame(eval, M61000[[2]])
names(eval)[names(eval)== "nn.index"] <- "t6b_cell_q75"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6b_q75"
eval <- data.frame(eval, M61000[[1]])
names(eval)[names(eval)== "nn.index"] <- "t6b_cell_q25"
names(eval)[names(eval)== "nn.dist"] <- "dis_t6b_q25"
# B
# ============================================================
eval <- data.frame(eval, dis_b_q100)
names(eval)[names(eval)== "nn.index"] <- "b_cell_q100"
names(eval)[names(eval)== "nn.dist"] <- "dis_b_q100"

#================================================================
## 3. Calculating Median of distances for all theoretical modells
#================================================================
eval <- rbind(eval, c(1)) # adding additional row
eval[length(eval[,1]), c(4:length(eval[1,]))] <- apply(eval[,c(4:length(eval[1,]))], 2, median, na.rm=FALSE)

#================================================================
## 4.Create Results DataFrame and Plot
#================================================================

teval <- t(eval)
median_quan <- data.frame(median=teval[c(4:length(eval[1,])),length(eval[,1])])
no <- 11 # number of models

df <- data.frame("Quantile"=c(rep("4th quantile", no)),
                 "Model"= c("M1", "M2", "M3", "M3*1000", "M4", "M4*1000", "M5", "M5*1000", "M6", "M6*1000", "B"),
                 "Distance" =c(median_quan[c(2,10,18,26,34,42,50,58,66,74,82),])
)

ggplot(df, aes(x=Model, y=Distance, group=df[,1], shape=df[,1], colour=df[,1])) + 
  geom_line(aes(size=df[,1]))+ scale_size_manual(values = c(0.5,0.75,1)) + 
  scale_colour_manual(values =c("#669900", "#E69F00", "#56B4E9"))+
  geom_point()+ 
  #scale_colour_manual(values =c("#669900", "#E69F00", "#56B4E9"))+
  geom_hline(yintercept=df[11,3], colour="#D55E00", linetype = "dotted", size=1)+
  xlab("Theoretische Modelle")+
  ylab("Distanz (Median in m)")+
  ggtitle("Vergleich der theoretischen Modelle")+
  theme(legend.title=element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey"),
        panel.grid.minor = element_line(colour = "grey", linetype = "dotted"),
        panel.background = element_rect(fill = "white"))

#================================================================
## Report
#================================================================
savehistory(file="6report/eval_pathe01.Rhistory")
