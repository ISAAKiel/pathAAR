################################################################################
## Collection of Functions needed to evaluate parameters of paths
## =============================================================================
## Author of Functions: Franziska Faupel
## Version: 04
## Date of last changes: 18.09.2018
## Licence Script: 
################################################################################

#===============================================================================
# Deviding Passage Values by Quantile
#===============================================================================

quan_theo <- function(ras_M1, sp_g){
  
  grid_t1 <- ras_M1
  grid_t1@data@values[which(grid_t1@data@values <=0)] <- NA
  quan_t1 <- quantile(grid_t1, na.rm = TRUE)
  # Using quantile values to subdevide rWalk values
    grid_t1[which(grid_t1@data@values <= quan_t1[1])] <- 1000
    grid_t1[which(grid_t1@data@values <= quan_t1[2])] <- 2000
    grid_t1[which(grid_t1@data@values <= quan_t1[3])] <- 3000
    grid_t1[which(grid_t1@data@values <= quan_t1[4])] <- 4000
    grid_t1[which(grid_t1@data@values <= quan_t1[5])] <- 5000
  # extracting coordinates of each quantile
    t1_q25 <- rbind(coordinates(grid_t1)[which(grid_t1@data@values == 2000 ),]) 
    t1_q50 <- rbind(coordinates(grid_t1)[which(grid_t1@data@values == 3000 ),]) 
    t1_q75 <- rbind(coordinates(grid_t1)[which(grid_t1@data@values == 4000 ),]) 
    t1_q100 <- rbind(coordinates(grid_t1)[which(grid_t1@data@values == 5000 ),]) 
  #creating a SPatial Points object
    poi_t1_q25 <- data.frame(t1_q25)
    sp::coordinates(poi_t1_q25)=~x+y 
    poi_t1_q50 <- data.frame(t1_q50)
    sp::coordinates(poi_t1_q50)=~x+y 
    poi_t1_q75 <- data.frame(t1_q75)
    sp::coordinates(poi_t1_q75)=~x+y 
    poi_t1_q100 <- data.frame(t1_q100)
    sp::coordinates(poi_t1_q100)=~x+y 
  #Calculating (and storing) shortest distance
    #calculates shortest distance between SpatialPoints Layers (second is the one, used in the record as response)
    dis_t1_q25 <- FNN::get.knnx(sp::coordinates(poi_t1_q25), sp::coordinates(sp_g),k=1) 
    dis_t1_q50 <- FNN::get.knnx(sp::coordinates(poi_t1_q50), sp::coordinates(sp_g),k=1) 
    dis_t1_q75 <- FNN::get.knnx(sp::coordinates(poi_t1_q75), sp::coordinates(sp_g),k=1) 
    dis_t1_q100 <- FNN::get.knnx(sp::coordinates(poi_t1_q100), sp::coordinates(sp_g),k=1) 
    # nn.index uses index numbers from poi_t... not from grave mounds, BUT distances are ordered, so first row, is first grave mound!!!!
  #Exporting Distribution of Quantile values
    ddf <- raster::rasterToPoints(grid_t1); ddf <- data.frame(ddf)
    colnames(ddf) <- c("X","Y","DEM")
    ddf$DEM[which(ddf$DEM <= 2000)] <- NA
    ddf$DEM[which(ddf$DEM == 3000)] <- 2
    ddf$DEM[which(ddf$DEM == 5000)] <- 4
    ddf$DEM[which(ddf$DEM == 4000)] <- 3
    ddf$DEM<- as.factor(ddf$DEM)

  plot_quant <- ggplot2::ggplot(ddf, ggplot2::aes(X, Y)) +
    ggplot2::geom_raster(aes(fill = DEM))+
    ggplot2::scale_fill_manual(values= c ("#CCCCCC", "#999999", "#666666"))+
    ggplot2::xlab("X (in m)")+
    ggplot2::ylab("Y (in m)")+
    ggplot2::ggtitle("Spread of Passpoints (in Quantiles)")+
    ggplot2::theme(legend.title=ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "white"),
          panel.grid.major = ggplot2::element_line(colour = "grey"),
          panel.grid.minor = ggplot2::element_line(colour = "grey", linetype = "dotted"),
          panel.background = ggplot2::element_rect(fill = "white"))

  dis_t1 <- list(q25 = dis_t1_q25, q50 = dis_t1_q50, q75 = dis_t1_q75, q100 = dis_t1_q100, plot_quant)
  
return(dis_t1)
}
