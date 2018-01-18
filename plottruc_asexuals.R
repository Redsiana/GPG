
############################# PLOT MEAN HETEROZYGOSITY

## input: truc= popR, popB, popRB, object, plotname, Blueworld, Redworld, s


plottruc <- function( popsize, asexpopsize, plotname, t. = t, popR. = popR, popB. = popB, popRB. = popRB, Blueworld. = Blueworld, Redworld. = Redworld, s. = s){
  
  y = 3/2 * s * Blueworld # get cartesian coordinates of centers on axis y
  x = sqrt(3) * s * ( Blueworld/2 + Redworld) # get cartesian coordinates of centers on axis x
  S <- SpatialPoints(cbind(x,y))
  hex_grid <- HexPoints2SpatialPolygons(S)
  plot(hex_grid, main=paste(plotname, "t =", t), border="white")
  abline(0,-0.5773503, col="red")
  abline(v=0, col="blue")
  

    for( i in 1:length( popsize )){
      trucR <- unique( popR[ popRB == names(popsize)[i] ] )
      trucB <- unique( popB[ popRB == names(popsize)[i] ] )
      ytruc <- 3/2 * s * trucB
      xtruc <- sqrt(3) * s * ( trucB/2 + trucR)
      colr <- 'olivedrab2'
      
      dxf <- sqrt( (popsize[i] / K) )
      # dxf[dxf>1] <- 1
      
      S <- SpatialPoints(cbind(xtruc,ytruc))
      hex_grid <- HexPoints2SpatialPolygons(S, dx=dxf*sqrt(3)*2/3)
      plot(hex_grid, col=colr, add=T)
    }
    
  
  
  if ( sum( asexpopsize > 0 )){ 
    for( i in 1:length( asexpopsize )){
      trucR <- unique( popR[ popRB == names(asexpopsize)[i] ] )
      trucB <- unique( popB[ popRB == names(asexpopsize)[i] ] )
      ytruc <- 3/2 * s * trucB
      xtruc <- sqrt(3) * s * ( trucB/2 + trucR)
      colr <- "black"
      
      dxf <- sqrt( (asexpopsize[i] / K) )
      # dxf[dxf>1] <- 1
      
      S <- SpatialPoints(cbind(xtruc,ytruc))
      hex_grid <- HexPoints2SpatialPolygons(S, dx=dxf*sqrt(3)*2/3)
      plot(hex_grid, col=colr, add=T)
    }
    
 }
  
  
 
  
  
  
#   legend(12, 10, seq(1,0,-0.1)  )
#   colorbar.plot(13, 0.5, strip=seq(0,1,0.1), strip.width = 0.05, strip.length = 0.48, 
#                 adj.x = 0.5, adj.y = 0.5, col = gray(seq(1,0,-0.1)), horizontal = F)
 
}




