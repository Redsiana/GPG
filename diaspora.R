### when the juveniles go see further if grass is greener there

# input: newbabyX, newbabyY, dispersal kernel, sex, repro
# output: newpopX, newpopY, newpopXY, coordalive

diaspora <- function( newbabyX, 
                      newbabyY, 
                      dispkernel,
                      sex, 
                      repro,
                      Xdim,
                      Ydim) {
  
 
  distance <- sample( x = dispkernel, size = length(newbabyX), replace = T )
  
  angle <- runif( n = length(newbabyX), min=0, max = 2*pi )
  
  # translate movement on X and Y grid axes
  
  movX <- round( distance * cos( angle ) )
  movY <- round( distance * sin( angle ) )
  
  newpopX <- newbabyX + movX
  newpopY <- newbabyY + movY
  
  # bouncing and wrapping boundaries
  newpopX[ newpopX < 0] <- abs( newpopX[ newpopX < 0] ) -1 # bouncing boundary at X=0
  newpopX[ newpopX > (Xdim-1) ] <- 2*(Xdim-1) - newpopX[ newpopX > (Xdim-1) ] # bouncing at (Xdim-1)
  
  newpopY[ newpopY > (Ydim-1) ] <- newpopY[ newpopY > (Ydim-1) ] - Ydim # if too far North, comes back South
  newpopY[ newpopY < 0 ] <- newpopY[ newpopY < 0 ] + Ydim # if too far South, comes back North
  
  ## do they end up in a patch without potential mating partners?
  newpopXY <- paste( newpopX,newpopY,sep="," )
  
  # create vector for females / males saying if there are mating partners in their patch 
  # (aim: save time as no mating partner = as good as dead, so no need to calculate further survival)
  
  haspartner <- vector( length = length( newpopXY ))
  haspartner[ sex == "fem" & repro == "s" ] <- newpopXY[ sex == "fem" & repro == "s" ] %in% newpopXY[ sex == "mal" ] 
  haspartner[ sex == "fem" & repro == "a" ] <- TRUE
  haspartner[ sex == "mal" ] <- newpopXY[ sex =="mal" ] %in% newpopXY[ sex =="fem" & repro == "s" ] 
  
  coordalive <- newpopXY
  coordalive[!haspartner] <- NA # coordinates of individuals, lonely and doomed ones are marked NA
  
  return( list( newpopX = newpopX, newpopY = newpopY, newpopXY = newpopXY, coordalive = coordalive, haspartner = haspartner ))
  
  ################### TO TEST - PLOT IT ALL
  # 
  # ############# Plot Hexagonia
  # y = 3/2 * s * Blueworld # get cartesian coordinates of centers on axis y
  # x = sqrt(3) * s * ( Blueworld/2 + Redworld) # get cartesian coordinates of centers on axis x
  # S <- SpatialPoints(cbind(x,y))
  # hex_grid <- HexPoints2SpatialPolygons(S)
  # plot(hex_grid, main="Hexatopia", border="white")
  # abline(0,-0.5773503, col="red")
  # abline(v=0, col="blue")
  
  # ############ Plot where the babies start
  # y = 3/2 * s * newbabyY # get cartesian coordinates of centers on axis y
  # x = sqrt(3) * s * ( newbabyY/2 + newbabyX) # get cartesian coordinates of centers on axis x
  # S <- SpatialPoints(cbind(x,y))
  # hex_grid <- HexPoints2SpatialPolygons(S, dx=1)
  # plot(hex_grid, main="Hexatopia", col="cadetblue3", add=T)
  
  
  
  ############# Plot trajectories on Hexagonia
  # 
  # initX <- sqrt(3) * s * ( newbabyY/2 + newbabyX)
  # initY <- 3/2 * s * newbabyY
  # finX <- sqrt(3) * s * ( newpopY/2 + newpopX)
  # finY <- 3/2 * s * newpopY
  # arrows(x0=initX, y0=initY, x1=finX, y1=finY, length=0.1, col="red")
  # 
  # nameplot <- paste('hexatopia', 2*t-1,'.png', sep="")
  # dev.copy(png, nameplot)
  # dev.off()
  
  # ############# Highlight newly colonized patches
  # y = 3/2 * s * newpopY # get cartesian coordinates of centers on axis y
  # x = sqrt(3) * s * ( newpopY/2 + newpopX) # get cartesian coordinates of centers on axis x
  # S <- SpatialPoints(cbind(x,y))
  # hex_grid <- HexPoints2SpatialPolygons(S, dx=sqrt(3)*2/3) # if don't specify a dx, plot will screw up when only one patch is plotted
  # plot(hex_grid, main= paste("Hexatopia", t), col="red", add=T)
  # # arrows(x0=initX, y0=initY, x1=finX, y1=finY, length=0.1, col="red")
  # 
  # nameplot <- paste('hexatopia', 2*t,'.png', sep="")
  # dev.copy(png, nameplot)
  # dev.off()
  
}


