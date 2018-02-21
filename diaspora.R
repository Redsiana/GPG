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
  newpopY[ newpopY > (Ydim-1) ] <- newpopY[ newpopY > (Ydim-1) ] - Ydim # x2 in case loong dispersal, and goes twice around
  newpopY[ newpopY < 0 ] <- newpopY[ newpopY < 0 ] + Ydim # if too far South, comes back North
  newpopY[ newpopY < 0 ] <- newpopY[ newpopY < 0 ] + Ydim # x2 in case looong dispersal
  
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

  
}

library(compiler)
diaspora <- cmpfun(diaspora)
