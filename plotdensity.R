## plot density of all (green) and asexuals (black, on top) individuals per patch

plotdensity <- function( K, Xdim, Ydim, Xinit, popsize, asexpopsize, plotname, t. = t, popX = popX, popY = popY, popXY = popXY){

  xtot <-  popX[ match(names(popsize), popXY) ]
  ytot <- popY[ match(names(popsize), popXY) ]
  
  xasex <- popX[ match(names(asexpopsize), popXY) ]
  yasex <- popY[ match(names(asexpopsize), popXY) ]
  
  side.tot <- sqrt( popsize / K )
  side.asex <- sqrt( asexpopsize / K )
  x1 <- c( (xtot - side.tot / 2), (xasex - side.asex / 2) )
  x2 <- c( (xtot + side.tot / 2), (xasex + side.asex / 2) )
  y1 <- c( (ytot - side.tot / 2), (yasex - side.asex / 2) )
  y2 <- c( (ytot + side.tot / 2), (yasex + side.asex / 2) )
  
  categ <- c( rep("total", length(xtot)), rep("asex", length(xasex)) ) 
  d <- data.frame(x1, x2, y1, y2, categ)
  
  ggplot() +
    scale_x_continuous(name="x", limits = c(0, Xdim/2)) +
    scale_y_continuous(name="y", limits = c(0, Ydim)) +
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = categ), color = "black") +
    scale_fill_manual(values=c("black", "deepskyblue1")) +
    theme(panel.background = element_rect( fill = 'burlywood3') ) +
    geom_vline(xintercept = Xinit, color = "white") +
    ggtitle( paste( plotname, "at t =", t. ) ) +
    coord_fixed(ratio = 1)
  
  # "olivedrab2"
}




