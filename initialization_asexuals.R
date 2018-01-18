##### creating the initial population kernel

# INPUT: K, Ninit, Xinit, Yini, fec, G
# 
# OUTPUT: sex, popgenome, newbabyX, newbabyY, repro


genesis <- function(K, 
                    fec, 
                    Ninit,
                    Xinit,
                    Yinit,
                    G){
  
  # surface of the starting band

  Sinit <- Xinit*Yinit
  
  # pop starts half female, half male, all sexuals, all on patch (0;0)
  sex <- c( rep( "fem", each = ceiling( Sinit * (Ninit/2) )),
            rep( "mal", each= floor( Sinit * (Ninit/2) )))   
  sex <- sample(sex)
  repro <- rep("s", Sinit*Ninit)
  newbabyX <- rep( 0:(Xinit-1), (Yinit*Ninit) )
  newbabyY <- rep( 0:(Yinit-1), each = (Xinit*Ninit) )
  
  # creating a collection of Ninit genotypes for one locus at equilibrium proportions (1/4 00, 1/4 11, 1/2 01)
  col1 <- head( sample( c( rep(1, ceiling( (Sinit*Ninit)/2 )), 
                           rep(0, ceiling( (Sinit*Ninit)/4 )), 
                           rep(2, ceiling( (Sinit*Ninit)/4 ))) ), (Sinit*Ninit) )
  
  # matrix with population genotype for 1 gene
  popgenome <- as.matrix( col1, drop = FALSE )
 
  # extended to G genes if G > 1
  if( G != 1 ){
    for ( g in 1:(G-1) ){
      popgenome <- cbind( popgenome, sample( col1, Ninit ) )
    }}
  
  popcloneline <- rep(0, length(repro))
  popclonalorigin <- matrix(0, nrow=length(repro), ncol=2)
  
  return( list( sex = sex, popgenome = popgenome, newbabyX = newbabyX, 
                newbabyY = newbabyY, repro = repro, popcloneline = popcloneline,
                popclonalorigin = popclonalorigin ))
}





# to homogeneize initial conditions between trials, pop always start with exactly half males and half females, all heterozygous


