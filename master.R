master <- function( s, 
                    K, 
                    fec,
                    fecasex,
                    G,
                    probamating,
                    Xdim,
                    Ydim,
                    Xinit,
                    Yinit,
                    dispkernel,
                    bsline,
                    Ka,
                    B,
                    M,
                    pmut,
                    tps, 
                    run,
                    compet
                    )
  {

  
  namerun <- paste( .compet, "K", K, '_fs', fec, '_fa',
                    fecasex, '_', probamating,'_G', G, '_b',
                    bsline, '_Ka', Ka, '_B', B,'M_', M,
                    '_pm', pmut, '_c', .c, '_',run,'.RData',
                    sep = "" )
  
  # gc() # clean unused memory
  # ani.record(reset = TRUE) # reset plot memory
  # dev.control('enable')
  
  Ninit <- K * fec
  

  
  # preparing data gathering
  
  INVASION_tot <- numeric( length = tps ) # pop size gathered once adults have settled and survived competition, before reproduction
  INVASION_asex <- numeric( length = tps )
  EDGE_t <- numeric( length = tps ) # position of the pop edge in time (defined as the X coordinates where all Y patches contain individuals)
  
  
  # Initialization
  
  # input: K, Ninit, G
  # output: popgenome, sex, newbabyX, newbabyY
  
  TEMP <- genesis( K = K, fec = fec, Xinit, Yinit, Ninit = Ninit, G = G) 
  
  popgenome <- TEMP$popgenome
  sex <- TEMP$sex
  newbabyX <- TEMP$newbabyX
  newbabyY <- TEMP$newbabyY
  repro  <- TEMP$repro
  popcloneline <- TEMP$popcloneline
  popclonalorigin <- TEMP$popclonalorigin
  

  
  #   # saves the reproductive population size per patch, after the first ever round of dispersion
  #   popsize_past <- tapply( !is.na(coordalive), coordalive, sum )
  #   popXY_past <- popXY
  #   popX_past <- popX
  #   popY_past <- popY
  
  
  
  
  # number of asexuals ever to arise through mutation
  count <- 0
  
  ######### --- enter loop --- ########
  
  for ( t in 1: tps ){
    
    # 0. Initial dispersal
    # 
    # input: newbabyX, newbabyY, 
    # output: newpopX, newpopY, newpopXY, coordalive, haspartner
    
    TEMP <- diaspora( newbabyX = newbabyX, newbabyY = newbabyY, 
                      dispkernel = dispkernel, sex = sex, repro = repro,
                      Xdim = Xdim, Ydim = Ydim)
    
    popXY <- TEMP$newpopXY
    popY <- TEMP$newpopY
    popX <- TEMP$newpopX
    coordalive <- TEMP$coordalive
    haspartner <- TEMP$haspartner
    
    
    # 1. Competition
    # 
    # input: haspartner, coordalive, popgenome, K, G
    # output: popsurvival
    
    ## to switch on genotype competition
    popsurvival <- competition( compet = compet,
                                haspartner = haspartner, 
                                coordalive = coordalive, 
                                popgenome = popgenome, K = K, G = G)
    

    # numer of patches with surviving individuals; and with more than half of their population asexual
    INVASION_tot[ t ] <- length( unique( coordalive[ popsurvival == 1 ] ) )
    temptable <- table( coordalive[ popsurvival == 1 ], 
                        factor( repro[ popsurvival == 1 ], levels = c('a', 's') ) )
    INVASION_asex[ t ] <- sum( temptable[,1]/temptable[,2] > .5 )
    
    xmax_t <- max( popX[ popsurvival == 1 ] ) # the furthest point the pop reaches
    while( EDGE_t[t] == 0 ){
      if( length( unique( popY[ popsurvival == 1 ][ popX[ popsurvival == 1] == xmax_t ] )) == Ydim ){
        EDGE_t[t] <- xmax_t
      } 
    xmax_t <- xmax_t - 1
    }
    
    
    # 2. Reproduction: Allee effects and offspring production - meetic
    # 
    # input: sex, coordalive, parameters Allee effect, popgenome, popsurvival
    # output: babyX, babyY, babysex, babygenome
    
    TEMP <- meetic( sex = sex, coordalive = coordalive, K = K, 
                    probamating = probamating, fec = fec, 
                    fecasex = fecasex, popgenome = popgenome, 
                    popsurvival = popsurvival, repro = repro, popXY = popXY, 
                    popX = popX, popY = popY, G = G, popcloneline = popcloneline,
                    popclonalorigin = popclonalorigin )
    
    if( !is.list( TEMP )) break
    
    babysex <- TEMP$babysex
    babyX <- TEMP$babyX
    babyY <- TEMP$babyY
    babygenome <- TEMP$babygenome
    babyrepro <- TEMP$babyrepro
    babycloneline <- TEMP$babycloneline
    babyclonalorigin <- TEMP$babyclonalorigin
    
    # 4. Survival of inbreeding to juvenile stage - stillbornn
    # 
    # input: babygenome, babysex, babyX, babyY
    # output: newbabygenome, newbabysex, newbabyX, newbabyY
    # 
    
    TEMP <- stillborn( babygenome = babygenome, babysex = babysex, 
                       babyX = babyX, babyY = babyY, B = B, M = M, Ka = Ka, 
                       bsline = bsline, pmut = pmut, babyrepro = babyrepro, 
                       babycloneline = babycloneline, babyclonalorigin = babyclonalorigin,
                       count = count, G = G, t = t)
    
    sex <- TEMP$newbabysex
    popgenome <- as.matrix( TEMP$newbabygenome, drop = FALSE )
    newbabyX <- TEMP$newbabyX
    newbabyY <- TEMP$newbabyY
    popcloneline <- TEMP$newbabycloneline
    popclonalorigin <- TEMP$newbabyclonalorigin
    count <- TEMP$count
    repro <- TEMP$newbabyrepro
    
    

    
    
    # Plot density of total, and asexual populations.
    
    if ( round( t/10 ) == t/10 ){
    
    # popsize <- tapply( popsurvival, coordalive, sum )
    # asexpopsize <- tapply( popsurvival[ repro == "a"], coordalive[ repro == "a"], sum )
    # dev.control('enable')
    # print( plotdensity( K = K, Xinit = Xinit, Xdim = Xdim, Ydim = Ydim,
    #                     popsize = popsize, asexpopsize = asexpopsize,
    #                     popX = popX, popY = popY, popXY = popXY,
    #                     plotname = 'Density per patch' ) )
    # ani.record()
    # 
    
      source("analysis_cloneposition.R", local = TRUE)
    
    }
    
    

    
    # Limit cases
    if( length( popXY ) == 0 | sum ( is.na(coordalive) ) == length(coordalive) ){
      popsurvival <- NA
      break
    } 
    if( sum( repro == "s" ) / sum( repro == "a") < 0.05 ){
      popsurvival <- competition( compet = compet, 
                                  haspartner = haspartner, 
                                  coordalive = coordalive, popgenome = popgenome,
                                  K = K, G = G)
      break
    } 
    
    ######### --- loop --- ########
    
    print( paste( run, ':', t ) )
    
    
    
    
  }
  
  
  
  
  # namerun <- paste( .vcompet[run], "K", .vK[run], '_fs', .vfec[run], '_fa',
  #                   .vfecasex[run], '_', .vprobamating[run],'_G', .vG[run], '_b',
  #                   .vbsline[run], '_Ka', .vKa[run], '_B', .vB[run],'M_', .vM[run],
  #                   '_pm', .vpmut[run], '_c', .vc[run], '_a', .va[run], '.RData',
  #                   sep = "" )
 # backup <- paste( 'C:/Users/Anais Tilquin/Desktop/SimulLaCiotat/', namerun, sep = "" )

 save( list = ls( all.names = TRUE ), file = namerun, envir = environment() )
 # save( list = ls( all.names = TRUE ), file = backup, envir = environment() )
  
 
 
 
 source("analysis_cloneposition.R")
 

 
  rm( list = setdiff( ls(), lsf.str() ) ) # removes all variables but functions
}

library(compiler)
cmpfun(master)