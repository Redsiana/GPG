master <- function( c,
                    mean_distance,
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
                    compet,
                    plot,
                    control,
                    autonomous,
                    asexmixis
                    )
  {

  
  namerun1 <- paste( compet, "K", K, '_fs', fec, '_fa',
                    fecasex, '_', probamating,'_G', G, '_b',
                    bsline, '_pm', pmut, '_d', mean_distance, '_c', c, '_',run,'_T1.RData',
                    sep = "" )
  namerun2 <- paste( compet, "K", K, '_fs', fec, '_fa',
                     fecasex, '_', probamating,'_G', G, '_b',
                     bsline, '_pm', pmut, '_d', mean_distance, '_c', c, '_',run,'_T2.RData',
                     sep = "" )
  
  # gc() # clean unused memory
  # ani.record(reset = TRUE) # reset plot memory
  # dev.control('enable')
  
  Ninit <- K * fec
  

  
  # preparing data gathering
  
  INVASION_tot <- numeric( length = 1 ) # pop size gathered once adults have settled and survived competition, before reproduction
  INVASION_sex <- numeric( length = 1 )
  EDGE_t <- NA # position of the pop edge in time (defined as the X coordinates where all Y patches contain individuals)
  EDGE_sex_t <- numeric( length = 1)
  pureness_sex_t <- numeric( length = 1)
  
  t_full_corridor <- numeric() # the time at which the entire pop reaches the end of the corridor (all patches full)
  notfull_boolean <- T # to detect only the first generation at which it gets full
  
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
  
  # count to attribute colors to clones in the plotting of cloneposition
  count_col <- 1
  old_clones_colors = NA
  old_clones_names = NA
  
  # give some random super high limit, but will be lowered when T2 is calculated (as twice the time it took to fill the corridor)
  T2 <- 10000
  t <- 1
  ######### --- enter loop --- ########
  
  while ( t < T2 ){
    bug = 0
    
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
    
    bug = 1
    
    # 1. Competition
    # 
    # input: haspartner, coordalive, popgenome, K, G
    # output: popsurvival
    
    ## to switch on genotype competition
    popsurvival <- competition( compet = compet,
                                haspartner = haspartner, 
                                coordalive = coordalive, 
                                popgenome = popgenome, K = K, G = G)
    
    if(is.character(popsurvival)) break

    
    
    
    bug = 2
    
    # 2. Reproduction: Allee effects and offspring production - meetic
    # 
    # input: sex, coordalive, parameters Allee effect, popgenome, popsurvival
    # output: babyX, babyY, babysex, babygenome
    
    TEMP <- meetic( sex = sex, coordalive = coordalive, K = K, 
                    probamating = probamating, fec = fec, 
                    fecasex = fecasex, popgenome = popgenome, 
                    popsurvival = popsurvival, repro = repro, popXY = popXY, 
                    popX = popX, popY = popY, G = G, popcloneline = popcloneline,
                    popclonalorigin = popclonalorigin,
                    control = control, 
                    asexmixis = asexmixis,
                    autonomous = autonomous)
    
    if( !is.list( TEMP )) break
    
    babysex <- TEMP$babysex
    babyX <- TEMP$babyX
    babyY <- TEMP$babyY
    babygenome <- TEMP$babygenome
    babyrepro <- TEMP$babyrepro
    babycloneline <- TEMP$babycloneline
    babyclonalorigin <- TEMP$babyclonalorigin
    
    bug = 3
    
    # 4. Survival of inbreeding to juvenile stage - stillbornn
    # 
    # input: babygenome, babysex, babyX, babyY
    # output: newbabygenome, newbabysex, newbabyX, newbabyY
    # 
    
    TEMP <- stillborn( babygenome = babygenome, babysex = babysex, 
                       babyX = babyX, babyY = babyY,  
                       bsline = bsline, pmut = pmut, babyrepro = babyrepro, 
                       babycloneline = babycloneline, babyclonalorigin = babyclonalorigin,
                       count = count, G = G, t = t, control = control)
    
    sex <- TEMP$newbabysex
    popgenome <- as.matrix( TEMP$newbabygenome, drop = FALSE )
    newbabyX <- TEMP$newbabyX
    newbabyY <- TEMP$newbabyY
    popcloneline <- TEMP$newbabycloneline
    popclonalorigin <- TEMP$newbabyclonalorigin
    count <- TEMP$count
    repro <- TEMP$newbabyrepro
    
    
    bug = 4
    
    

    
    
    # numer of individuals of either types
    INVASION_tot[ t ] <- length( popcloneline )
    INVASION_sex[ t ] <- sum( popcloneline == 0 )
    
    xmax_t <- max( newbabyX ) # the furthest point the pop reaches
    while( is.na( EDGE_t[t]) & xmax_t > 0){ # selects as EDGE the farthest point when every one of the Y patches is colonized
      if( length( unique( popY[ popsurvival == 1 ][ popX[ popsurvival == 1] == xmax_t ] )) == Ydim ){
        EDGE_t[t] <- xmax_t
      } 
      xmax_t <- xmax_t - 1
    }
    
    nb_sexuals <- sum(popcloneline == 0 )
    # X coordinate before which 95% of the sexual population is
    EDGE_sex_t[t] <- sort( newbabyX[ popcloneline == 0 ] )[ round(nb_sexuals*0.95) ]
    # fraction of sexuals within this area (to know if there are clusters of asex)
    pureness_sex_t[t] <- sum( newbabyX[popcloneline == 0 ] <= EDGE_sex_t[t] ) / 
      sum( newbabyX <= EDGE_sex_t[t] )

    
    # Plot density of total, and asexual populations.
    
    if ( plot == TRUE & round( t/10 ) == t/10 ){
      
      # popsize <- tapply( popsurvival, coordalive, sum )
      # asexpopsize <- tapply( popsurvival[ repro == "a"], coordalive[ repro == "a"], sum )
      # dev.control('enable')
      # print( plotdensity( K = K, Xinit = Xinit, Xdim = Xdim, Ydim = Ydim,
      #                     popsize = popsize, asexpopsize = asexpopsize,
      #                     popX = popX, popY = popY, popXY = popXY,
      #                     plotname = 'Density per patch' ) )
      # ani.record()
      # 
      TEMP <- analysis_cloneposition( popcloneline = popcloneline, 
                                      popclonalorigin = popclonalorigin, 
                                      newbabyX = newbabyX, 
                                      namerun = namerun2, 
                                      t = t,
                                      old_clones_names = old_clones_names,
                                      old_clones_colors = old_clones_colors,
                                      count_col = count_col, 
                                      EDGE_t = EDGE_t)
      
      count_col <- TEMP$count_col
      old_clones_colors <- TEMP$main_clones_colors
      old_clones_names <- TEMP$main_clones_names
    }
    
    
    # when the population reaches the end of the corridor
    if( notfull_boolean ){
      if( EDGE_t[t] == (Xdim-1) ){
      T1 <- t
      save( list = ls( all.names = TRUE ), file = namerun1, envir = environment() )
      T2 <- 2*T1
      notfull_boolean <- FALSE
      
      analysis_cloneposition( popcloneline = popcloneline, 
                              popclonalorigin = popclonalorigin, 
                              newbabyX = newbabyX, 
                              namerun = namerun1, 
                              t = t,
                              old_clones_names = old_clones_names,
                              old_clones_colors = old_clones_colors,
                              count_col = count_col,
                              EDGE_t = EDGE_t)
      
      }
    }
    
    t <- t+1
    
    # Limit cases
    if( length( popXY ) == 0 | sum ( is.na(coordalive) ) == length(coordalive) ){
      popsurvival <- NA
      break
    } 
    if( sum( repro == "s" ) / sum( repro == "a") < 0.05 ){
      break
    }
    # if( EDGE_sex_t[t] == Xdim ) break
    
    ######### --- loop --- ########
    
    bug = 5
    
    print( paste( run, ':', t ) )
    
  }
  
  res <- c(max(EDGE_sex_t), #                                 max X position attained by sexuals
           max( which(EDGE_sex_t==max(EDGE_sex_t))), #        time when they attained it
           pureness_sex_t[ mean(c( max( which(EDGE_sex_t==max(EDGE_sex_t))), t)) ], #  pureness halfway between recession point and end simulation
           max(INVASION_sex), #                               max number of sexuals attained in the pop
           max( which(INVASION_sex==max(INVASION_sex))), #    time when they attain it
           INVASION_tot[t-1], #                                 how many individuals when simulation ended
           INVASION_sex[t-1], #                                 how many sexuals when simulation ended
           EDGE_t[t-1], #                                       how far the entire pop was when simulation ended
           EDGE_sex_t[t-1], #                                   how far the sexual pop was when simulation ended
           (t-1), #                                             when the simulation ended
           namerun2
           )
  
  
  # namerun <- paste( .vcompet[run], "K", .vK[run], '_fs', .vfec[run], '_fa',
  #                   .vfecasex[run], '_', .vprobamating[run],'_G', .vG[run], '_b',
  #                   .vbsline[run], '_Ka', .vKa[run], '_B', .vB[run],'M_', .vM[run],
  #                   '_pm', .vpmut[run], '_c', .vc[run], '_a', .va[run], '.RData',
  #                   sep = "" )
 # backup <- paste( 'C:/Users/Anais Tilquin/Desktop/SimulLaCiotat/', namerun, sep = "" )

 save( list = ls( all.names = TRUE ), file = namerun2, envir = environment() )
 # save( list = ls( all.names = TRUE ), file = backup, envir = environment() )
  
 
 analysis_cloneposition( popcloneline = popcloneline, 
                         popclonalorigin = popclonalorigin, 
                         newbabyX = newbabyX, 
                         namerun = namerun2, 
                         t = (t-1),
                         old_clones_names = old_clones_names,
                         old_clones_colors = old_clones_colors,
                         count_col = count_col,
                         EDGE_t = EDGE_t)
 

 return(res = res)
 
  # rm( list = setdiff( ls(), lsf.str() ) ) # removes all variables but functions
  
  
}

library(compiler)
master <- cmpfun(master)