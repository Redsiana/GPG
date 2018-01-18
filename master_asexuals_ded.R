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
  
#   # saves the reproductive population size per patch, after the first ever round of dispersion
#   popsize_past <- tapply( !is.na(coordalive), coordalive, sum )
#   popXY_past <- popXY
#   popX_past <- popX
#   popY_past <- popY
  
  # preparing data gathering
  
  INVASION_tot <- numeric( length = tps ) # pop size gathered once adults have settled and survived competition, before reproduction
  INVASION_asex <- numeric( length = tps )
  EDGE_t <- numeric( length = tps ) # position of the pop edge in time (defined as the X coordinates where all Y patches contain individuals)
  
  # number of asexuals ever to arise through mutation
  count <- 0
  
  ######### --- enter loop --- ########
  
  for ( t in 1: tps ){
    
    
    # 1. Competition
    # 
    # input: haspartner, coordalive, popgenome, K, G
    # output: popsurvival
    
    ## to switch on genotype competition
    popsurvival <- competition( haspartner = haspartner, 
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
    
      source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/analysis_cloneposition.R", local = TRUE)
    
    }
    
    
    
    
    
    # 5. Dispersal - voyage
    # 
    # input: newbabyX, newbabyY, dispersal kernel, sex, repro
    # output: newpopX, newpopY, newpopXY, coordalive, haspartner
    
    # TEMP <- diaspora( newbabyX = newbabyX, newbabyY = newbabyY, short_dist_sd = short_dist_sd, plongdist = plongdist, long_dist_sd = long_dist_sd, long_dist_mean = long_dist_mean)
    
    TEMP <- diaspora( newbabyX = newbabyX, newbabyY = newbabyY, 
                      dispkernel = dispkernel, sex = sex, repro = repro,
                      Xdim = Xdim, Ydim = Ydim)
    
    popXY <- TEMP$newpopXY
    popY <- TEMP$newpopY
    popX <- TEMP$newpopX
    coordalive <- TEMP$coordalive
    haspartner <- TEMP$haspartner
    
    # Limit cases
    if( length( popXY ) == 0 | sum ( is.na(coordalive) ) == length(coordalive) ){
      popsurvival <- NA
      break
    } 
    if( sum( repro == "s" ) / sum( repro == "a") < 0.1 ){
      popsurvival <- competition( haspartner = haspartner, 
                                  coordalive = coordalive, popgenome = popgenome,
                                  K = K, G = G)
      break
    } 
    
    ######### --- loop --- ########
    
    print( paste( run, ':', t ) )
    
  }
  
  if( t == tps ) popsurvival <- competition( haspartner = haspartner, 
                                             coordalive = coordalive, 
                                             popgenome = popgenome, K = K, G = G)
  


  
  
  
  # namerun <- paste( .vcompet[run], "K", .vK[run], '_fs', .vfec[run], '_fa',
  #                   .vfecasex[run], '_', .vprobamating[run],'_G', .vG[run], '_b',
  #                   .vbsline[run], '_Ka', .vKa[run], '_B', .vB[run],'M_', .vM[run],
  #                   '_pm', .vpmut[run], '_c', .vc[run], '_a', .va[run], '.RData',
  #                   sep = "" )
 # backup <- paste( 'C:/Users/Anais Tilquin/Desktop/SimulLaCiotat/', namerun, sep = "" )

 save( list = ls( all.names = TRUE ), file = namerun, envir = environment() )
 # save( list = ls( all.names = TRUE ), file = backup, envir = environment() )
  
 
 
 
 
 
 
 
 sorted_clone_prevalence <- sort( table(popcloneline)[-1])# gives table with clone number and number of individuals
 main_clones_names <- as.numeric( names( tail(sorted_clone_prevalence, 5) ))
 
 clone_bars <- as.numeric( table(popX[ !(popcloneline %in% c(0,main_clones_names)) ]) )
 clone_patches <- as.numeric( names( table(popX[ !(popcloneline %in% c(0,main_clones_names)) ])) )
 clone_idcol <- as.factor( rep( "other cl", length(clone_bars)) )
 df <- data.frame(x = clone_patches, y = clone_bars, fac = clone_idcol )
 
 x_emergence <- numeric(length = length(main_clones_names))
 g_emergence <- numeric(length = length(main_clones_names))
 g_edge <- numeric(length = length(main_clones_names))
 label_clone <- numeric(length = length(main_clones_names))
 
 for(clone in 1:length(main_clones_names)){
   clone_bars <- as.numeric( table(popX[ popcloneline == main_clones_names[clone]]) )
   clone_patches <- as.numeric( names( table(popX[ popcloneline == main_clones_names[clone]])) )
   clone_idcol <- as.factor( rep( main_clones_names[clone], length(clone_bars)) )
   df <- rbind(df, data.frame(x = clone_patches, y = clone_bars, fac = clone_idcol ))
   
   x_emergence[clone] <- unique( popclonalorigin[popcloneline %in% main_clones_names[clone],1] )
   g_emergence[clone] <- unique( popclonalorigin[popcloneline %in% main_clones_names[clone],2] )
   g_edge[clone] <- EDGE_t[g_emergence[clone]]
   label_clone[clone] <- paste(g_emergence[clone], ":", x_emergence[clone], '-', g_edge[clone])
 }
 
 sexual_bars <- as.numeric( table(popX[ popcloneline == 0]) )
 sexual_patches <- as.numeric( names( table(popX[ popcloneline == 0])) )
 sexual_idcol <- rep("sex", length(sexual_bars))
 df <- rbind( df, data.frame(x = sexual_patches, y = sexual_bars, fac = sexual_idcol ) )
 
 
 
 # this df has 3 columns, the X-coordinates, the number of individuals, the identity (5 main clones, other cl or sex)
 mypalette <- c( brewer.pal((length(main_clones_names)+1),  "YlOrRd"), "black" )
 
 labels_plot <- c("other cl", label_clone, "sex")
 
 ggplot(data=df, aes(x=x, y=y, fill=fac, order=fac)) +
   geom_bar(stat="identity" ) +   theme_minimal() +
   scale_fill_manual(values= mypalette, labels = labels_plot) +
   geom_vline(xintercept = jitter(x_emergence), col=mypalette[2:(length(main_clones_names)+1)], size = .8 )+
   ggtitle(paste("pop composition after", t, "generations")) +
   theme(aspect.ratio=5/10)
 
 # this prints, for each of the 5 main clones, the generation at which it appeared, the x coordinates, and the coordinates of the edge back then
 
 ggsave(filename = paste(namerun, ".jpg", sep=""))
 
 
 

 
  rm( list = setdiff( ls(), lsf.str() ) ) # removes all variables but functions
}

library(compiler)
cmpfun(master)