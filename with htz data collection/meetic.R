### where they find a soul mate

# INPUT: sex, coordalive, K, parameters Allee effect: probamating, theta, fecundity f, popXY, popX, popY, G
# OUTPUT: babysex, babyR, babyB, babygenome, babyrepro
# 

## NEW ARGUMENT: asexmixis = apo, auto, self
## CHECK LINE 79 - 87


meetic <- function( sex, 
                    K, 
                    probamating, 
                    theta, 
                    fec, 
                    fecasex, 
                    popgenome, 
                    popsurvival, 
                    repro, 
                    popXY, 
                    popX, 
                    popY,
                    G,
                    popcloneline,
                    popclonalorigin,
                    control,
                    asexmixis,
                    autonomous,
                    hermaphrodite,
                    resident_selfer,
                    facultative_parthenogen ){
  
  if(hermaphrodite){
    
    
    # table of number of mates per patch, EXCLUDING self as mate. For facultative selfer, they'll find out later
    m_per_patch <- table( repro[popsurvival==1], popXY[popsurvival==1])['s', , drop=F] - 1 # removing self as mate
    m_per_patch[ m_per_patch < 0] <- 0
    if ( length( m_per_patch ) == 1) names(m_per_patch) <- unique( popXY[popsurvival==1 & sex=='mal'] )
    
    # if( resident_selfer ){
    #   
    #   # table of number of mates per patch, counting self as mate
    #   m_per_patch <- table( repro[popsurvival==1], popXY[popsurvival==1])['s', , drop=F] 
    #   if ( length( m_per_patch ) == 1) names(m_per_patch) <- unique( popXY[popsurvival==1 & sex=='mal'] )
    #   
    # } else {
    #   
    #   # table of number of mates per patch, EXCLUDING self as mate
    #   m_per_patch <- table( repro[popsurvival==1], popXY[popsurvival==1])['s', , drop=F] - 1 # removing self as mate
    #   m_per_patch[ m_per_patch < 0] <- 0
    #   if ( length( m_per_patch ) == 1) names(m_per_patch) <- unique( popXY[popsurvival==1 & sex=='mal'] )
    #   
    # }
    
  } else {
    
    ## mate-finding Allee-effect
    if( sum( sex[popsurvival==1]=='mal' ) == 0 ) return( 'nomales' )
    
    # table of number of males per patch
    m_per_patch <- table( sex[popsurvival==1], popXY[popsurvival==1])['mal', , drop=F] 
    if ( length( m_per_patch ) == 1) names(m_per_patch) <- unique( popXY[popsurvival==1 & sex=='mal'] )
    
  }
  

  pmating <- probamating * popsurvival
  pmating[ popXY %in% colnames(m_per_patch)[m_per_patch==0] ] <- 0 # only females with mates on their patch have pmating, otherwise 0

  # case where we want sexuals to be facultative selfers, if they don't have a SEXUAL partner on the patch
  realized_selfing <- 0
  
  if( hermaphrodite & resident_selfer ){
    
    lonely_sexual <- table(popXY[ popsurvival==1 & repro=="s" ] ) == 1 # patches with only one sexual individual
    selfing_sexual <- names( popsurvival[popsurvival==1 & repro=="s"] ) %in% 
      names( lonely_sexual[lonely_sexual] ) # vector of SURVIVING indiv, with T if is a lonely sexual on a patch
    names(selfing_sexual) <- names( popsurvival[popsurvival==1 & repro=="s"] )
  #NB: selfing_sexual is vector of size of surviving sexuals, giving T for isolated sexuals
    realized_selfing <- selfing_sexual
    if(probamating!=1) realized_selfing[selfing_sexual] <- mapply( 
      FUN = rbinom, prob = probamating, size = 1, n = sum(selfing_sexual) )
  
  }
  
  
  # for all other sexual cases (won't have to self)
  # Only females get a 1 or 0 (so the vector is the size of the female pop), and males are picked by females randomly; saves some random draw for the males
  if(probamating==1){
    realized_mating <- pmating[ sex =="fem" & repro == "s" ]
  } else {
    realized_mating <- mapply( FUN = rbinom, prob = pmating[ sex =="fem" & repro == "s" ], size = 1, n = 1 ) 
  }
  
  # for asexual females that need the presence of a dude around
  
  ## this first version is if asexuals can only accept pollen from sexuals - not from fellow asexual hermaphrodites
  # if( autonomous == F & hermaphrodite == F ){
  #   if(probamating==1){
  #     realized_parthenogenesis <- pmating[ sex =="fem" & repro == "a" ] # pmating takes into account the presence of males in a patch
  #   } else {
  #     realized_parthenogenesis <- mapply( FUN = rbinom, prob = pmating[ sex =="fem" & repro == "a" ], size = 1, n = 1 ) 
  #   } 
  # }
  
  ## asexual hermaphrodites can accept pollen from fellow asexuals, or from sexuals
  if( autonomous == F & hermaphrodite == T ){
  
    if(probamating == 1){
      realized_parthenogenesis <- popsurvival[ sex =="fem" & repro == "a" ] # baseline: everyone who survived reproduces
      lonely <- table(popXY[popsurvival==1]) == 1 # patches with only one individual
      
      if( sum( lonely ) != 0 ){ # BUT if some individuals are lonely, they won't
        realized_parthenogenesis[ 
          names(realized_parthenogenesis) %in% names( lonely[lonely] ) ] <- 0
      }
      
    } else {
      realized_parthenogenesis <- popsurvival[ sex =="fem" & repro == "a" ] # baseline: everyone reproduces
      realized_parthenogenesis <- mapply( FUN = rbinom, prob = realized_parthenogenesis, size = 1, n = 1)
      
      if( sum( lonely ) != 0 ){ # BUT if some individuals are lonely, they won't
        realized_parthenogenesis[ 
          names(realized_parthenogenesis) %in% names( lonely[lonely] ) ] <- 0
      }
    }
  }
  
  if( autonomous == T){
    realized_parthenogenesis <- popsurvival[ sex =="fem" & repro == "a" ] # if you survive, you reproduce
  }
  
  
  if( sum(realized_mating) == 0 ) return( 'nomating')
  
  # evaluates condition early to speed up code
  reproductive.females <- which(sex == "fem" & repro == "s")[ realized_mating == 1 ]
  clonal.females <- which(repro == 'a')[realized_parthenogenesis == 1]
  
  # coordinates of the sexual mothers
  mum_patch <- popXY[ reproductive.females ]
  
  # total number of sexually + asexually produced offspring
  newpopsize <- fec * sum( realized_mating ) + 
    fecasex * length(clonal.females) + 
    fec * sum( realized_selfing )
  
  
  # for each sexual mother, picks a partner randomly among those on her patch
  # first: what partners are there on the patch?
  if( hermaphrodite ){
    
    male_patch <- popXY[popsurvival == 1 & repro == "s" ]
    
  } else {
    
    male_patch <- popXY[popsurvival == 1 & sex == "mal" ]
    
  }
  
  # second: which one will each female pick?
  if( hermaphrodite ){
    
    # we want the dad's id to be different from the mother's, because the residents can't self (desperation selfing treated separately)
    dad_id <- sapply( X = reproductive.females,
                      FUN = function(x) sample( 
                        as.character( # because if only one mate on patch, sample behaves pathologically
                          which( male_patch == popXY[x] )[ 
                            which( male_patch == popXY[x] ) !=x ]), # this excludes the mother's id
                        1 ))
    dad_id <- as.numeric( dad_id )
    
  } else { 
    
    dad_id <- sapply( X = mum_patch,
                      FUN = function(x) sample( which( male_patch == x ), 1 ))
  }
  
  # genotypes of babies produced sexually
  mumgenes <- popgenome[ reproductive.females , , drop = F]
  dadgenes <- popgenome[ dad_id, ]
  mendel_parents <- mumgenes + dadgenes 
  mendel <- mendel_parents[rep(1:nrow(mendel_parents), times = fec), ]
  
  mendel0 <- mendel == 0
  mendel4 <- mendel == 4
  mendel1 <- mendel == 1
  mendel3 <- mendel == 3
  mendel2 <- mendel == 2
  mendelmom <- mumgenes[rep(1:nrow(mumgenes), times = fec), ] == 1
  
  
  
  # REMARK: ( 1 baby per female ), fec times ; NOT ( fec baby ), nb of female times
  sex_babygenome <- matrix( nrow = nrow(mendel) , ncol = (G ))
  sex_babygenome[ mendel0 ] <- 0
  sex_babygenome[ mendel4 ] <- 2
  sex_babygenome[ mendel1 ] <- sample(0:1, size = sum( mendel1 ), replace = TRUE)
  sex_babygenome[ mendel3 ] <- sample(1:2, size = sum( mendel3 ), replace = TRUE)
  sex_babygenome[ mendel2 & mendelmom ] <- sample( c(0,1,1,2) , size = sum( mendel2 & mendelmom ), replace = TRUE )
  sex_babygenome[ mendel2 & !mendelmom ] <- 1
  
  # give the neutral mutation line number for control
  if(control == TRUE){
    inherited_sexcloneline <- c( 
      rep( popcloneline[ reproductive.females ], times = fec ),
      rep( popcloneline[ popsurvival == 1 & repro == "s" ][realized_selfing], times = fec ) )
    inherited_sexX <- c( 
      rep( popclonalorigin[ reproductive.females, 1 ], times = fec ),
      rep( popclonalorigin[ popsurvival == 1 & repro == "s" ,1][realized_selfing], times = fec ))
    inherited_sext <- c(
      rep( popclonalorigin[ reproductive.females, 2 ], times = fec ),
      rep( popclonalorigin[ popsurvival == 1 & repro == "s" ,2][realized_selfing], times = fec ))
  } else {
    inherited_sexcloneline <- rep(0, fec*( sum( realized_mating) + sum(realized_selfing)) )
    inherited_sexX <- rep(0, fec*( sum( realized_mating) + sum(realized_selfing)) )
    inherited_sext <- rep(0, fec*( sum( realized_mating) + sum(realized_selfing)) )
  }
  
  ## Genotypes of babies produced by selfing of lonely sexual facultative selfers
  selfed_babygenome <-  matrix( nrow = 0  , ncol = G )
  
  if( resident_selfer ){
    selfed_babygenome <- matrix( nrow = sum( realized_selfing ) * fec , ncol = G )
    
    if( sum( realized_selfing ) > 0 ){
      selfed_babygenome_unique <- popgenome[ popsurvival == 1 & repro == "s",][ realized_selfing, , drop = F ]
      selfed_babygenome <- selfed_babygenome_unique[rep(1:nrow(selfed_babygenome_unique), times = fec), ] # this line repeats the ENTIRE matrix fecasex fois
      htz_loci <- sum( selfed_babygenome == 1 )
      selfmix <- sample(c(0,1,1,2), htz_loci, replace=T) # this transforms htz loci into 00, 11, 01 or 10
      selfed_babygenome[selfed_babygenome == 1] <- selfmix 
    }
  }
  
  ## Asexual reproduction, if the asexuals are pure asexuals
  
  # if( facultative_parthenogen == F ){ ## let's forget about that for now!!
    if(asexmixis == "apo"){
      # genotypes of babies produced asexually
      asex_babygenome <- matrix( nrow = length( clonal.females ) * fecasex , ncol = G )
      if( length( clonal.females ) > 0 ){
        asex_babygenome_unique <- popgenome[ clonal.females, , drop = F ]
        asex_babygenome <- asex_babygenome_unique[rep(1:nrow(asex_babygenome_unique), times = fecasex), ] # this line repeats the ENTIRE matrix fecasex fois
      }
    }
    if(asexmixis == "auto"){
      # genotypes of babies produced asexually
      asex_babygenome <- matrix( nrow = length( clonal.females ) * fecasex , ncol = G )
      if( length( clonal.females ) > 0 ){
        asex_babygenome_unique <- popgenome[ clonal.females, , drop = F ]
        asex_babygenome <- asex_babygenome_unique[rep(1:nrow(asex_babygenome_unique), times = fecasex), ] # this line repeats the ENTIRE matrix fecasex fois
        htz_loci <- sum( asex_babygenome == 1 )
        automix <- sample(c(0,2), htz_loci, replace=T) # this transforms htz loci into 00 or 11
        asex_babygenome[asex_babygenome == 1] <- automix
      }
    }
    if(asexmixis == "self"){
      # genotypes of babies produced asexually
      asex_babygenome <- matrix( nrow = length( clonal.females ) * fecasex , ncol = G )
      if( length( clonal.females ) > 0 ){
        asex_babygenome_unique <- popgenome[ clonal.females, , drop = F ]
        asex_babygenome <- asex_babygenome_unique[rep(1:nrow(asex_babygenome_unique), times = fecasex), ] # this line repeats the ENTIRE matrix fecasex fois
        htz_loci <- sum( asex_babygenome == 1 )
        selfmix <- sample(c(0,1,1,2), htz_loci, replace=T) # this transforms htz loci into 00, 11, 01 or 10
        asex_babygenome[asex_babygenome == 1] <- selfmix
      }
    }
  # }
 #  if( facultative_parthenogen == T ){ # note that the asexuals have fecasex regardless of mating status
 #   
 #    # in pmating, survival is already taken into account
 #   realized_mating_asexuals <- pmating[ sex =="fem" & repro == "a" ]
 #   realized_parthenogenesis_asexuals <- (!realized_mating_asexuals) * popsurvival[sex=="fem"&repro=="a"]
 #   
 #   # evaluates condition early to speed up code
 #   reproductive.asexual.females <- which(sex == "fem" & repro == "a")[ realized_mating_asexuals == 1 ]
 #   clonal.asexual.females <- which(repro == 'a')[realized_parthenogenesis_asexuals == 1]
 #   
 #   # coordinates of the sexual mothers
 #   asexual_mum_patch <- popXY[ reproductive.asexual.females ]
 #   
 #   # total number of sexually + asexually produced offspring
 #   sexual.prod.asexuals <- fecasex * length(reproductive.asexual.females)
 #   asexual.prod.asexuals <- fecasex * length(clonal.asexual.females)
 #   
 #   male_patch <- popXY[popsurvival == 1 & sex == "mal" ]
 #   
 #   
 #   
 #   
 #   
 #   
 #   ################"
 #   ################"""
 #   ############"" THAT'S WHERE I HAVE A PROBLEM
 #   ##############################################
 #   ########################"
 #   
 #   dad_id <- sapply( X = asexual_mum_patch,
 #                     FUN = function(x) sample( which( male_patch == x ), 1 ))
 # 
 #  
 #   # genotypes of babies produced sexually
 #   mumgenes <- popgenome[ reproductive.females , , drop = F]
 #   dadgenes <- popgenome[ dad_id, ]
 #   mendel_parents <- mumgenes + dadgenes 
 #   mendel <- mendel_parents[rep(1:nrow(mendel_parents), times = fec), ]
 #   
 #   mendel0 <- mendel == 0
 #   mendel4 <- mendel == 4
 #   mendel1 <- mendel == 1
 #   mendel3 <- mendel == 3
 #   mendel2 <- mendel == 2
 #   mendelmom <- mumgenes[rep(1:nrow(mumgenes), times = fec), ] == 1
 #   
 #   
 #   
 #   # REMARK: ( 1 baby per female ), fec times ; NOT ( fec baby ), nb of female times
 #   sex_babygenome <- matrix( nrow = nrow(mendel) , ncol = (G ))
 #   sex_babygenome[ mendel0 ] <- 0
 #   sex_babygenome[ mendel4 ] <- 2
 #   sex_babygenome[ mendel1 ] <- sample(0:1, size = sum( mendel1 ), replace = TRUE)
 #   sex_babygenome[ mendel3 ] <- sample(1:2, size = sum( mendel3 ), replace = TRUE)
 #   sex_babygenome[ mendel2 & mendelmom ] <- sample( c(0,1,1,2) , size = sum( mendel2 & mendelmom ), replace = TRUE )
 #   sex_babygenome[ mendel2 & !mendelmom ] <- 1
 #    
 #    
 # }
  
  
  
  
  
  
  
  inherited_cloneline <- NULL
  inherited_X <- NULL
  inherited_t <- NULL
  if( length( clonal.females ) > 0 ){
    inherited_cloneline <- rep( popcloneline[ clonal.females ], times = fecasex )
    inherited_X <- rep( popclonalorigin[ clonal.females, 1 ], times = fecasex )
    inherited_t <- rep( popclonalorigin[ clonal.females, 2 ], times = fecasex )
  }
  
  babygenome <- rbind(sex_babygenome, selfed_babygenome, asex_babygenome)
  
  # sex of babies, sexuals then asexuals, in one vector
  if( hermaphrodite ){
    babysex <- rep("fem", fec * ( sum( realized_mating ) + sum( realized_selfing ) ) )
  } else {
    babysex <- sample( c("fem", "mal"), replace = T, size = fec * sum( realized_mating ) )
  }
  
  babysex <- c( babysex, rep("fem", fecasex * length( clonal.females )) ) 
  
  # reproductive mode of babies, sexuals then asexuals, in one vector
  babyrepro <- c( rep("s", fec*( sum( realized_mating ) + sum( realized_selfing ) ) ),
                  rep( "a", fecasex * length( clonal.females )) )
  babycloneline <- c(inherited_sexcloneline, inherited_cloneline )
  babyclonalorigin <- matrix(0, nrow = length(babyrepro), ncol = 2)
  babyclonalorigin[,1] <-  c(inherited_sexX, inherited_X )
  babyclonalorigin[,2] <-  c(inherited_sext, inherited_t )
  
  # babies birthplace, sexuals then asexuals, in one vector
  babyX <- c( rep( popX[  reproductive.females ], fec ),
              rep( popX[ popsurvival ==1 & repro == "s"][realized_selfing], fec),
              rep( popX[ clonal.females ], fecasex ) )
  
  babyY <- c( rep( popY[ reproductive.females ], fec ), 
              rep( popY[ popsurvival ==1 & repro == "s"][realized_selfing], fec),
              rep( popY[ clonal.females ], fecasex ) )
  
  
  
  
  return( list( babysex = babysex, babyX = babyX, babyY = babyY, 
                babygenome = babygenome, babyrepro = babyrepro, 
                babycloneline = babycloneline, babyclonalorigin = babyclonalorigin ) )
  
}

library(compiler)
meetic <- cmpfun(meetic)
