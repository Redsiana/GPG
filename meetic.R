### where they find a soul mate

# INPUT: sex, coordalive, K, parameters Allee effect: probamating, theta, fecundity f, popXY, popX, popY, G
# OUTPUT: babysex, babyR, babyB, babygenome, babyrepro
# 

## NEW ARGUMENT: asexmixis = apo, auto, self
## CHECK LINE 79 - 87


meetic <- function( sex, 
                    coordalive, 
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
                    autonomous){
  
  ## mate-finding Allee-effect
  
  if( sum( sex[popsurvival==1]=='mal' ) == 0 ) return( 'nomales' )
  
  # table of number of males per patch
  m_per_patch <- table( sex[popsurvival==1], coordalive[popsurvival==1])['mal', , drop=F] 
  if ( length( m_per_patch ) == 1) names(m_per_patch) <- unique( coordalive[popsurvival==1 & sex=='mal'] )
  
  
  # if ( is.numeric(probamating) ) {
  pmating <- probamating * popsurvival
  pmating[ coordalive %in% colnames(m_per_patch)[m_per_patch==0] ] <- 0 # only females with males on their patch have pmating, otherwise 0
  # }
  # if ( probamating == "allee") {
  #   if ( missing('theta') ) {
  #     theta <- 3 / ( 0.5 * K ) # this theta gives such a curve for mating probability:
  #     # curve(expr= 1 - exp ( - 6 * x ), from=0, to=1, xlab='percent of K that are males', ylab='female mating proba')
  #   }
  #   pmating <-  (1 - exp ( - theta * m_per_patch[popXY] )) 
  #   pmating <- pmating * popsurvival # removing the females that didn't survive the competition earlier
  #   pmating[is.na(pmating)] <- 0 # otherwise rbinom errors over the NA, see if can't code them as 0 from the start
  # }
  
  
  # selfer <- coordalive %in% colnames(m_per_patch)[m_per_patch==0]
  # selfer <- selfer
  
  # Only females get a 1 or 0 (so the vector is the size of the female pop), and males are picked by females randomly; saves some random draw for the males
  if(probamating==1){
    realized_mating <- pmating[ sex =="fem" & repro == "s" ]
  } else {
    realized_mating <- mapply( FUN = rbinom, prob = pmating[ sex =="fem" & repro == "s" ], size = 1, n = 1 ) 
  }
  
  # for asexual females that need the presence of a dude around
  if( autonomous == F){
    if(probamating==1){
      realized_parthenogenesis <- pmating[ sex =="fem" & repro == "a" ]
    } else {
      realized_parthenogenesis <- mapply( FUN = rbinom, prob = pmating[ sex =="fem" & repro == "a" ], size = 1, n = 1 ) 
    } 
  } else {
    realized_parthenogenesis <- popsurvival[ sex =="fem" & repro == "a" ] # if you survive, you reproduce
  }
  
  
  if( sum(realized_mating) == 0 ) return( 'nomating')
  
  # evaluates condition early to speed up code
  reproductive.females <- which(sex == "fem" & repro == "s")[ realized_mating == 1 ]
  clonal.females <- which(repro == 'a')[realized_parthenogenesis == 1]
  
  # coordinates of the sexual mothers
  mum_patch <- coordalive[ reproductive.females ]
  mum_patch2 =  coordalive[ reproductive.females ]
  
  # total number of sexually + asexually produced offspring
  newpopsize <- fec * sum( realized_mating ) + fecasex * length(clonal.females)
  
  
  # for each sexual mother, picks a partner randomly among those on her patch
  male_patch <- popXY[popsurvival == 1 & sex == "mal" ]
  dad_id <- sapply( X = mum_patch,
                    FUN = function(x) sample( which( male_patch == x ), 1 ))
  
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
    inherited_sexcloneline <- rep( popcloneline[ reproductive.females ], times = fec )
    inherited_sexX <- rep( popclonalorigin[ reproductive.females, 1 ], times = fec )
    inherited_sext <- rep( popclonalorigin[ reproductive.females, 2 ], times = fec )
  } else {
    inherited_sexcloneline <- rep(0, fec*sum( realized_mating) )
    inherited_sexX <- rep(0, fec*sum( realized_mating) )
    inherited_sext <- rep(0, fec*sum( realized_mating) )
  }
  
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
  
  inherited_cloneline <- NULL
  inherited_X <- NULL
  inherited_t <- NULL
  if( length( clonal.females ) > 0 ){
    inherited_cloneline <- rep( popcloneline[ clonal.females ], times = fecasex )
    inherited_X <- rep( popclonalorigin[ clonal.females, 1 ], times = fecasex )
    inherited_t <- rep( popclonalorigin[ clonal.females, 2 ], times = fecasex )
  }
  
  babygenome <- rbind(sex_babygenome, asex_babygenome)
  
  # sex of babies, sexuals then asexuals, in one vector
  babysex <- sample( c("fem", "mal"), replace = T, size = fec * sum( realized_mating ) )
  babysex <- c( babysex, rep("fem", fecasex * length( clonal.females )) ) 
  
  # reproductive mode of babies, sexuals then asexuals, in one vector
  babyrepro <- c( rep("s", fec*sum( realized_mating ) ), rep( "a", fecasex * length( clonal.females )) )
  babycloneline <- c(inherited_sexcloneline, inherited_cloneline )
  babyclonalorigin <- matrix(0, nrow = length(babyrepro), ncol = 2)
  babyclonalorigin[,1] <-  c(inherited_sexX, inherited_X )
  babyclonalorigin[,2] <-  c(inherited_sext, inherited_t )
  
  # babies birthplace, sexuals then asexuals, in one vector
  babyX <- c( rep( popX[  reproductive.females ], fec ), rep( popX[ clonal.females ], fecasex ) )
  babyY <- c( rep( popY[ reproductive.females ], fec ), rep( popY[ clonal.females ], fecasex ) )
  
  
  
  
  return( list( babysex = babysex, babyX = babyX, babyY = babyY, 
                babygenome = babygenome, babyrepro = babyrepro, 
                babycloneline = babycloneline, babyclonalorigin = babyclonalorigin ) )
  
}

library(compiler)
meetic <- cmpfun(meetic)
