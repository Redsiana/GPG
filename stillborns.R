### where viability selection occurs depending on the levels of homozygosity

# INPUT: babygenome, babysex, babyX, babyY, G
# OUTPUT: newbabygenome, newbabysex, newbabyX, newbabyY 

# the relationship between heterozygosity and survival is:
# curve(expr = Ka /( 1 + 1*exp( -B*(x-M) ) ), from=0, to=1)
# 
# #### TO-DO : Still need to add the popclonal origin to the sexuals with neutral mutation
# so that analysis cloneposition works

stillborn <- function( babygenome, 
                       babysex, 
                       babyX, 
                       babyY, 
                       bsline, 
                       pmut, 
                       babyrepro,
                       babycloneline,
                       babyclonalorigin,
                       count,
                       G,
                       t,
                       control){
  
  
  # each juvenile's proportion of heterozygous loci
  heteroz <- apply( babygenome == 1, 1 , sum ) / G
  
  # each juvenile's survival probability to inbreeding depression
  # inbreeding <- Ka /( 1 + 1*exp( -B*(heteroz-M) )) + bsline
  inbreeding <- 1-(1-bsline)*((heteroz-1))^2
  
  
  # each juvenile's actual survival
  babysurvival <- as.logical( mapply( FUN = rbinom, prob = inbreeding, size = 1, n = 1 ) )
  
  newbabygenome <- babygenome[ babysurvival, ]
  newbabysex <- babysex[ babysurvival ]
  newbabyX <- babyX[ babysurvival ]
  newbabyY <- babyY[ babysurvival ]
  newbabyrepro <- babyrepro[ babysurvival ]
  newbabycloneline <- babycloneline[ babysurvival ]
  newbabyclonalorigin <- babyclonalorigin[ babysurvival, ]
  
  
  # instead of newbabyrepro == "s", we have newbabycloneline == 0 for control, 
  # and !control, since the ones with 0 are non mutated lineages
  # evaluating condition here speeds up code tremendously
  non.mutated.females = which(newbabysex == "fem" & newbabycloneline == 0)
  
  
  # instead of newbabyrepro == "s", we have newbabycloneline == 0 for control, 
  # since they're all sexuals. The ones with 0 are non mutated lineages
  nmutants <- rbinom(n = 1, size = length( non.mutated.females ), prob = pmut)
  if(nmutants > 0){
    # chose which female babies from non-mutated lineages will mutated
    idmutants <- sample( 1:length( non.mutated.females ), nmutants)
    # gives a new number to each new mutated lineage appeared this generation
    idline <- seq((count + 1), (count + nmutants))
    # the X coordinates on which those babies mutate
    idXorigin <- newbabyX[ non.mutated.females][ idmutants ]
    # time at which those babies mutate
    idtime <- rep(t, length(idmutants))
  
    newbabyclonalorigin[ non.mutated.females ,][ idmutants, 1 ] <- idXorigin
    newbabyclonalorigin[ non.mutated.females ,][ idmutants, 2 ] <- idtime
    # finally, change their cloneline from 0 to their lineage id
    newbabycloneline[ non.mutated.females][ idmutants ] <- idline
    # in the control the mutated baby stays sexual, if not control (!control), they become asexual
    if(!control) newbabyrepro[ non.mutated.females ][ idmutants ] <- "a"
  
  }
  
  count <- count + nmutants
  
  
  
  
  return( list( newbabygenome = newbabygenome, newbabysex = newbabysex, 
                newbabyX = newbabyX, newbabyY = newbabyY, 
                newbabyrepro = newbabyrepro, newbabycloneline = newbabycloneline,
                newbabyclonalorigin = newbabyclonalorigin, count = count))
  
}

library(compiler)
stillborn <- cmpfun(stillborn)

