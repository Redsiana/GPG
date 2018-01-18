### where viability selection occurs depending on the levels of homozygosity

# INPUT: babygenome, babysex, babyX, babyY, G
# OUTPUT: newbabygenome, newbabysex, newbabyX, newbabyY 

# the relationship between heterozygosity and survival is:
# curve(expr = Ka /( 1 + 1*exp( -B*(x-M) ) ), from=0, to=1)

stillborn <- function( babygenome, 
                       babysex, 
                       babyX, 
                       babyY, 
                       B, 
                       M, 
                       Ka, 
                       bsline, 
                       pmut, 
                       babyrepro,
                       babycloneline,
                       babyclonalorigin,
                       count,
                       G,
                       t){
  
  # each juvenile's proportion of heterozygous loci
  heteroz <- apply( babygenome == 1, 1 , sum ) / G
  
  # each juvenile's survival probability to inbreeding depression
  inbreeding <- Ka /( 1 + 1*exp( -B*(heteroz-M) )) + bsline
  
  # each juvenile's actual survival
  babysurvival <- as.logical( mapply( FUN = rbinom, prob = inbreeding, size = 1, n = 1 ) )
  
  newbabygenome <- babygenome[ babysurvival, ]
  newbabysex <- babysex[ babysurvival ]
  newbabyX <- babyX[ babysurvival ]
  newbabyY <- babyY[ babysurvival ]
  newbabyrepro <- babyrepro[ babysurvival ]
  newbabycloneline <- babycloneline[ babysurvival ]
  newbabyclonalorigin <- babyclonalorigin[ babysurvival, ]
  
  nmutants <- rbinom(n = 1, size = sum( newbabysex == "fem" & newbabyrepro == "s" ), prob = pmut)
  idmutants <- sample( 1: sum( newbabysex == "fem" & newbabyrepro == "s" ), nmutants)
  idline <- seq((count + 1), (count + nmutants))
  idXorigin <- newbabyX[ newbabysex == 'fem' & newbabyrepro == "s"][ idmutants ]
  idtime <- rep(t, length(idmutants))
  count <- count + nmutants
  
  newbabycloneline[ newbabysex == 'fem' & newbabyrepro == "s"][ idmutants ] <- idline
  newbabyclonalorigin[ newbabysex == 'fem' & newbabyrepro == "s" ,][ idmutants, 1 ] <- idXorigin
  newbabyclonalorigin[ newbabysex == 'fem' & newbabyrepro == "s" ,][ idmutants, 2 ] <- idtime
  newbabyrepro[ newbabysex == 'fem' & newbabyrepro == "s"][ idmutants ] <- "a"
  
  return( list( newbabygenome = newbabygenome, newbabysex = newbabysex, 
                newbabyX = newbabyX, newbabyY = newbabyY, 
                newbabyrepro = newbabyrepro, newbabycloneline = newbabycloneline,
                newbabyclonalorigin = newbabyclonalorigin, count = count))
  
}

