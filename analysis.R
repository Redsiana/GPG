
plot(EDGE_t, type='l', col="green", main = "temporal dynamics")
abline(h=10, col="grey")
lines(pureness_sex_t*10, col="red")
points(EDGE_sex_t)
abline(v=max( which(EDGE_sex_t==max(EDGE_sex_t))) )





















############# OLD ANALYSIS - LOTS OF STUFF TO RECOVER

.K = 20
.fec = 6
.fecasex = 3
.G = 3
.probamating = 1 
.compet = '_ded_'
.mean_distance = 1
.c = 1
.bsline = 0.6 
.Ka = 0.4 
.B = 20
.M = 0.3 
.pmut = 0.005
.tps = 20
.s = 2/3 


subDir = paste("K", .K, '_',.compet, '_fs', .fec, '_fa', .fecasex, '_', .probamating,
               '_G', .G, '_bKaBM', .bsline, '_', .Ka, '_', .B,'_', .M, 
               '_pm', .pmut, '_mdist', .mean_distance, '_c', .c, sep = "" )
mainDir = "E:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/simulations"
setwd( paste( mainDir, '/', subDir, sep = "" ) )
runs <- list.files(pattern = ".RData")



# for PLOTTING
Dworld = 61 # diameter of the world in units (odd number) for PLOTTING

Sworld <- (Dworld + 1) / 2 # number of cells composing a side
Segmentworld <- c( seq(Sworld,Dworld), seq(Dworld -1,Sworld) ) # number of cells to have each value on the Red axis
range <- seq( -Sworld+1 , Sworld-1 ) # range of coordinates on both red and blue axis
Redworld <- rep( range, Segmentworld) # generate the Red coordinates of all the cells in the world

temp1 <- rep( list( -range ), Sworld )
Blueworld <- append( mapply(FUN = head, temp1, Segmentworld[1: Sworld ] ) , mapply(FUN = tail, temp1[-1], Segmentworld[ (Sworld+1):Dworld ] ))
Blueworld <- unlist( Blueworld )









########### TEMPORAL DYNAMICS OF INVASION
invasion_tot <- matrix( nrow = length(runs), ncol = .tps  )
invasion_asex <- matrix( nrow = length(runs), ncol = .tps  )

########### MARGINALITY
MEANasexdist <- numeric( length = length(runs) )
MEANsexdist <- numeric( length = length(runs) )
QUANTasexdist <- matrix( nrow = length(runs), ncol = 4 )
QUANTsexdist <- matrix( nrow = length(runs), ncol = 4 )


########### HETEROZYGOSITY
QHTZasex <- matrix( nrow = length(runs), ncol = 4) 
QHTZsex <- matrix( nrow = length(runs), ncol = 4) 

########### POPULATION DIVERSITY IN RESOURCE USE
QSHANNONasex <- matrix( nrow = length(runs), ncol = 4) 
QSHANNONsex <- matrix( nrow = length(runs), ncol = 4) 


for ( i in 1 : length( runs ) )
  {
  load( runs[i] )
  invasion_tot[i,] <- INVASION_tot
  invasion_asex[i,] <- INVASION_asex
  
  # -----------
  popy <- 3/2 * s * popB[ popsurvival==1 ]
  popx <- sqrt(3) * s * ( popB[ popsurvival==1 ]/2 + popR[ popsurvival==1 ])
  popdist <- sqrt( popx^2 + popy^2 ) 
  
  MEANasexdist[i] <- mean( popdist[ repro[ popsurvival == 1 ] == 'a' ] )
  MEANsexdist[i] <- mean( popdist[ repro[ popsurvival == 1 ] == 's' ] )
  
  QUANTasexdist[i,] <- quantile( popdist[ repro[ popsurvival == 1 ] == 'a' ], probs = c(.25, .5, .75, 1) )
  QUANTsexdist[i,] <- quantile( popdist[ repro[ popsurvival == 1 ] == 's' ], probs = c(.25, .5, .75, 1) )
  
  # -----------
  
  # -----------
  
  #### PHENOTYPIC DIVERSITY -> Shannon index based on RESOURCE CONSUMPTION
  # some kind of shannon index based on the resource consumed. Of the total resource consumed, 
  # how much was eaten of each category of the 2*G. Pop of total heteroz would consume 1/2G of each.
  
  occupiedpatches <- sort( unique( coordalive ) )
  patch_compet <- matrix( nrow = length(occupiedpatches) , ncol= 2*G )
  rownames( patch_compet ) = occupiedpatches
  for (g in 1:G){
    Tgenepatch <- table( factor( popgenome[,g], levels = 0:2), coordalive )
    #indivperpatch <- colSums( Tgenepatch )
    patch_compet[,(g*2-1)] <- Tgenepatch[1,]*2 + Tgenepatch[2,] # vector with number of 0 alleles on the patch
    patch_compet[,(g*2)] <- Tgenepatch[3,]*2 + Tgenepatch[2,] # vector with number of 1 alleles on the patch
  }
  freq_sp <- ( patch_compet / apply(patch_compet, 1, sum)) 
  log_freq_sp <- log( freq_sp )
  log_freq_sp[ is.infinite( log_freq_sp ) ] <- 0
  shannon <- - apply( freq_sp * log_freq_sp, 1, sum )
  std_shannon <- shannon / ( - log ((1/(2*G))) )
  
  ind_shannon_a <- std_shannon[ popRB[ repro == 'a' & popsurvival == 1 ] ] 
  ind_shannon_s <- std_shannon[ popRB[ repro == 's' & popsurvival == 1 ] ] 
  
  # ---------
  
  for( j in 1:4)
    {
    Ahtz <- popgenome[ popsurvival == 1,][ repro[ popsurvival == 1 ] == 'a' & popdist < QUANTasexdist[i,j], ] ==1 
    Shtz <- popgenome[ popsurvival == 1,][ repro[ popsurvival == 1 ] == 's' & popdist < QUANTasexdist[i,j], ] ==1 
    QHTZasex[i,j] <- sum( Ahtz ) / ( G*length( Ahtz ) )
    QHTZsex[i,j] <- sum( Shtz ) / ( G*length( Shtz ) )
    QSHANNONasex[i,j] <- mean( ind_shannon_a[ popdist[ repro[popsurvival==1] == 'a'] < QUANTasexdist[i,j] ] )
    QSHANNONsex[i,j] <- mean( ind_shannon_s[ popdist[ repro[popsurvival==1] == 's'] < QUANTasexdist[i,j] ] )
    }
  
 
  
  
  
  }




## TEMPORAL DYNAMICS OF INVASION

par(mfrow = c(3,4) )

rad_invasion_tot <-  sqrt( invasion_tot / ( 2 * pi ) )
rad_invasion_asex <- sqrt( invasion_asex / ( 2 * pi ) )
rad_invasion_sex <-  sqrt( ( invasion_tot - invasion_asex ) / ( 2 * pi ) )

MAX <- max(rad_invasion_tot)
MIN <- min(rad_invasion_asex)

plot( rad_invasion_tot[1,], type='l', main = '', col = 'red', ylim = c(MIN, MAX) )
points( rad_invasion_sex[1,], type = 'l', col = 'green')
points( rad_invasion_asex[1,], pch = 4, type = 'l', col = 'black')

for ( i in 2:length(runs)){
  points( rad_invasion_tot[i,], type='l', col = 'red', ylim = c(MIN, MAX)  )
  points( rad_invasion_sex[i,], type = 'l', col = 'green')
  points( rad_invasion_asex[i,], pch = 4, type = 'l', col = 'black')
}

# averaged over all the runs:
points(  apply( rad_invasion_tot, 2, mean), type='l', col = 'blue', ylim = c(MIN, MAX), main = 'radius AVERAGE', lwd=5 )
points( apply( rad_invasion_sex, 2, mean), type = 'l', col = 'blue', lwd=5)
points( apply( rad_invasion_asex, 2, mean), pch = 4, type = 'l', col = 'blue', lwd=5)



##################### 

## MARGINALITY

boxplot( cbind( QUANTasexdist, QUANTsexdist), 
         names = rep(c("Q1", "Q2", "Q3", "Q4"), 2),
         col = c(rep('grey', 4), rep('green', 4)), 
         main = "Quantiles of distance to center of asexuals and sexuals" )

## HETEROZYGOSITY BY DISTANCE AT t = tps

boxplot( cbind( QHTZasex, QHTZsex), 
         names = rep(c("Q1", "Q2", "Q3", "Q4"), 2),
         col = c(rep('grey', 4), rep('green', 4)), 
         main = "Mean htz for asexual and sexual distance quantiles" )

## PATCH DIVERSITY IN RESOURCE USE experienced by different class at different distances from the center, t = tps

boxplot( cbind( QHTZasex, QHTZsex), 
         names = rep(c("Q1", "Q2", "Q3", "Q4"), 2),
         col = c(rep('grey', 4), rep('green', 4)), 
         main = "Mean resource use index experienced by asexuals and sexuals, as a function of distance" )





