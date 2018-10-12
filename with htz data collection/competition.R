### where different genotypes can exploid different (limited) resources

# input: compet, popXY, haspartner,  popgenome, G, K
# output : popsurvival

# each patch can hold at max K individuals, as it has K units of each resource type (ie for each gene)

# competition <- function( compet, popXY, haspartner, popgenome, G , K){
competition <- function( compet, popXY, popgenome, G , K, niches){  
  
  if(compet=="_ded_"){
    
    occupiedpatches <- sort( unique( popXY ) )
    
    if( niches==3 ){
      # prepare for computing the number of copies of each allele present in each patch
      patch_compet <- matrix( nrow = length(occupiedpatches) , ncol= 3*G )
      rownames( patch_compet ) = occupiedpatches # /!/ names are sorted
      
      for (g in 1:G){
        Tgenepatch <- table( factor( popgenome[,g], levels = 0:2), popXY ) # table with allelic sums per patch for gene g
        patch_compet[,(g*3-2)] <- Tgenepatch[1,] # vector with number of 00 genotypes for gene g on the patch
        patch_compet[,(g*3-1)] <- Tgenepatch[2,] # vector with number of 01 genotypes for gene g on the patch
        patch_compet[,(g*3)] <- Tgenepatch[3,] # vector with number of 11 genotypes for gene g on the patch
      }
      
      # table with reward of possessing each allele (col: al0Gen1 al1Gen1 al0Gen2 al1Gen2...) on each patch (row)
      patch_compet <- (K/3) / patch_compet # total resource of each type available in the patch, divided by number of corresponding alleles in the patch
      patch_compet[ patch_compet > 1 ] <- 1 # on patches with excess food, individual share is bounded at 1; Inf values also become 1 but that doesn't matter as they stand for alleles no-one has
      
      reward <-  patch_compet[ popXY, ] # matrix of the reward got by each individual (row) from each of his alleles (g1a0, g1a1, g2a0 etc)
      
      popfitness <- vector( length = length(popXY) )
      for( g in 1:G ){
        fitness_g <- (popgenome[,g] == 0) * reward[, g*3-2] + 
          ( popgenome[,g] == 1 ) * reward[, 3*g-1] +
          ( popgenome[,g] == 2 ) * reward[, 3*g]
        popfitness <- popfitness + fitness_g
      }
      popfitness <- popfitness / G # the survival probability resulting from resource acquisition, of each individual with reproductive prospects
      
    } else {
      
      # prepare for computing the number of copies of each allele present in each patch
      patch_compet <- matrix( nrow = length(occupiedpatches) , ncol= 2*G )
      rownames( patch_compet ) = occupiedpatches # /!/ names are sorted
      
      for (g in 1:G){
        Tgenepatch <- table( factor( popgenome[,g], levels = 0:2), popXY ) # table with allelic sums per patch for gene g
        patch_compet[,(g*2-1)] <- Tgenepatch[1,]*2 + Tgenepatch[2,] # vector with number of copies of the 0 allele of gene g on the patch
        patch_compet[,(g*2)] <- Tgenepatch[3,]*2 + Tgenepatch[2,] # vector with number of copies of the 1 allele of gene g on the patch
      }
      
      # table with reward of possessing each allele (col: al0Gen1 al1Gen1 al0Gen2 al1Gen2...) on each patch (row)
      patch_compet <- K / patch_compet # total resource of each type available in the patch, divided by number of corresponding alleles in the patch
      patch_compet[ patch_compet > 1 ] <- 1 # on patches with excess food, individual share is bounded at 1; Inf values also become 1 but that doesn't matter as they stand for alleles no-one has
      
      # popgenome_haspartner <- popgenome[haspartner==T,]
      reward <-  patch_compet[ popXY, ] # matrix of the reward got by each individual (row) from each of his alleles (g1a0, g1a1, g2a0 etc)
      
      popfitness <- vector( length = length(popXY) ) 
      for( g in 1:G ){
        fitness_g <- popgenome[,g] * reward[, 2*g] + ( 2- popgenome[,g] ) * reward[, 2*g-1]
        popfitness <- popfitness + fitness_g
      }
      popfitness <- popfitness / (2*G) # the survival probability resulting from resource acquisition, of each individual with reproductive prospects
      
    }
    
    
    # each individual of the population now stochastically survives (1) or not (0). Without reproductive prospects, the individual dies.
    

    popsurvival <- numeric( length(popXY) )
    popsurvival <- mapply( FUN = rbinom, prob = popfitness, size = 1, n = 1)

    
    return(popsurvival)
  }
  
  if(compet=='_fc_'){
    

    pop_patch <- table( popXY )
    popfitness <- K / pop_patch[ popXY ]
    popfitness[ popfitness > 1 ] <- 1
    condition <- !is.na( popfitness )
    popsurvival <- rep( 0, length(popXY) )
    popsurvival[ condition ] <- mapply( FUN = rbinom, prob = popfitness[ condition ], size = 1, n = 1)
    
    
    
    return( popsurvival )
  }
  
  
}


library(compiler)
competition <- cmpfun(competition)



