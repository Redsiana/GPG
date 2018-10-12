
  # -----------
  analysis_shannon <- function( newbabyX = newbabyX,
                    newbabyY = newbabyY,
                    popgenome = popgenome,
                    G = G, Xdim = Xdim ){
    #### PHENOTYPIC DIVERSITY -> Shannon index based on RESOURCE CONSUMPTION
    # some kind of shannon index based on the resource consumed. Of the total resource consumed, 
    # how much was eaten of each category of the 2*G. Pop of total heteroz would consume 1/2G of each.
    newbabyXY = paste(newbabyX, ",", newbabyY, sep="")
    occupiedpatches <- sort( unique( newbabyXY ) )
    X_position <- as.numeric( unlist( strsplit(occupiedpatches, ",") )[ c(TRUE,FALSE) ] )
    # patch_compet <- matrix( nrow = length(occupiedpatches) , ncol= 3*G )
    # rownames( patch_compet ) = occupiedpatches
    # for (g in 1:G){
    #   Tgenepatch <- table( factor( popgenome[,g], levels = 0:2), newbabyXY )
    #   #indivperpatch <- colSums( Tgenepatch )
    #   patch_compet[,(g*3-2)] <- Tgenepatch[1,] # vector with number of 00 genotypes on the patch
    #   patch_compet[,(g*3-1)] <- Tgenepatch[2,] # vector with number of 01 genotypes on the patch
    #   patch_compet[,(g*3)] <- Tgenepatch[3,]
    # }
    # freq_sp <- ( patch_compet / (apply(patch_compet, 1, sum)/G) ) 
    # log_freq_sp <- log2( freq_sp )
    # log_freq_sp[ is.infinite( log_freq_sp ) ] <- 0
    # 
    # shannon <- - apply( freq_sp * log_freq_sp, 1, sum )
    # std_shannon <- shannon / ( - log2 (3) )
    # 
    # X_position <- as.numeric( unlist( strsplit(occupiedpatches, ",") )[ c(TRUE,FALSE) ] )
    # Y_position <- as.numeric( unlist( strsplit(occupiedpatches, ",") )[ c(FALSE, TRUE) ] )
    # 
    # shannon_gradient <- tapply(std_shannon, X_position, mean) 
    # res_shannon <- rep(NA, Xdim)
    # res_shannon[ as.numeric( names(shannon_gradient)) + 1] <- shannon_gradient
    # 
    # 
    
    shannon_function = function(x) {
      ## BEWARE this is not a true Shannon Index. The normalizing factor is if
      ## all individuals on the patch were different, or if all possible genotypes are represented
      ( - sum( ( table(x) / length(x) ) * log2 ( table(x) / length(x) ) ) ) / min( log2( length(x) ), log2( 3^G ) )
    }
    
    genome <- apply(popgenome, 1, paste, collapse="") # collapses loci in 1 string
    shannon_per_patch <- tapply( genome, newbabyXY, shannon_function )
    shannon_per_X <- tapply( shannon_per_patch, X_position, mean )
    res_shannon <- rep(NA, Xdim)
    res_shannon[ as.numeric( names(shannon_per_X)) + 1] <- shannon_per_X

    return(res_shannon)
  }
  
  
  # ind_shannon_a <- std_shannon[ popRB[ repro == 'a' & popsurvival == 1 ] ] 
  # ind_shannon_s <- std_shannon[ popRB[ repro == 's' & popsurvival == 1 ] ] 
  
  # ---------

  analysis_htz <- function( popgenome = popgenome,
                            G = G,
                            newbabyX = newbabyX, Xdim = Xdim){
    # Ahtz <- rowSums( popgenome[ popcloneline == 0, ] ==1 )/G
    # Shtz <- rowSums( popgenome[ popcloneline != 0, ] ==1 )/G
    Thtz <- rowSums( popgenome ==1 )/G
    # 
    # Thtz_gradient <- tapply(Thtz, newbabyX, mean)
    # Ahtz_gradient <- tapply(Ahtz, newbabyX, mean)
    Thtz_gradient <- tapply(Thtz, newbabyX, mean)
    
    res_htz <- rep(NA, Xdim)
    res_htz[ as.numeric( names(Thtz_gradient)) + 1] <- Thtz_gradient
    return(res_htz)
  }

  
  library(compiler)
  analysis_shannon <- cmpfun(analysis_shannon)
  analysis_htz <- cmpfun(analysis_htz)
