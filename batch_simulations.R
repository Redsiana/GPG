rm( list = ls(all.names=TRUE) )

require(fields)
require(sp)
require(animation)
require(ggplot2)
require(RColorBrewer)

source("parameters.R")
source("master.R")
source("genesis.R")
source("diaspora.R")
source("competition.R")
source("meetic.R")
source("stillborns.R")
source("analysis_cloneposition.R")
# 
# ############ generate the distribution of dispersal distances
# a <- .mean_distance / ( gamma( 2/.c ) / gamma( 1/.c ) )
# 
# nsamples <- 1E5
# x <- runif(nsamples, min = 0.5, max = 1)
# Gamma <- function(a, z) pgamma(z, a, lower = FALSE) * gamma(a) # define the same incomplete gamma function as Mathematica, which is gamma(a) from z to infinity
# f2 <- function(x, u) 1 - ( x* (abs(x/a)^.c) ^(-1/.c) * Gamma( 1/.c , abs(x/a)^.c ) ) / ( 2*a*gamma(1/.c) ) - u
# my.uniroot2 <- function(x) uniroot(f2, interval = c(0.000001, 500), u = x, tol = 0.0001)$root
# .dispkernel <- dispkernel <- vapply(x, my.uniroot2, numeric(1))
# 
# kernel0.6 = .dispkernel
# kernel0.5 = .dispkernel
# kernek1.5 = .dispkernel
# # kernel0.8 = .dispkernel
# # kernel1 = .dispkernel
# # kernel1.2 = .dispkernel
# 
# par(mfrow = c(3,2))
# plot(table(round(kernel0.5)), xlim=c(0,18), main="0.5", ylim=c(0,50000))
# plot(table(round(kernel0.6)), xlim=c(0,18), main="0.6", ylim=c(0,50000))
# plot(table(round(kernel0.8)), xlim=c(0,18), main="0.8", ylim=c(0,50000))
# plot(table(round(kernel1)), xlim=c(0,18), main="1", ylim=c(0,50000))
# plot(table(round(kernel1.2)), xlim=c(0,18), main="1.2", ylim=c(0,50000))
# plot(table(round(kernel5)), xlim=c(0,18), main="5", ylim=c(0,50000))
# 
# hist(kernel0.8, xlim=c(0,15), breaks=15)
# hist(kernel1, xlim=c(0,15),  breaks=15)
# hist(truc, xlim=c(0,15),  breaks=15)

################################################################################
##########################    running simulations    ###########################
################################################################################




if ( exists('.compet') == F ) print('how is compet?')
subDir = paste("K", .K, '_',.compet, '_fs', .fec, '_fa', .fecasex, '_', .probamating,
               '_G', .G, '_bKaBM', .bsline, '_', .Ka, '_', .B,'_', .M, 
               '_pm', .pmut, '_mdist', .mean_distance, '_c', .c, '_Xdim', .Xdim, sep = "" )

mainDir = paste( getwd(), "/Simulations", sep="")
# mainDir = 'C:/Users/Anais Tilquin/Dropbox/Simulations'
# mainDir = "E:/Kokkonuts/Where sex/modelling/hexagonal grid/WITH ASEXUALS/simulations"
dir.create( file.path( mainDir, subDir ), showWarnings = FALSE)

setwd( file.path( mainDir, subDir ) )

.start = 1

for ( .run in .start:.n ) {
  
  compet <- .results$compet[.run]
  K  <- .results$K[.run]
  fec <- .results$fs[.run]
  fecasex <- .results$fa[.run]
  probamating <- .results$probamating[.run]
  G <- .results$G[.run]
  bsline <- .results$b[.run]
  pmut <- .results$pmut[.run]
  mean_distance <- .results$mean_dist[.run]
  c <- .results$c[.run]

  
  if(.run == .start){
    a <- mean_distance / ( gamma( 2/c ) / gamma( 1/c ) )
    
    nsamples <- 1E5
    x <- runif(nsamples, min = 0.5, max = 1)
    Gamma <- function(a, z) pgamma(z, a, lower = FALSE) * gamma(a) # define the same incomplete gamma function as Mathematica, which is gamma(a) from z to infinity
    f2 <- function(x, u) 1 - ( x* (abs(x/a)^c) ^(-1/c) * Gamma( 1/c , abs(x/a)^c ) ) / ( 2*a*gamma(1/c) ) - u
    my.uniroot2 <- function(x) tryCatch( 
      uniroot(f2, 
              interval = c(0.000001, 500),
              u = x, 
              tol = 0.0001, 
              extendInt = "yes" )$root,
      error=function(e) NA)
    # .dispkernel <- dispkernel <- vapply(x, my.uniroot2, numeric(1))
    .dispkernel <- dispkernel <- na.omit( sapply(x, my.uniroot2) )
    
  } else {
    rm(dispkernel)
    if(  mean_distance != .results$mean_dist[.run-1] | c != .results$c[.run-1]  ){
      a <- mean_distance / ( gamma( 2/c ) / gamma( 1/c ) )
      
      nsamples <- 1E5
      x <- runif(nsamples, min = 0.5, max = 1)
      Gamma <- function(a, z) pgamma(z, a, lower = FALSE) * gamma(a) # define the same incomplete gamma function as Mathematica, which is gamma(a) from z to infinity
      f2 <- function(x, u) 1 - ( x* (abs(x/a)^c) ^(-1/c) * Gamma( 1/c , abs(x/a)^c ) ) / ( 2*a*gamma(1/c) ) - u
      my.uniroot2 <- function(x) tryCatch( 
        uniroot(f2, 
                interval = c(0.000001, 500),
                u = x, 
                tol = 0.0001, 
                extendInt = "yes" )$root,
        error=function(e) NA)
      # .dispkernel <- dispkernel <- vapply(x, my.uniroot2, numeric(1))
      .dispkernel <- dispkernel <- na.omit( sapply(x, my.uniroot2) )
      
      }
  }
  
  plot = .plot
  seed = sample(1:1000, 1)
  set.seed(seed)
  res <- master(c = c,
          mean_distance = mean_distance,
          K = K, 
          fec = fec, 
          fecasex = fecasex,
          G = G, 
          probamating = probamating,
          dispkernel = .dispkernel,
          Xinit = .Xinit,
          Yinit = .Yinit,
          Xdim = .Xdim,
          Ydim = .Ydim,
          bsline = bsline, 
          pmut = pmut,
          tps = .tps,
          run = .run,
          compet = compet,
          plot = .plot)
  
  .results[.run, 12:21] <- res
  .results[.run, 22] <- seed
  
  rm(list = setdiff(ls(), lsf.str()))
  write.table(.results, "results_1-200.txt")
  
  }

