rm(list=ls())
require(fields)
require(sp)
require(animation)
require(ggplot2)
require(RColorBrewer)

source("parameters.R")


############ generate the distribution of dispersal distances
a <- .mean_distance / ( gamma( 2/.c ) / gamma( 1/.c ) )

nsamples <- 1E5
x <- runif(nsamples, min = 0.5, max = 1)
Gamma <- function(a, z) pgamma(z, a, lower = FALSE) * gamma(a) # define the same incomplete gamma function as Mathematica, which is gamma(a) from z to infinity
f2 <- function(x, u) 1 - ( x* (abs(x/a)^.c) ^(-1/.c) * Gamma( 1/.c , abs(x/a)^.c ) ) / ( 2*a*gamma(1/.c) ) - u
my.uniroot2 <- function(x) uniroot(f2, interval = c(0.000001, 500), u = x, tol = 0.0001)$root
.dispkernel <- dispkernel <- vapply(x, my.uniroot2, numeric(1))



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

results <- data.frame(matrix(NA, nrow = 10, ncol = 21))
colnames(results) <- c('compet',
                       'K',
                       'fs',
                       'fa',
                       'probamating',
                       'G',
                       'b',
                       'pmut',
                       'mean_dist',
                       'c',
                       'run',
                       'max.X.sex', 
                       'max.X.sex.t',
                       'mix',
                       'max.nb.sex',
                       'max.nb.sex.t', 
                       'nb.all.end',
                       'nb.sex.end',
                       'X.all.end',
                       'X.sex.end',
                       'end.t')


.run=1

for ( .run in 1:10 )
{
  
  results[i,1:11] <- c(.compet, .K, .fec, .fecasex, .probamating, .G, .bsline, .pmut, .mean_distance, .c, .run) 

  
  res <- master(K = .K, 
          fec = .fec, 
          fecasex = .fecasex,
          G = .G, 
          probamating = .probamating,
          dispkernel = .dispkernel,
          Xinit = .Xinit,
          Yinit = .Yinit,
          Xdim = .Xdim,
          Ydim = .Ydim,
          bsline = .bsline, 
          pmut = .pmut,
          tps = .tps,
          run = .run,
          compet = .compet,
          plot = .plot)
  
  results[i, 12:21] <- res
  
  }

