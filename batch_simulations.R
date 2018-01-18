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






for ( .run in 1:10 )
{
  master( s = .s, 
          K = .K, 
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
          Ka = .Ka,
          B = .B,
          M = .M,
          pmut = .pmut,
          tps = .tps,
          run = .run,
          compet = .compet)
  }

saveVideo(expr=master( s = .s,
                       K = .K, 
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
                       Ka = .Ka,
                       B = .B,
                       M = .M,
                       pmut = .pmut,
                       tps = .tps,
                       run = run ), video.name = "animation.mp4", img.name = "Rplot",
          ffmpeg = ani.options("ffmpeg"), other.opts = if (grepl("[.]mp4$",
                                                                 video.name)) "-pix_fmt yuv420p")
