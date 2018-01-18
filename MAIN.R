rm(list=ls())
require(fields)
require(sp)
require(animation)
require(ggplot2)
require(RColorBrewer)

################################################################################
########################   Parameter definition    #############################
################################################################################


########################################
## Population parameters
.K = 10 # 'carrying capacity' of the patches
.fec = 5 # fecundity per female that reaches reproduction
.fecasex = 2
.G = 3 # number of genes modelled (each has two alleles, 0 and 1)
.probamating = 1 # can be set to a constant probability of mating, or "allee" (male density-dependent)
.compet = '_ded_' # '_fc_' or '_ded_'
.mean_distance = 1
.Xinit = 10
.Yinit = 10
.Xdim = 200
.Ydim = 10
.c = 1
.bsline = 0.8 # survival of total homozygotes
.Ka = 1-.bsline # how much more to total heterozygotes survive
.B = 20 # steepness of the sigmoid
.M = 0.3 # heterozygosity at the middle of the sigmoid
.pmut = 0.005
.tps = 500




## Population parameters
K = 10 # 'carrying capacity' of the patches
fec = 5 # fecundity per female that reaches reproduction
fecasex = 2
G = 3 # number of genes modelled (each has two alleles, 0 and 1)
probamating = 1 # can be set to a constant probability of mating, or "allee" (male density-dependent)
compet = '_ded_' # '_fc_' or '_ded_'
mean_distance = 1
Xinit = 5
Yinit = 10
Xdim = 200
Ydim = 10
c = 1
bsline = 0.6 # survival of total homozygotes
Ka = 0.4 # how much more to total heterozygotes survive
B = 20 # steepness of the sigmoid
M = 0.3 # heterozygosity at the middle of the sigmoid
pmut = 0.005
tps = 500
s = 2/3 # this way the distance between two hexagon centers on the R and B axes is 1



.Ninit <- .K * .fec


########################################
## type of competition (for directory name, remember to use the right function)


########################################
## DISPERSAL KERNEL (no d-dep)
# short_dist_sd <- 0.5 # sd of the normal law centered on 0 for short-distance dispersers
# plongdist <- 0.1 # probability of being a long-distance disperser
# long_dist_mean <- 6
# long_dist_sd <- 3


########################################
## Dispersal kernel, exponential power law as found in KLEIN06

# thin-tailed for c > 1, fat-tailed for c < 1, exponentiel for c = 1
# mean distance is ( a * gamma( 2/c ) ) / gamma( 1/c )



a <- .mean_distance / ( gamma( 2/.c ) / gamma( 1/.c ) )

nsamples <- 1E5
x <- runif(nsamples, min = 0.5, max = 1)
Gamma <- function(a, z) pgamma(z, a, lower = FALSE) * gamma(a) # define the same incomplete gamma function as Mathematica, which is gamma(a) from z to infinity
f2 <- function(x, u) 1 - ( x* (abs(x/a)^.c) ^(-1/.c) * Gamma( 1/.c , abs(x/a)^.c ) ) / ( 2*a*gamma(1/.c) ) - u
my.uniroot2 <- function(x) uniroot(f2, interval = c(0.000001, 500), u = x, tol = 0.0001)$root
.dispkernel <- vapply(x, my.uniroot2, numeric(1))



########################################
## Homozygosity depression


# curve( Ka /( 1 + 1*exp( -B*(x-M) )) +bsline , from = 0, to = 1)

########################################
## Probability to switch to asexuality


########################################
## number of generations


########################################
## GEOGRAPHY OF HEXATOPIA


# # for PLOTTING
# Dworld = 61 # diameter of the world in units (odd number) for PLOTTING
# 
# Sworld <- (Dworld + 1) / 2 # number of cells composing a side
# Segmentworld <- c( seq(Sworld,Dworld), seq(Dworld -1,Sworld) ) # number of cells to have each value on the Red axis
# range <- seq( -Sworld+1 , Sworld-1 ) # range of coordinates on both red and blue axis
# Redworld <- rep( range, Segmentworld) # generate the Red coordinates of all the cells in the world
# 
# temp1 <- rep( list( -range ), Sworld )
# Blueworld <- append( mapply(FUN = head, temp1, Segmentworld[1: Sworld ] ) , mapply(FUN = tail, temp1[-1], Segmentworld[ (Sworld+1):Dworld ] ))
# Blueworld <- unlist( Blueworld )




################################################################################
##########################    running simulations    ###########################
################################################################################




if ( exists('.compet') == F ) print('how is compet?')
subDir = paste("K", .K, '_',.compet, '_fs', .fec, '_fa', .fecasex, '_', .probamating,
               '_G', .G, '_bKaBM', .bsline, '_', .Ka, '_', .B,'_', .M, 
               '_pm', .pmut, '_mdist', .mean_distance, '_c', .c, '_Xdim', .Xdim, sep = "" )

mainDir = 'C:/Users/Anais Tilquin/Dropbox/Simulations'
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
