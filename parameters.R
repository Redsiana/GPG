
########################################
## Population parameters
.K <- K <- 10 # 'carrying capacity' of the patches
.fec <- fec <- 5 # fecundity per female that reaches reproduction
.fecasex <- fec <- 2
.G <- G <- 3 # number of genes modelled (each has two alleles, 0 and 1)
.probamating <- probamating <- 1 # can be set to a constant probability of mating, or "allee" (male density-dependent)


########################################
## Homozygosity depression
.bsline <- bsline <- 0.8 # survival of total homozygotes
Ka <- .Ka <- 1-.bsline # how much more to total heterozygotes survive
.B <- B <- 20 # steepness of the sigmoid
.M <- M <- 0.3 # heterozygosity at the middle of the sigmoid

# curve( Ka /( 1 + 1*exp( -B*(x-M) )) +bsline , from = 0, to = 1)

########################################
## Probability to switch to asexuality
.pmut <- pmut <- 0.005

########################################
## number of generations
.tps <- tps <- 500

########################################
## type of competition
.compet <- compet <- '_ded_' # '_fc_' or '_ded_'

########################################
## Dispersal kernel, exponential power law as found in KLEIN06

# thin-tailed for c > 1, fat-tailed for c < 1, exponentiel for c = 1
# mean distance is ( a * gamma( 2/c ) ) / gamma( 1/c )
.mean_distance <- mean_distance <- 1
.c <- c <- 1

########################################
## geography world
.Xinit <- Xinit <- 10
.Yinit <- Yinit <- 10
.Xdim <- Xdim <- 200
.Ydim <- Ydim <- 10


.Ninit <- Ninit <- .K * .fec