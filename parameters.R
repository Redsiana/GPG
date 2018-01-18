
########################################
## Population parameters
.K <- K <- 10 # 'carrying capacity' of the patches
.fec <- fec <- 5 # fecundity per female that reaches reproduction
.fecasex <- fec <- 2
.G <- G <- 3 # number of genes modelled (each has two alleles, 0 and 1)
.probamating <- probamating <- 1 # can be set to a constant probability of mating, or "allee" (male density-dependent)
.compet <- compet <- '_ded_' # '_fc_' or '_ded_'
.mean_distance <- mean_distance <- 1
.Xinit <- Xinit <- 10
.Yinit <- Yinit <- 10
.Xdim <- Xdim <- 200
.Ydim <- Ydim <- 10
.c <- c <- 1
.bsline <- bsline <- 0.8 # survival of total homozygotes
Ka <- .Ka <- 1-.bsline # how much more to total heterozygotes survive
.B <- B <- 20 # steepness of the sigmoid
.M <- M <- 0.3 # heterozygosity at the middle of the sigmoid
.pmut <- pmut <- 0.005
.tps <- tps <- 500

.Ninit <- Ninit <- .K * .fec