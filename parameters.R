# .plot <- plot <- T
# 
# 
# ########################################
# ## Population parameters
# .K <- K <- 10 # 'carrying capacity' of the patches
# .fec <- fec <- 5 # fecundity per female that reaches reproduction
# .fecasex <- fecasex <- 2
# .G <- G <- 5 # number of genes modelled (each has two alleles, 0 and 1)
# .probamating <- probamating <- 1 # can be set to a constant probability of mating, or "allee" (male density-dependent)
# 
# ########################################
# ## Homozygosity depression
# .bsline <- bsline <- 0.8 # survival of total homozygotes
# # Ka <- .Ka <- 1-.bsline # how much more to total heterozygotes survive
# # .B <- B <- 20 # steepness of the sigmoid
# # .M <- M <- 0.3 # heterozygosity at the middle of the sigmoid
# 
# # curve( Ka /( 1 + 1*exp( -B*(x-M) )) +bsline , from = 0, to = 1)
# 
# ########################################
# ## Probability to switch to asexuality
# .pmut <- pmut <- 0.005
# 
# ########################################
# ## number of generations
# .tps <- tps <- 500
# 
# ########################################
# ## type of competition
# .compet <- compet <- '_ded_' # '_fc_' or '_ded_'
# 
# ########################################
# ## Dispersal kernel, exponential power law as found in KLEIN06
# 
# # thin-tailed for c > 1, fat-tailed for c < 1, exponentiel for c = 1
# # mean distance is ( a * gamma( 2/c ) ) / gamma( 1/c )
# .mean_distance <- mean_distance <- 1
# .c <- c <- 1
# 
# ########################################
# ## geography world
# .Xinit <- Xinit <- 10
# .Yinit <- Yinit <- 10
# .Xdim <- Xdim <- 200
# .Ydim <- Ydim <- 10
# 
# 
# .Ninit <- Ninit <- .K * .fec


####################### for batch
.plot <- F

.tps = 1000
.Xinit <- Xinit <- 10
.Yinit <- Yinit <- 10
.Xdim <- Xdim <- 200
.Ydim <- Ydim <- 10

.seed = 42
set.seed(.seed)
n <- 500
vK <- sample( 5:40, n, replace=T)
vfec <- sample( 3:10, n, replace=T)
vdisadvantage <- sample( seq( 0.2, 0.8, 0.1 ), n, replace = T)
vfecasex <- round(vdisadvantage * vfec)
vG <- sample( 2:20, n, replace=T)
vprobamating <- rep( 1, n )
vbsline <- sample( seq( 0.2,1, 0.1 ), n, replace=T)
vpmut <- sample( c( 0.005, 0.0005, 0.00005 ), n, replace=T)
vcompet <- sample(c( "_ded_", "_fc_" ), n, replace = T)
vmean_distance <- sample( seq(0.5,3,0.1), n, replace = T) 
vc <- sample( c(0.5, 0.8, 1, 1.2, 5), n, replace = T)

results <- data.frame(matrix(NA, nrow = n, ncol = 22))
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
                       'end.t',
                       'seed')

results[,2:11] <- cbind( vK, vfec, vfecasex, vprobamating, vG, vbsline,
                        vpmut, vmean_distance, vc, 1:n)
results[,1] <- vcompet
.results <- results[
  with(results, order(vmean_distance, vc)),
  ]

