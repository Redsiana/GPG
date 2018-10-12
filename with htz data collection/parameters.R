.plot <- F

.control <- F
.autonomous <- T
.asexmixis <- "apo"
.hermaphrodite <- T
.resident_selfer <- F
.niches <- 3
.facultative_parthenogen <- F


.tps = 1000
.Xinit <- Xinit <- 10
.Yinit <- Yinit <- 10
.Xdim <- Xdim <- 200
.Ydim <- Ydim <- 10

# .seed = 42
.seed = 66
set.seed(.seed)
.n <- 200
vbsline <- c( rep( 1, (2*.n) ), rep( 0.5, (2*.n) ) )
vcompet <- c( rep('_ded_', .n), rep('_fc_', .n), rep('_ded_', .n), rep('_fc_', .n))

vK <- sample( 5:40, .n, replace=T)
vK <- rep(vK, 4)
vfec <- sample( c(4,6,8,10,12), .n, replace=T)
vfec <- rep(vfec, 4)
if( .hermaphrodite ){ 
  vfecasex <- vfec
  } else{ 
    vdisadvantage <- rep(.5, .n)
    vdisadvantage <- rep(vdisadvantage, 4)
    vfecasex <- round(vdisadvantage * vfec)
    vfecasex[vfecasex==1] <- 2
  }

vG <- sample( 2:20, .n, replace=T)
vG <- rep(vG, 4)
vprobamating <- rep( 1, (4*.n) )

# bsline = 0.2
# curve(expr=1-(1-bsline)*((x-1))^2, from=0, to=1)
vpmut <- sample( c( 0.005, 0.0005, 0.00005 ), .n, replace=T)
vpmut <- rep(vpmut, 4)
# vcompet <- sample(c( "_ded_", "_fc_" ), .n, replace = T)
vmean_distance <- sample( seq(0.5,3,0.1), .n, replace = T) 
vmean_distance <- rep( vmean_distance, 4)
vc <- sample( c(0.5, 0.8, 1, 1.2, 5), .n, replace = T)
vc <- rep(vc, 4)

results <- data.frame(matrix(NA, nrow = (4*.n), ncol = 23))
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
                       'namerun',
                       'seed'
                       )

results[,2:11] <- cbind( vK, vfec, vfecasex, vprobamating, vG, vbsline,
                        vpmut, vmean_distance, vc, 1:(.n*4))
results[,1] <- vcompet
.results <- results[
  with(results, order(vmean_distance, vc)),
  ]

write.table(results, "parameters.txt", col.names = T)
