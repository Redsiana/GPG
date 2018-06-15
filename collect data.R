setwd("C:/Users/Anais Tilquin/Desktop/Simulations February/seed66")
SEED = 66

files <- list.files(pattern = "\\.RData$")
length(files)

pureness = vector(length = length(files))

results <- data.frame(matrix(NA, nrow = length(files), ncol = 26))
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
                       'mean.X.sex',
                       'mean.X.asex',
                       'sd.X.sex',
                       'sd.X.asex',
                       'seed')


results$seed <- SEED

reaches_end = NA
inflection = NA
inflection2 = NA

results[is.na(results$nb.all.end),]

for(id in 1:length(files)){
  
  load(files[id])
  
  # if(max(EDGE_t) < (Xdim-1)){
  #   reaches_end[id] <- NA
  # } else {
  #   reaches_end[id] <- min( which( EDGE_t == (Xdim-1)))
  # }
  
  inflection[id] <- max( which( EDGE_sex_t == (max(EDGE_sex_t)) ))
  inflection2[id] <- max( which( EDGE_sex_t >= ((max(EDGE_sex_t)-2)) ))


  
  
  # results[id,1:11] <- cbind( compet, K, fec, fecasex, probamating, G, bsline,
  #                          pmut, mean_distance, c, id)
  # results[id,12:21] <- res
  # nums = na.omit(as.numeric( str_extract_all(files[id],"\\(?[0-9,.]+\\)?")[[1]] ))
  # results$mean_dist[id] <- nums[length(nums)-2]
  # results$c[id] <- nums[length(nums)-1]
  # results$run[id] <- nums[length(nums)]
  # results[results$seed==SEED,22][id] <- mean(newbabyX[repro == "s"])
  # results[results$seed==SEED,23][id] <- mean(newbabyX[repro == "a"])
  # results[results$seed==SEED,24][id] <- sd(newbabyX[repro == "s"])
  # results[results$seed==SEED,25][id] <- sd(newbabyX[repro == "a"])
  # results$mix[id] = pureness_sex_t[ mean(c( max( which(EDGE_sex_t==max(EDGE_sex_t))), t)) ]
  print(id)
}


plot(inflection2~reaches_end)
abline(0,2)


results = na.omit(results)
results[results[,1]==1,1] <- "_ded_"
results[results[,1]==2, 1] <- "_fc_"

results42 = results
write.table(results, "results_seed42.txt", col.names=T)

results = rbind(results24, results42)
write.table(results, "results_24_42_extended.txt", col.names=T)

# check if any simulation was done with the same parameters
sum(duplicated(results[,1:11]))



# results = results[with(results, order(compet, K, fs, fa,G, b, pmut, mean_dist, c )), ]
parameters = data.frame(vcompet, vK, vfec, vfecasex, vG, vbsline, vpmut, vmean_distance, vc)
# parameters = parameters[with(parameters, order(vcompet, vK, vfec, vfecasex,vG, vbsline, vpmut, vmean_distance, vc )), ]


# putting back the correct dispersal parameters in dataframe results


for(i in 1:500){files_theoric[i] = namerun <- paste( parameters$vcompet[i], "K", parameters$vK[i], '_fs', parameters$vfec[i], '_fa',
                                                     parameters$vfecasex[i], '_', parameters$vprobamating[i],'_G', parameters$vG[i], '_b',
                                                     parameters$vbsline[i], '_pm', parameters$vpmut[i], '_d', parameters$vmean_distance[i], '_c', parameters$vc[i], '_','.RData',
                                                     sep = "" )}


parameters = cbind(parameters, files_theoric)
parameters = parameters[order(files_theoric),]
results[results$seed == 42,]$mean_dist = parameters$vmean_distance[-c(1,2,3)] # for some reason the first simulations in the folder don't give results
results[results$seed == 42,]$c = parameters$vc[-c(1,2,3)]

results$compet = droplevels(results$compet)
write.table(results, "results_24_42.txt", col.names = T)




## order the table and the vector with the names of the datafiles (T2) the same
results_alphab <- (results[order(results$compet,as.character(results$K), as.character(results$fs), as.character(results$G), as.character(results$b) ),])
filesT2 <- list.files(pattern = "\\T2.RData$")


holes_df = results_alphab[is.na(results_alphab$nb.all.end),]
files_holes = filesT2[is.na(results_alphab$nb.all.end) ]

for(id in 1:length(files_holes)){
  
  load(files_holes[id])

  holes_df[id, 17:20] <-  c( INVASION_tot[t-1], #                                 how many individuals when simulation ended
  INVASION_sex[t-1], #                                 how many sexuals when simulation ended
  EDGE_t[t-1], #                                       how far the entire pop was when simulation ended
  EDGE_sex_t[t-1])
  
}
results_alphab[is.na(results_alphab$nb.all.end),] <- holes_df

results <- results_alphab
write.table(results, "results_filled_gaps.txt", col.names = T)











#### EXTRACT COMPOSITION OF STARTING AND ENDING PATCHES

results = read.table("results.txt", h=T)

# to erase from envir
ghost = c("B","babyclonalorigin" ,"babycloneline"    ,
"babygenome"     ,   "babyrepro"   ,      "babysex"      ,    
 "babyX"       ,      "babyY"    ,         "bsline"      ,     
 "bug"        ,       "c"        ,         "compet"     ,      
 "coordalive"     ,   "count"          ,   "count_col"    ,    
"dispkernel" ,       "EDGE_sex_t"      ,  "EDGE_t"        ,   
 "fec"         ,      "fecasex"          , "G"              ,  
 "haspartner"   ,     "INVASION_sex"      ,"INVASION_tot"    , 
 "K"             ,    "Ka"               , "M"                ,
 "mean_distance"  ,   "namerun1"   ,       "namerun2"         ,
"nb_sexuals"      ,  "newbabyX"    ,      "newbabyY"         ,
 "Ninit"            , "notfull_boolean"  , "old_clones_colors",
 "old_clones_names" , "plot"       ,       "pmut"             ,
 "popclonalorigin"  , "popcloneline" ,     "popgenome"        ,
 "popsurvival",       "popX"          ,    "popXY"            ,
 "popY"        ,      "probamating"    ,   "pureness_sex_t"   ,
 "repro"        ,       "run"              ,
 "sex"           ,    "t"                , "t_full_corridor"  ,
                "T2"           ,     "TEMP"             ,
 "tps"             ,  "Xdim"          ,    "Xinit"            ,
 "xmax_t"           , "Ydim"           ,   "Yinit"     )


files <- list.files(pattern = "\\.RData$")
results_alphab <- (results[order(results$compet,as.character(results$K), as.character(results$fs), as.character(results$G), as.character(results$b) ),])
filesT1 <- list.files(pattern = "\\T1.RData$")
filesT2 <- list.files(pattern = "\\T2.RData$")
filesT1theoretical <- gsub('T2', 'T1', filesT2)

compo_farT1 <- rep(NA, length(filesT2))
compo_farT2 <- NA
compo_closeT1 <- rep(NA, length(filesT2))
compo_closeT2 <- NA
extinct_before_T1 <- NA
extinct_before_T2 <- NA

for(id in 1:length(filesT2)){
  
  if( filesT1theoretical[id] %in% files ){
    load(filesT1theoretical[id])
    
    # percentage of sexuals in the last 10x10
    compo_farT1[id] <- sum( newbabyX[popcloneline== 0] > (Xdim- Xinit)) / 
      (sum( newbabyX[popcloneline==0] > (Xdim-Xinit)) + sum( newbabyX[popcloneline!= 0] > (Xdim-Xinit)))
    
    compo_closeT1[id] <- sum( newbabyX[popcloneline==0] < Xinit ) / 
      (sum( newbabyX[popcloneline==0] < Xinit) + sum( newbabyX[popcloneline!=0] < Xinit))
    
    extinct_before_T1[id] <- FALSE
  } else {
    extinct_before_T1[id] <- TRUE
  }
  
  rm(list=ghost)
  
  load(filesT2[id])
  
  if(exists("T1")) {
    if( t == 2*T1 ){
      extinct_before_T2[id] <- FALSE
    } else { extinct_before_T2[id] <- TRUE }
  } else {
    extinct_before_T2[id] <- TRUE
  }
  
  # percentage of sexuals in the last 10x10
  compo_farT2[id] <- sum( newbabyX[popcloneline==0] > (Xdim- Xinit)) / 
    (sum( newbabyX[popcloneline==0] > (Xdim-Xinit)) + sum( newbabyX[popcloneline!=0] > (Xdim-Xinit)))
  
  compo_closeT2[id] <- sum( newbabyX[popcloneline==0] < Xinit ) / 
    (sum( newbabyX[popcloneline==0] < Xinit) + sum( newbabyX[popcloneline!=0] < Xinit))

  rm(list=c(ghost, 'T1'))
  
  print(id)
  }

write.table(data.frame(results_alphab, compo_closeT1 = compo_closeT1,
                       compo_farT1 = compo_farT1,
                       compo_closeT2 = compo_closeT2,
                       compo_farT2 = compo_farT2,
                       extinct_before_T1 = extinct_before_T1,
                       extinct_before_T2 = extinct_before_T2,
                       namerun = filesT2), "df_compo_both_ends.txt")
DF = read.table("df_compo_both_ends.txt", h=T)


## collect cumulative distribution of sex and asex (as a percentage) for each regime
T1_cumul_ded1 = T1_cumul_ded0.5 = T1_cumul_fc1 = T1_cumul_fc0.5 = numeric(200)
T2_cumul_ded1 = T2_cumul_ded0.5 = T2_cumul_fc1 = T2_cumul_fc0.5 = numeric(200)
T1_matrix_ded1 = T1_matrix_ded0.5 = T1_matrix_fc1 = T1_matrix_fc0.5 = matrix(NA, 800,200) 
T2_matrix_ded1 = T2_matrix_ded0.5 = T2_matrix_fc1 = T2_matrix_fc0.5 = matrix(NA, 800,200)

for(id in 1:length(filesT2)){
  
  if( filesT1theoretical[id] %in% files ){
    load(filesT1theoretical[id])

    tab = table(newbabyX, popcloneline==0)
    distri = tab[,1]/(rowSums(tab))

    # if(compet=="_ded_" & bsline==1){ T1_cumul_ded1 <- T1_cumul_ded1 + distri}
    # if(compet=="_ded_" & bsline==0.5){ T1_cumul_ded0.5 <- T1_cumul_ded0.5 + distri}
    # if(compet=="_fc_" & bsline==1){ T1_cumul_fc1 <- T1_cumul_fc1 + distri}
    # if(compet=="_fc_" & bsline==0.5){ T1_cumul_fc0.5 <- T1_cumul_fc0.5 + distri}
    
    if(compet=="_ded_" & bsline==1){ T1_matrix_ded1[id,1:length(distri)] <- distri }
    if(compet=="_ded_" & bsline==0.5){ T1_matrix_ded0.5[id,1:length(distri)] <- distri }
    if(compet=="_fc_" & bsline==1){ T1_matrix_fc1[id,1:length(distri)] <- distri }
    if(compet=="_fc_" & bsline==0.5){ T1_matrix_fc0.5[id,1:length(distri)] <- distri }

    # ONLY TAKE THE DISTRIBUTION DATA IF THE SIMULATION REACHED THE END OF THE CORRIDOR

  }

  rm(list=c(ghost, "tab", "distri"))
  
  load(filesT2[id])
  
  if(exists("T1")) {  # means it takes data for runs when not extinct before T1, but perhaps before
    tab = table(newbabyX, popcloneline==0)
    distri = tab[,1]/(rowSums(tab))
    
    # if(compet=="_ded_" & bsline==1){ T2_cumul_ded1 <- T2_cumul_ded1 + distri}
    # if(compet=="_ded_" & bsline==0.5){ T2_cumul_ded0.5 <- T2_cumul_ded0.5 + distri}
    # if(compet=="_fc_" & bsline==1){ T2_cumul_fc1 <- T2_cumul_fc1 + distri}
    # if(compet=="_fc_" & bsline==0.5){ T2_cumul_fc0.5 <- T2_cumul_fc0.5 + distri}
    
    if(compet=="_ded_" & bsline==1){ T2_matrix_ded1[id,1:length(distri)] <- distri }
    if(compet=="_ded_" & bsline==0.5){ T2_matrix_ded0.5[id,1:length(distri)] <- distri }
    if(compet=="_fc_" & bsline==1){ T2_matrix_fc1[id,1:length(distri)] <- distri }
    if(compet=="_fc_" & bsline==0.5){ T2_matrix_fc0.5[id,1:length(distri)] <- distri }
    
  }
  
  rm(list=c(ghost, 'T1', 'tab', 'distri'))
  
  print(id)
}

T1_matrix_ded1 = na.omit(T1_matrix_ded1)
T1_matrix_ded0.5 = na.omit(T1_matrix_ded0.5) 
T1_matrix_fc1 = na.omit(T1_matrix_fc1)
T1_matrix_fc0.5 = na.omit(T1_matrix_fc0.5) 
T2_matrix_ded1 = na.omit(T2_matrix_ded1)
T2_matrix_ded0.5 = na.omit(T2_matrix_ded0.5)
T2_matrix_fc1 = na.omit(T2_matrix_fc1)
T2_matrix_fc0.5 = na.omit(T2_matrix_fc0.5)


par(mfrow = c(4,2))
matplot(t(T1_matrix_ded1), type='l')
points(colSums(T1_matrix_ded1)/dim(T1_matrix_ded1)[1], type="l", lwd=10)

matplot(t(T2_matrix_ded1), type='l')
points(colSums(T2_matrix_ded1)/dim(T2_matrix_ded1)[1], type="l", lwd=10)

matplot(t(T1_matrix_fc1), type='l')
points(colSums(T1_matrix_fc1)/dim(T1_matrix_fc1)[1], type="l", lwd=10)

matplot(t(T2_matrix_fc1), type='l')
points(colSums(T2_matrix_fc1)/dim(T2_matrix_fc1)[1], type="l", lwd=10)

matplot(t(T1_matrix_ded0.5), type='l')
points(colSums(T1_matrix_ded0.5)/dim(T1_matrix_ded0.5)[1], type="l", lwd=10)

matplot(t(T2_matrix_ded0.5), type='l')
points(colSums(T2_matrix_ded0.5)/dim(T2_matrix_ded0.5)[1], type="l", lwd=10)

matplot(t(T1_matrix_fc0.5), type='l')
points(colSums(T1_matrix_fc1)/dim(T1_matrix_fc0.5)[1], type="l", lwd=10)

matplot(t(T2_matrix_fc0.5), type='l')
points(colSums(T2_matrix_fc0.5)/dim(T2_matrix_fc0.5)[1], type="l", lwd=10)


maxlim = round( max(c(T1_cumul_ded0.5,T1_cumul_ded1, T1_cumul_fc0.5, T1_cumul_fc0.5, T1_cumul_fc1,
      T2_cumul_ded0.5, T2_cumul_ded1, T2_cumul_fc0.5, T2_cumul_fc1))/id, 2 ) +0.01
par(mfrow = c(4,2))
barplot(as.matrix(T1_cumul_ded1/id), ylim = c(0,maxlim), beside=T, main='T1_ded1')
abline(v=weighted.mean(0:199, T1_cumul_ded1/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T2_cumul_ded1/id), ylim = c(0,maxlim), beside=T, main='T2_ded1')
abline(v=weighted.mean(0:199, T2_cumul_ded1/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T1_cumul_fc1/id), ylim = c(0,maxlim), beside=T, main='T1_fc1')
abline(v=weighted.mean(0:199, T1_cumul_fc1/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T2_cumul_fc1/id), ylim = c(0,maxlim), beside=T, main='T2_fc1')
abline(v=weighted.mean(0:199, T2_cumul_fc1/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T1_cumul_ded0.5/id), ylim = c(0,maxlim), beside=T, main='T1_ded0.5')
abline(v=weighted.mean(0:199, T1_cumul_ded0.5/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T2_cumul_ded0.5/id), ylim = c(0,maxlim), beside=T, main='T2_ded0.5')
abline(v=weighted.mean(0:199, T2_cumul_ded0.5/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T1_cumul_fc0.5/id), ylim = c(0,maxlim), beside=T, main='T1_fc0.5')
abline(v=weighted.mean(0:199, T1_cumul_fc0.5/id), col="red")
abline(v=99.5, lty=2)

barplot(as.matrix(T2_cumul_fc0.5/id), ylim = c(0,maxlim), beside=T, main='T2_fc0.5')
abline(v=weighted.mean(0:199, T2_cumul_fc0.5/id), col="red")
abline(v=99.5, lty=2)

# 
# ## with bugged simulations replaced by NA:
# missing.lines = which(is.na(results_alphab$X.sex.end))
# df1 = data.frame(compo_closeT1 = compo_closeT1,
#            compo_farT1 = compo_farT1,
#            compo_closeT2 = compo_closeT2,
#            compo_farT2 = compo_farT2,
#            extinct_before_T1 = extinct_before_T1,
#            extinct_before_T2 = extinct_before_T2,
#            namerun = filesT2)
# df2 <- rbind( df1[1:( missing.lines[1]-1 ),], 
#               rep(NA,7),
#               df1[( missing.lines[1] ):( missing.lines[2] +1),],
#               rep(NA, 7),
#               df1[( missing.lines[2]+2 ):dim(df1)[1], ])

# select subset where both times are reached
DF2 = DF[(DF$extinct_before_T1 == FALSE &DF$extinct_before_T2 == FALSE),]


dT1 <- data.frame(y = c(compo_closeT1, compo_farT1),
                group = as.factor(rep(c('close', 'far'), each = length(filesT2))),
                id = rep(1:length(filesT2), 2))
dT2 <- data.frame(y = c(compo_closeT2, compo_farT2),
                  group = as.factor(rep(c('close', 'far'), each = length(filesT2))),
                  id = rep(1:length(filesT2), 2))
dclose <- data.frame(y = c(compo_closeT1, compo_closeT2),
                     group = as.factor(rep(c('close', 'far'), each = length(filesT2))),
                     id = rep(1:length(filesT2), 2))
dfar <- data.frame(y = c(compo_farT1, compo_farT2),
                   group = as.factor(rep(c('close', 'far'), each = length(filesT2))),
                   id = rep(1:length(filesT2), 2))

dfarT1 <- data.frame(y = c(DF$compo_closeT1[DF$compet == "_ded_"],
                           DF$compo_closeT1[DF$compet == "_fc_"]),
                  group = as.factor(rep(c('ded', 'fc'), each = 400)),
                  id = rep(1:400, 2))

dinbfarT1 <- data.frame(y = c(DF$compo_farT2[DF$b == 1],
                           DF$compo_farT2[DF$b == 0.5]),
                     group = as.factor(rep(c('no inb', 'inb'), each = 400)),
                     id = rep(1:400, 2))

d = dinbfarT1
ggplot(d, aes(y = y)) +
  geom_boxplot(aes(x = rep(c(-3, 3), each = (dim(d)[1])/2), group = group), fill = 'steelblue') +
  geom_point(aes(x = rep(c(-1, 1), each = (dim(d)[1])/2)), size = 5) +
  geom_line(aes(x = rep(c(-1, 1), each = (dim(d)[1])/2), group = id))
