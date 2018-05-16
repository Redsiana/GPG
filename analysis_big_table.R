source('display plots function.R') # defines function showme

results = na.omit( read.table('results_24_42_extended.txt', h=T))


plot(results)

fec_ratio = results$fa / results$fs
results = cbind(results, fec_ratio)
results$max.nb.sex = results$max.nb.sex / results$K
results$nb.all.end = results$nb.all.end / results$K
results$nb.sex.end = results$nb.sex.end / results$K

X.frac.sexuals = results$X.sex.end / results$X.all.end
results = cbind(results, X.frac.sexuals, results$mean.X.sex/results$mean.X.asex )
names(results)[30] = "mean.X.frac"
results_full = results

par(mfrow=c(2,2))
# mean X of asexuals and sexuals
matplot(t(cbind(results$mean.X.asex/results$X.all.end, results$mean.X.sex/results$X.all.end)), 
        type="b", pch=19, col=1, lty=1,
        xlab = "1=Asex, 2 = Sex",
        ylab = "mean X coordinate (end)", 
        main = 'mean X (end)') 
matplot(t(cbind(results$sd.X.asex, results$sd.X.sex)), 
        type="b", pch=19, col=1, lty=1, 
        xlab = "1=Asex, 2 = Sex",
        ylab = "sd X coordinate (end)",
        main = "sd X (end)") 
hist(results$mean.X.asex/results$X.all.end, main = "asexuals mean position", xlim=c(0,1))
hist(results$mean.X.sex/results$X.all.end, main = "sexuals mean position", xlim=c(0,1))


## if I want to remove the simulations that ended in 1000 generations
# and where asexuals don't have an a leat 2 fold cost
results = results_full[results_full$end.t<1000,]
results = results_full[results_full$end.t==1000,]
results = results[results$fec_ratio<=0.5,]

## comparer ded and fc
results_ded = results_full[results_full$compet == '_ded_',]
results_fc = results_full[results_full$compet == '_fc_',]
par(mfrow=c(3,2))
# mean X of asexuals and sexuals
matplot(t(cbind(results_ded$mean.X.asex/results_ded$X.all.end, results_ded$mean.X.sex/results_ded$X.all.end)), type="b", pch=19, col=1, lty=1, main = "ded mean") 
matplot(t(cbind(results_fc$mean.X.asex/results_fc$X.all.end, results_fc$mean.X.sex/results_fc$X.all.end)), type="b", pch=19, col=1, lty=1, main = "fc mean")
hist(results_ded$mean.X.asex/results_ded$X.all.end, main = "ded asex")
hist(results_fc$mean.X.asex/results_fc$X.all.end, main = "fc asex")
hist(results_ded$mean.X.sex/results_ded$X.all.end, main = "ded sex")
hist(results_fc$mean.X.sex/results_fc$X.all.end, main = "fc sex")
par(mfrow = c(1,1))
boxplot(results_ded$mean.X.asex/results_ded$X.all.end, results_fc$mean.X.asex/results_fc$X.all.end )
t.test(results_ded$mean.X.asex/results_ded$X.all.end, results_fc$mean.X.asex/results_fc$X.all.end )


results = results[results$fec_ratio<0.5,]
results2 = results[results$end.t==1000,]
results3 = results[results$end.t<1000,]
matplot(t(cbind(results2$mean.X.asex/results2$X.all.end, results2$mean.X.sex/results2$X.all.end)), type="b", pch=19, col=1, lty=1) 
matplot(t(cbind(results2$sd.X.asex, results2$sd.X.sex)), type="b", pch=19, col=1, lty=1) 
hist(results2$mean.X.sex/results2$X.all.end)
hist(results2$mean.X.asex/results2$X.all.end)

matplot(t(cbind(results3$mean.X.asex/results3$X.all.end, results3$mean.X.sex/results3$X.all.end)), type="b", pch=19, col=1, lty=1) 
matplot(t(cbind(results3$sd.X.asex, results3$sd.X.sex)), type="b", pch=19, col=1, lty=1) 
hist(results3$mean.X.sex/results3$X.all.end)
hist(results3$mean.X.asex/results3$X.all.end)


## fraction of asexuals in those runs that lasted 1000 generations
frac.asex.1000 <- (results2$nb.all.end - results2$nb.sex.end)/results2$nb.all.end

# when no asex establishes, they're distributed randomly centered around the middle
# when a sizeable amount establishes, they tend to be on the margin, 
# but the more the more centered -> they've actually taken over
plot(results2$mean.X.asex ~ frac.asex.1000)

# in what simulations do handicapped asexuals reach such high fractions?? -> a bit of everything...
results2[frac.asex.1000>0.8,]
showme( results2[frac.asex.1000>0.2,] )


plot(results2$mean.X.asex ~ results2$fec_ratio)
plot(frac.asex.1000 ~ results2$fec_ratio)
plot( results2$mean.X.asex ~ results2$fec_ratio )
abline( lm( frac.asex.1000 ~ results2$fec_ratio ) )
plot(lm( results2$mean.X.asex ~ results2$fec_ratio )$residuals ~ frac.asex.1000)

results3 = results2[frac.asex.1000>0.1,]


results.ded = results[results$compet == "_ded_",]
results.fc = results[results$compet == "_fc_",]



for(jj in c(22:26, 30)){
  par(mfrow=c(5,5))
  for(ii in c(1:7, 9, 10, 12:21)){
    plot( results[,jj] ~ results[,ii], xlab=names(results[ii]), ylab= names(results[jj]))
    # if(ii > 1) abline(lm( results[,jj] ~ results[,ii] ))
  }
  boxplot(results[,jj]~results$c, xlab="c")
  boxplot(results[,jj]~results$pmut, xlab="pmut")
}

for(jj in 12:21){
  par(mfrow=c(3,3))
  for(ii in c(1:3,23, 6,7,9,10)){
    plot( results.ded[,jj] ~ results.ded[,ii], xlab=names(results.ded[ii]), ylab= names(results.ded[jj]))
    if(ii > 1)  abline(lm( results.ded[,jj] ~ results.ded[,ii] ))
  }
  boxplot(results.ded[,jj]~results.ded$pmut, xlab="pmut")
}

for(jj in 12:21){
  par(mfrow=c(3,3))
  for(ii in c(1:3,23, 6,7,9,10)){
    plot( results.fc[,jj] ~ results.fc[,ii], xlab=names(results.fc[ii]), ylab= names(results.fc[jj]))
    if(ii > 1) abline(lm( results.fc[,jj] ~ results.fc[,ii] ))
  }
  boxplot(results.fc[,jj]~results.fc$pmut, xlab="pmut")
}


# for pmut = 5e-04, sexuals never reached the 200 mark - but does not seem to be because of
# randomly not having cases with low relative fec
sum( results$fec_ratio[results$pmut==0.0005] < 0.5 ) / length( results$fec_ratio[results$pmut==0.0005] < 0.5  )
sum( results$fec_ratio[results$pmut==0.00005] < 0.5 ) / length( results$fec_ratio[results$pmut==0.00005] < 0.5  )
sum( results$fec_ratio[results$pmut==0.005] < 0.5 ) / length( results$fec_ratio[results$pmut==0.005] < 0.5  )

# do the mutation rate correlate with the percentage of asexuals at the middle and end?
boxplot(results$mix~results$pmut, main = "btw max sex and end")
boxplot( (results$nb.sex.end / results$nb.all.end) ~ results$pmut, main = "end sim")
#### OH MY GOD WHAT'S GOING ON WITH INTERMEDIATE MUTATION LEVELS

plot(results$nb.sex.end[results$pmut==0.0005])

par(mfrow=c(1,1))
plot(results$X.sex.end~results$G)
identify((results$X.sex.end~results$G))


## what are those points where the sexuals reach just below 200?
results[results$X.sex.end>150,c(12,20)]



results$mean.X.asex = results$mean.X.asex/results$X.all.end
results$mean.X.sex = results$mean.X.sex/results$X.all.end

results = results_full

results = results_full[results_full$end.t<1000,]
results = results_full[results_full$end.t==1000,]
results = results[results$fec_ratio<=0.5,]
# switch on/off homozygosity
homodep = results$b
homodep[homodep < .9] <- "ON"
homodep[homodep >= .9] <- "OFF"
homodep = as.factor(homodep)
results = results_full
results = results[results$b <= .9,]

par(mfrow=c(3,2))
contrast <- results$compet
contrast <- homodep
for(ii in c(1)){
  # par(mfrow=c(3,4))
  # for(jj in 12:21){
  #   plot( results[,jj] ~ results[,ii], xlab=names(results[ii]), ylab= names(results[jj]))
  #   if(ii > 1) abline(lm( results[,jj] ~ results[,ii] ))
  # }
  plot( results$mean.X.sex ~ contrast )
  plot( results$mean.X.asex ~ contrast )
  plot( results$sd.X.sex ~ contrast )
  plot( results$sd.X.asex ~ contrast )
  plot( (results$mean.X.asex - results$mean.X.sex) ~ contrast )
  # boxplot(results[,jj]~results$c, xlab="c")
  # boxplot(results[,jj]~results$pmut, xlab="pmut")
}
t.test((results$mean.X.asex - results$mean.X.sex) ~ contrast )

par(mfrow=c(2,1))
hist(results$mean.X.asex[results$compet == "_ded_"] - results$mean.X.sex[results$compet == "_ded_"])
hist(results$mean.X.asex[results$compet == "_fc_"] - results$mean.X.sex[results$compet == "_fc_"])

library(vioplot)
vioplot(na.omit(results$mean.X.asex[results$compet == "_ded_"] - results$mean.X.sex[results$compet == "_ded_"]),
        na.omit(results$mean.X.asex[results$compet == "_fc_"] - results$mean.X.sex[results$compet == "_fc_"]))
