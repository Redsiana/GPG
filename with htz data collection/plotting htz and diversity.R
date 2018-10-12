setwd("C:/Users/Anais Tilquin/Desktop/Simulations February/hermaphro/b0.8 3niches/auto apo auto control/auto_apo")

library(lattice)
levelplot(t(htz_matrix[1:t,]), main = "average heterozygosity",col.regions = gray(95:0/100) )
levelplot(t(diversity_matrix[1:t,]), main = "genetic diversity", col.regions = grey(95:0/100))
levelplot(t(sex_matrix[1:t,]), main = "frequency of sexuals", col.regions = gray(95:0/100))
plot(INVASION_tot)
points(INVASION_tot-INVASION_sex, col="blue")
points(INVASION_sex, col="red")


require(gridExtra)
require(grid)

files <- list.files(pattern = "\\T2.RData$")

for(idx in 1:length(files)){
  load(files[idx])
  plot1 <- levelplot(t(sex_matrix[1:t,]), main = "frequency of sexuals", 
                     col.regions = gray(95:0/100), xlab = "Space", ylab = "Time", panel=function(...) { 
                       grid.rect(gp=gpar(col=NA, fill="bisque")) 
                       panel.levelplot(...) 
                     } #,
                     #scales=list(x=x.scale, y=y.scale)
                     )
  plot2 <- levelplot(t(htz_matrix[1:t,]), main = "average heterozygosity",
                     col.regions = gray(95:0/100), xlab = "Space", ylab = "Time", panel=function(...) { 
                       grid.rect(gp=gpar(col=NA, fill="bisque")) 
                       panel.levelplot(...) 
                     } )
  plot3 <- levelplot(t(diversity_matrix[1:t,]), main = "genetic diversity", 
                     col.regions = grey(95:0/100), xlab = "Space", ylab = "Time", panel=function(...) { 
                       grid.rect(gp=gpar(col=NA, fill="bisque")) 
                       panel.levelplot(...) 
                     } )
  
  
  png(filename=paste("plot", namerun2, ".png", sep=""), 
      # type="cairo",
      units="in", 
      width=40,
      height=20,
      pointsize=40,
      res=96)
  grid.arrange(plot1,plot2, plot3, ncol=3, top=namerun2)
  dev.off()
  
print(idx)
}

