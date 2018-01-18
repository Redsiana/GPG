#### invectigate clone identity, who is the most common, when did she originate
## and where
## watch out: popX is taken now after dispersal, before survival

library(ggplot2)
library(RColorBrewer)

# load(file = '_ded_K10_fs4_fa2_1_G3_b0.6_Ka0.4_B20M_0.3_pm0.005_c1.RData')



sorted_clone_prevalence <- sort( table(popcloneline)[-1])# gives table with clone number and number of individuals
main_clones_names <- as.numeric( names( tail(sorted_clone_prevalence, 5) ))

clone_bars <- as.numeric( table(popX[ !(popcloneline %in% c(0,main_clones_names)) ]) )
clone_patches <- as.numeric( names( table(popX[ !(popcloneline %in% c(0,main_clones_names)) ])) )
clone_idcol <- as.factor( rep( "other cl", length(clone_bars)) )
df <- data.frame(x = clone_patches, y = clone_bars, fac = clone_idcol )

x_emergence <- numeric(length = length(main_clones_names))
g_emergence <- numeric(length = length(main_clones_names))
g_edge <- numeric(length = length(main_clones_names))
label_clone <- numeric(length = length(main_clones_names))

for(clone in 1:length(main_clones_names)){
  clone_bars <- as.numeric( table(popX[ popcloneline == main_clones_names[clone]]) )
  clone_patches <- as.numeric( names( table(popX[ popcloneline == main_clones_names[clone]])) )
  clone_idcol <- as.factor( rep( main_clones_names[clone], length(clone_bars)) )
  df <- rbind(df, data.frame(x = clone_patches, y = clone_bars, fac = clone_idcol ))

  x_emergence[clone] <- unique( popclonalorigin[popcloneline %in% main_clones_names[clone],1] )
  g_emergence[clone] <- unique( popclonalorigin[popcloneline %in% main_clones_names[clone],2] )
  g_edge[clone] <- EDGE_t[g_emergence[clone]]
  label_clone[clone] <- paste(g_emergence[clone], ":", x_emergence[clone], '-', g_edge[clone])
}

sexual_bars <- as.numeric( table(popX[ popcloneline == 0]) )
sexual_patches <- as.numeric( names( table(popX[ popcloneline == 0])) )
sexual_idcol <- rep("sex", length(sexual_bars))
df <- rbind( df, data.frame(x = sexual_patches, y = sexual_bars, fac = sexual_idcol ) )



# this df has 3 columns, the X-coordinates, the number of individuals, the identity (5 main clones, other cl or sex)
mypalette <- c( brewer.pal((length(main_clones_names)+1),  "YlOrRd"), "black" )

labels_plot <- c("other cl", label_clone, "sex")

ggplot(data=df, aes(x=x, y=y, fill=fac, order=fac)) +
  geom_bar(stat="identity" ) +   theme_minimal() +
  scale_x_continuous(limits=c(0,200)) +
  scale_y_continuous(limits=c(0,600)) +
  scale_fill_manual(values= mypalette, labels = labels_plot) +
  geom_vline(xintercept = jitter(x_emergence), col=mypalette[2:(length(main_clones_names)+1)], size = .8 )+
  ggtitle(paste("pop composition after", t, "generations")) +
  theme(aspect.ratio = 5/10)

# this prints, for each of the 5 main clones, the generation at which it appeared, the x coordinates, and the coordinates of the edge back then

ggsave(filename = paste(namerun,t, ".jpg", sep=""))





# table(popX[ popcloneline == 41])
# 
# unique( popclonalorigin[popcloneline==126,] ) # emerged on 1. km 15 at 2. generation 9
# EDGE_t[9] # it's 14, so the clone arose beyond the edge
# 
# hist( popX[ popcloneline == 126], breaks=100 )
# abline(v = 15, col="red") # it ran off! but didn't go back
# 
# 
# unique( popclonalorigin[popcloneline==92,] ) # emerged on km 9 at generaion 6
# EDGE_t[6] # it's 11, so the clone arose beyond the edge
# 
# hist( popX[ popcloneline == 92], breaks=100 )
# abline(v = 16) # it progressed, and went back!
# 
# hist( popX[ popcloneline == 0] , xlim=c(0,75))
