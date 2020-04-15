#R code for multivariate analysis

install.packages("vegan") #vegetation analysis
library(vegan)
setwd("C:/lab/")
biomes<-read.table("biomes.csv", header=T, sep=",")
head(biomes)

#make use of multivarite analysis
#DEtrended CORrespondance ANAlysis---to reduce the dimenctions
#function: decorana

multivar<-decorana(biomes)
plot(multivar)
plot(multivar, cex=1.2)

 #               DCA1   DCA2    DCA3    DCA4
Eigenvalues     0.5117 0.3036 0.12125 0.14267
Decorana values 0.5360 0.2869 0.08136 0.04814
Axis lengths    3.7004 3.1166 1.30055 1.47888

# 0.5117 52% of percent..+DCA2%. arrived at 80%, thus we lost 20% ( DCA3,DCA4).. the species we don't see inside the convex shape are out side because we lost 20%
#the points into the graph are the plots, so the 20 dimentions

#put the plots in their own biomes
biomes_types <- read.table("biomes_types.csv", header=T, sep=",")
#to use the columns we need to attach
attach(biomes_types)


# to link the biome types with the plots X species
ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3) #col for each biomes, kind for the type of grafh, lwd for the dimention of the colored line

ordispider(multivar, type, col=1:4, label = T) # type is the column we want to see, label is refered the the name of type column we want to see in the graph


