#R code for multivariate analysis

install.packages("vegan") #vegetation and diversity analysis
library(vegan)
setwd("C:/lab/")
biomes<-read.table("biomes.csv", header=T, sep=",") # sep="," in cvs file the variables are separeted by comma, biomes.cvs is the data frame
head(biomes)

## Multivariate analysis
#DEtrended CORrespondance ANAlysis---to reduce the dimenctions (variables)
#the function is decorana

multivar<-decorana(biomes)
plot(multivar)
plot(multivar, cex=1.2)

 #               DCA1   DCA2    DCA3    DCA4
Eigenvalues     0.5117 0.3036 0.12125 0.14267
Decorana values 0.5360 0.2869 0.08136 0.04814
Axis lengths    3.7004 3.1166 1.30055 1.47888

# DCA1 describe the 0.5117 52% of percent of the variables..+DCA2%: DEtrended CORrespondance ANAlysis descrbes the 80%, thus the 20% is lost( DCA3,DCA4).. the species taht are not inside the convex line are outside because there is an inaccurancy 20%
#the points into the graph are the plots, so the 20 dimentions

#how referring the plots (species belonging to the plot) wiyh their own biomes
biomes_types <- read.table("biomes_types.csv", header=T, sep=",")
#to use the columns we need to attach
attach(biomes_types)


# linking the biome types with the plots X species
ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3) #col for each biomes, kind for the type of grafh, lwd for the dimention of the colored line

ordispider(multivar, type, col=1:4, label = T) # type is the column we want to see, label is refered to the name of type column we want to see in the graph


