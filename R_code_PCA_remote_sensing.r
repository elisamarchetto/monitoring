
setwd("C:/lab/")
library(raster)
library(RStoolbox) # for performing PCA

p224r63_2011 <- brick("p224r63_2011_masked.grd") 
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")


## An other way to visualize RGB image
library(ggplot2) 
ggRGB(p224r63_2011,5,4,3)

p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
ggRGB(p224r63_1988,5,4,3)

par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

#reducing the nlayers (bands): PCA
#PricipalComponentAnalysis

names(p224r63_2011) #to see the names of variables
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre) #to see if the variables bands are correlated, better with pairs function
par(mfrow=c(3,1)) #for exemple to see the pattern correlation 
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)
plot(p224r63_2011$B1_sre, p224r63_2011$B4_sre)
plot(p224r63_2011$B1_sre, p224r63_2011$B2_sre)
#or
pairs(p224r63_2011)

p224r63_2011_res <- aggregate(p224r63_2011, fact=10) # accurancy is lower but R calculation is faster
#library needed RStoolbox
p224r63_2011_pca <- rasterPCA(p224r63_2011_res)
p224r63_2011_pca
plot(p224r63_2011_pca$map)#plot all of PC. PC1 is accounting for most of the variations 
#in $model show the PC in proportion
summary(p224r63_2011_pca$model) # PC weight in percentage, how a PC summarizes the variables
pairs(p224r63_2011)
#plot PC1,PC2,PC3
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin") # 1,2,3 related to the number associated to the PC

p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res)
p224r63_1988_pca
plotRGB(p224r63_1988_pca$map, r=1, g=2, b=3, stretch="Lin")
#difference in PCA in time
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
cldif <- colorRampPalette(c('blue','black','yellow'))(100)
plot(difpca, col=cldif)
plot(difpca$PC1,col=cldif)
plotRGB(difpca, r=1, g=2, b=3, stretch="Lin") #highest possible variation and where




