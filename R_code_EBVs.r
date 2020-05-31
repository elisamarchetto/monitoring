### dealing with Essential Biodiversity Variables
# understanding heterogeneity

setwd("C:/lab/")
library(raster)
library(RStoolbox) #for PCA
## raster function import only a single layer, brick multiple layers
snt <- brick("snt_r10.tif")
plot(snt)
plotRGB(snt,3,2,1, stretch="lin")
plotRGB(snt,4,2,1, stretch="lin")
# how the different layers(bands) are related... for a better understanding of PCA analysis
pairs(snt)
sntpca <- rasterPCA(snt)
#information about the output of the model, in other words the percentage of variance related to the components
summary(sntpca$model) # proportion of variance: 0.7015076 for the PC1(good approximation)
plotRGB(sntpca$map, 1, 2, 3, stretch="lin")
#calculte the standard deviation using a moving window (5x5)
# create the moving window: it is a matrix
window <- matrix(1, nrow = 5, ncol = 5) # all the values are set to 1, empty window
# focal function for sd, in this case, works only for RasterLayer
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)
cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) 
plot(sd_snt, col=cl)

##folcal function can be applied also to a image directly taken in field: cladonia.jpg
#library(raster), to open the image with several layers: brick function
clad <- brick("cladonia_stellaris_calaita.JPG")
plotRGB(clad, 1,2,3, stretch="Lin")
window <- matrix(1, nrow = 3, ncol = 3) # 3x3 pixel is the extent of moving window, 1 is a number that not impact the calculation
# recall library RStoolbox for performing PCA analysis. The goal is to apply standard deviation to the principal component
cladpca <- rasterPCA(clad) #the bands are correlated each other hence the PC describe 98%
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd) # show the structural complexity of cladonia, variability of a single organism
PC1_agg <- aggregate(cladpca$map$PC1, fact=10) # resempling the reslution to make focal caculation faster
sd_cladagg <- focal(PC1_agg, w=window, fun=sd)
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
par(mfrow=c(1,2))
plot(sd_clad, col=cl)
plot(sd_cladagg, col=cl) # less accurancy
