### dealing with Essential Biodiversity Variables
# understanding heterogeneity

setwd("C:/lab/")
library(raster)
library(RStoolbox)
## raster function import only a single layer, brick multiple layers
snt <- brick("snt_r10.tif")
plot(snt)
plotRGB(snt,3,2,1, stretch="lin")
plotRGB(snt,4,2,1, stretch="lin")
# how the different layers(bands) are related... for PCA analysis
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
