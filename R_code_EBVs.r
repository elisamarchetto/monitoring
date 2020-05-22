### dealing with Essential Biodiversity Variables
# understanding heterogeneity

setwd("C:/lab/")
library(raster)
## raster function import only a single layer, brick multiple layers
snt <- brick("snt_r10.tif")
