# R_code_snow.r

setwd("C:/lab/")
install.packages("ncdf4") # in order to import ncdf file
library(ncdf4)
library(raster)
snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc") # only one layer
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)

setwd("C:/lab/snow/")
snow2000 <- raster("snow2000r.tif")
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

# Faster function to import many data: lapply
# first define the list of files to aggregate together 
rlist <- list.files(pattern="snow20") #all the files have in common "snow20"
# lapply import all of the raster
import <- lapply(rlist, raster)
#  create one raster from several RasterLayers 
snow.multitemp <- stack(import)# MULTITEMPORAL ANALYSISSSSSS, import the RasterStack
plot(snow.multitemp, col=cl)

# let's make a prediction
source("prediction.r") #read r code
