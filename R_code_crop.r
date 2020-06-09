setwd("C:/lab/")
library(ncdf4)
library(raster)

snow <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
ext <- c(0, 20, 35, 50) # set the extention desired
zoom(snow, ext=ext)
snowitaly <- crop(snow, ext)
zoom(snow, ext=drawExtent()) # draw directly the extent from the graph
