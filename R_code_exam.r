# R code for exam project

library(rgdal)
library(gdalUtils)
library(raster)
library(rasterVis)
library(RStoolbox)

setwd("C/lab/exam/")

# Sentinel-2 satellite data

gdal_translate("T30SWG_20191016T110041_B02.jp2", "B02.tif")
gdal_translate("T30SWG_20191016T110041_B03.jp2", "B03.tif")
gdal_translate("T30SWG_20191016T110041_B04.jp2", "B04.tif")
gdal_translate("T30SWG_20191016T110041_B08.jp2", "B08.tif")

octB02 <- raster("B02.tif")
octB03 <- raster("B03.tif")
octB04 <- raster("B04.tif")
octB08 <- raster("B08.tif")

#or directly in SNAP
rlist <- list.files(pattern="sub")
import1 <- lapply(rlist, raster)
biocrast <- stack(import1)
agg_biocrast <- aggregate(biocrast, fact=10)
dvi <- agg_biocrast$subset_3_of_T30SWG_20191016T110041_B08 - agg_biocrast$subset_2_of_T30SWG_20191016T110041_B04
cl <-  colorRampPalette(c("black", "green", "red"))(100)
clb <-  colorRampPalette(c("black", "gold", "blue"))(100)
plot(dvi, col=cl)
writeRaster(biocrast, "biocrast.tif")
writeRaster(agg_biocrast, "agg_biocrast.tif")

ndvi <- (biocrast$biocrast.4 - biocrast$biocrast.3) / (biocrast$biocrast.4 + biocrast$biocrast.3)
ci <- 1 - (biocrast$biocrast.3 - biocrast$biocrast.1) / (biocrast$biocrast.3 - biocrast$biocrast.1) # it doesn't work
plot(ndvi, col=cl)
plot(ci, col=cl)




mean <- calc(STACK1, fun = mean)


