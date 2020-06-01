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




