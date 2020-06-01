# R code for exam project

library(rgdal)
library(gdalUtils)
library(raster)
library(rasterVis)
library(RStoolbox)

setwd("C/lab/")

gdal_translate("T33XVJ_20170803T125711_B01.jp2", "B02.tif")
gdal_translate("T33XVJ_20170803T125711_B01.jp2", "B03.tif")
gdal_translate("T33XVJ_20170803T125711_B01.jp2", "B04.tif")
gdal_translate("T33XVJ_20170803T125711_B01.jp2", "B08.tif")

