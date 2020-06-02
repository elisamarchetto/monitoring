# R code for exam project

library(rgdal)
library(gdalUtils)
library(raster)
library(rasterVis)
library(RStoolbox)

setwd("C/lab/exam/")

# Sentinel-2 satellite data
# extrapolation directly in R to optain GeoTIFF
gdal_translate("T30SWG_20191016T110041_B02.jp2", "B02.tif")
gdal_translate("T30SWG_20191016T110041_B03.jp2", "B03.tif")
gdal_translate("T30SWG_20191016T110041_B04.jp2", "B04.tif")
gdal_translate("T30SWG_20191016T110041_B08.jp2", "B08.tif")

octB02 <- raster("B02.tif")
octB03 <- raster("B03.tif")
octB04 <- raster("B04.tif")
octB08 <- raster("B08.tif")

## or directly in SNAP
# Creating a RasterStack

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
ci <- 1 - (biocrast$biocrast.3 - biocrast$biocrast.1) / (biocrast$biocrast.3 - biocrast$biocrast.1) # Crust Index: it doesn't work for such a big extention
plot(ndvi, col=cl)
plot(ci, col=cl)
### let's calulate BSCI: Biological Soil Crust Index 
# mean of Bgreen, Bred, B NIR
b8  <- raster("subset_3_of_T30SWG_20191016T110041_B08.tif")
b4 <- raster("subset_2_of_T30SWG_20191016T110041_B04.tif")
b3 <- raster("subset_1_of_T30SWG_20191016T110041_B03.tif")
r_brick <- brick(b3, b4, b8)
r_brickPCA <- rasterPCA(r_brick)

mean <- calc(r_brickPCA$map$PC1, fun = mean)

# Abs of Bgreen, Bred
b4 <- raster("subset_2_of_T30SWG_20191016T110041_B04.tif")
b3 <- raster("subset_1_of_T30SWG_20191016T110041_B03.tif")
r_brick1 <- brick(b3, b4)
r_brick1PCA <- rasterPCA(r_brick1)
abs <- calc(r_brick1PCA$map$PC1, fun = abs)
# BSCI = (1-L*|B4-B3|) / (meanB8, B4, B3), L=2

bsci <- (1 -2*(abs)) / (mean)
plot(bsci, col=cl)







