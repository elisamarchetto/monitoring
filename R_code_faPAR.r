#how to look at chemical cycling from sitellite

library(raster)
library(rasterVis) #for levelplot
library(rasterdiv)

setwd("C:/lab/")
plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA))
levelplot(copNDVI)
faPAR10 <- raster("faPAR10.tif")
levelplot(faPAR10)

pdf("copNDVI.pdf")
levelplot(copNDVI)
dev.off()

 

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

ls() #lookin for faPAR10
faPAR10
#let's see how much space need for 8-bits images
writeRaster(copNDVI, "copNDVI.tif") #write the data copNDVI in .tif, 5.5MB 
#faPAR10 in in bits from 0 to 0.93. Change from 0 to 255
faPAR <- stretch(faPAR10,minv=0,maxv=250)
writeRaster(faPAR, "faPAR.tif")
faPAR
