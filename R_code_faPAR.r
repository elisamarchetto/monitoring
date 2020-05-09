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
