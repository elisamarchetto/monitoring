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

### regression model between faPAR and NDVI

#exemple:
erosion <- c(12, 14, 16, 24, 26, 40, 55, 67)
hm <- c(30, 100, 150, 200, 260, 340, 460, 600)
plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals")
model1 <- lm(hm ~ erosion) # hm is y ax, and erosion x
summary(model1)
# hm= 9.2752erosion-26.9888, R-squared:  0.9747, random situation??  p-value: 5.127e-06
# line described by slope b and intersection a: abline function
abline(model1)
faPAR10 <- raster("faPAR10.tif") #library(raster)
#library(rasterdiv)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

install.packages("sf")
library(sf) # to call st_* functions
random.points <- function(x,n)
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}
pts <- random.points(faPAR10,1000) # in this case we use 1000 points in stead of all the values in faPAR10

copNDVIp <- extract(copNDVI, pts) # correlate the values of copNDVI with each random points choosen
faPAR10p <- extract(faPAR10,pts)
#build the linear model between copNDVIp and faPAR10 (copNDVIp because the calculation is faster with less values)
# the line is calculated by reducing the distance between the points (x;y) in the graph
# phothosythesis vs biomass
model2 <- lm(faPAR10p ~ copNDVIp) # R = 04 because in conifer forest biomass is high and phothosynthesis but not that high. p = 2 ^-16 two variables are related each other
plot(faPAR10p, copNDVIp)
abline(model2, col="red")
