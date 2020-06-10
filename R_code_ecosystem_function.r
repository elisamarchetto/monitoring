## Global biomass and the effect of its changes in ecosystem functions
# Data copernicus with Sentinel-2
# Chemical cycling

install.packages("rasterdiv") #functions to calculate indices of diversity
library(rasterdiv)
install.packages("rasterVis")#raster visualisation, es. levelplot function
library(rasterVis)
data(copNDVI) # from rasterdiv library. It is a RasterLayer sets at 8-bits. The dataset is the Copernicus Long-term (1999-2017) average Normalise Difference Vegetation Index
plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253, 255, NA), right=TRUE) #removing water pixels 253, 255, NA using cbind argument
copNDVI10 <- aggregate(copNDVI, fact=10)
levelplot(copNDVI10)

library(ggplot2)
myPalette <- colorRampPalette(c('white','green','dark green'))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

 
# to plot NDVI score on the map
ggR(copNDVI, geom_raster = TRUE) +
scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
labs(x="Longitude",y="Latitude", fill="")+
#   theme(legend.position = "bottom") +
  NULL
# +
# ggtitle("NDVI")


setwd("C:/lab/")
defor1 <- brick("defor1_.jpg")
defor2 <- brick("defor2_.jpg")
#band1=NIR, Band2=red, band3=green defor1_.1 or defor2_.1 (band1)

par(mfrow=c(2,1))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

#calculate the DVI for both images

dvi1 <- defor1$defor1_.1 - defor1$defor1_.2 # $ symbol to link the layers, each layer rappresents a band
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100)
par(mfrow=c(1,2))
plot(dvi1, col=cl)
plot(dvi2, col=cl)

difdvi <- dvi1 - dvi2 # diffecence of DVI between the 2 images. Observation at the rate of loss of vagetation between 2011 and 1988
dev.off()
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col=cld) #loss of ecosystem services
hist(difdvi) #high loss in biomass ecosystem services
