# Remote Sensing

install.packages("RStoolbox")# toolbox for remote sensing image processing and analysis
install.packages("raster")

setwd("C:/lab/")
library(raster)

#import images
p224r63_2011 <- brick("p224r63_2011_masked.grd") #import RasterBrick: nlayers
plot(p224r63_2011)
cl <- colorRampPalette(c('black','grey','light grey'))(100)
plot(p224r63_2011, col=cl)
#Landsat satellite: image with resolution of 30m (each pixel)...B1 blu, B2 green, B3 red, B4 NIR ecc

#multiframe of different plots
par(mfrow=c(2,2)) # mf for multiframe with a graph 2 x 2 to visualize 4 separated ghaphs of 4 bands( B1,B2,B3,B4)

clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) #B1
plot(p224r63_2011$B1_sre, col=clb)

clg <- colorRampPalette(c('dark green','green','light green'))(100)#B2 
plot(p224r63_2011$B2_sre, col=clg)

clr <- colorRampPalette(c('dark red','red','pink'))(100)#B3
plot(p224r63_2011$B1_sre, col=clr)

cln <- colorRampPalette(c('red','orange','yellow'))(100) #B4
plot(p224r63_2011$B4_sre, col=cln)

# an other way to visualize the 4 bands in a page
par(mfrow=c(4,1))

#plot as the human eyes see the image using RGB components 
dev.off() # close the previous work
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") # it is possible to use only 3 bands at the time
# stretch: stretching improves the appearance of the data by spreading the pixel values along a histogra
plotRGB(p224r63_2011, r=4, g=2, b=1, stretch="Lin") # r for the B4: NIR, the plants reflect much more in NIR

setwd("C:/lab/")
load("remote_sensing.RData")
library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plot(p224r63_1988)

par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

par(mfrow=c(2,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
# to see the noise (clouds/humidity) in the images. Enhance the noise with stretch="hist"
par(mfrow=c(2,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist") # hist: calcutation of the area under the integral that describes the shock of enhancing the color
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist") # amount of humudity was high because of evapotrasnpirantion

# to calculate DVI of 2011
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre #BA: NIR - B3: red.  p224r63_2011$B4_sre linking the layer of B4
cl <- colorRampPalette(c('yellow','light blue','lightpink4'))(100)
plot(dvi2011, col=cl)

# to calculate DVI of 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre
cl <- colorRampPalette(c('yellow','light blue','red'))(100)
plot(dvi1988, col=cl)

# consider the difference between DVI
diff <- dvi2011 - dvi1988
plot(diff)

#change the grain= the dimention of pixels. To see for exemple the corridors. In RS grain=resolution
# function is aggregate()
p224r63_2011res1 <- aggregate(p224r63_2011, fact=10) #fact: factor is the amount of time increase the pixels
p224r63_2011res2 <- aggregate(p224r63_2011, fact=100)

par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res1, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res2, r=4, g=3, b=2, stretch="Lin")



