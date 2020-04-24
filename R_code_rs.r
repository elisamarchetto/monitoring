# remote sensing

install.packages("RStoolbox")
install.packages("raster")

setwd("C:/lab/")
library(raster)

#import images
p224r63_2011 <- brick("p224r63_2011_masked.grd")
plot(p224r63_2011)
cl <- colorRampPalette(c('black','grey','light grey'))(100)
plot(p224r63_2011, col=cl)
#landsat resolution of 30m (each pixel)...B1 blu, B2 green, B3 red, B4 NIR ecc

#multiframe of different plots
par(mfrow=c(2,2)) # mf for multiframe with a graph 2 x 2 because we are using 4 bands( B1,B2,B3,B4)

clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) #B1
plot(p224r63_2011$B1_sre, col=clb)

clg <- colorRampPalette(c('dark green','green','light green'))(100)#B2 
plot(p224r63_2011$B2_sre, col=clg)

clr <- colorRampPalette(c('dark red','red','pink'))(100)#B3
plot(p224r63_2011$B1_sre, col=clr)

cln <- colorRampPalette(c('red','orange','yellow'))(100) #B4
plot(p224r63_2011$B4_sre, col=cln)

#or I can change with 
par(mfrow=c(4,1))

#plot as the human eyes see the image using RGB components related to the bands
dev.off() #previous work
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") #we can use only 3 bands at the time

plotRGB(p224r63_2011, r=4, g=2, b=1, stretch="Lin") #baceuse the plants are reflecting a lot the NIR





