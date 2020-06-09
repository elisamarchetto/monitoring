# R_code_no2.r
library(raster)

setwd("C:/lab/no2/")
rlist <- list.files(pattern="EN")
import <- lapply(rlist, raster)
EN <- stack(import)
cl <- colorRampPalette(c('red','orange','yellow'))(100) #
plot(EN, col=cl)

par(mfrow=c(1,2)) # to see the change through time
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)
# 3 layers RGB image from EN. See where and when Europe was more affected by NO2 pollution
plotRGB(EN, r=1, g=7, b=13, stretch="lin")
#difference map between 2 situations
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) # 
plot(dif, col=cld) 

# Quantitative decrease of no2 
boxplot(EN) # five-number summary is the minimum, first quartile, median, third quartile, and maximum. The first quartile is the median of the data points to the left of the median.
boxplot(EN,outline=F)
boxplot(EN,outline=F, horizontal=T)
boxplot(EN,outline=F, horizontal=T, axes=T) 

plot(EN$EN_0001, EN$EN_0013)
abline(0,1,col="red") # to see if the values are under the line, this means decrease!! 

# fast version of import and plot of many data for lazy people!
rlist <- list.files(pattern="snow")
import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1) # most of the value under the curve
plot(snow.multitemp$snow2000r, snow.multitemp$snow2020r) # better change in time
abline(0,1,col="red")
