# R code for exam 

#####
#1. R_code_first.r
#2. R_code_multipanel.r
#3. R_code_Spatial.r
#4. R_code_point_pattern_analysis.r
#5. R_code_multivariate_analysis.r
#6. R_code_remote_sensing.r
#7. R_code_ecosystem_functions.r
#8. R_code_multivariate_analysis_RS_data.r
#9. R_code_ecosystem's_reflectance.r
#10. R_code_faPAR.r
#11. R_code_EBVs.r
#12. R_code_snow_cover.r
#13. R_code_monitoring_air_pollution_no2.r
#14. R_code_crop.r
#15. R_code_interpolation.r
#16. R_code_sdm.r
#17. R_code_myproject_exam.r
######


######1. R_code_first.r
# test for the first time functions and codes in R software

install.packages("sp") # function to install a package. sp is a Classes and methods for spatial data
library(sp)# function to recall the package installed

data(meuse) #recalling the dataset

# have a look at the structure of dataset:
meuse

#having a look at the first part of dataset
head(meuse)

#correlete two varialables
attach(meuse) # attach variables for using plot function
plot(zinc,copper) #displaying graph
plot(zinc,copper, col="green") # col means color 
plot(zinc,copper, col="green",pch=19) # pch is a command to change the symbol in the graph
plot(zinc,copper, col="green",pch=19,cex=2) # cex is the command for the dimention of pch in this case

######2. R_code_multipanel.r

#multipanel in R: monitoring ecosystems

install.packages("GGally") #'GGally' extends 'ggplot2' by adding several functions to reduce the complexity of combining geometric objects
library(sp) # recall the package sp
#require(sp) is the same of library(sp)
data(meuse) #recalling the dataset
attach(meuse) #make use of data set of meuse

meuse #to see all the variables
plot(cadmium) # plot a variable

#make pairs vairables in plot
pairs(meuse) #plot matrix, consisting of scatterplots for each variable-combination of a data frame. pairs is a function stored in GGally (in this case)

#to correlate only some variables in the plot. It needs to be indicated the dataset that is going to be used
pairs(~ cadmium+copper+lead+zinc,data=meuse)
#or
pairs(meuse[,3:6]) # [,3:6] from 3(cadmium) to 6(zinc), it is the condition

pairs(meuse[,3:6],pch=19)

library(GGally)
ggpairs(meuse[,3:6]) # to prettify the scatterplots and having more information about the dispersion of the values for each variables

#######3. R_code_Spatial.r

#see how to use coordinates in space related to the variables.....using the function: coordinates 

library(sp)

data(meuse) #to recall the data
head(meuse) # head function shows first rows of the dataset

#coordinates is the functionn to visualize the variables in the space 
coordinates(meuse)=~x+y #thinking spatialy
plot(meuse) # plot all the variables in meuse
spplot(meuse, "zinc") #for a spatial variable (it dipends on coordinates function), spplot dislpays the amount of a variable for each location

#exercise: spatial amount of copper
spplot(meuse, "copper") # asign the dataset and the variable

#change the title in the graph using main= "....
spplot(meuse, "copper", main="Copper concentration")

# displaying "zinc" spatialy (plot) duplicating the size of the symbol according with the dimention of the value
bubble(meuse, "zinc")
#exercise: change the color of zinc in red
bubble(meuse, "zinc", col="red")

#download covid_aff in lab into C:
#inport data base covid in R
#setting the working directory: lab
setwd("C:/lab")
covid<-read.table("covid_agg.csv", head=T) # the file contains the names of the variables as its first line, header of the colums true
head(covid)

attach(covid) # make available the variables of the dataset
plot(country,cases)
plot(country,cases, las=0) #paralel labels
plot(country,cases, las=1) #horizontal labels
plot(country,cases, las=2) #perperdincolar labels
plot(country,cases, las=3) #vertical labels
plot(country,cases, las=3, cex.axis=0.5) #cex for the size of labels axis
install.packages("ggplot2") # A system for "declaratively" creating graphics, based on "The Grammar of Graphics"

# to use the function in ggplot2 require the package
library(ggplot2)

#to open the saved work
setwd("C:/lab") # set the working directory in which the data are located
load("R_code_spatial")
ls() #list of objects used

library(ggplot2)
#dataset to use
data(mpg)
head(mpg)
#function we are going to use ggplot.. components: data, aes, geometry
ggplot(mpg,aes(x=displ,y=hwy)) + geom_point()#aes: aesthetic mapping choosing the variables of the dataset to visualize, geom_ :Specifies the geometric objects that define the graph type
ggplot(mpg,aes(x=displ,y=hwy)) + geom_line() # values linked as line
ggplot(mpg,aes(x=displ,y=hwy)) + geom_polygon() # as polygon

head(covid)
ggplot(covid,aes(x=lon,y=lat, size=cases)) + geom_point() #lat and lon dipend on size, the size is the cases variable

######4. R_code_point_pattern_analysis.r

#Point pattern analysis: density map

install.packages("spatstat") # function to install a package. spatstat is a toolbox for analysing Spatial Point Patterns
library(spatstat) #recalling the package
attach(covid) #make use of variables
head(covid) # seeing the first rows

covids <- ppp(lon, lat, c(-180, 180), c(-90, 90)) #ppp means panel point pattern. The general form is ppp(x.coordinates, y.coordinates, x.range, y.range) it creates a point pattern dataset in the two-dimensional plane for the range of values of lat and long. 
#to create a density map
d<- density(covids) 
plot(d) 
#to see the point of covids object on the map
points(covids)

#open the last work: load()
setwd("C:/lab/")
load("point_pattern_analysis.RData") #loading the working directory saved
library(spatstat)

#to use vector format in coastline
install.packages("rgdal") #Geospatial Data Abstraction. Data georeferenced
library(rgdal)

coastlines <- readOGR("ne_10m_coastline.shp") #readORG is rgdal function to read shapefile goreferenced
plot(d)
points(covids)
plot(coastlines, add=T)
#or
install.packages("rnaturalearth")
library(rnaturalearth)
coastlines <- rnaturalearth::ne_download(scale = 10, type = 'coastline', category = 'physical')

#let's change the color of the graph
cl <- colorRampPalette(c("yellow","orange","red")) (100)
plot(d, col=cl)
points(covids)
plot(coastlines, add=T)

#export in Pdf
pdf("covid_density.pdf")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off() # to close the operation
#or export in png
png("covid_density.png")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()

## 5. R_code_multivariate_analysis.r

#R code for multivariate analysis

install.packages("vegan") #vegetation and diversity analysis
library(vegan) #recalling the package
setwd("C:/lab/")
biomes<-read.table("biomes.csv", header=T, sep=",") #header of the colums, sep="," in cvs file the variables are separeted by comma, biomes.cvs is the data frame
head(biomes) # seeing first rows

## Multivariate analysis
#DEtrended CORrespondance ANAlysis---to reduce the dimenctions (variables)
#the function is decorana

multivar<-decorana(biomes) # decorana is the function to reduce the dimentions
plot(multivar)
plot(multivar, cex=1.2)

 #               DCA1   DCA2    DCA3    DCA4
Eigenvalues     0.5117 0.3036 0.12125 0.14267
Decorana values 0.5360 0.2869 0.08136 0.04814
Axis lengths    3.7004 3.1166 1.30055 1.47888

# DCA1 describe the 0.5117 52% of percent of the variables..+DCA2%: DEtrended CORrespondance ANAlysis descrbes the 80%, thus the 20% is lost( DCA3,DCA4).. the species taht are not inside the convex line are outside because there is an inaccurancy 20%
#the points into the graph are the plots, so the 20 dimentions

#how referring the plots (species belonging to the plot) wiyh their own biomes
biomes_types <- read.table("biomes_types.csv", header=T, sep=",")
#to use the columns we need to attach
attach(biomes_types)


# linking the biome types with the plots X species
ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3) # function to connect the plots by a ellipse with their own biome. col for each biomes, kind for the type of grafh, lwd for the dimention of the colored line

ordispider(multivar, type, col=1:4, label = T) #  return the coordinates to which each point is connected; type is the column we want to see, label is refered to the name of type column we want to see in the graph

######6. R_code_remote_sensing.r

# Remote Sensing

install.packages("RStoolbox")# toolbox for remote sensing image processing and analysis
install.packages("raster") # Analyzing and modeling of gridded spatial data

setwd("C:/lab/")
library(raster) # recalling the package to make use of the functions

#import images
p224r63_2011 <- brick("p224r63_2011_masked.grd") #import a file as RasterBrick: nlayers
plot(p224r63_2011)
cl <- colorRampPalette(c('black','grey','light grey'))(100) 
plot(p224r63_2011, col=cl)
#Landsat satellite: image with resolution of 30m (each pixel)...B1 blu, B2 green, B3 red, B4 NIR for this RasterBrick

#multiframe of different plots
par(mfrow=c(2,2)) # multiple graphs in a single plot: 2rows x 2 columns to visualize 4 separated ghaphs of 4 bands( B1,B2,B3,B4)

clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) #B1
plot(p224r63_2011$B1_sre, col=clb)# $B1_sre linking layer 1 of the RasterBrick

clg <- colorRampPalette(c('dark green','green','light green'))(100)#B2 
plot(p224r63_2011$B2_sre, col=clg) # linking B2

clr <- colorRampPalette(c('dark red','red','pink'))(100)#B3
plot(p224r63_2011$B1_sre, col=clr)

cln <- colorRampPalette(c('red','orange','yellow'))(100) #B4
plot(p224r63_2011$B4_sre, col=cln)

# an other way to visualize the 4 bands in a page
par(mfrow=c(4,1))

#plot as the human eyes see the image using RGB components. Graph with colors that our eyes perceive
dev.off() # close the previous work
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin") # it is possible to use only 3 bands at the time; red component "r" links layer 3
# stretch: stretching improves the appearance of the data by spreading the pixel values along a histogram; improve clarity and contrast 
plotRGB(p224r63_2011, r=4, g=2, b=1, stretch="Lin") # r for the B4: NIR, the plants reflect much more in NIR

setwd("C:/lab/")
load("remote_sensing.RData")
library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd") #import rasterBrick: nlayers
plot(p224r63_1988) # displaying nlayers on the map

par(mfrow=c(2,1)) # plot RGB images to look at the differences among the years
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

par(mfrow=c(2,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin") # r=4 refered to Band NIR; highest reflectance of NIR for vegetation
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin")
# to see the noise (clouds/humidity) in the images. Enhance the noise with stretch="hist"
par(mfrow=c(2,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist") # hist: calcutation of the area under the integral that describes the shock of enhancing the color
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist") # amount of humudity was high because of evapotrasnpirantion

# to calculate DVI of 2011
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre #BA: NIR - B3: red. DVI is used to quantify vegetation greenness(how higher is NIR than red, more healthy is the plant). p224r63_2011$B4_sre linking the layer of B4
cl <- colorRampPalette(c('yellow','light blue','lightpink4'))(100)
plot(dvi2011, col=cl)

# to calculate DVI of 1988
dvi1988 <- p224r63_1988$B4_sre - p224r63_1988$B3_sre 
cl <- colorRampPalette(c('yellow','light blue','red'))(100)
plot(dvi1988, col=cl)

# consider the difference between DVI
diff <- dvi2011 - dvi1988
plot(diff) # loss and gain of vegetation

#change the grain= the dimention of pixels. In RS grain=resolution
# function is aggregate(), obtaining lower resolution and bigger cell pixel.
p224r63_2011res1 <- aggregate(p224r63_2011, fact=10) #fact: horizontal and vertical aggregation factor of the pixel
p224r63_2011res2 <- aggregate(p224r63_2011, fact=100)

par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res1, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res2, r=4, g=3, b=2, stretch="Lin")

######7. R_code_ecosystem_functions.r

## Global biomass and the effect of its changes in ecosystem functions
# Data copernicus with Sentinel-2
# Chemical cycling

install.packages("rasterdiv") #functions to calculate indices of diversity
library(rasterdiv) # recalling package
install.packages("rasterVis")#raster visualisation, es. levelplot function
library(rasterVis)
data(copNDVI) # from rasterdiv library. It is a RasterLayer sets at 8-bits. The dataset is the Copernicus Long-term (1999-2017) average Normalise Difference Vegetation Index
plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253, 255, NA), right=TRUE) #removing water pixels 253, 255, NA using cbind argument; cbind(combines vector, matrix or data frame by columns): combine from 253 to 255 and NA pixel values
copNDVI10 <- aggregate(copNDVI, fact=10) # fact of 10 is a horizoltal and vertical aggregation of pixel by a value of 10 pixel
levelplot(copNDVI10) #plot a raster object and displaying the level of the values along x and y axis

#library(ggplot2)
#myPalette <- colorRampPalette(c('white','green','dark green'))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

 
## to plot NDVI score on the map
#ggR(copNDVI, geom_raster = TRUE) +
#scale_fill_gradientn(name = "NDVI", colours = myPalette(100))+
#labs(x="Longitude",y="Latitude", fill="")+
#   theme(legend.position = "bottom") +
  #NULL
# +
# ggtitle("NDVI")


setwd("C:/lab/") # setting the working directory
defor1 <- brick("defor1_.jpg") # import file as a raster with multiple layers
defor2 <- brick("defor2_.jpg")
#band1=NIR, Band2=red, band3=green defor1_.1 or defor2_.1 (band1)

par(mfrow=c(2,1)) # compare the RGB images
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
dev.off() # close the graph
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col=cld) #loss of ecosystem services
hist(difdvi) #histogram of difdvi: high loss in biomass ecosystem services

######8. R_code_multivariate_analysis_RS_data.r

setwd("C:/lab/")
library(raster) #recalling the package
library(RStoolbox) # for performing PCA

p224r63_2011 <- brick("p224r63_2011_masked.grd") #import file as rasterBrick
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin") # RGB image, r=5 for band NIR


## An other way to visualize RGB image
library(ggplot2) 
ggRGB(p224r63_2011,5,4,3)

p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
ggRGB(p224r63_1988,5,4,3)

par(mfrow=c(1,2)) # see the difference between the images in one map
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

#reducing the nlayers (bands): PCA
#PricipalComponentAnalysis

names(p224r63_2011) #to see the names of variables
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre) #to see if the variables bands are correlated, better with pairs function
par(mfrow=c(3,1)) #for exemple to see the pattern correlation 
plot(p224r63_2011$B1_sre, p224r63_2011$B3_sre)
plot(p224r63_2011$B1_sre, p224r63_2011$B4_sre)
plot(p224r63_2011$B1_sre, p224r63_2011$B2_sre)
#or
pairs(p224r63_2011) # easier and in one map

p224r63_2011_res <- aggregate(p224r63_2011, fact=10) # reducing the resolution of a factor of 10, fact: horizontal and vertical aggregation factor of the pixel, accurancy is lower but R calculation is faster
#library needed RStoolbox
p224r63_2011_pca <- rasterPCA(p224r63_2011_res) # function for applying PCA
p224r63_2011_pca
plot(p224r63_2011_pca$map)#plot all of PC. PC1 is accounting for most of the variations 
#in $model show the PC in proportion
summary(p224r63_2011_pca$model) # PC weight in percentage, how a PC summarizes the variables
pairs(p224r63_2011)
#plot PC1,PC2,PC3
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin") # 1,2,3 related to the number associated to the PC

p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res)
p224r63_1988_pca
plotRGB(p224r63_1988_pca$map, r=1, g=2, b=3, stretch="Lin")
#difference in PCA in time
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
cldif <- colorRampPalette(c('blue','black','yellow'))(100)
plot(difpca, col=cldif)
plot(difpca$PC1,col=cldif)
plotRGB(difpca, r=1, g=2, b=3, stretch="Lin") #highest possible variation and where

######9. R_code_ecosystem's_reflectance.r
#R_code_ecosystem's_reflectance.r

library(raster) # recalling the package
toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2) # for creating the matrix as a layer
values(toy) <- c(1.13,1.44,1.55,3.4) #asign the values to each cells
plot(toy) # displaying on a map
text(toy, digits=2) # max of digits after comma 2

#bits: possible combinations of pixels based on NDs
toy2bits <- stretch(toy,minv=0,maxv=3) # matrix with 2^2 pixel (2 bits), so 4 values. Stretch the pixel values from a min of 0 to max 3
storage.mode(toy2bits[]) = "integer" # integer values
plot(toy2bits) 
text(toy2bits, digits=2) #lower is the amount of bits lower is the "diversity"
toy4bits <- stretch(toy,minv=0,maxv=15) #2^4 combinations
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2)
toy8bits <- stretch(toy,minv=0,maxv=255) #better discrimination from one objects to an other. 2^8 combinations
storage.mode(toy8bits[]) = "integer"
plot(toy8bits)
text(toy8bits, digits=2) 

par(mfrow=c(1,4)) # see on 1 map the difference of the pixel after having asigning bits values

plot(toy)
text(toy, digits=2)
plot(toy2bits)
text(toy2bits, digits=2)
plot(toy4bits)
text(toy4bits, digits=2)
plot(toy8bits)
text(toy8bits, digits=2)

######10. R_code_faPAR.r

#faPAR: Essential Climate Variable
# Chemical cycling from sitellite: carbon cycle
# Regression Analysis

library(raster) # recalling the package
library(rasterVis) #for levelplot
library(rasterdiv)

setwd("C:/lab/") # set the working directory
plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # remove pixels related to water 253, 255, NA using cbind argument
levelplot(copNDVI) #plot a raster object and displaying the level of the values along x and y axis
faPAR10 <- raster("faPAR10.tif") # import file as one layer; file already aggregated of fact 10
levelplot(faPAR10) #### faPAR: Fraction of Absorbed Photosynthetically Active Radiation. Proxy of carbon dioxide assimilation and canopy's energy absorption capacity
#creating a pdf
pdf("copNDVI.pdf")
levelplot(copNDVI)
dev.off()

 

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

ls() #lists all the objects in the working environment
faPAR10
#let's see how much space need for 8-bits images
writeRaster(copNDVI, "copNDVI.tif") #write the data copNDVI in .tif, 5.5MB 
#function to pass at 8-bits(from 0 to 255 combinations of DNs)
faPAR <- stretch(faPAR10,minv=0,maxv=255)
writeRaster(faPAR, "faPAR.tif")
faPAR

### regression model between faPAR and NDVI

#exemple:
erosion <- c(12, 14, 16, 24, 26, 40, 55, 67) # vector of values
hm <- c(30, 100, 150, 200, 260, 340, 460, 600)
plot(erosion, hm, col="red", pch=19, xlab="erosion", ylab="heavy metals")
model1 <- lm(hm ~ erosion) # lm is the function to set the linear model; hm is y axis, and erosion x axis
summary(model1)
# equation: hm= 9.2752erosion-26.9888, R-squared:  0.9747(prediction for the correctness of the model), p-value: 5.127e-06 (the result of the correlation is not a casuality)
# line described by slope b and intersection a: abline function
abline(model1)
faPAR10 <- raster("faPAR10.tif") #library(raster)
#library(rasterdiv)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

install.packages("sf")
library(sf) # to call st_* functions; to encode spatial vector data
random.points <- function(x,n)
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}
pts <- random.points(faPAR10,1000) # in this case we use 1000 points in stead of all the values in faPAR10

copNDVIp <- extract(copNDVI, pts) # extract 1000 random points from copNDVI
faPAR10p <- extract(faPAR10,pts)
#build the linear model between copNDVIp and faPAR10 (copNDVIp because the calculation is faster with less values)
# the line is calculated by reducing the distance between the points (x;y) in the graph
# phothosythesis vs biomass
model2 <- lm(faPAR10p ~ copNDVIp) # R = 0.4 because in conifer forest biomass is high and phothosynthesis but not that high. p = 2 ^-16 two variables are related each other
plot(faPAR10p, copNDVIp)
abline(model2, col="red")

######11. R_code_EBVs.r

### dealing with Essential Biodiversity Variables
# understanding ecosystem structure, the heterogeneity with standard deviation

setwd("C:/lab/") # set working directory
library(raster) # recall package
library(RStoolbox) #for PCA
## raster function import only a single layer, brick multiple layers
snt <- brick("snt_r10.tif")
plot(snt)
plotRGB(snt,3,2,1, stretch="lin")# r=3, g=2, b=1; plot RGB image
plotRGB(snt,4,2,1, stretch="lin") # r=4, band 4 is NIR

pairs(snt)# how the different layers(bands) are related... for a better understanding of PCA analysis
sntpca <- rasterPCA(snt)
#information about the output of the model, in other words the percentage of variance related to the components
summary(sntpca$model) # proportion of variance: 0.7015076 for the PC1(good approximation), almost the 70%
plotRGB(sntpca$map, 1, 2, 3, stretch="lin")
#calculte the standard deviation using a moving window (5x5)
# create the moving window: it is a matrix that moves by 5x5 pixel and the result is 1 final pixel. It reduces the calculation of focal function
window <- matrix(1, nrow = 5, ncol = 5) # all the values are set to 1, empty window
# focal function for sd, standard deviation, it works only for RasterLayer. In raster package
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd)# fun is the function to be calculated by focal, sd is the standard deviation:  measure of the amount of variation of a set of a values, w is the moving window
cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) 
plot(sd_snt, col=cl)

##folcal function can be applied also to an image directly taken in field: cladonia.jpg
# detecting the organism structure by standard deviation

#library(raster), to open the image with several layers: brick function
clad <- brick("cladonia_stellaris_calaita.JPG")
plotRGB(clad, 1,2,3, stretch="Lin")
window <- matrix(1, nrow = 3, ncol = 3) # 3x3 pixel is the extent of moving window (moving through the matrix by 3x3 pixels, obtainig 1 pixel). 1 is a number that not impact the calculation
# recall library RStoolbox for performing PCA analysis. The goal is to apply standard deviation to the principal component
cladpca <- rasterPCA(clad) #the bands are correlated each other hence the PC describe 98%
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd) # show the structural complexity of cladonia, variability of a single organism
PC1_agg <- aggregate(cladpca$map$PC1, fact=10) # resempling the reslution to make focal caculation faster
sd_cladagg <- focal(PC1_agg, w=window, fun=sd)
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
par(mfrow=c(1,2))
plot(sd_clad, col=cl)
plot(sd_cladagg, col=cl) # less accurancy

######12. R_code_snow_cover.r

# R_code_snow_cover.r

setwd("C:/lab/") #set working directory
install.packages("ncdf4") # in order to import ncdf file
library(ncdf4) # recall package
library(raster)
snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc") # import the file: only one layer
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)

setwd("C:/lab/snow/") #creating snow folder is easier to apply lapply function
snow2000 <- raster("snow2000r.tif") #import file as 1 layer
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3)) # see the graphs related to the each year on 1 map. Multitemporal change in snow cover
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

# Faster function to import many data: lapply
# first define the list of files to aggregate together 
rlist <- list.files(pattern="snow20") #all the files have in common "snow20"
# lapply import all of the raster
import <- lapply(rlist, raster)
#  create one raster from several RasterLayers 
snow.multitemp <- stack(import)# MULTITEMPORAL ANALYSISSSSSS, import the RasterStack
plot(snow.multitemp, col=cl)

# let's make a prediction
## prediction

#require(raster)
#require(rgdal)

# define the extent
#ext <- c(-180, 180, -90, 90)
#extension <- crop(snow.multitemp, ext)
    
# make a time variable (to be used in regression)
#time <- 1:nlayers(snow.multitemp)

# run the regression
#fun <- function(x) {if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }} 
#predicted.snow.2025 <- calc(extension, fun) # time consuming: make a pause!
#predicted.snow.2025.norm <- predicted.snow.2025*255/53.90828

source("prediction.r") #read r code

#setwd("C:/lab/snow/") and library(raster)

load("snow_.RData")
prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl)
# how to export the output
writeRaster(prediction, "final.tif")
#how to make the pdf of the graph, ex. stack of all of the outputs
final.stack <- stack(snow.multitemp, prediction)
plot(final.stack, col=cl)
# export the graph in PDF
pdf("my_final_exciting_graph.pdf")
plot(final.stack, col=cl)
dev.off()
# export in png
png("my_final_exciting_graph.png")
plot(final.stack, col=cl)
dev.off()

######13. R_code_monitoring_air_pollution_no2.r

# R_code_monitoring_air_pollution_no2.r
library(raster) # recalling the package

setwd("C:/lab/no2/")
# create RasterStack
rlist <- list.files(pattern="EN") #list the the files (from the folder) that have in common part of the title
import <- lapply(rlist, raster) # import all the files as raster
EN <- stack(import) # stack them together in one raster object
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
boxplot(EN,outline=F) # not draw the outliner values
boxplot(EN,outline=F, horizontal=T) # horizonal boxplot
boxplot(EN,outline=F, horizontal=T, axes=T) #displaying the axis

plot(EN$EN_0001, EN$EN_0013)
abline(0,1,col="red") # to see if the values are under the line, this means a decrease in NO2 concentration!! 45 degree line dividing x and y plan in equal parts 

# fast version of import and plot of many data for lazy people!
rlist <- list.files(pattern="snow")
import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1) # most of the value under the curve
plot(snow.multitemp$snow2000r, snow.multitemp$snow2020r) # better change in time
abline(0,1,col="red")

######14. R_code_crop.r

# change the extention of a raster image by using zoom or crop functions

setwd("C:/lab/")
library(ncdf4)
library(raster)

snow <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
ext <- c(0, 20, 35, 50) # set the extention desired
zoom(snow, ext=ext) # ext means extent
snowitaly <- crop(snow, ext) # write directly the object name for the extent
zoom(snow, ext=drawExtent()) # draw directly the extent from the graph

######15. R_code_interpolation.r

# R_code_interpolation.r
# interpolation field data

setwd("C:/lab/")
library(spatstat) #analysing Spatial Point Patterns
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T) # separator is ; and header of title is true
head(inp) # have a look to the first rows of the dataset
attach(inp) #attach the dataset, making use of the variables

#estimate canopy cover

plot(X,Y) # Y west coordinates
summary(inp) #let's see the minumum and maximum of X and Y, in order to give an extent to spatstat
#convert the data to a point pattern object
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000)) # range of min and max of X and Y, asign the coordinates to spatstat
# giving information about the variable:
#names(inp) #names of the variables
marks(inppp) <- Canopy.cov # add variable to inpp object
canopy <- Smooth(inppp)# visualize the data were they are not been measured. Smooth canopy.cov as interpolation of the surface estimatating where the dat have not been recorded
plot(canopy)# plot density of canopy
points(inppp, col="green")

#lichens for detecting the quality of air
marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)

par(mfrow=c(1,2)) # is there any correlation between canopy and lichens density?... no congruence with canopy cover and lichens
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)

par(mfrow=c(1,3))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)
plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2) # to see the negative correlation

# Psammofile

inp.psam <- read.table("dati_psammofile.csv", sep=";", head=T)
attach(inp.psam)
head(inp.psam)
 
plot(E,N) #clumped distribution
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150)) # coordinates
marks(inp.psam.ppp) <- C_org # ecological data
C <- Smooth(inp.psam.ppp) # warning message:lower amount of data (or no point) for some part because of clumped set
plot(C)
points(inp.psam.ppp) # solution: mean value for each clumped zone or select the zone of the graph and zoom on the top of them: separation of the main graphs in several graphs
 
######16. R_code_sdm.r

# Model species distribution
# Probability of distribution of the species

install.packages("sdm") # function to installing the package
#dataset in library(sdm)

library(sdm) # recall package
library(raster) # for environmental variables
library(rgdal) #input vector layers(for species)
file <- system.file("external/species.shp", package="sdm") # system.file is the function to import the file stored in sdm package
species <- shapefile(file) # to read the .shp file
plot(species)
plot(species[species$Occurrence == 1,],col='blue',pch=16) # condition with [] , == equal to , comma for ending the condition
points(species[species$Occurrence == 0,],col='red',pch=16) # add the occurrence to the plot as point symbol

# use of ecological variables
path <- system.file("external", package="sdm")
lst <- list.files(path=path,pattern='asc$',full.names = T) # to list the predictors' files
preds <- stack(lst) # import the file in one raster with multiple layers, the predictors
cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)

# look at the correlation with predictors and species
plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16) # species preference: low elevation
plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16) #high temperature
plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16) #medium precipitation
plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16) #medium vegetation

## Prediction model
#set the data
d <- sdmData(train=species, predictors=preds)#assign the data
d

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm")# y axis is occurence, x the predictors: y= a +bx1 + cx2 +..
# logistic curve is better (asintot) but is possible to use linear model: methods="...." is the type of model
p1 <- predict(m1, newdata=preds) # let's predict the distribution of the species. Predicting sdm using the method set m1 using the object preds

###probability of distribution of species in space###
plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)

s1 <- stack(preds,p1)
plot(s1, col=cl) # to see final relationship with predictors and final prediction!!!!!

## R_code_MyProject_exam.r

# R code for exam project

#library(rgdal)
#library(gdalUtils)
library(raster)
library(RStoolbox)

#set the working directory
setwd("C:/lab/exam/")

# Sentinel-2 L1C. Date 16-10-2019. Localisation: Catalonia, datum=WGS84, coordinates in UTM
# optain GeoTIFF bands using rgdal and gdalUtils packages

#gdal_translate("T30SWG_20191016T110041_B02.jp2", "B02.tif")
#gdal_translate("T30SWG_20191016T110041_B03.jp2", "B03.tif")
#gdal_translate("T30SWG_20191016T110041_B04.jp2", "B04.tif")
#gdal_translate("T30SWG_20191016T110041_B08.jp2", "B08.tif")
#gdal_translate("T30SWG_20191016T110041_B011.jp2", "B011.tif")

#octB02 <- raster("B02.tif")
#octB03 <- raster("B03.tif")
#octB04 <- raster("B04.tif")
#octB08 <- raster("B08.tif")
#octB11 <- raster("B11.tif")

## or elaborate them directly in SNAP 

## next ##
# Creating a RasterStack

# first of all, resampling B11:
b11 <- raster("subset_0_of_T30SWG_20191016T110041_B11.tif") # move in C:/lab/ to not be in conflict with the new file
b11_dis <- disaggregate(b11, fact=2) # B SWIR2 from resolution 20m to 10m as the other bands
writeRaster(b11_dis, "sub_11.tif") # into C:/lab/exam/

# rasterStack of bands 2,3,4,8,11
rlist <- list.files(pattern="sub") # or pattern="B" using gdalUtils
import <- lapply(rlist, raster)
dry_land_2019 <- stack(import)
plotRGB(dry_land_2019, 4,3,2, stretch="lin")

# vegetation indicies: NDVI and NDWI2
# require RStoolbox

ndvi <- spectralIndices(dry_land_2019 , red = "subset_2_of_T30SWG_20191016T110041_B04", nir = "subset_3_of_T30SWG_20191016T110041_B08", indices = "NDVI")
ndwi2 <- spectralIndices(dry_land_2019 , nir = "subset_3_of_T30SWG_20191016T110041_B08", swir2 = "sub_11", indices = "NDWI2")
ndvi_bit <- stretch(ndvi, minv=0, maxv=255) # to compare ndvi and ndwi2 with the same range of pixel values
ndwi2_bit <- stretch(ndwi2, minv=0, maxv=255)

cl <-  colorRampPalette(c("black", "green", "red"))(100)
clb <-  colorRampPalette(c("black", "gold", "blue"))(100)

# displaying on one map the indicies
par(mfrow=c(2,1))
plot(ndvi_bit, col=cl, main="NDVI")
plot(ndwi2_bit, col=clb, main="NDWI2")

# albedo: Black Sky Albedo, satellite MCD43A3, 500m resolution, bands: 1(red),2(NIR),3(blue),4(green), Date: 16-10-2019, coordinates in decimal degrees lon -2.941667, -2.083333, lat 36.95417, 37.85417 
b1 <- raster("MCD43A3.006_Albedo_BSA_Band1_doy2019289_aid0001.tif")
b2 <- raster("MCD43A3.006_Albedo_BSA_Band2_doy2019289_aid0001.tif")
b3 <- raster("MCD43A3.006_Albedo_BSA_Band3_doy2019289_aid0001.tif")
b4 <- raster("MCD43A3.006_Albedo_BSA_Band4_doy2019289_aid0001.tif")

albedo <- stack(b1, b2, b3, b4)
plot(albedo)

# reducing the dimentions for having a rappresentation of all of the bands together
pairs(albedo)# have a fast look at the correlation among bands
albedoPCA2019 <- rasterPCA(albedo)
summary(albedoPCA$model)# PC1 describes 94%
plot(albedoPCA$map)

# checking and displaying on a map NDVI and Albedo
par(mfrow=c(1,2))
plot(albedoPCA2019$map$PC1)
plot(ndvi_bit)

# Applying Linear Regression for albedo and NDVI

# align coordinates systems, extent and ncell
pr <- projectRaster(ndvi, albedoPCA2019$map$PC1)
rm <- stack(pr, albedoPCA2019$map$PC1)
names(rm) <- c("NDVI", "AlbedoPCA")
cla <-  colorRampPalette(c("white", "red", "blue"))(100)
plot(rm, col=cla) #displaying the result and comparing ndvi and albedo

ndvi2019_ex <- extract(rm$NDVI, c(1:1000, 1000, 1000))
albedo2019_ex <- extract(rm$AlbedoPCA, c(1:1000, 1000, 1000))
model <- lm(ndvi2019_ex ~ albedo2019_ex) # y is ndvi, dipendent variable; x is albedo, indipendent variable
plot(albedo2019_ex, ndvi2019_ex)
abline(model, col="red")
summary(model) #R-squared:  0.4976, p-value: < 2.2e-16

# Applying Linear Regression for albedo and NDWI2: as a proxy of soil moisture content
prw <- projectRaster(ndwi2, albedoPCA2019$map$PC1)
rmw <- stack(prw, albedoPCA2019$map$PC1)
names(rm) <- c("NDWI2", "AlbedoPCA")
plot(rmw)

ndwi22019_ex <- extract(rmw$NDWI2, c(1:1000, 1000, 1000))
albedo2019_ex <- extract(rmw$AlbedoPCA, c(1:1000, 1000, 1000))
modelw <- lm( ndwi22019_ex ~ albedo2019_ex)
plot(albedo2019_ex, ndwi22019_ex)
abline(modelw, col="green")
summary(modelw) # R-squared:  0.3429,  p-value: < 2.2e-16

par(mfrow=c(1,2))
plot(albedo2019_ex, ndvi2019_ex, xlab="Albedo 2019", ylab="NDVI 2019", main="NDVI ~ Albedo" )
abline(model, col="red")
plot(albedo2019_ex, ndwi22019_ex, xlab="Albedo 2019", ylab="NDWI2 2019", main="NDWI2 ~ Albedo")
abline(modelw, col="green")

##The Terra Moderate Resolution Imaging Spectroradiometer (MODIS) Vegetation Indices (MOD13Q1) Version 6, 16 days at 250 meter (m) spatial resolution as a Level 3 product
# SDS name: 250m 16 days NDVI
ndvi_2019 <- raster("MOD13Q1_NDVI_2019.tif")
ndvi_2015 <- raster("MOD13Q1_NDVI_2015.tif")
ndvi_2010 <- raster("MOD13Q1_NDVI_2010.tif")

#ndvi layers need to have the same extent of albedo2019 to compare them with albedo layers
extent <- c(-2.941667, -2.083333, 36.95417, 37.85417)
ndvi_2019ext <- crop(ndvi_2019, extent)
ndvi_2015ext <- crop(ndvi_2015, extent)
ndvi_2010ext <- crop(ndvi_2010, extent)

#Difference in biomass 
dif_ndvi <- ndvi_2019ext- ndvi_2010ext
clc <- colorRampPalette(c('red','green','yellow')) (100)
plot(dif_ndvi, col= clc)
hist(dif_ndvi)

# rise or decrease in biomass
plot(ndvi_2010ext, ndvi_2019ext, xlab="NDVI 2010", ylab="NDVI 2019", main="Trend of NDVI")
abline(0,1,col="red") # 45° line that divides the cartesian plane in 2 equal parts

NDVI <- stack(ndvi_2010ext, ndvi_2015ext, ndvi_2019ext)
boxplot(NDVI,outline=F, horizontal=T, axes=T, names=c("NDVI 2010", "NDVI 2015", "NDVI 2019"), main="Boxplot NDVI") 

## The Moderate Resolution Imaging Spectroradiometer (MODIS) MCD43A3 Version 6 Albedo Model dataset, 16 days of Terra and Aqua MODIS data at 500 meter (m) resolution
# Date 16 October 2015 and 2010

# 2015
alb_2015b1 <- raster("MCD43A3.006_Albedo_BSA_Band1_doy2015.tif")
alb_2015b2 <- raster("MCD43A3.006_Albedo_BSA_Band2_doy2015.tif")
alb_2015b3 <- raster("MCD43A3.006_Albedo_BSA_Band3_doy2015.tif")
alb_2015b4 <- raster("MCD43A3.006_Albedo_BSA_Band4_doy2015.tif")

# albedo layers need to have the same extent of albedo2019
extent <- c(-2.941667, -2.083333, 36.95417, 37.85417)
alb_2015b1ex <- crop(alb_2015b1, extent)
alb_2015b2ex <- crop(alb_2015b2, extent)
alb_2015b3ex <- crop(alb_2015b3, extent)
alb_2015b4ex <- crop(alb_2015b4, extent)
albedo_2015 <- stack(alb_2015b1ex, alb_2015b2ex, alb_2015b3ex, alb_2015b4ex)
albedo_2015PCA <- rasterPCA(albedo_2015)
summary(albedoPCA2015$model)

# 2010
alb_2010b1 <- raster("MCD43A3.006_Albedo_BSA_Band1_doy2010.tif")
alb_2010b2 <- raster("MCD43A3.006_Albedo_BSA_Band2_doy2010.tif")
alb_2010b3 <- raster("MCD43A3.006_Albedo_BSA_Band3_doy2010.tif")
alb_2010b4 <- raster("MCD43A3.006_Albedo_BSA_Band4_doy2010.tif")

extent <- c(-2.941667, -2.083333, 36.95417, 37.85417)
alb_2010b1ex <- crop(alb_2010b1, extent)
alb_2010b2ex <- crop(alb_2010b2, extent)
alb_2010b3ex <- crop(alb_2010b3, extent)
alb_2010b4ex <- crop(alb_2010b4, extent)
albedo_2010 <- stack(alb_2010b1ex, alb_2010b2ex, alb_2010b3ex, alb_2010b4ex)
albedo_2010PCA <- rasterPCA(albedo_2010)
summary(albedoPCA2010$model)

# in my opinion faster passage
writeRaster(albedoPCA2019$map$PC1, "albedoPCA2019.tif")
writeRaster(albedo_2010PCA$map$PC1, "albedoPCA2010.tif")
writeRaster(albedo_2015PCA$map$PC1, "albedoPCA2015.tif")

albedoPCA2010 <- raster("albedoPCA2010.tif")
albedoPCA2015 <- raster("albedoPCA2015.tif")
albedoPCA2019 <- raster("albedoPCA2019.tif")

# Difference in albedo 
albedoPCA2019b <- stretch(albedoPCA2019, minv=0, maxv=255)
albedoPCA2015b <- stretch(albedoPCA2015, minv=0, maxv=255)
albedoPCA2010b <- stretch(albedoPCA2010, minv=0, maxv=255)

dif_albedo <- albedoPCA2019b - albedoPCA2010b
plot(dif_albedo)
plot(albedoPCA2010b, albedoPCA2019b, xlab="Albedo 2010", ylab="Albedo 2019", main="Trend of Albedo")
abline(0,1,col="red")

albedobox <- stack(albedoPCA2010b, albedoPCA2015b, albedoPCA2019b)
boxplot(albedobox, outline=F, horizontal=T, axes=T, names=c("Albedo 2010", "Albedo 2015", "Albedo 2019"), main="Boxplot Albedo")



# let's check if with higher amount of data multitemporal analysis of NDVI changes
# satellite data day by day for the month of October

library(MODISTools)
VI <- mt_subset(product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                lat = 37,
                lon = -2,
                start = "2019-10-01",
                end = "2019-10-31",
                km_lr = 100,
                km_ab = 100,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

V2 <- mt_subset(product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                lat = 37,
                lon = -2,
                start = "2015-10-01",
                end = "2015-10-31",
                km_lr = 100,
                km_ab = 100,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

V3 <- mt_subset(product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                lat = 37,
                lon = -2,
                start = "2010-10-01",
                end = "2010-10-31",
                km_lr = 100,
                km_ab = 100,
                site_name = "testsite",
                internal = TRUE,
                progress = FALSE)

VI_r <- mt_to_raster(df = VI)
V2_r <- mt_to_raster(df = V2)
V3_r <- mt_to_raster(df = V3)

multit_NDVI <- stack(VI_r, V2_r, V3_r)
# multitemporal NDVI needs to have the same coordinates system and extent of previous NDVI data to compare the boxplot results
multit_NDVIex <- projectRaster(multit_NDVI,alb_2015b1ex)
boxplot(multit_NDVIex,outline=F, horizontal=T, axes=T, names=c("NDVI 2010", "NDVI 2015", "NDVI 2019"), main="Boxplot multitemporal NDVI")





