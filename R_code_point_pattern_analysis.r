#Point pattern analysis: density map

#install library packages for density analysis
install.packages("spatstat") # toolbox for analysing Spatial Point Patterns
library(spatstat)
attach(covid)
head(covid)

covids <- ppp(lon, lat, c(-180, 180), c(-90, 90)) #ppp means is panel point pattern, it creates a point pattern dataset in the two-dimensional plane. for the range of values of lat and long. 
#to build the density map
d<- density(covids) 
plot(d) 
#to see the point of covids object on the map
points(covids)

#to open last work
setwd("C:/lab/")
load("point_pattern_analysis.RData")
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

#export your R in Pdf
pdf("covid_density.pdf")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()
#or export in png
png("covid_density.png")
clr <- colorRampPalette(c("light green", "yellow","orange","violet")) (100)
plot(d, col=clr, main="Densities of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()

