# Model species distribution
# Probability of distribution of the species

install.packages("sdm")
#dataset in library(sdm)

library(sdm)
library(raster) # for environmental variables
library(rgdal) #input vector layers(for species)
file <- system.file("external/species.shp", package="sdm") # system.file is the function to import the file into de sdm package
species <- shapefile(file)
plot(species)
plot(species[species$Occurrence == 1,],col='blue',pch=16) # condition with [] , == equal to , comma for end the condition
points(species[species$Occurrence == 0,],col='red',pch=16) # points to add the occurrence to the plot

# use of ecological variables
path <- system.file("external", package="sdm")
lst <- list.files(path=path,pattern='asc$',full.names = T) 
preds <- stack(lst)
cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)

# look at the correlation with predictors and species
plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16) #low elevation
plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16) #high temperature
plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16) #medium precipitation
plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16) #medium vegetation

## Prediction model
#set the data
d <- sdmData(train=species, predictors=preds)
d

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm")# y axis is occurence, x the predictors: y= a +bx1 + cx2 +..
# logistic curve is better (asintot) but is possible to use linear model: methods="...." is the type of model
p1 <- predict(m1, newdata=preds)

###probability of distribution of species in space###
plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)

s1 <- stack(preds,p1)
plot(s1, col=cl) # to see final relationship with predictors and final prediction!!!!!


