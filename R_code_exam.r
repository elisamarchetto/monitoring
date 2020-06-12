# R code for exam 

##1. R_code_first.r
# test for the first time functions and codes in R software

install.packages("sp") # Classes and methods for spatial data
library(sp)

data(meuse)

# have a look at the structure of dataset:
meuse

#having a look at the first part of dataset
head(meuse)

#correlete two varialables
attach(meuse) # attach values for using plot function
plot(zinc,copper)
plot(zinc,copper, col="green")
plot(zinc,copper, col="green",pch=19) # pch is a command to change the symbol in the graph
plot(zinc,copper, col="green",pch=19,cex=2) # cex is the command for the dimention of pch in this case

##2. R_code_Spatial.r

#see how to use coordinates in space related to the variables.....using the function: coordinates 

library(sp)

data(meuse) #to recall the data
head(meuse)

#coordinates is the functionn to visualize the variables in the space 
coordinates(meuse)=~x+y #thinking spatialy
plot(meuse) # plot all the variables in meuse
spplot(meuse, "zinc") #for a spatial variable (dipens on coordinates function)

#exercise: spatial amount of copper
spplot(meuse, "copper")

#to change the title in the graph using main= "....
spplot(meuse, "copper", main="Copper concentration")

#for duplicating the size of a variable
bubble(meuse, "zinc")
#exercise: change the color of zinc in red
bubble(meuse, "zinc", col="red")

#download covid_aff in lab into C:
#inport data base covid in R
#setting the working directory: lab
setwd("C:/lab")
covid<-read.table("covid_agg.csv", head=T) # the file contains the names of the variables as its first line
head(covid)

attach(covid)
plot(country,cases)
plot(country,cases, las=0) #paralel labels
plot(country,cases, las=1) #horizontal labels
plot(country,cases, las=2) #perperdincolar labels
plot(country,cases, las=3) #vertical labels
plot(country,cases, las=3, cex.axis=0.5) #cex for the size of labels axis
install.packages("ggplot2") # A system for 'declaratively' creating graphics, based on "The Grammar of Graphics"

#use the function
library(ggplot2)

#to open the saved work
setwd("C:/lab")
load("R_code_spatial")
ls() #list of objects i've done, es. the functions of my work space
# I'll obtain "covid" "meuse"

library(ggplot2)
#data set to use
data(mpg)
head(mpg)
#function we are going to use ggplot.. components: data, aes, geometry
ggplot(mpg,aes(x=displ,y=hwy)) + geom_point()#aes: aesthetic mapping choosing the variables of the dataset to visualize, geom_ :Specifies the geometric objects that define the graph type
ggplot(mpg,aes(x=displ,y=hwy)) + geom_line()
ggplot(mpg,aes(x=displ,y=hwy)) + geom_polygon()

head(covid)
ggplot(covid,aes(x=lon,y=lat, size=cases)) + geom_point() #lat and lon dipend on size

## 3. R_code_multivar.r

#R code for multivariate analysis

install.packages("vegan") #vegetation and diversity analysis
library(vegan)
setwd("C:/lab/")
biomes<-read.table("biomes.csv", header=T, sep=",") # sep="," in cvs file the variables are separeted by comma, biomes.cvs is the data frame
head(biomes)

## Multivariate analysis
#DEtrended CORrespondance ANAlysis---to reduce the dimenctions (variables)
#the function is decorana

multivar<-decorana(biomes)
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
ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3) #col for each biomes, kind for the type of grafh, lwd for the dimention of the colored line

ordispider(multivar, type, col=1:4, label = T) # type is the column we want to see, label is refered to the name of type column we want to see in the graph








