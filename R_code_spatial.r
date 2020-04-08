#see how to use coordinates in space related to the variables.....using the function: coordinates 

library(sp)

data(meuse) #to recall the data
head(meuse)

#to use the coordinates in group
coordinates(meuse)=~x+y #thinking spatialy
plot(meuse) #the whole variables in meuse
spplot(meuse, "zinc") #for a spatial variable

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
covid<-read.table("covid_agg.csv", head=T)
head(covid)

attach(covid)
plot(country,cases)
plot(country,cases, las=0) #paralel labels
plot(country,cases, las=1) #horizontal labels
plot(country,cases, las=2) #perperdincolar labels
plot(country,cases, las=3) #vertical labels
plot(country,cases, las=3, cex.axis=0.5) #cex for the size of labels axis
install.packages("ggplot2")

#use the function
library(ggplot2)

#to open the saved work
setwd("C:/lab")
load("R_code_spatial")
ls() #list of objects i've done
# I'll obtain "covid" "meuse"

library(ggplot2)
#data set to use
data(mpg)
head(mpg)
#function we are going to use ggplot.. components: data, aes, geometry
ggplot(mpg,aes(x=displ,y=hwy)) + geom_point()
ggplot(mpg,aes(x=displ,y=hwy)) + geom_line()
ggplot(mpg,aes(x=displ,y=hwy)) + geom_polygon()

head(covid)
ggplot(covid,aes(x=lon,y=lat, size=cases)) + geom_point() #lat and lon dipend on size

