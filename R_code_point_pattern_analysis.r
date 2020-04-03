#Point pattern analysis: density map

#install library packages for density analysis
install.packages("spatstat")
library(spatstat)
attach(covid)
head(covid)

covids <- ppp(lon, lat, c(-180, 180), c(-90, 90)) #c for the range of lat and long, ppp is panel point pattern
#to build the density map
d<- density(covids) #covids in other to write the object of ppp
plot(d) #to see the map
#to see the point of covids on the map
points(covids)
