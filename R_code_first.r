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
