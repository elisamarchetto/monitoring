#multipanel in R: monitoring ecosystems

install.packages("GGally") #'GGally' extends 'ggplot2' by adding several functions to reduce the complexity of combining geometric objects
library(sp)
#require(sp) is the same of library(sp)
data(meuse)
attach(meuse) #make use of data set of meuse

meuse #to see all the variables
plot(cadmium)

#make che pairs vairables in plot
pairs(meuse) #intersections of variables in a kind of matrix. pairs is a function stored in GGally (in this case)

#to correlate only some variables in the plot. It needs to be indicated the dataset that is going to be used
pairs(~ cadmium+copper+lead+zinc,data=meuse)
#or
pairs(meuse[,3:6]) # [,3:6] from 3(cadmium) to 6(zinc)
#16 and 15 are the same function!

pairs(meuse[,3:6],pch=19)

library(GGally)
ggpairs(meuse[,3:6])
