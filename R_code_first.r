install.packages("sp")
library(sp)

data(meuse)

#structure of data set:
meuse

# look of first part of data set
head(meuse)

#correlete two varialables
attach(meuse)
plot(zinc,copper)
plot(zinc,copper, col="green")
plot(zinc,copper, col="green",pch=19)
plot(zinc,copper, col="green",pch=19,cex=2)
