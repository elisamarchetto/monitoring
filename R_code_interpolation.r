# R_code_interpolation.r
# interpolation field data

setwd("C:/lab/")
library(spatstat)
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T) # separator is ; and header of title is true
head(inp)
attach(inp) #attach the dataset

#estimate canopy cover

plot(X,Y) # Y west coordinates
summary(inp) #let's see the minumum and maximum of X and Y, in order to give an extent to spatstat
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000)) # range of min and max of X and Y, assign the coordinates to spatstat
# giving information about the variable: lables
#names(inp) #names of the variables
marks(inppp) <- Canopy.cov
canopy <- Smooth(inppp)# visualize the data were they are not been measured in pixel
##list validation distance of value from the line that record the values: means measured the amount of error
plot(canopy)#density
points(inppp, col="green")

#lichens for detecting the quality of air
marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp) # no congruence with canopy cover and lichens

par(mfrow=c(1,2))
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
points(inp.psam.ppp) # solution: mean value for each clumped zone or select the zone of the graph and zoom on the top of them: separation of the main graph in several graphs
 


 
