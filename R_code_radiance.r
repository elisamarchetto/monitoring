#R_code_radiance.r

library(raster)
toy <- raster(ncol=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2) # for creating the matrix as a layer
values(toy) <- c(1.13,1.44,1.55,3.4)
plot(toy)
text(toy, digits=2)
toy2bits <- stretch(toy,minv=0,maxv=3) # reshape the matrix with 2^2 pixel (2 bits), so 4 values
storage.mode(toy2bits[]) = "integer" # integer values
plot(toy2bits) #bits: possible combinations of pixels based on NDs
text(toy2bits, digits=2) #lower is the amount of bits lower is the diversity
toy4bits <- stretch(toy,minv=0,maxv=15)
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2)
toy8bits <- stretch(toy,minv=0,maxv=255) #better discrimination from one objects to an other
storage.mode(toy8bits[]) = "integer"
plot(toy8bits)
text(toy8bits, digits=2) 

par(mfrow=c(1,4))

plot(toy)
text(toy, digits=2)
plot(toy2bits)
text(toy2bits, digits=2)
plot(toy4bits)
text(toy4bits, digits=2)
plot(toy8bits)
text(toy8bits, digits=2)



