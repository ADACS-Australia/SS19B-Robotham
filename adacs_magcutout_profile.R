library(ProFound)
library(pryr)
library(profvis)
image_big=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
for(i in 1:4){image_big=rbind(image_big,image_big)}
for(i in 1:4){image_big=cbind(image_big,image_big)}
box=c(101,101)
profvis({
n=10000L
ooselect = matrix(0.0,box[1],box[2])
randoms1 = runif(n,min=1,max=dim(image_big)[1])
randoms2 = runif(n,min=1,max=dim(image_big)[2])
for (i in 1:n) {
  loc=c(floor(randoms1[i])+0.5,floor(randoms2[i])+0.5)
  select=magcutout(image_big, loc=loc, box=box, shiftloc=FALSE, paddim=TRUE)$image
}
gc()
for (i in 1:n) {
  loc=c(floor(randoms1[i])+0.5,floor(randoms2[i])+0.5)
  magcutoutADACSInPlace(image_big, loc=loc, box=box, ooselect)
}
})

