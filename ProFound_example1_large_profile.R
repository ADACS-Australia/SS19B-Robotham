library(ProFound)
library(profvis)
image_big=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
profound=profoundProFoundADACSInPlace(image_big, skycut=1.5, magzero=30, verbose=TRUE, plot=FALSE)
for(i in 1:4){image_big=rbind(image_big,image_big)}
for(i in 1:4){image_big=cbind(image_big,image_big)}

gc()
profvis({
  profound=profoundProFoundADACSInPlace(image_big, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE,doclip=FALSE)
  gc()
  profound=profoundProFound(image_big, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE,doclip=FALSE)
})
