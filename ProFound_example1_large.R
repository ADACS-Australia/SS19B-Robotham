library(ProFound)
image_big=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
for(i in 1:4){image_big=rbind(image_big,image_big)}
for(i in 1:4){image_big=cbind(image_big,image_big)}
profound=profoundProFound(image_big, skycut=1.5, magzero=30, verbose=TRUE, plot=FALSE)
