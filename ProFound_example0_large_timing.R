library(ProFound)
library(profvis)
library(tictoc)
image_big=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat

profound=1
for (i in 1:6) {
  print("double grid size")
  image_big=rbind(image_big,image_big)
  image_big=cbind(image_big,image_big)
  gc()
  tic(paste("Done pass::",i," dim:",dim(image_big)[1],dim(image_big)[2]," input memory:",object.size(image_big), " ",object.size(profound)))
  profound=profoundProFoundADACSInPlace(image_big, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE,doclip=TRUE)
  elapsed = toc()
  print(paste("sleeping 10 seconds. speed per pixel was:",(elapsed$toc-elapsed$tic)/prod(dim(image_big))))
  Sys.sleep(10)
}

