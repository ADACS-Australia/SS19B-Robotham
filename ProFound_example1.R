library(ProFound)
#image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
image=readFITS("../testdata/TAO_mockimages/image.1154.fits")$imDat
print(paste("input size", dim(image)[1], dim(image)[2]))
profound=profoundProFoundADACSInPlace(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
#profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
