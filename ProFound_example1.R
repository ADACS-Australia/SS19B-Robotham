library(ProFound)
library(stringi)
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
#image=readFITS("../testdata/TAO_mockimages/image.1154.fits")$imDat
#image=readFITS("../testdata/TAO_mockimages/image.2906.fits")$imDat
image[1,1] = NaN
print(paste("input size", dim(image)[1], dim(image)[2]))
profound=profoundProFoundADACSInPlace(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE, type="bicubic", skytype="median", skyRMStype = "quanlo")
#profound=profoundProFoundADACSInPlace(image, verbose=TRUE, box=c(39, 39))
#profound=profoundProFound(image, verbose=TRUE, box=c(39, 39))
#profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)