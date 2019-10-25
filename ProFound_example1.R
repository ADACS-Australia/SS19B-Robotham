library(ProFound)
library(stringi)
library(bitmatrix)
b<-new(BitMatrix, 10, 10)
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
#image=readFITS("../testdata/TAO_mockimages/image.1154.fits")$imDat
#image=readFITS("../testdata/TAO_mockimages/image.2906.fits")$imDat
#image[1,1] = NaN
bmask=new(BitMatrix, dim(image)[1],dim(image)[2])
bmask$settrue(1,1)
profound=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE, type="bicubic", skytype="median", skyRMStype = "quanlo", bmask=bmask)
#print(paste("Nseg=", profound$Nseg))
#profound=adacs_ProFound(image, verbose=TRUE, box=c(20, 20), bmask=bmask)
#profound=profoundProFound(image, verbose=TRUE, box=c(39, 39))
#profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=FALSE, type="bicubic", skytype="median", skyRMStype = "quanlo")
