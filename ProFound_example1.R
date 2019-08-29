library(ProFound)
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFoundADACSInPlace(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
