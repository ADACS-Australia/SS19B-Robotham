library(ProFound)
library(RUnit)

image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
print("############# START OF TESTS: adacs_test_1.R ##########")
result=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE, type="bicubic", skytype="median", skyRMStype = "quanlo")
checkEquals(56, result$Nseg, "call with bicubic, median, quanlo")
result=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE, type="bicubic", skytype="adacs_median", skyRMStype = "adacs_quanlo")
checkEquals(56, result$Nseg, "call with bicubic, adacs_median, adacs_quanlo")
result=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE, type="bilinear", skytype="adacs_median", skyRMStype = "adacs_quanlo")
checkEquals(66, result$Nseg, "call with bilinear, adacs_median, adacs_quanlo")
result=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE, type="bilinear", skytype="mean", skyRMStype = "quanlo")
checkEquals(62, result$Nseg, "call with bilinear, mean, quanlo")
result=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE, type="bilinear", skytype="adacs_mean", skyRMStype = "adacs_quanlo")
checkEquals(62, result$Nseg, "call with bilinear, adacs_mean, adacs_quanlo")

image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
result=adacs_ProFound(image, skycut=1.5, magzero=30, verbose=FALSE, plot=FALSE, type="bicubic", skytype="median", skyRMStype = "quanlo")
checkEquals(56, result$Nseg, "call with bicubic, median, quanlo")

print("############# END OF TESTS: adacs_test_1.R ##########")
