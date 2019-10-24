library(ProFound)
library(stringi)
library(RUnit)

initialiseGlobals(TRUE)

# keywords to enum check (The C/C++ uses integers to represent enums)
print("############# START OF TESTS #########")

# checking enumForKeyword vs Globals

checkEquals(doRMNA, FALSE, "doRMNA")
checkIdentical(MEDIAN, enumForKeyword("median"), "median")
checkIdentical(MEAN, enumForKeyword("mean"), "mean")
checkIdentical(MODE, enumForKeyword("mode"), "mode")
checkIdentical(RMEDIAN, enumForKeyword("rmedian"), "R version of 'median'")
checkIdentical(RMEAN, enumForKeyword("rmean"), "R version of 'mean'")
checkIdentical(RMODE, enumForKeyword("rmode"), "R version of 'mode'")
checkIdentical(BOTH, enumForKeyword("quanboth"), "quanboth")
checkIdentical(LO, enumForKeyword("quanlo"), "quanlo")
checkIdentical(HI, enumForKeyword("quanhi"), "quanhi")
checkIdentical(SD, enumForKeyword("sd"), "sd")
checkIdentical(RBOTH, enumForKeyword("rquanboth"), "R version of quanboth")
checkIdentical(RLO, enumForKeyword("rquanlo"), "R version of quanlo")
checkIdentical(RHI, enumForKeyword("rquanhi"), "R version of quanhi")
checkIdentical(RSD, enumForKeyword("rsd"), "R version of sd")
checkIdentical(AUTO, enumForKeyword("auto"), "auto")
checkIdentical(SET, enumForKeyword("set"), "set")
checkIdentical(CLASSIC_BILINEAR, enumForKeyword("bilinear"), "Classic Bilinear")
checkIdentical(AKIMA_BICUBIC, enumForKeyword("bicubic"), "Akima Bicubic")

# check compareNA
checkIdentical(c(TRUE, TRUE, FALSE), compareNA(c(10, 10, 10), c(10, 10, 20)))
checkIdentical(is.na(NULL), compareNA(NULL, NULL))


print("############# END OF TESTS ##########")
