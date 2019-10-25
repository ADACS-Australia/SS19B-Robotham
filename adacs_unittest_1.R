library(ProFound)
library(RUnit)

initialiseGlobals(TRUE)

# keywords to enum check (The C/C++ uses integers to represent enums)
print("############# START OF TESTS #########")

# checking enumForKeyword vs Globals

checkEquals(doRMNA, FALSE, "doRMNA")
checkIdentical(MEDIAN, enumForKeyword("adacs_median"), "median")
checkIdentical(MEAN, enumForKeyword("adacs_mean"), "mean")
checkIdentical(MODE, enumForKeyword("adacs_mode"), "mode")
checkIdentical(RMEDIAN, enumForKeyword("median"), "R version of 'median'")
checkIdentical(RMEAN, enumForKeyword("mean"), "R version of 'mean'")
checkIdentical(RMODE, enumForKeyword("mode"), "R version of 'mode'")
checkIdentical(BOTH, enumForKeyword("adacs_quanboth"), "quanboth")
checkIdentical(LO, enumForKeyword("adacs_quanlo"), "quanlo")
checkIdentical(HI, enumForKeyword("adacs_quanhi"), "quanhi")
checkIdentical(SD, enumForKeyword("adacs_sd"), "sd")
checkIdentical(RBOTH, enumForKeyword("quanboth"), "R version of quanboth")
checkIdentical(RLO, enumForKeyword("quanlo"), "R version of quanlo")
checkIdentical(RHI, enumForKeyword("quanhi"), "R version of quanhi")
checkIdentical(RSD, enumForKeyword("sd"), "R version of sd")
checkIdentical(AUTO, enumForKeyword("auto"), "auto")
checkIdentical(SET, enumForKeyword("set"), "set")
checkIdentical(CLASSIC_BILINEAR, enumForKeyword("bilinear"), "Classic Bilinear")
checkIdentical(AKIMA_BICUBIC, enumForKeyword("bicubic"), "Akima Bicubic")

# check compareNA
checkIdentical(c(TRUE, TRUE, FALSE), compareNA(c(10, 10, 10), c(10, 10, 20)))
checkIdentical(is.na(NULL), compareNA(NULL, NULL))


print("############# END OF TESTS ##########")
