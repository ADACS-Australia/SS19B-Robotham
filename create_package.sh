#!/bin/bash
cat > tmp_create_package.R << !
library(Rcpp)
Rcpp.package.skeleton("bitmatrix",
    author="Ray Seikel",
    email="rseikel@bigpond.com",
    example_code=FALSE,
    module=TRUE,
    cpp_files=c("bitmatrix.h", "rcpp_module.cpp", "rcpp_module.hpp"),
    attributes=TRUE
    )
!
./remove_installed_package.sh
RScript tmp_create_package.R

rm -f bitmatrix/src/Num.cpp
rm -f bitmatrix/src/stdVector.cpp
cat > bitmatrix/R/zzz.R << !
## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}


## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
loadModule("yada", TRUE)
!

cat > tmp_patch_package.R << !
library(Rcpp)
compileAttributes(pkgdir = 'bitmatrix', verbose = TRUE)
!
RScript tmp_patch_package.R
R CMD INSTALL --build bitmatrix

cat > tmp_test_package.R << "!"
library(bitmatrix)
x=c(1,2,3)
y=c(4,5,6)
w<-new(BitMatrix, 10, 10)
w$settrue(5,5)
w$istrue(5,5)
w$trues()
a<-new(Adacs)
mask=matrix(0L,20,20)
mask[1,1] = 1;
bmask=new(BitMatrix, mask)
bmask$trues()
!
#cp bitmatrix.h.installed /Library/Frameworks/R.framework/Resources/library/bitmatrix/include/bitmatrix.h
Rscript tmp_test_package.R
rm tmp_*
