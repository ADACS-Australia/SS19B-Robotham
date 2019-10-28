#!/bin/bash
cat > tmp_remove_installed_package.R << !
remove.packages("bitmatrix","/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
!
rm -rf bitmatrix
Rscript tmp_remove_installed_package.R
