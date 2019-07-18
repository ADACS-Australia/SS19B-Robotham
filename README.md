# SS19B-Robotham

Use R version 3.5.1 as a minimum

First time install, build, and run
----------------------------------

On MacOS assuming R is already installed

1.  Build straight from repo using: ./run_Build_from_scratch.sh
2.  Run the example: ./run_example1_non_interactive.sh
3.  To run the example interactively: Cut and paste the following into an R session (started with command R, ended with CTRL D)
```
library(ProFound)
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
```

On OzSTAR assuming R 3.5.1 is available as a module

1.  module load r/3.5.1
2.  From R run install.packages("devtools") and say you want a private library when it asks
3.  From R run install.packages("BiocManager")
4.  From R run BiocManager::install("EBImage")
5.  ./run_Build_from_scratch.sh
6.  ./run_example1_non_interactive.sh

On Windows using R version 3.6.x:

From R (RStudio) -

1. install.packages("devtools")
2. install.packages("provis")
3. install.packages("akima")
4. install.packages("imager")
5. install.packages("BiocManager")
6. BiocManager::install("EBIMage")
7. library(devtools)
8. install_github("asgr/magicaxis")
9. install_github("asgr/Profound")

To run an example (on windows):

1.library(ProFound)
2. image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
3. profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
             
Building from cloned repo
-------------------------

Something like this (not sure of the minimum steps):

On MacOS assuming RStudio is already installed

1.  Launch RStudio
2.  Use File/Open Project... and navigate to the "ProFound" sub directory of the repo clone
3.  If required, "Build/Configure Build Tools..." and under "Build Tools" select "Package" for "Project build tools" ...
4.  Select "Use devtools package if available" and OK
5.  To build select "Build/Clean and Rebuild"

