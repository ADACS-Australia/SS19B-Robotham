# SS19B-Robotham

Use R version 3.5.1 as a minimum

First time install, build, and run
----------------------------------

## On MacOS assuming R (>=3.5.0) is already installed ##

*  Build straight from repo using: ./run_Build_from_scratch.sh
*  Run the example: ./run_example1_non_interactive.sh
*  To run the example interactively: Cut and paste the following into an R session (started with command R, ended with CTRL D)

## On OzSTAR assuming R 3.5.1 is available as a module ##

*  `module load r/3.5.1`
*  From R run `install.packages("devtools")` and say you want a private library when it asks
*  From R run `install.packages("BiocManager")`
*  From R run `BiocManager::install("EBImage")`
*  `./run_Build_from_scratch.sh`
*  `./run_example1_non_interactive.sh`

## On Windows using R version 3.6.x: ##

From R session (or from RStudio console) -

* `install.packages("devtools")`
* `install.packages("provis")`
* `install.packages("akima")`
* `install.packages("imager")`
* `install.packages("BiocManager")`
* `BiocManager::install("EBIMage")`
* `library(devtools)`
* `install_github("asgr/magicaxis")`
* `install_github("asgr/Profound")`

## To run an example: ##

```
library(ProFound)
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
```
             
Building from cloned repo
-------------------------

Something like this (not sure of the minimum steps):

On MacOS assuming RStudio is already installed

1.  Launch RStudio
2.  Use File/Open Project... and navigate to the "ProFound" sub directory of the repo clone
3.  If required, "Build/Configure Build Tools..." and under "Build Tools" select "Package" for "Project build tools" ...
4.  Select "Use devtools package if available" and OK
5.  To build select "Build/Clean and Rebuild"

