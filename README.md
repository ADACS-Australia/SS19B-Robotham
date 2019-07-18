# SS19B-Robotham
On MacOS assuming R is already installed

1.  Build straight from repo using: ./run_Build_from_scratch.sh
2.  Run the example: ./run_example1_non_interactive.sh
3.  To run the example interactively: Cut and paste the following into an R session (started with command R, ended with CTRL D)
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
