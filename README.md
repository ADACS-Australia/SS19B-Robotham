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
             
