library(ProFound)
library(stringi)
doclip = TRUE
initialiseGlobals(doclip)
args = commandArgs(TRUE)
box_size = as.integer(args[1])
image_resize_steps = as.integer(args[2])

image = readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
for (i in seq_along(seq_len(image_resize_steps))) {
	if (i %% 2 == 0) {
		image = cbind(image, image)
	}
	else {
		image = rbind(image, image)
	}
}
# prepare
# Ensure dimensions of box are odd numbers
box = c(box_size,box_size)
if (box[1]%%2 == 0) {
  box = c(box[1]+1, box[2])
}
if (box[2]%%2 == 0) {
  box = c(box[1], box[2]+1)
}
# Ensure dimensions of boxadd are even numbers (so that box + boxadd stays odd)
scratch <- list(scratchN1=matrix(0.0,box[1],box[2]), scratchN2=matrix(0.0,box[1],box[2]),scratchI1=matrix(FALSE,box[1],box[2]), scratchI2=matrix(FALSE,box[1],box[2]), scratchSKY=matrix(0.0,dim(image)[1],dim(image)[2]),scratchSKYRMS=matrix(0.0,dim(image)[1],dim(image)[2]))

#
result = profoundMakeSkyGridADACSInPlace(image, box=c(box_size, box_size),scratch=scratch)

image_width = dim(image)[1]
image_height = dim(image)[2]
image_size = image_width * image_height
image_kbytes = image_size * 8 / 1024
msg = paste(image_width, image_height, image_size, image_kbytes, box_size, sep=', ')
print(msg)
