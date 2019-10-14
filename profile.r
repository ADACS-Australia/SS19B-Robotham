library(ProFound)
library(stringi)
library(bitmatrix)

args = commandArgs(TRUE)
box_size = as.integer(args[1])
image_resize_steps = as.integer(args[2])
what = args[3]

# Prepare image of requested size
image = readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
for (i in seq_along(seq_len(image_resize_steps))) {
	if (i %% 2 == 0) {
		image = cbind(image, image)
	}
	else {
		image = rbind(image, image)
	}
}
image[1,1] = NaN

# Go, go, go!
if (what == 'profound-original') {
	result = profoundProFound(image, box=c(box_size, box_size))
} else if (what == 'profound-adacs') {
	result = profoundProFoundADACSInPlace(image, box=c(box_size, box_size))
} else if (what == 'skygrid-original') {
	result = profoundMakeSkyGrid(image, box=c(box_size, box_size))
} else if (what == 'skygrid-adacs') {
	# Ensure dimensions of box are odd numbers
	box = c(box_size,box_size)
	if (box[1]%%2 == 0) {
	  box = c(box[1]+1, box[2])
	}
	if (box[2]%%2 == 0) {
	  box = c(box[1], box[2]+1)
	}
	scratch <- list(
			scratchN1=matrix(0.0,box[1],box[2]),
			scratchN2=matrix(0.0,box[1],box[2]),
			scratchI1=matrix(FALSE,box[1],box[2]),
			scratchI2=matrix(FALSE,box[1],box[2]),
			scratchSKY=matrix(0.0,dim(image)[1],dim(image)[2]),
			scratchSKYRMS=matrix(0.0,dim(image)[1],dim(image)[2])
	)
	initialiseGlobals(TRUE)
	result = profoundMakeSkyGridADACSInPlace(image, box=box, scratch=scratch)
} else {
        print(paste("Option not recognised: ",what))
}

# Report and good bye
image_width = dim(image)[1]
image_height = dim(image)[2]
image_size = image_width * image_height
image_kbytes = image_size * 8 / 1024
msg = paste(image_width, image_height, image_size, image_kbytes, box_size, sep=', ')
print(msg)
