library(ProFound)
library(stringi)

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

# Go, go, go!
if (what == 'profound-original') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
	mask[1,1] = 1L
	result = profoundProFound(image, box=c(box_size, box_size), mask=mask)
} else if (what == 'profound-adacs') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
	mask[1,1] = 1L
	result = profoundProFoundADACSInPlace(image, box=c(box_size, box_size), mask=mask)
} else if (what == 'skygrid-original') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
	mask[1,1] = 1L
	result = profoundMakeSkyGrid(image, box=c(box_size, box_size), mask=mask)
} else if (what == 'skygrid-adacs') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
	mask[1,1] = 1L
	# Ensure dimensions of box are odd numbers
	box = c(box_size,box_size)
	if (box[1]%%2 == 0) {
	  box = c(box[1]+1, box[2])
	}
	if (box[2]%%2 == 0) {
	  box = c(box[1], box[2]+1)
	}
	scratch <- list(
			scratchSKY=matrix(0.0,dim(image)[1],dim(image)[2]),
			scratchSKYRMS=matrix(0.0,dim(image)[1],dim(image)[2])
	)
	initialiseGlobals(TRUE)
	result = adacs_MakeSkyGrid(image, box=box, mask=mask, scratch=scratch)
}

# Report and good bye
image_width = dim(image)[1]
image_height = dim(image)[2]
image_size = image_width * image_height
image_kbytes = image_size * 8 / 1024
msg = paste(image_width, image_height, image_size, image_kbytes, box_size, sep=', ')
print(msg)
