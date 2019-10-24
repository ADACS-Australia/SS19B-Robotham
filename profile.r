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
result=gc()

# Go, go, go!
if (what == 'profound-original') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
        mask[1,1] = 1L
	result = profoundProFound(image, box=c(box_size, box_size), mask=mask)
} else if (what == 'profound-adacs') {
	bmask = new(BitMatrix,dim(image)[1],dim(image)[2])
        bmask$settrue(1,1)
	result = adacs_ProFound(image, box=c(box_size, box_size), bmask=bmask)
} else if (what == 'skygrid-original') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
        mask[1,1] = 1L
	result = profoundMakeSkyGrid(image, box=c(box_size, box_size), mask=mask)
} else if (what == 'skygrid-adacs') {
	bmask = new(BitMatrix,dim(image)[1],dim(image)[2])
        bmask$settrue(1,1)
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
	result = adacs_MakeSkyGrid(image, box=box, scratch=scratch, bmask=bmask)
} else if (what == 'baseline1') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
        mask[dim(image)[1],dim(image)[2]] = 1L
} else if (what == 'baseline2') {
	bmask = new(BitMatrix,dim(image)[1],dim(image)[2])
        bmask$settrue(dim(image)[1],dim(image)[2])
} else if (what == 'baseline3') {
	mask = matrix(0L,dim(image)[1],dim(image)[2])
        for (i in 1:dim(image)[1]) {
            mask[i,i] = 1L
        }
        result = which(mask>0)
} else if (what == 'baseline4') {
	bmask = new(BitMatrix,dim(image)[1],dim(image)[2])
        bmask$settrue(dim(image)[1],dim(image)[2])
        for (i in 1:dim(image)[1]) {
            bmask$settrue(i,i)
        }
        result = bmask$trues()
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
