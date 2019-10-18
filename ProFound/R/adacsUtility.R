compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}
initialiseGlobals = function(doclip)
{
  if (doclip) {
    doRMNA <<- FALSE
  } else {
    doRMNA <<- TRUE
  }
  # some enums
  MEDIAN <<- 1  # 'median'
  MEAN <<- 2    # 'mean'
  MODE <<- 3    # 'mode'
  RMEDIAN <<- 4  # R version of 'median'
  RMEAN <<- 5    # R version of 'mean'
  RMODE <<- 6    # R version of 'mode'
  
  BOTH <<- 1  # 'quanboth'
  LO <<- 2    # 'quanlo'
  HI <<- 3    # 'quanhi'
  SD <<- 4    # 'sd'
  RBOTH <<- 5  # R version of 'quanboth'
  RLO <<- 6    # R version of 'quanlo'
  RHI <<- 7    # R version of 'quanhi'
  RSD <<- 8    # R version of 'sd'
  
  AUTO <<- 1
  SET <<- 2
  
  CLASSIC_BILINEAR <<- 1
  AKIMA_BICUBIC <<- 2
}
enumForKeyword = function(keyword)
{
  result = NULL
  if (stri_detect_fixed(keyword,"median",case_insensitive=TRUE)) {
    result = MEDIAN
  }
  if (stri_detect_fixed(keyword,"mean",case_insensitive=TRUE)) {
    result = MEAN
  }
  if (stri_detect_fixed(keyword,"mode",case_insensitive=TRUE)) {
    result = MODE
  }
  if (stri_detect_fixed(keyword,"rmedian",case_insensitive=TRUE)) {
    result = RMEDIAN
  }
  if (stri_detect_fixed(keyword,"rmean",case_insensitive=TRUE)) {
    result = RMEAN
  }
  if (stri_detect_fixed(keyword,"rmode",case_insensitive=TRUE)) {
    result = RMODE
  }
  if (stri_detect_fixed(keyword,"quanboth",case_insensitive=TRUE)) {
    result = BOTH
  }
  if (stri_detect_fixed(keyword,"quanlo",case_insensitive=TRUE)) {
    result = LO
  }
  if (stri_detect_fixed(keyword,"quanhi",case_insensitive=TRUE)) {
    result = HI
  }
  if (stri_detect_fixed(keyword,"sd",case_insensitive=TRUE)) {
    result = SD
  }
  if (stri_detect_fixed(keyword,"rquanboth",case_insensitive=TRUE)) {
    result = RBOTH
  }
  if (stri_detect_fixed(keyword,"rquanlo",case_insensitive=TRUE)) {
    result = RLO
  }
  if (stri_detect_fixed(keyword,"rquanhi",case_insensitive=TRUE)) {
    result = RHI
  }
  if (stri_detect_fixed(keyword,"rsd",case_insensitive=TRUE)) {
    result = RSD
  }
  if (stri_detect_fixed(keyword,"auto",case_insensitive=TRUE)) {
    result = AUTO
  }
  if (stri_detect_fixed(keyword,"set",case_insensitive=TRUE)) {
    result = SET
  }
  if (stri_detect_fixed(keyword,"bilinear",case_insensitive=TRUE)) {
    result = CLASSIC_BILINEAR
  }
  if (stri_detect_fixed(keyword,"bicubic",case_insensitive=TRUE)) {
    result = AKIMA_BICUBIC
  }
  invisible(result)
}
# Helper function for Rcpp access
adacs_mode = function(clip)
{
  temp=density(clip, na.rm=TRUE)
  result=temp$x[which.max(temp$y)]
  invisible(result)
}
# This version assumes the input image is of type double
magcutoutADACSInPlace = function (image, loc = dim(image)/2, box = c(101, 101), oimage=NULL) 
{
  # Assume box is an odd size
  # Assume loc and box are integers
  shiftloc = FALSE
  paddim = TRUE
  expand = FALSE
  loc = as.numeric(loc)
  xcen = loc[1]
  ycen = loc[2]
  loc = ceiling(loc)
  boxhalfx = floor(box[1]/2)
  boxhalfy = floor(box[2]/2)
  if(length(box)==1){box=rep(box,2)}
  xlo = ceiling(loc[1] - (box[1]/2 - 0.5))
  xhi = ceiling(loc[1] + (box[1]/2 - 0.5))
  ylo = ceiling(loc[2] - (box[2]/2 - 0.5))
  yhi = ceiling(loc[2] + (box[2]/2 - 0.5))
  
  loc.diff = c(x=xlo-1, y=ylo-1)
  
  expand = paddim && shiftloc
  diffxlo = xlo - 1
  if (diffxlo < 0) {
    xlo = 1
    if(expand) xhi = xlo + (box[1] - 1)
  }
  diffxhi = xhi - dim(image)[1]
  if (diffxhi > 0) {
    xhi = dim(image)[1]
    if(expand) {
      xlo = xlo - diffxhi
      if(xlo < 1) xlo = 1
    }
  }
  diffylo = ylo - 1
  if (diffylo < 0) {
    ylo = 1
    if(expand) yhi = ylo + (box[2] - 1)
  }
  diffyhi = yhi - dim(image)[2]
  if (diffyhi > 0) {
    yhi = dim(image)[2]
    if(expand) {
      ylo = ylo - diffyhi
      if(ylo < 1) ylo = 1
    }
  }
  xoffset = floor(xcen-(box[1]/2-0.5))
  yoffset = floor(ycen-(box[2]/2-0.5))
  if (yhi-ylo < 0 | xhi-xlo < 0) {
    #image = matrix(NA, box[1], box[2])
    # TODO: test
    subset_cpp_inplace(image, ylo, yhi, xlo, xhi, yoffset, xoffset, oimage)
    image = oimage
  }
  else {
    subset_cpp_inplace(image, ylo, yhi, xlo, xhi, yoffset, xoffset, oimage)
    image = oimage
  }
  output = list(image = image)
  invisible(output)
}
# This version assumes the input image is of type int
magcutoutADACSInPlaceI = function (image, loc = dim(image)/2, box = c(100, 100), oimage=NULL) 
{
  # Assume box is an odd size
  # Assume loc and box are integers
  shiftloc = FALSE
  paddim = TRUE
  expand = FALSE
  loc = as.numeric(loc)
  xcen = loc[1]
  ycen = loc[2]
  loc = ceiling(loc)
  boxhalfx = floor(box[1]/2)
  boxhalfy = floor(box[2]/2)
  if(length(box)==1){box=rep(box,2)}
  xlo = ceiling(loc[1] - (box[1]/2 - 0.5))
  xhi = ceiling(loc[1] + (box[1]/2 - 0.5))
  ylo = ceiling(loc[2] - (box[2]/2 - 0.5))
  yhi = ceiling(loc[2] + (box[2]/2 - 0.5))
  
  loc.diff = c(x=xlo-1, y=ylo-1)
  
  expand = paddim && shiftloc
  diffxlo = xlo - 1
  if (diffxlo < 0) {
    xlo = 1
    if(expand) xhi = xlo + (box[1] - 1)
  }
  diffxhi = xhi - dim(image)[1]
  if (diffxhi > 0) {
    xhi = dim(image)[1]
    if(expand) {
      xlo = xlo - diffxhi
      if(xlo < 1) xlo = 1
    }
  }
  diffylo = ylo - 1
  if (diffylo < 0) {
    ylo = 1
    if(expand) yhi = ylo + (box[2] - 1)
  }
  diffyhi = yhi - dim(image)[2]
  if (diffyhi > 0) {
    yhi = dim(image)[2]
    if(expand) {
      ylo = ylo - diffyhi
      if(ylo < 1) ylo = 1
    }
  }
  xoffset = floor(xcen-(box[1]/2-0.5))
  yoffset = floor(ycen-(box[2]/2-0.5))
  if (yhi-ylo < 0 | xhi-xlo < 0) {
    #image = matrix(NA, box[1], box[2])
    # TODO: test
    subset_cpp_inplaceI(image, ylo, yhi, xlo, xhi, yoffset, xoffset, oimage)
    image = oimage
  }
  else {
    subset_cpp_inplaceI(image, ylo, yhi, xlo, xhi, yoffset, xoffset, oimage)
    image = oimage
  }
  output = list(image = image)
  invisible(output)
}
FindSkyCellValues=function(image=NULL, objects=NULL, mask=NULL, loc=dim(image)/2, box=c(100,100),
                           skypixmin=prod(box)/2, boxadd=box/2, boxiters=0,shiftloc=FALSE,paddim=TRUE){
  # Work out if we need to expand the box (and keep indices of NON source pixels)
  skyN=0
  iterN=0
  tempcomb={}
  while(skyN<skypixmin & iterN<=boxiters){
    if(!is.null(objects)){
      tempcomb=magcutout(image=objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
      if(!is.null(mask)){
        tempcomb=tempcomb & (magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0)
      }
    }else{
      tempcomb=magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
    }
    tempcomb[is.na(tempcomb)]=FALSE
    if(!is.null(tempcomb)){
      tempcomb=which(tempcomb)
      skyN=length(tempcomb)
    }else{
      skyN=0
    }
    box=box+boxadd
    iterN=iterN+1
  }
  box=box-boxadd #since one too many boxadds will have occurred when it terminates
  
  if(skyN>0){
    select=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image[tempcomb]
  }else{
    select=NA
  }
  invisible(list(select=select, box=box, skyN=skyN))
}

adacs_MakeSkyGrid=function(image=NULL, bobjects=NULL, bmask=NULL, box=c(100,100), grid=box, type='bicubic', skytype='median', skyRMStype='quanlo', sigmasel=1,
                                               skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, cores=1,
                                               scratch=NULL) {
  temp_bi_sky = scratch[['scratchSKY']]
  temp_bi_skyRMS = scratch[['scratchSKYRMS']]
  adacs<-new(Adacs)
  if (is.null(bmask)) {
    bmask = new(BitMatrix, dim(image)[1], dim(image)[2])
  }
  adacs$Cadacs_MakeSkyGrid(image, bobjects, bmask, 
                     box[1], box[2],
                     grid[1], grid[2],
                     boxadd[1], boxadd[2],
                     enumForKeyword(type), skypixmin, boxiters,
                     doclip, enumForKeyword(skytype), enumForKeyword(skyRMStype), sigmasel,
                     temp_bi_sky, temp_bi_skyRMS)
  
  invisible(list(sky=temp_bi_sky, skyRMS=temp_bi_skyRMS))
}

adacs_MakeMask=function(image, maskflag=NULL) {
  bmask=new(BitMatrix, dim(image)[1], dim(image)[2])
  if (!is.null(maskflag)) {
    bmask$maskValue(maskflag)
  } else {
    bmask$maskNaN(image)
  }
}

mlook=function(H,r)
{
  # Once the input values
  # are added to H, the median value lies in the first index for which
  # the sum of values to that index reaches 2r2 + 2r + 1.
  sum = 0
  rmax = 2*r*r+2*r+1
  for (i in 1:256) {
    sum = sum + H[i]
    if (sum>=rmax) {
      break;
    }
  }
  invisible(i)
}
median2DR=function(input,r)
{
  S = dim(input)[1]-2*r
  xsel = as.integer((r+1):(r+S))
  ysel = as.integer((r+1):(r+S))
  output = input[xsel,ysel]
  
  H = replicate(256,0)
  
  startRow = 2
  endRow = S
  
  for (col in 1:S) {
    # Note we are 1 relative (The fastmedian_5506.pdf is 0 relative)
    # initialize H to I[0 .. 2r][0 .. 2r]. // yellow region
    for (j in 1:256) {
      H[j] = 0
    }
    for (row in 1:(2*r+1)) {
      for (j in 1:(2*r+1)) {
        index = input[row, j]
        H[index] = H[index] + 1
      }
    }
    # find median value m in H, write m to O[0][0].
    output[1,col] = mlook(H,r)
  
    # for row = 1 to S - 1:
    for (row in startRow:endRow) {
      # add values I[2r + row][0 .. 2r] to H.
      for (j in col:(col+2*r)) {
        index = input[2*r+row, j]
        H[index] = H[index] + 1
      }
      # subtract values I[row - 1][0 .. 2r] from H.
      for (j in col:(col+2*r)) {
        index = input[row-1, j]
        H[index] = H[index] - 1
      }
      # find median value m in H; write m to O[row][0].
      output[row,col] = mlook(H,r)
    }
  }
  
  invisible(output)
}

