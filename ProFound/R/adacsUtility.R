compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}
magcutoutADACSInPlace = function (image, loc = dim(image)/2, box = c(100, 100), oimage=NULL) 
{
  # Assume box is an even size
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
  #print(paste("xoffset ",xoffset," xlo",xlo))
  #print(paste("yoffset ",yoffset," ylo",ylo))
  if (yhi-ylo < 0 | xhi-xlo < 0) {
    image = matrix(NA, box[1], box[2])
  }
  else {
    # TODO: Implement this as a C/C++ method
    #print(paste("loc ",loc))
    subset_cpp_inplace(image, ylo, yhi, xlo, xhi, yoffset, xoffset, oimage)
    image = oimage
  }
  output = list(image = image)
  invisible(output)
}
