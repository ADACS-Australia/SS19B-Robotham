compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}
magcutoutADACS = function (image, loc = dim(image)/2, box = c(100, 100)) 
{
  # Assume box is an even size
  # Assume loc and box are integers
  shiftloc = FALSE
  paddim = TRUE
  expand = FALSE
  loc = as.numeric(loc)
  xcen = ceiling(loc[1])
  ycen = ceiling(loc[2])
  if (length(box) == 1) {
    box = rep(box, 2)
  }
  boxhalfx = floor(box[1]/2)
  boxhalfy = floor(box[2]/2)
  
  
  xloo = xcen - boxhalfx
  xhio = xcen + boxhalfx
  yloo = ycen - boxhalfy
  yhio = ycen + boxhalfy
  xlo = xloo
  xhi = xhio
  ylo = yloo
  yhi = yhio
  needsPadding = FALSE
  if (xlo <= 0) {
    xlo = 1
    needsPadding = TRUE
  }
  if (xhi > dim(image)[1]) {
    xhi = dim(image)[1]
    needsPadding = TRUE
  }
  if (ylo <= 0) {
    ylo = 1
    needsPadding = TRUE
  }
  if (yhi > dim(image)[2]) {
    yhi = dim(image)[2]
    needsPadding = TRUE
  }
  xoffset = xcen-(boxhalfx+1)
  yoffset = ycen-(boxhalfy+1)
  xsel = as.integer(xlo:xhi)
  ysel = as.integer(ylo:yhi)
  if (length(xsel) == 0 | length(ysel) == 0) {
    image = matrix(NA, box[1], box[2])
  }
  else {
    # TODO: Implement this as a C/C++ method
    image <- subset_cpp(image, ylo, yhi, xlo, xhi)
    #image = image[xsel, ysel]
    #if (all(image == yyy)) {
    #  print("ok")
    #} else {
    #  print("$$$$$$$$$$$$$$ NOT OK $$$$$$$$$$$$$")
    #}
    
    if (needsPadding) {
      padded = matrix(NA, box[1], box[2])
      #print(paste("pad ",xlo," ",xoffset," ",ylo," ",yoffset))
      padded[xsel - xoffset, ysel - yoffset] = image
      image = padded
    }
  }
  output = list(image = image)
  invisible(output)
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
magcutoutADACSInPlaceI = function (image, loc = dim(image)/2, box = c(100, 100), oimage=NULL) 
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
    print(paste("subset_cpp_inplaceI ",ylo, " ", yhi, " ", xlo, " ", xhi))
    subset_cpp_inplaceI(image, ylo, yhi, xlo, xhi, yoffset, xoffset, oimage)
    image = oimage
  }
  output = list(image = image)
  invisible(output)
}
magclipADACS=function(x, sigma='auto', clipiters=5, sigmasel=1, estimate='both'){
  sel = !is.na(x) & !is.nan(x) & !is.null(x) & is.finite(x)
  clipx=sort(x[sel])
  #oclipx=order(x)
  
  if(clipiters>0 & length(x)>0){
    newlen=length(clipx)
    sigcut=pnorm(sigmasel)
    
    for(i in 1:clipiters){
      if(newlen<=1){break}
      oldlen=newlen
      roughmed=clipx[newlen/2]
      if(sigma=='auto'){
        clipsigma=qnorm(1-2/max(newlen,2,na.rm=TRUE))
      }else{
        clipsigma=sigma
      }
      if(estimate=='both'){
        #vallims=clipsigma*diff(quantile(clipx,c(1-sigcut,sigcut)))/2/sigmasel
        vallims=clipsigma*(clipx[sigcut*newlen]-clipx[(1-sigcut)*newlen])/2/sigmasel
      }
      if(estimate=='lo'){
        #vallims=clipsigma*(roughmed-quantile(clipx,1-sigcut))/sigmasel
        vallims=clipsigma*(roughmed-clipx[(1-sigcut)*newlen])/sigmasel
      }
      if(estimate=='hi'){
        #vallims=clipsigma*(quantile(clipx,sigcut)-roughmed)/sigmasel
        vallims=clipsigma*(clipx[sigcut*newlen]-roughmed)/sigmasel
      }
      
      clipx=clipx[clipx>=(roughmed-vallims) & clipx<=(roughmed+vallims)]
      newlen=length(clipx)
      
      if(oldlen==newlen){break}
    }
  }else{
    clipx=x
  }
  cliplogic=NA
  range=NA
  invisible(list(x=clipx, clip=cliplogic, range=range, clipiters=i))
}

magclipADACSOrder=function(x, sigma='auto', clipiters=5, sigmasel=1, estimate='both'){
  sel = !is.na(x) & !is.nan(x) & !is.null(x) & is.finite(x)
  clipx=sort(x[sel])
  
  if(clipiters>0 & length(x)>0){
    newlen=length(clipx)
    sigcut=pnorm(sigmasel)
    
    for(i in 1:clipiters){
      if(newlen<=1){break}
      oldlen=newlen
      roughmed=clipx[newlen/2]
      if(sigma=='auto'){
        clipsigma=qnorm(1-2/max(newlen,2,na.rm=TRUE))
      }else{
        clipsigma=sigma
      }
      if(estimate=='both'){
        #vallims=clipsigma*diff(quantile(clipx,c(1-sigcut,sigcut)))/2/sigmasel
        vallims=clipsigma*(clipx[sigcut*newlen]-clipx[(1-sigcut)*newlen])/2/sigmasel
      }
      if(estimate=='lo'){
        #vallims=clipsigma*(roughmed-quantile(clipx,1-sigcut))/sigmasel
        vallims=clipsigma*(roughmed-clipx[(1-sigcut)*newlen])/sigmasel
      }
      if(estimate=='hi'){
        #vallims=clipsigma*(quantile(clipx,sigcut)-roughmed)/sigmasel
        vallims=clipsigma*(clipx[sigcut*newlen]-roughmed)/sigmasel
      }
      
      clipx=clipx[clipx>=(roughmed-vallims) & clipx<=(roughmed+vallims)]
      newlen=length(clipx)
      
      if(oldlen==newlen){break}
    }
  }else{
    clipx=x
  }
  cliplogic=NA
  range=NA
  invisible(list(x=clipx, clip=cliplogic, range=range, clipiters=i))
}
adacsFindSkyCellValues=function(image=NULL, objects=NULL, mask=NULL, loc=dim(image)/2, box=c(100,100),
                                         skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, scratch=NULL){
      # Work out if we need to expand the box (and keep indices of NON source pixels)
      skyN=0
      iterN=0
      tempcomb={}
      while(skyN<skypixmin & iterN<=boxiters){
        #print("SKY ITERATION INVOLVES LOGICALS")
        #rm(tempcomb)
        if(!is.null(objects)){
          #print("SUBSET BOOLEAN MATRIX WITH OBJECTS")
          ootempcomb = scratch[['scratchI1']]
          #tempcomb=magcutout(image=objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
          magcutoutADACSInPlaceI(image=objects, loc=loc, box=box, oimage=ootempcomb)
          #DO above ootempcomb = ootempcomb==0
          tempcomb = ootempcomb
          # assert ootempcomb == select
          #if(!compareNA(ootempcomb,tempcomb)){
          #  print("Problem at C")
          #}
          if(!is.null(mask)){
            print("ADD MASKING!!")
            # HANG-ON! what are you doing?
            #tempcomb=tempcomb & (magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0)
            oootempcomb = matrix(0.0,box[1],box[2])
            magcutoutADACSInPlaceI(image=mask, loc=loc, box=box, oimage=oootempcomb)
            #DO above oootempcomb = oootempcomb==0
            tempcomb = tempcomb & oootempcomb
          }
        }else{
          #print("SUBSET BOOLEAN MATRIX WITHOUT OBJECTS - MUST BE MASK ONLY?")
          tempcomb = scratch[['scratchI1']]
          #tempcomb=magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
          magcutoutADACSInPlaceI(image=mask, loc=loc, box=box, oimage=tempcomb)
          # Invert the booleans
          # DO above tempcomb = tempcomb==0
        }
        tempcomb[is.na(tempcomb)]=FALSE # is this even possible?
        if(!is.null(tempcomb)){
          # this is counting how many pixels are NOT source pixels (because we inverted the booleans)
          tempcomb=which(tempcomb) # keep the indices
          skyN=length(tempcomb)
        }else{
          skyN=0
        }
        box=box+boxadd
        iterN=iterN+1
      }
      box=box-boxadd #since one too many boxadds will have occurred when it terminates
      
      if(skyN>0){
        select = scratch[['scratchN2']]
        #select=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image[tempcomb]
        magcutoutADACSInPlace(image,loc=loc,box=box,oimage=select)
        select = select[tempcomb] # this changes the matrix to a vector
        
        #rm(tempcomb)
      }else{
        select=NA
      }
      invisible(list(select=select, box=box, skyN=skyN))
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

