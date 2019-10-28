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
  if (tolower(keyword) == 'adacs_median') {
    result = MEDIAN
  }
  if (tolower(keyword) == 'adacs_mean') {
    result = MEAN
  }
  if (tolower(keyword) == 'adacs_mode') {
    result = MODE
  }
  if (tolower(keyword) == 'median') {
    result = RMEDIAN
  }
  if (tolower(keyword) == 'mean') {
    result = RMEAN
  }
  if (tolower(keyword) == 'mode') {
    result = RMODE
  }
  if (tolower(keyword) == 'adacs_quanboth') {
    result = BOTH
  }
  if (tolower(keyword) == 'adacs_quanlo') {
    result = LO
  }
  if (tolower(keyword) == 'adacs_quanhi') {
    result = HI
  }
  if (tolower(keyword) == 'adacs_sd') {
    result = SD
  }
  if (tolower(keyword) == 'quanboth') {
    result = RBOTH
  }
  if (tolower(keyword) == 'quanlo') {
    result = RLO
  }
  if (tolower(keyword) == 'quanhi') {
    result = RHI
  }
  if (tolower(keyword) == 'sd') {
    result = RSD
  }
  if (tolower(keyword) == 'auto') {
    result = AUTO
  }
  if (tolower(keyword) == 'set') {
    result = SET
  }
  if (tolower(keyword) == 'bilinear') {
    result = CLASSIC_BILINEAR
  }
  if (tolower(keyword) == 'bicubic') {
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
# Helper function for users wanting to create BitMatrix masks
adacs_MakeMask=function(image, maskflag=NULL) {
  bmask=new(BitMatrix, dim(image)[1], dim(image)[2])
  if (!is.null(maskflag)) {
    bmask$maskValue(maskflag)
  } else {
    bmask$maskNaN(image)
  }
}

adacs_MakeSkyGrid=function(image=NULL, bobjects=NULL, bmask=NULL, box=c(100,100), grid=box, type='bicubic', skytype='median', skyRMStype='quanlo', sigmasel=1,
                                               skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, cores=1,
                                               scratch=NULL) {
  temp_bi_sky = scratch[['scratchSKY']]
  temp_bi_skyRMS = scratch[['scratchSKYRMS']]
  adacs<-new(Adacs)
  if (is.null(bobjects)) {
    bobjects = new(BitMatrix, dim(image)[1],dim(image)[2])
  }
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

