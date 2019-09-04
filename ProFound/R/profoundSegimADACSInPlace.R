profoundMakeSegimADACSInPlace=function(image=NULL, mask=NULL, objects=NULL, skycut=1, pixcut=3, tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE, SBlim=NULL, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, watershed = 'ProFound', ...){
  
  call=match.call()
  if(verbose){message(' - Running MakeSegim:')}
  timestart = proc.time()[3]

  #Treat image NAs as masked regions:
  
  if(!is.null(mask)){
    if(anyNA(image)){
      mask[is.na(image)]=1L
    }
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  hassky=!is.null(sky)
  hasskyRMS=!is.null(skyRMS)
  
  image_orig=image
  
  image_sky=image-sky
  
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
  
  if(smooth){
    if(verbose){message(paste(" - Smoothing the image -", round(proc.time()[3]-timestart,3), "sec"))}
    if(requireNamespace("imager", quietly = TRUE)){
      image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
    }else{
      if(!requireNamespace("EBImage", quietly = TRUE)){
        stop('The imager or EBImage package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
      }
      message(" - WARNING: imager package not installed, using EBImage gblur smoothing!")
      image=as.matrix(EBImage::gblur(image,sigma))
    }
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!is.null(SBlim) & !missing(magzero)){
    #image[image<skycut | image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
    image[image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
  }
  if(!is.null(mask)){
    image[mask>0]=0
  }
  if(verbose){message(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  if(any(image>0)){
    if(watershed=='EBImage'){
      if(!requireNamespace("EBImage", quietly = TRUE)){
        stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
      }
      image[image<skycut]=0
      segim=EBImage::imageData(EBImage::watershed(image,tolerance=tolerance,ext=ext))
      segtab=tabulate(segim)
      segim[segim %in% which(segtab<pixcut)]=0L
      mode(segim)='integer'
    }else if(watershed=='ProFound'){
      segim=water_cpp(image=image, nx=dim(image)[1], ny=dim(image)[2], abstol=tolerance, reltol=reltol, cliptol=cliptol, ext=ext, skycut=skycut, pixcut=pixcut, verbose=verbose)
    }else if(watershed=='ProFound-old'){
      segim=water_cpp_old(image=image, nx=dim(image)[1], ny=dim(image)[2], abstol=tolerance, reltol=reltol, cliptol=cliptol, ext=ext, skycut=skycut, pixcut=pixcut, verbose=verbose)
    }else{
      stop('watershed option must either be EBImage/ProFound/ProFound-old!')
    }
  }else{
    segim=image
  }
  
  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profoundSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  objects=matrix(0L,dim(segim)[1],dim(segim)[2])
  objects[]=as.logical(segim)
  
  if(stats & any(image>0)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE or no segments")}
    segstats=NULL
  }
  
  #if(!is.null(SBlim) & !missing(magzero)){
  #  SBlim=min(SBlim, profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=doRMNA)
  #}else if(is.null(SBlim) & !missing(magzero) & skycut>0){
  #  SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  #}else{
  #  SBlim=NULL
  #}
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - MakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  invisible(list(segim=segim, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, call=call))
}
profoundMakeSegimDilateADACSInPlace=function(image=NULL, segim=NULL, mask=NULL, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
  if(verbose){message(' - Running MakeSegimDilate:')}
  timestart = proc.time()[3]
  
  call=match.call()
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  #Treat image NAs as masked regions:
  
  if(!is.null(mask) & !is.null(image)){
    mask[is.na(image)]=1L
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(verbose){message(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(is.null(expand) | length(expand)==0){
    objects=matrix(0L,dim(segim)[1],dim(segim)[2])
    objects[]=as.logical(segim)
    
    if(stats){
      if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
      segstats=profoundSegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
    }else{
      if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    return(invisible(list(segim=segim, objects=objects, segstats=segstats, header=header, call=call)))
  }
  
  if(expand[1]=='all'){
    segim_new=segim
    maxorig=max(segim_new, na.rm=doRMNA)+1L
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    segim_new=EBImage::imageData(EBImage::dilate(segim_new, kern)) #Run Dilate
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    replace=which(segim!=0) #put back non-dilated segments
    segim_new[replace]=segim[replace] #put back non-dilated segments
  }else{
    segim_new=segim
    #segim_new[!(segim_new %in% expand)]=0L #remove things that will not be dilated
    if('fastmatch' %in% .packages()){ #remove things that will not be dilated
      segim_new[fastmatch::fmatch(segim_new, expand, nomatch = 0L) == 0L] = 0L
    }else{
      segim_new[!(segim_new %in% expand)] = 0L
    }
    maxorig=max(segim_new, na.rm=doRMNA)+1L
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    segim_new=EBImage::imageData(EBImage::dilate(segim_new, kern)) #Run Dilate
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    replace=which(segim!=0) #put back non-dilated segments
    segim_new[replace]=segim[replace] #put back non-dilated segments
    rm(replace)
  }
  mode(segim_new)='integer'
  
  rm(segim)
  
  if(!is.null(mask) & !is.null(image)){
    segim_new[mask!=0]=0
    image[mask!=0]=NA
  }
  
  if(stats & !is.null(image)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  objects=matrix(0L,dim(segim_new)[1],dim(segim_new)[2])
  objects[]=as.logical(segim_new)
  
  if(plot & !is.null(image)){
    profoundSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - profoundMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return(invisible(list(segim=segim_new, objects=objects, segstats=segstats, header=header, call=call)))
}