adacs_MakeSegim=function(image=NULL, mask=NULL, bmask=NULL, objects=NULL, skycut=1, pixcut=3, tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE, SBlim=NULL, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, watershed = 'ProFound', ...){
  
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
    segstats=adacs_SegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
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
adacs_MakeSegimDilate=function(image=NULL, segim=NULL, mask=NULL, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
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
      segstats=adacs_SegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
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
    segstats=adacs_SegimStats(image=image, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
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

adacs_SegimStats=function(image=NULL, segim=NULL, mask=NULL, sky=NULL, skyRMS=NULL, magzero=0, gain=NULL, pixscale=1, header=NULL, sortcol='segID', decreasing=FALSE, rotstats=FALSE, boundstats=FALSE, offset=1, cor_err_func=NULL, app_diam=1){
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
  }
  
  Napp=ceiling(pi*(app_diam/2/pixscale)^2)
  
  if(!is.null(sky)){
    hassky=any(is.finite(sky))
    if(hassky & length(sky)==1){
      sky=rep(sky,length(image))
    }
  }else{
    hassky=FALSE
  }
  if(hassky){
    image=image-sky
  }
  if(!is.null(skyRMS)){
    hasskyRMS=any(is.finite(skyRMS))
    if(hasskyRMS & length(skyRMS)==1){
      skyRMS=rep(skyRMS,length(image))
    }
  }else{
    hasskyRMS=FALSE
  }
  
  #Treat image NAs as masked regions:
  
  if(!is.null(mask)){
    mask[is.na(image)]=1L
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  #Set masked things to NA, to be safe:
  
  if(!is.null(mask)){
    image[mask!=0]=NA
    #segim[mask!=0]=NA
  }
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  #segvec=which(tabulate(segim)>0)
  #segvec=segvec[segvec>0]
  
  segsel=which(segim>0)
  
  xloc = rep(1:xlen, times = ylen)[segsel]
  yloc = rep(1:ylen, each = xlen)[segsel]
  
  if(hassky & hasskyRMS){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), sky=as.numeric(sky[segsel]), skyRMS=as.numeric(skyRMS[segsel]))
    rm(sky)
    rm(skyRMS)
  }
  if(hassky & hasskyRMS==FALSE){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), sky=as.numeric(sky[segsel]))
    rm(sky)
  }
  if(hassky==FALSE & hasskyRMS){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), skyRMS=as.numeric(skyRMS[segsel]))
    rm(skyRMS)
  }
  if(hassky==FALSE & hasskyRMS==FALSE){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]))
  }
  
  setkey(tempDT, segID, flux)
  
  rm(xloc)
  rm(yloc)
  rm(image)
  
  #tempDT[is.na(tempDT)]=0
  segID=tempDT[,.BY,by=segID]$segID
  
  x=NULL; y=NULL; flux=NULL; sky=NULL; skyRMS=NULL
  
  fluxout=tempDT[,.fluxcalc(flux,Napp=Napp), by=segID]
  fluxout$flux_app[which(fluxout$flux_app>fluxout$flux)]=fluxout$flux[which(fluxout$flux_app>fluxout$flux)]
  mag=profoundFlux2Mag(flux=fluxout$flux, magzero=magzero)
  mag_app=profoundFlux2Mag(flux=fluxout$flux_app, magzero=magzero)
  
  if(any(fluxout$flux==0, na.rm=doRMNA)){
    fluxout$N50seg[fluxout$flux==0]=fluxout$N100seg[fluxout$flux==0]
    fluxout$N90seg[fluxout$flux==0]=fluxout$N100seg[fluxout$flux==0]
  }
  
  if(hassky){
    #With one version of data.table sd doesn't work when all numbers are identical (complains about gsd and negative length vectors). Fixed with explicit sd caclulation until this gets fixed.
    flux_err_sky=tempDT[,sd(sky, na.rm=doRMNA)*1, by=segID]$V1*fluxout$N100seg
    #flux_err_sky=tempDT[,sqrt(sum((sky-mean(sky, na.rm=doRMNA))^2)/(.N-1)), by=segID]$V1*fluxout$N100seg
  }else{
    flux_err_sky=0
  }
  
  if(hasskyRMS){
    flux_err_skyRMS=tempDT[,sqrt(sum(skyRMS^2, na.rm=doRMNA)), by=segID]$V1
    pchi=pchisq(tempDT[,sum((flux/skyRMS)^2, na.rm=doRMNA), by=segID]$V1, df=fluxout$N100seg, log.p=TRUE)
    signif=qnorm(pchi, log.p=TRUE)
    FPlim=qnorm(1-fluxout$N100seg/(xlen*ylen))
  }else{
    flux_err_skyRMS=0
    signif=NA
    FPlim=NA
  }
  
  if(!is.null(gain)){
    flux_err_shot=sqrt(fluxout$flux)/gain
  }else{
    flux_err_shot=0
  }
  
  if(!is.null(cor_err_func)){
    cor_seg=cor_err_func(fluxout$N100seg)
    flux_err_cor=sqrt((flux_err_sky^2)/(1-cor_seg)-flux_err_sky^2)
  }else{
    cor_seg=0
    flux_err_cor=0
  }
  
  flux_err_sky[!is.finite(flux_err_sky)]=0
  flux_err_skyRMS[!is.finite(flux_err_skyRMS)]=0
  flux_err_shot[!is.finite(flux_err_shot)]=0
  flux_err_cor[!is.finite(flux_err_cor)]=0
  
  flux_err=sqrt(flux_err_sky^2+flux_err_skyRMS^2+flux_err_shot^2+flux_err_cor^2)
  mag_err=(2.5/log(10))*abs(flux_err/fluxout$flux)
  
  if(hassky){
    sky_mean=tempDT[,mean(sky, na.rm=doRMNA), by=segID]$V1
  }else{
    sky_mean=0
  }
  
  if(hasskyRMS){
    skyRMS_mean=tempDT[,mean(skyRMS, na.rm=doRMNA), by=segID]$V1
  }else{
    skyRMS_mean=0
  }
  
  xcen=tempDT[,.meanwt(x-0.5, flux),by=segID]$V1
  ycen=tempDT[,.meanwt(y-0.5, flux),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x-0.5,flux)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y-0.5,flux)),by=segID]$V1
  covxy=tempDT[,.covarwt(x-0.5,y-0.5,flux),by=segID]$V1
  
  xmax=xcen
  ymax=ycen
  xmax[!is.na(fluxout$flux)]=tempDT[,x[which.max(flux)]-0.5,by=segID]$V1
  ymax[!is.na(fluxout$flux)]=tempDT[,y[which.max(flux)]-0.5,by=segID]$V1
  
  sep=sqrt((xcen-xmax)^2+(ycen-ymax)^2)*pixscale
  
  pad=10^ceiling(log10(ylen+1))
  uniqueID=ceiling(xmax)*pad+ceiling(ymax)
  
  if(rotstats){
    asymm=tempDT[,.asymm(x-0.5,y-0.5,flux),by=segID]$V1
    flux_reflect=tempDT[,.reflect(x-0.5,y-0.5,flux),by=segID]$V1
    mag_reflect=profoundFlux2Mag(flux=flux_reflect, magzero=magzero)
  }else{
    asymm=NA
    flux_reflect=NA
    mag_reflect=NA
  }
  
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(abs(rad$hi)+0.08333333) #Added variance of uniform in quadrature (prevents zeros)
  rad$lo=sqrt(abs(rad$lo)+0.08333333) #Added variance of uniform in quadrature (prevents zeros)
  axrat=rad$lo/rad$hi
  eigvec=.cov2eigvec(xsd, ysd, covxy)
  ang=.eigvec2ang(eigvec)
  
  R50seg=sqrt(fluxout$N50seg/(axrat*pi))*pixscale
  R90seg=sqrt(fluxout$N90seg/(axrat*pi))*pixscale
  R100seg=sqrt(fluxout$N100seg/(axrat*pi))*pixscale
  
  con=R50seg/R90seg
  con[R90seg==0]=NA
  
  SB_N50=profoundFlux2SB(flux=fluxout$flux*0.5/fluxout$N50seg, magzero=magzero, pixscale=pixscale)
  SB_N90=profoundFlux2SB(flux=fluxout$flux*0.9/fluxout$N90seg, magzero=magzero, pixscale=pixscale)
  SB_N100=profoundFlux2SB(flux=fluxout$flux/fluxout$N100seg, magzero=magzero, pixscale=pixscale)
  
  if(!is.null(header)){
    coord=magWCSxy2radec(xcen, ycen, header=header)
    RAcen=coord[,1]
    Deccen=coord[,2]
    coord=magWCSxy2radec(xmax, ymax, header=header)
    RAmax=coord[,1]
    Decmax=coord[,2]
  }else{
    RAcen=NA
    Deccen=NA
    RAmax=NA
    Decmax=NA
  }
  
  if(boundstats){
    segim_inner=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)]
    off_down=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)-offset]
    off_left=segim[(offset+1):(xlen-offset)-offset,(offset+1):(ylen-offset)]
    off_up=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)+offset]
    off_right=segim[(offset+1):(xlen-offset)+offset,(offset+1):(ylen-offset)]
    
    inner_segim=segim_inner>0 & off_down==segim_inner & off_left==segim_inner & off_up==segim_inner & off_right==segim_inner
    
    segim_edge=segim_inner
    segim_edge[inner_segim==1]=0
    tab_edge=tabulate(segim_edge)
    tab_edge=c(tab_edge,rep(0,max(segID)-length(tab_edge)))
    tab_edge=cbind(1:max(segID),tab_edge)
    Nedge=tab_edge[match(segID,tab_edge[,1]),2]
    
    outer_sky=segim_inner>0 & (off_down==0 | off_left==0 | off_up==0 | off_right==0)
    
    segim_sky=segim_inner
    segim_sky[outer_sky==0]=0
    tab_sky=tabulate(segim_sky)
    tab_sky=c(tab_sky,rep(0,max(segID)-length(tab_sky)))
    tab_sky=cbind(1:max(segID),tab_sky)
    Nsky=tab_sky[match(segID,tab_sky[,1]),2]
    
    rm(off_down)
    rm(off_left)
    rm(off_up)
    rm(off_right)
    rm(tab_edge)
    rm(tab_sky)
    
    BorderBottom=segim[segim[,1]>0,1]
    BorderLeft=segim[1,segim[1,]>0]
    BorderTop=segim[segim[,ylen]>0,ylen]
    BorderRight=segim[xlen,segim[xlen,]>0]
    tab_bottom=tabulate(BorderBottom)
    tab_left=tabulate(BorderLeft)
    tab_top=tabulate(BorderTop)
    tab_right=tabulate(BorderRight)
    
    tab_border=cbind(1:max(segID),0,0,0,0)
    tab_border[1:length(tab_bottom),2]=tab_border[1:length(tab_bottom),2]+tab_bottom
    tab_border[1:length(tab_left),3]=tab_border[1:length(tab_left),3]+tab_left
    tab_border[1:length(tab_top),4]=tab_border[1:length(tab_top),4]+tab_top
    tab_border[1:length(tab_right),5]=tab_border[1:length(tab_right),5]+tab_right
    bordersel=match(segID,tab_border[,1])
    Nborder=tab_border[bordersel,2]+tab_border[bordersel,3]+tab_border[bordersel,4]+tab_border[bordersel,5]
    flag_border=1*(tab_border[bordersel,2]>0)+2*(tab_border[bordersel,3]>0)+4*(tab_border[bordersel,4]>0)+8*(tab_border[bordersel,5]>0)
    
    if(!is.null(mask)){
      outer_mask=segim_inner>0 & (mask[2:(xlen-1)+1,2:(ylen-1)]==1 | mask[2:(xlen-1)-1,2:(ylen-1)]==1 | mask[2:(xlen-1),2:(ylen-1)+1]==1 | mask[2:(xlen-1),2:(ylen-1)-1]==1)
      segim_mask=segim_edge
      segim_mask[outer_mask==0]=0
      tab_mask=tabulate(segim_mask)
      tab_mask=c(tab_mask,rep(0,max(segID)-length(tab_mask)))
      tab_mask=cbind(1:max(segID),tab_mask)
      Nmask=tab_mask[match(segID,tab_mask[,1]),2]
      
      rm(segim_inner)
      rm(tab_mask)
      
    }else{
      Nmask=0
    }
    
    Nedge=Nedge+Nborder
    #Nsky=Nsky-Nmask #Raw Nsky-Nmask-Nborder, to correct for masked pixels
    Nobject=Nedge-Nsky-Nborder # Nedge-Nsky
    edge_frac=Nsky/Nedge
    
    #Using Ramanujan approximation from Wikipedia:
    
    A=R100seg/pixscale
    B=R100seg*axrat/pixscale
    h=(A-B)^2/(A+B)^2
    C=pi*(A+B)*(1+(3*h)/(10+sqrt(4-3*h)))
    edge_excess=Nedge/C
    
  }else{
    Nedge=NA
    Nsky=NA
    Nobject=NA
    Nborder=NA
    Nmask=NA
    flag_border=NA
    edge_frac=NA
    edge_excess=NA
  }
  
  if(anyNA(fluxout$flux)){
    bad=is.na(fluxout$flux)
    asymm[bad]=NA
    flux_reflect[bad]=NA
    mag_reflect[bad]=NA
    signif[bad]=NA
    FPlim[bad]=NA
    edge_excess[bad]=NA
  }
  
  segstats=data.table(segID=segID, uniqueID=uniqueID, xcen=xcen, ycen=ycen, xmax=xmax, ymax=ymax, RAcen=RAcen, Deccen=Deccen, RAmax=RAmax, Decmax=Decmax, sep=sep, flux=fluxout$flux, mag=mag, flux_app=fluxout$flux_app, mag_app=mag_app, cenfrac=fluxout$cenfrac, N50=fluxout$N50seg, N90=fluxout$N90seg, N100=fluxout$N100seg, R50=R50seg, R90=R90seg, R100=R100seg, SB_N50=SB_N50, SB_N90=SB_N90, SB_N100=SB_N100, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, semimaj=rad$hi, semimin=rad$lo, axrat=axrat, ang=ang, signif=signif, FPlim=FPlim, flux_err=flux_err, mag_err=mag_err, flux_err_sky=flux_err_sky, flux_err_skyRMS=flux_err_skyRMS, flux_err_shot=flux_err_shot, flux_err_cor=flux_err_cor, cor_seg=cor_seg, sky_mean=sky_mean, sky_sum=sky_mean*fluxout$N100seg, skyRMS_mean=skyRMS_mean, Nedge=Nedge, Nsky=Nsky, Nobject=Nobject, Nborder=Nborder, Nmask=Nmask, edge_frac=edge_frac, edge_excess=edge_excess, flag_border=flag_border)
  invisible(as.data.frame(segstats[order(segstats[[sortcol]], decreasing=decreasing),]))
}