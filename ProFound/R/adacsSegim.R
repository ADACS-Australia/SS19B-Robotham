.adacs_FluxCalcMin=function(image=NULL, segim=NULL, bmask=NULL){
  
  #Set masked things to NA, to be safe:
  image[bmask$trues()]=NA
  
  segsel=which(segim>0)
  segID=flux=NULL
  tempDT=data.table(segID=as.integer(segim[segsel]),flux=as.numeric(image[segsel]))
  
  output=tempDT[,.fluxcalcmin(flux), by=segID]
  setkey(output, segID)
  
  return(as.data.frame(output))
}
adacs_MakeSegim=function(image=NULL, bmask=NULL, bobjects=NULL, skycut=1, pixcut=3, tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE, SBlim=NULL, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, watershed = 'ProFound', ...){
  
  call=match.call()
  if(verbose){message(' - Running MakeSegim:')}
  timestart = proc.time()[3]
  
  #Treat image NAs as masked regions:
  
  bmask$maskNaN(image)
  
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
      stop('The imager package is needed for smoothing to work.', call. = FALSE)
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
  image[bmask$trues()]=0
  if(verbose){message(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  if(any(image>0)){
    if(watershed=='ProFound'){
      segim=water_cpp(image=image, nx=dim(image)[1], ny=dim(image)[2], abstol=tolerance, reltol=reltol, cliptol=cliptol, ext=ext, skycut=skycut, pixcut=pixcut, verbose=verbose)
    }else if(watershed=='ProFound-old'){
      segim=water_cpp_old(image=image, nx=dim(image)[1], ny=dim(image)[2], abstol=tolerance, reltol=reltol, cliptol=cliptol, ext=ext, skycut=skycut, pixcut=pixcut, verbose=verbose)
    }else{
      stop('watershed option must either be ProFound/ProFound-old!')
    }
  }else{
    segim=image
  }
  
  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    adacs_SegimPlot(image=image_orig, segim=segim, bmask=bmask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  bobjects=new(BitMatrix, segim)
  
  if(stats & any(image>0)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=adacs_SegimStats(image=image_orig, segim=segim, bmask=bmask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE or no segments")}
    segstats=NULL
  }
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - MakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  invisible(list(segim=segim, bobjects=bobjects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, call=call))
}
adacs_MakeSegimDilate=function(image=NULL, segim=NULL, bmask=NULL, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
  if(verbose){message(' - Running MakeSegimDilate:')}
  timestart = proc.time()[3]
  
  
  
  call=match.call()
  
  #Treat image NAs as masked regions:
  
  if (!is.null(image)) {
    if (is.null(bmask)) {
      bmask=new(BitMatrix, dim(image)[1],dim(image)[2])
    }
    bmask$maskNaN(image)
  }
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  kern = .makeBrush(size, shape=shape)
  
  if(verbose){message(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(is.null(expand) | length(expand)==0){
    
    if(stats){
      if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
      segstats=adacs_SegimStats(image=image, segim=segim, bmask=bmask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
    }else{
      if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    return(invisible(list(segim=segim, segstats=segstats, header=header, call=call)))
  }
  if(expand[1]=='all'){
    segim_new = .dilate_cpp(segim, kern)

  }else{
    segim_new=segim
    if('fastmatch' %in% .packages()){ #remove things that will not be dilated
      segim_new[fastmatch::fmatch(segim_new, expand, nomatch = 0L) == 0L] = 0L
    }else{
      segim_new[!(segim_new %in% expand)] = 0L
    }
    segim_new = .dilate_cpp(segim_new, kern)
    replace=which(segim!=0) #put back non-dilated segments
    segim_new[replace]=segim[replace] #put back non-dilated segments
    rm(replace)
  }
  mode(segim_new)='integer'
  
  rm(segim)
  
  if(!is.null(image)){
    segim_new[bmask$trues()]=0
    image[bmask$trues()]=NA
  }
  
  if(stats & !is.null(image)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=adacs_SegimStats(image=image, segim=segim_new, bmask=bmask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  if(plot & !is.null(image)){
    adacs_SegimPlot(image=image, segim=segim_new, bmask=bmask, sky=sky, ...)
  }
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - profoundMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return(invisible(list(segim=segim_new, segstats=segstats, header=header, call=call)))
}
adacs_MakeSegimDilateBitMatrix=function(bobjects=bobjects, bmask=NULL, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){

  if(verbose){message(' - Running MakeSegimDilate:')}
  timestart = proc.time()[3]
  
  call=match.call()
  
  kern = .makeBrush(size, shape=shape)
  
  if(verbose){message(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(is.null(expand) | length(expand)==0){
    
    segstats=NULL
    
    return(invisible(list(bobjects=bobjects, call=call)))
  }
  
  if(expand[1]=='all'){
    bobjects$dilatesparse(kern)
  }else{
    stop('expand must be \'all\' for MakeSegimDilateBitMatrix', call. = FALSE)
  }
  
  if(verbose){message(paste(" - profoundMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return(invisible(list(bobjects=bobjects, call=call)))
}

adacs_SegimStats=function(image=NULL, segim=NULL, bmask=NULL, sky=NULL, skyRMS=NULL, magzero=0, gain=NULL, pixscale=1, header=NULL, sortcol='segID', decreasing=FALSE, rotstats=FALSE, boundstats=FALSE, offset=1, cor_err_func=NULL, app_diam=1){
  
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
  
  bmask$maskNaN(image)
  
  #Set masked things to NA, to be safe:
  
  image[bmask$trues()]=NA
  
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
    
    mask = matrix(0L,xlen,ylen)
    mask[bmask$trues()]=1L
      outer_mask=segim_inner>0 & (mask[2:(xlen-1)+1,2:(ylen-1)]==1 | mask[2:(xlen-1)-1,2:(ylen-1)]==1 | mask[2:(xlen-1),2:(ylen-1)+1]==1 | mask[2:(xlen-1),2:(ylen-1)-1]==1)
      segim_mask=segim_edge
      segim_mask[outer_mask==0]=0
      tab_mask=tabulate(segim_mask)
      tab_mask=c(tab_mask,rep(0,max(segID)-length(tab_mask)))
      tab_mask=cbind(1:max(segID),tab_mask)
      Nmask=tab_mask[match(segID,tab_mask[,1]),2]
      
      rm(segim_inner)
      rm(tab_mask)
      
    
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

adacs_SegimPlot=function(image=NULL, segim=NULL,  bmask=NULL, sky=NULL, header=NULL, col=rainbow(max(segim), end=2/3), profound=NULL, ...){
  if(!is.null(image)){
    if(class(image)=='profound'){
      if(is.null(segim)){segim=image$segim}
      if(is.null(bmask)){bmask=image$bmask}
      if(is.null(sky)){sky=image$sky}
      if(is.null(header)){header=image$header}
      image=image$image
      if(is.null(image)){stop('Need image in profound object to be non-Null')}
    }
  }
  if(!is.null(profound)){
    if(class(profound) != 'profound'){
      stop('Class of profound input must be of type \'profound\'')
    }
    if(is.null(image)){image=profound$image}
    if(is.null(image)){stop('Need image in profound object to be non-Null')}
    if(is.null(segim)){segim=profound$segim}
    if(is.null(bmask)){bmask=profound$bmask}
    if(is.null(sky)){sky=profound$sky}
    if(is.null(header)){header=profound$header}
  }
  if(!is.null(image)){
    if(any(names(image)=='imDat') & is.null(header)){
      header=image$hdr
      image=image$imDat
    }else if(any(names(image)=='imDat') & !is.null(header)){
      image=image$imDat
    }
    if(any(names(image)=='dat') & is.null(header)){
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }else if(any(names(image)=='dat') & !is.null(header)){
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & is.null(header)){
      header=image$header
      image=image$image
    }else if(any(names(image)=='image') & !is.null(header)){
      image=image$image
    }
  }
  
  if(!is.null(sky)){
    image=image-sky
  }
  
  segim[is.na(segim)]=0L
  
  if(is.null(header)){header=NULL}
  if(is.null(header)){
    temp=magimage(image, ...)
  }else{
    temp=magimageWCS(image, header=header, ...)
  }
  if(min(segim,na.rm=doRMNA)!=0){segim=segim-min(segim,na.rm=doRMNA)}
  segvec=which(tabulate(segim)>0)
  for(i in segvec){
    z=segim==i
    z=z[ceiling(temp$x), ceiling(temp$y)]
    contour(temp$x,temp$y,z,add=T,col=col[i],zlim=c(0,1),drawlabels=FALSE,nlevels=1)
  }
  if(!is.null(bmask)){
    mask=matrix(0L,bmask$nrow(), bmask$ncol())
    bmask$copyTo(mask)
    magimage(mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))
    rm(mask)
  }
}