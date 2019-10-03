profoundProFoundADACSInPlace=function(image=NULL, segim=NULL, objects=NULL, mask=NULL, skycut=1, pixcut=3, tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE, SBlim=NULL, size=5, shape='disc', iters=6,
                                      threshold=1.05, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL, redosegim=FALSE, redosky=TRUE, redoskysize=21, box=c(101,101), grid=box,
                                      type='bicubic', skytype='median', skyRMStype='quanlo', roughpedestal=FALSE, sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, iterskyloc=TRUE, deblend=FALSE, df=3, radtrunc=2, iterative=FALSE, doclip=TRUE,
                                      shiftloc = FALSE, paddim = TRUE, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, nearstats=boundstats, groupstats=boundstats, group=NULL, groupby='segim_orig',
                                      offset=1, haralickstats=FALSE, sortcol="segID", decreasing=FALSE, lowmemory=FALSE, keepim=TRUE, watershed='ProFound', pixelcov=FALSE, deblendtype='fit', psf=NULL, fluxweight='sum',
                                      convtype = 'brute', convmode = 'extended', fluxtype='Raw', app_diam=1, Ndeblendlim=Inf, ...){
  initialiseGlobals(doclip)
  if(verbose){message('Running ProFound:')}
  timestart=proc.time()[3]
  
  call=match.call()
  
  if(length(box)==1){
    box=rep(box,2)
    if(missing(grid)){grid=box}
    if(missing(boxadd)){boxadd=box/2}
    if(missing(skypixmin)){skypixmin=prod(box)/2}
  }
  if(length(grid)==1){
    grid=rep(grid,2)
  }
  if(length(boxadd)==1){
    boxadd=rep(boxadd,2)
  }
  
  # Ensure dimensions of box are odd numbers
  if (box[1]%%2 == 0) {
    box = c(box[1]+1, box[2])
  }
  if (box[2]%%2 == 0) {
    box = c(box[1], box[2]+1)
  }
  # Ensure dimensions of boxadd are even numbers (so that box + boxadd stays odd)
  if (boxadd[1]%%2 == 1) {
    boxadd = c(boxadd[1]+1, boxadd[2])
  }
  if (box[2]%%2 == 1) {
    boxadd = c(boxadd[1], boxadd[2]+1)
  }
  
  fluxtype=tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    if(verbose){message('Using raw flux units')}
    fluxscale=1
  }else if (fluxtype=='jansky'){
    if(verbose){message('Using Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero-8.9))
  }else{
    stop('fluxtype must be Jansky / Raw!')
  }
  
  #Split out image and header parts of input:
  
  if(!is.null(image)){
    if(any(names(image)=='imDat') & is.null(header)){
      if(verbose){message('Supplied image contains image and header components')}
      header=image$hdr
      image=image$imDat
    }else if(any(names(image)=='imDat') & !is.null(header)){
      if(verbose){message('Supplied image contains image and header but using specified header')}
      image=image$imDat
    }
    if(any(names(image)=='dat') & is.null(header)){
      if(verbose){message('Supplied image contains image and header components')}
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }else if(any(names(image)=='dat') & !is.null(header)){
      if(verbose){message('Supplied image contains image and header but using specified header')}
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & is.null(header)){
      if(verbose){message('Supplied image contains image and header components')}
      header=image$header
      image=image$image
    }else if(any(names(image)=='image') & !is.null(header)){
      if(verbose){message('Supplied image contains image and header but using specified header')}
      image=image$image
    }
  }else{
    stop('Missing image - this is a required input!')
  }
  
  if(verbose){message(paste('Supplied image is',dim(image)[1],'x',dim(image)[2],'pixels'))}
  
  #Treat image NAs as masked regions:
  
  badpix=NULL
  if(!is.null(mask)){
    mask=mask*1L #Looks silly, but this ensures a logical mask becomes integer.
    if(length(mask)==1 & !is.na(mask[1])){
      maskflag=mask
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[image==maskflag]=1L
    }
    if(anyNA(image)){
      badpix=which(is.na(image))
      mask[badpix]=1L
      image[badpix]=0
    }
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      badpix=which(is.na(image))
      mask[badpix]=1L
      image[badpix]=0
    }
  }
  
  #if(!is.null(segim) & !is.null(mask)){
  #  segim=segim*(1-mask) #I don't think we actually need this
  #}
  
  #Get the pixel scale, if possible and not provided:
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste('Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel'))}
  }else{
    if(verbose){message(paste('Using suggested pixel scale:',round(pixscale,3),'asec/pixel'))}
  }
  
  skyarea=prod(dim(image))*pixscale^2/(3600^2)
  if(verbose){message(paste('Supplied image is',round(dim(image)[1]*pixscale/60,3),'x',round(dim(image)[2]*pixscale/60,3),'amin, ', round(skyarea,3),'deg-sq'))}
  
  if(is.null(objects)){
    if(!is.null(segim)){
      objects=matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[]=as.logical(segim)
    }
  }else{
    objects=objects*1 #Looks silly, but this ensures a logical mask becomes integer.
  }
  
  # Create scratch matrices
  #scratch <- list(scratchN1=matrix(0.0,box[1],box[2]), scratchN2=matrix(0.0,box[1],box[2]),scratchI1=matrix(FALSE,box[1],box[2]), scratchI2=matrix(FALSE,box[1],box[2]))
  scratch <- list(scratchN1=matrix(0.0,box[1],box[2]), 
                  scratchN2=matrix(0.0,box[1],box[2]),
                  scratchI1=matrix(FALSE,box[1],box[2]),
                  scratchI2=matrix(FALSE,box[1],box[2]),
                  scratchSKY=matrix(0.0,dim(image)[1],dim(image)[2]),
                  scratchSKYRMS=matrix(0.0,dim(image)[1],dim(image)[2]))
  
  #Check for user provided sky, and compute if missing:
  
  hassky=!is.null(sky)
  hasskyRMS=!is.null(skyRMS)
  
  if((hassky==FALSE | hasskyRMS==FALSE) & is.null(segim)){
    if(verbose){message(paste('Making initial sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    roughsky=profoundMakeSkyGridADACSInPlace(image=image, objects=objects, mask=mask, box=box, grid=grid, boxadd=boxadd,
                                             type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, 
                                             skypixmin=skypixmin, boxiters=0, doclip=doclip, shiftloc=shiftloc, paddim=paddim,scratch=scratch,final=FALSE)
    if(roughpedestal){
      roughsky$sky=median(roughsky$sky,na.rm=doRMNA)
      roughsky$skyRMS=median(roughsky$skyRMS,na.rm=doRMNA)
    }
    if(hassky==FALSE){
      sky=roughsky$sky
      if(verbose){message(' - Sky statistics :')}
      if(verbose){print(summary(as.numeric(sky)))}
    }
    if(hasskyRMS==FALSE){
      skyRMS=roughsky$skyRMS
      if(verbose){message(' - Sky-RMS statistics :')}
      if(verbose){print(summary(as.numeric(skyRMS)))}
    }
      rm(roughsky)
  }else{
    if(verbose){message("Skipping making initial sky map - User provided sky and sky RMS, or user provided segim")}
  }
  
  #Make the initial segmentation map, if not provided.
  
  if(is.null(segim)){
    if(verbose){message(paste('Making initial segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
    segim=profoundMakeSegimADACSInPlace(image=image, objects=objects, mask=mask, sky=sky, skyRMS=skyRMS,
                            tolerance=tolerance, ext=ext, reltol=reltol, cliptol=cliptol, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,
                            magzero=magzero, pixscale=pixscale, verbose=verbose, watershed=watershed, plot=FALSE, stats=FALSE)
    objects=segim$objects
    segim=segim$segim
  }else{
    redosegim=FALSE
    if(verbose){message("Skipping making an initial segmentation image - User provided segim")}
  }
  if(any(segim>0)){
    if((hassky==FALSE | hasskyRMS==FALSE)){
      if (TRUE) {
      if(redosky){
        if(verbose){message(paste('Doing initial aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
        objects_redo=profoundMakeSegimDilateADACSInPlace(segim=objects, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      }else{
        objects_redo=objects
      }
      if(verbose){message(paste('Making better sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      bettersky=profoundMakeSkyGridADACSInPlace(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, boxadd=boxadd,
                                                type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                                                skypixmin=skypixmin, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim, scratch=scratch,final=TRUE)
      if(hassky==FALSE){
        sky=bettersky$sky
        if(verbose){message(' - Sky statistics :')}
        if(verbose){print(summary(as.numeric(sky)))}
      }
      if(hasskyRMS==FALSE){
        skyRMS=bettersky$skyRMS
        if(verbose){message(' - Sky-RMS statistics :')}
        if(verbose){print(summary(as.numeric(skyRMS)))}
      }
      if(redosegim){
        if(verbose){message(paste('Making better segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
        imagescale=(image-sky)/skyRMS
        imagescale[!is.finite(imagescale)]=0
        if(!is.null(SBlim) & !missing(magzero)){
          imagescale[imagescale<skycut | sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
        }else{
          imagescale[imagescale<skycut]=0
        }
        if(!is.null(mask)){
          imagescale[mask!=0]=0
        }
        segim[imagescale==0]=0
        objects[segim==0]=0
      }
      }
    }else{
      if(verbose){message("Skipping making better sky map - User provided sky and sky RMS")}
    }
    
      if (TRUE) {
    if(iters>0 | iterskyloc){
      if(verbose){message(paste('Calculating initial segstats -',round(proc.time()[3]-timestart,3),'sec'))}
      segstats=.profoundFluxCalcMin(image=image, segim=segim, mask=mask)
      skystats=.profoundFluxCalcMin(image=sky, segim=segim, mask=mask)
      skystats=skystats$flux/skystats$N100
      skymed=median(skystats, na.rm=doRMNA)
      origfrac=segstats$flux - (skystats*segstats$N100)
      
      if(iterskyloc){
        localadd=1
      }else{
        localadd=0
      }
      
      #compmat=matrix(0,nrow = dim(segstats)[1], ncol = iters+1+localadd)
      #Nmat=compmat
      #compmat[,1]=segstats[,'flux']
      #Nmat[,1]=segstats[,'N100']
      #flux_old=segstats[,'flux']
      #N100_old=segstats[,'N100']
      
      #segim_array=array(0L, dim=c(dim(segim),iters+1+localadd))
      #segim_array[,,1]=segim
      
      segim_orig=segim
      expand_segID=segstats[,'segID']
      SBlast=rep(Inf,length(expand_segID))
      selseg=rep(0,length(expand_segID))
      
      if(verbose){message('Doing dilations:')}
        
      for(i in 1:(iters)){
        if(verbose){message(paste('Iteration',i,'of',iters,'-',round(proc.time()[3]-timestart,3),'sec'))}
        segim_new=profoundMakeSegimDilateADACSInPlace(segim=segim, expand=expand_segID, size=size, shape=shape, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$segim
        segstats_new=.profoundFluxCalcMin(image=image, segim=segim_new, mask=mask)
        SBnew=(segstats_new$flux - segstats$flux) / (segstats_new$N100 - segstats$N100)
        fluxgrowth = (segstats_new$flux - skystats * segstats_new$N100) / (segstats$flux - skystats * segstats$N100) #account for flux growth
        skyfrac = abs(((skystats-skymed) * (segstats_new$N100-segstats$N100)) / (segstats_new$flux - segstats$flux)) #account for sky growth
        expand_segID=segstats[which(segstats_new$flux>0 & fluxgrowth > threshold & SBnew < (SBlast/threshold) & skyfrac < 0.5 & selseg==(i-1)),'segID']
        expand_segID=expand_segID[is.finite(expand_segID)]
        if(length(expand_segID)==0){break}
        updateID=which(segstats$segID %in% expand_segID)
        selseg[updateID] = i
        segstats[updateID,] = segstats_new[updateID,]
        SBlast = SBnew
        if('fastmatch' %in% .packages()){ #dilate segments that pass tests
          selpix = which(fastmatch::fmatch(segim_new, expand_segID, nomatch = 0L) > 0) 
        }else{
          selpix = which(segim_new %in% expand_segID)
        }
        segim[selpix]=segim_new[selpix]
      }
      
      if(iterskyloc){
        segim_skyloc=profoundMakeSegimDilateADACSInPlace(segim=segim, size=size, shape=shape, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$segim
        segstats_new=.profoundFluxCalcMin(image=image, segim=segim_skyloc, mask=mask)
        skyseg_mean=(segstats_new$flux-segstats$flux)/(segstats_new$N100-segstats$N100)
        skyseg_mean[!is.finite(skyseg_mean)]=0
      }else{
        skyseg_mean=NA
      }
      
      objects=matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[]=as.logical(segim)
      
      origfrac = origfrac / (segstats$flux - (skystats * segstats$N100))
    }else{
      if(verbose){message('Iters set to 0 - keeping segim un-dilated')}
      segim_orig=segim
      selseg=0
      origfrac=1
    }
    
    if(redosky){
      if(redoskysize %% 2 == 0){redoskysize=redoskysize+1}
      if(verbose){message(paste('Doing final aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      objects_redo=profoundMakeSegimDilateADACSInPlace(segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      if(verbose){message(paste('Making final sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      sky=profoundMakeSkyGridADACSInPlace(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, boxadd=boxadd,
                                          type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                                          skypixmin=skypixmin, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim, scratch=scratch,final=TRUE)
      skyRMS=sky$skyRMS
      sky=sky$sky
      if(verbose){message(' - Sky statistics :')}
      if(verbose){print(summary(as.numeric(sky)))}
      if(verbose){message(' - Sky-RMS statistics :')}
      if(verbose){print(summary(as.numeric(skyRMS)))}
    }else{
      if(verbose){message("Skipping making final sky map - redosky set to FALSE")}
      objects_redo=NULL
    }
    
    Norig=tabulate(segim_orig)
    
    if(pixelcov){
      if(verbose){message(paste('Calculating pixel covariance -',round(proc.time()[3]-timestart,3),'sec'))}
      cor_err_func=profoundPixelCorrelation(image=image, objects=objects, mask=mask, sky=sky, skyRMS=skyRMS, fft=FALSE, lag=apply(expand.grid(c(1,2,4),c(1,10,100,1000,1e4)),MARGIN=1,FUN=prod))$cor_err_func
    }else{
      cor_err_func=NULL
    }
    
    if(lowmemory){
      image=image-sky
      sky=0
      skyRMS=0
      segim_orig=NULL
      objects=NULL
      objects_redo=NULL
    }
    
    if(stats & !is.null(image)){
      if(verbose){message(paste('Calculating final segstats for',length(which(tabulate(segim)>0)),'objects -',round(proc.time()[3]-timestart,3),'sec'))}
      if(verbose){message(paste(' - magzero =', round(magzero,3)))}
      if(verbose){
        if(is.null(gain)){
          message(paste(' - gain = NULL (ignored)'))
        }else{
          message(paste(' - gain =', round(gain,3)))
        }
      }
      if(verbose){message(paste(' - pixscale =', round(pixscale,3)))}
      if(verbose){message(paste(' - rotstats =', rotstats))}
      if(verbose){message(paste(' - boundstats =', boundstats))}
      segstats=profoundSegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS,
                                  magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset, cor_err_func=cor_err_func, app_diam=app_diam)
      segstats=cbind(segstats, iter=selseg, origfrac=origfrac, Norig=Norig[segstats$segID], skyseg_mean=skyseg_mean)
      segstats=cbind(segstats, flag_keep=segstats$origfrac>= median(segstats$origfrac[segstats$iter==iters]) | segstats$iter<iters)
    }else{
      if(verbose){message("Skipping segmentation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    if(nearstats){
      near=profoundSegimNear(segim=segim, offset=offset)
    }else{
      near=NULL
    }
    
    if(deblend){
      groupstats=TRUE
    }
    
    if(groupstats){
      if(verbose){message(paste(' - groupstats = TRUE - ',round(proc.time()[3]-timestart,3),'sec'))}
      if(groupby=='segim'){
        if(is.null(group)){
          group=profoundSegimGroup(segim)
        }
      }else if(groupby=='segim_orig'){
        if(is.null(group)){
          #message(round(proc.time()[3]-timestart,3))
          group=profoundSegimGroup(segim_orig)
          if(any(group$groupsegID$Ngroup>1)){
            #message(round(proc.time()[3]-timestart,3))
            group$groupim=profoundSegimKeep(segim=segim, segID_merge=group$groupsegID[group$groupsegID$Ngroup>1,'segID'])
            #message(round(proc.time()[3]-timestart,3))
            group$groupsegID$Npix=tabulate(group$groupim)[group$groupsegID$groupID]
          }
        }
      }else{
        stop('Non legal groupby option, must be segim or segim_orig!')
      }

      if(stats & !is.null(image) & !is.null(group)){
        groupstats=profoundSegimStats(image=image, segim=group$groupim, mask=mask, sky=sky, skyRMS=skyRMS,
                                      magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset, cor_err_func=cor_err_func, app_diam=app_diam)
        colnames(groupstats)[1]='groupID'
      }else{
        groupstats=NULL
      }
    }else{
      if(verbose){message(' - groupstats = FALSE')}
      group=NULL
      groupstats=NULL
    }
    
    if(deblend & stats & !is.null(image) & any(group$groupsegID$Ngroup>1)){
      if(verbose){message(paste(' - deblend = TRUE - ',round(proc.time()[3]-timestart,3),'sec'))}
      tempblend=profoundFluxDeblend(image=image-sky, segim=segim, segstats=segstats,
                                    groupim=group$groupim, groupsegID=group$groupsegID, magzero=magzero, df=df, radtrunc=radtrunc, iterative=iterative, doallstats=TRUE, deblendtype=deblendtype, psf=psf,
                                    fluxweight=fluxweight, convtype=convtype, convmode=convmode, Ndeblendlim = Ndeblendlim)
      if(!is.null(tempblend)){
        segstats=cbind(segstats,tempblend[,-2])
      }
    }else{
      if(verbose){message(' - deblend = FALSE')}
    }
    
    if(haralickstats){
      if(requireNamespace("EBImage", quietly = TRUE)){
        scale=10^(0.4*(30-magzero))
        temphara=(image-sky)*scale
        if(!is.null(mask)){
          temphara[mask!=0]=0
        }
        temphara[!is.finite(temphara)]=0
        haralick=as.data.frame(EBImage::computeFeatures.haralick(segim,temphara))
        haralick=haralick[segstats$segID,]
      }else{
        if(verbose){
          message('The EBImage package is needed to compute Haralick statistics.')
          haralick=NULL
        }
      }
    }else{
      haralick=NULL
    }
    
    if(plot){
      if(verbose){message(paste('Plotting segments -',round(proc.time()[3]-timestart,3),'sec'))}
      if(any(is.finite(sky))){
        profoundSegimPlot(image=image-sky, segim=segim, mask=mask, header=header, ...)
      }else{
        profoundSegimPlot(image=image, segim=segim, mask=mask, header=header, ...)
      }
    }else{
      if(verbose){message("Skipping segmentation plot - plot = FALSE")}
    }
    
    if(is.null(SBlim)){
      SBlim=NULL
    }else if(is.numeric(SBlim)){
      SBlimtemp=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
      SBlimtemp=matrix(SBlimtemp,dim(skyRMS)[1],dim(skyRMS)[2])
      SBlimtemp[which(SBlimtemp>SBlim)]=SBlim
      SBlim=SBlimtemp
    }else if(SBlim[1]=='get' & skycut> -Inf){
      SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
    }
    
    if(is.null(header)){header=NULL}
    if(keepim==FALSE){image=NULL; mask=NULL}
    if(is.null(mask)){mask=NULL}
    if(!is.null(badpix)){image[badpix]=NA}
    row.names(segstats)=NULL
    
    segstats[,grep('flux',colnames(segstats))]=fluxscale*segstats[,grep('flux',colnames(segstats))]
    
    if(verbose){message(paste('ProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    output=list(segim=segim, segim_orig=segim_orig, objects=objects, objects_redo=objects_redo, sky=sky, skyRMS=skyRMS, image=image, mask=mask, segstats=segstats, Nseg=dim(segstats)[1], near=near, group=group, groupstats=groupstats, haralick=haralick, header=header, SBlim=SBlim, magzero=magzero, dim=dim(segim), pixscale=pixscale, skyarea=skyarea, gain=gain, call=call, date=date(), time=proc.time()[3]-timestart, ProFound.version=packageVersion('ProFound'), R.version=R.version)
      }  else {
        # just to keep it happy
      output = list(image=image)
    }
    }else{
    if(is.null(header)){header=NULL}
    if(keepim==FALSE){image=NULL; mask=NULL}
    if(is.null(mask)){mask=NULL}
    if(!is.null(badpix)){image[badpix]=NA}
    if(verbose){message('No objects in segmentation map - skipping dilations and CoG')}
    if(verbose){message(paste('ProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    output=list(segim=NULL, segim_orig=NULL, objects=NULL, objects_redo=NULL, sky=sky, skyRMS=skyRMS, image=image, mask=mask, segstats=NULL, Nseg=0, near=NULL, group=NULL, groupstats=NULL, haralick=NULL, header=header, SBlim=NULL,  magzero=magzero, dim=dim(segim), pixscale=pixscale, skyarea=skyarea, gain=gain, call=call, date=date(), time=proc.time()[3]-timestart, ProFound.version=packageVersion('ProFound'), R.version=R.version)
  }
  class(output)='profound'
  invisible(output)
}
