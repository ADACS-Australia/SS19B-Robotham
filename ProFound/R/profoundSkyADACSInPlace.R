profoundSkyEstLocADACSInPlace=function(image=NULL, objects=NULL, mask=NULL, loc=dim(image)/2, box=c(100,100), skytype='median', skyRMStype='quanlo', sigmasel=1,
                                       skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, plot=FALSE, scratch=NULL, final=FALSE,...){
  if(!is.null(objects) | !is.null(mask)){
    select = adacsFindSkyCellValuesC(image, objects, mask, loc[1], loc[2], box[1], box[2], boxadd[1], boxadd[2], skypixmin, boxiters)
    if (plot) {
      newbox = adacsFindSkyCellValuesBoxC(image, objects, mask, loc[1], loc[2], box[1], box[2], boxadd[1], boxadd[2], skypixmin, boxiters)
      box = c(newbox[1],newbox[2])
    }
  }else{
    #print("NO OBJECTS AT ALL")
    select = scratch[['scratchN1']]
    magcutoutADACSInPlace(image,loc=loc,box=box,oimage=select)
  }
  if(plot){
    image=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image
    imout=magimage(image, ...)
    if(!is.null(mask)){
      contour(x=imout$x, y=imout$y, magcutout(mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image, add=T, col='red', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
    if(!is.null(objects)){
      contour(x=imout$x, y=imout$y, magcutout(objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image, add=T, col='blue', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
  }
  # object "select" may be either a matrix or a vector at this point
  if(doclip){
    #suppressWarnings({
      #clip=magclipADACS(select, sigmasel=sigmasel, estimate = LO)$x
      if (final)
        clip = adacsmagclipV(select,AUTO,5,sigmasel,LO)
      else
        clip = adacsmagclip(select,AUTO,5,sigmasel,LO)
      #if (!compareNA(clip,clip2)) {
      #  print("problem in adacsmagclip")
      #}
    #})
  }else{
    clip=select
  }
  
  if(length(clip)==1){
    if(is.na(clip)){
      return(invisible(list(val=c(NA, NA), clip=NA)))
    }
  }
  
  if(skytype=='median'){
    if('Rfast' %in% .packages()){
      skyloc=try(Rfast::med(clip, na.rm=doRMNA), silent=TRUE)
      if(class(skyloc)=='try-error'){skyloc=NA}
    }else{
      # RS try calling this from C
      skyloc=stats::median(clip, na.rm=doRMNA)
      #skylocC=get_median(clip,stats::median)
    }
  }else if(skytype=='mean'){
    skyloc=mean(clip, na.rm=doRMNA)
  }else if(skytype=='mode'){
    temp=density(clip, na.rm=doRMNA)
    skyloc=temp$x[which.max(temp$y)]
  }
  
  if(skyRMStype=='quanlo'){
    temp=clip-skyloc
    temp=temp[temp<0]
    skyRMSloc=abs(quantile(temp, pnorm(-sigmasel)*2,na.rm=doRMNA))/sigmasel
  }else if(skyRMStype=='quanhi'){
    temp=clip-skyloc
    temp=temp[temp>0]
    skyRMSloc=abs(quantile(temp, (pnorm(sigmasel)-0.5)*2,na.rm=doRMNA))/sigmasel
  }else if(skyRMStype=='quanboth'){
    temp=clip-skyloc
    templo=temp[temp<0]
    temphi=temp[temp>0]
    skyRMSloclo=abs(quantile(templo, pnorm(-sigmasel)*2,na.rm=doRMNA))/sigmasel
    skyRMSlochi=abs(quantile(temphi, (pnorm(sigmasel)-0.5)*2,na.rm=doRMNA))/sigmasel
    skyRMSloc=(skyRMSloclo+skyRMSlochi)/2
  }else if(skyRMStype=='sd'){
    skyRMSloc=sqrt(.varwt(clip, wt=1, xcen=skyloc))
  }
  
  return(invisible(list(val=c(skyloc, skyRMSloc), clip=clip)))
}

profoundMakeSkyGridADACSInPlaceOldSpline=function(image=NULL, objects=NULL, mask=NULL, box=c(100,100), grid=box, type='bicubic', skytype='median', skyRMStype='quanlo', sigmasel=1,
                                         skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, cores=1,
                                         scratch=NULL, final=FALSE){
  if(!requireNamespace("akima", quietly = TRUE)){
    if(type=='bicubic'){
      stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
    }
    if(type=='bilinear'){
      useakima=FALSE
    }
  }else{
    useakima=TRUE
  }
  print(paste("type=",type," useakima=",useakima))
  # box MUST NOT be larger than the input image
  if(box[1]>dim(image)[1]){box[1]=dim(image)[1]}
  if(box[2]>dim(image)[2]){box[2]=dim(image)[2]}
  if(grid[1]>dim(image)[1]){grid[1]=dim(image)[1]}
  if(grid[2]>dim(image)[2]){grid[2]=dim(image)[2]}
  
  # tile over input image with tile size (grid) and no overlap
  # xseq,yseq give the centres of each tile
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)
  
  if(cores>1){
    registerDoParallel(cores=cores)
    i=NULL
    tempsky=foreach(i = 1:dim(tempgrid)[1], .combine='rbind')%dopar%{
      profoundSkyEstLocADACSInPlace(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                                    skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim, scratch=scratch, final=final)$val
    }
    tempsky=rbind(tempsky)
  }else{
    tempsky=matrix(0,dim(tempgrid)[1],2)
    for(i in 1:dim(tempgrid)[1]){
      tempsky[i,]=profoundSkyEstLocADACSInPlace(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                                                skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim, scratch=scratch, final=final)$val
    }
  }
  
  # Take these boxcar median values as anchors for akima splines to expand to cover all the original cells
  
  if (TRUE) {
  xseq=c(-grid[1]/2,xseq,max(xseq)+grid[1]/2)
  yseq=c(-grid[2]/2,yseq,max(yseq)+grid[2]/2)
  
  tempmat_sky=matrix(0,length(xseq),length(yseq))
  tempmat_sky[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,1]
  tempmat_sky[is.na(tempmat_sky)]= stats::median(tempmat_sky, na.rm = TRUE)
  
  tempmat_skyRMS=matrix(0,length(xseq),length(yseq))
  tempmat_skyRMS[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,2]
  tempmat_skyRMS[is.na(tempmat_skyRMS)]=stats::median(tempmat_skyRMS, na.rm = TRUE)
  
  xstart=min(3,dim(tempmat_sky)[1]-1)
  ystart=min(3,dim(tempmat_sky)[2]-1)
  xend=max(length(xseq)-2,2)
  yend=max(length(yseq)-2,2)
  
  tempmat_sky[1,]=tempmat_sky[2,]*2-tempmat_sky[xstart,]
  tempmat_sky[length(xseq),]=tempmat_sky[length(xseq)-1,]*2-tempmat_sky[xend,]
  tempmat_sky[,1]=tempmat_sky[,2]*2-tempmat_sky[,ystart]
  tempmat_sky[,length(yseq)]=tempmat_sky[,length(yseq)-1]*2-tempmat_sky[,yend]
  
  tempmat_skyRMS[1,]=tempmat_skyRMS[2,]*2-tempmat_skyRMS[xstart,]
  tempmat_skyRMS[length(xseq),]=tempmat_skyRMS[length(xseq)-1,]*2-tempmat_skyRMS[xend,]
  tempmat_skyRMS[,1]=tempmat_skyRMS[,2]*2-tempmat_skyRMS[,ystart]
  tempmat_skyRMS[,length(yseq)]=tempmat_skyRMS[,length(yseq)-1]*2-tempmat_skyRMS[,yend]
  
  if(dim(tempmat_sky)[1]>1){
    
    #expand out map here!! and then use akima::bilinear function
    
    bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
    bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
    
    if(type=='bilinear'){
      if(useakima){
        tempgrid=expand.grid(xseq, yseq)
        temp_bi_sky=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_sky),xo=bigridx, yo=bigridy)$z
        temp_bi_skyRMS=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_skyRMS),xo=bigridx, yo=bigridy)$z
      }else{
        temp_bi_sky=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_sky))
        temp_bi_skyRMS=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_skyRMS))
      }
    }else if(type=='bicubic'){
      xxx = dim(image)[1];
      yyy = dim(image)[2];
      cinterp = interpolateSkyGrid(xseq,yseq,tempmat_sky,xxx,yyy);
      temp_bi_sky=akima::bicubic(xseq, yseq, tempmat_sky, bigridx, bigridy)$z
      temp_bi_skyRMS=akima::bicubic(xseq, yseq, tempmat_skyRMS, bigridx, bigridy)$z
    }else{
      stop('type must be one of bilinear / bicubic !')
    }
    
    rm(bigridx)
    rm(bigridy)
  
    temp_bi_sky=matrix(temp_bi_sky, dim(image)[1], dim(image)[2])
    temp_bi_skyRMS=matrix(temp_bi_skyRMS, dim(image)[1], dim(image)[2])
  }else{
    temp_bi_sky=matrix(tempmat_sky[1,1], dim(image)[1], dim(image)[2])
    temp_bi_skyRMS=matrix(tempmat_skyRMS[1,1], dim(image)[1], dim(image)[2])
  }
  
  if(!is.null(mask)){
    temp_bi_sky[mask>0]=NA
    temp_bi_skyRMS[mask>0]=NA
  }
  
  invisible(list(sky=temp_bi_sky, skyRMS=temp_bi_skyRMS))
  }
}
profoundMakeSkyGridADACSInPlace=function(image=NULL, objects=NULL, mask=NULL, box=c(100,100), grid=box, type='bicubic', skytype='median', skyRMStype='quanlo', sigmasel=1,
                                         skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, cores=1,
                                         scratch=NULL, final=FALSE){
  if(!requireNamespace("akima", quietly = TRUE)){
    if(type=='bicubic'){
      stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
    }
    if(type=='bilinear'){
      useakima=FALSE
    }
  }else{
    useakima=TRUE
  }
  # box MUST NOT be larger than the input image
  if(box[1]>dim(image)[1]){box[1]=dim(image)[1]}
  if(box[2]>dim(image)[2]){box[2]=dim(image)[2]}
  if(grid[1]>dim(image)[1]){grid[1]=dim(image)[1]}
  if(grid[2]>dim(image)[2]){grid[2]=dim(image)[2]}
  
  # tile over input image with tile size (grid) and no overlap
  # xseq,yseq give the centres of each tile
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)
  
  if(cores>1){
    registerDoParallel(cores=cores)
    i=NULL
    tempsky=foreach(i = 1:dim(tempgrid)[1], .combine='rbind')%dopar%{
      profoundSkyEstLocADACSInPlace(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                                    skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim, scratch=scratch, final=final)$val
    }
    tempsky=rbind(tempsky)
  }else{
    tempsky=matrix(0,dim(tempgrid)[1],2)
    for(i in 1:dim(tempgrid)[1]){
      tempsky[i,]=profoundSkyEstLocADACSInPlace(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                                                skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim, scratch=scratch, final=final)$val
    }
  }
  
  # Take these boxcar median values as anchors for akima splines to expand to cover all the original cells
  
  if (TRUE) {
    xseq=c(-grid[1]/2,xseq,max(xseq)+grid[1]/2)
    yseq=c(-grid[2]/2,yseq,max(yseq)+grid[2]/2)
    
    tempmat_sky=matrix(0,length(xseq),length(yseq))
    tempmat_sky[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,1]
    tempmat_sky[is.na(tempmat_sky)]= stats::median(tempmat_sky, na.rm = TRUE)
    
    tempmat_skyRMS=matrix(0,length(xseq),length(yseq))
    tempmat_skyRMS[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,2]
    tempmat_skyRMS[is.na(tempmat_skyRMS)]=stats::median(tempmat_skyRMS, na.rm = TRUE)
    
    xstart=min(3,dim(tempmat_sky)[1]-1)
    ystart=min(3,dim(tempmat_sky)[2]-1)
    xend=max(length(xseq)-2,2)
    yend=max(length(yseq)-2,2)
    
    tempmat_sky[1,]=tempmat_sky[2,]*2-tempmat_sky[xstart,]
    tempmat_sky[length(xseq),]=tempmat_sky[length(xseq)-1,]*2-tempmat_sky[xend,]
    tempmat_sky[,1]=tempmat_sky[,2]*2-tempmat_sky[,ystart]
    tempmat_sky[,length(yseq)]=tempmat_sky[,length(yseq)-1]*2-tempmat_sky[,yend]
    
    tempmat_skyRMS[1,]=tempmat_skyRMS[2,]*2-tempmat_skyRMS[xstart,]
    tempmat_skyRMS[length(xseq),]=tempmat_skyRMS[length(xseq)-1,]*2-tempmat_skyRMS[xend,]
    tempmat_skyRMS[,1]=tempmat_skyRMS[,2]*2-tempmat_skyRMS[,ystart]
    tempmat_skyRMS[,length(yseq)]=tempmat_skyRMS[,length(yseq)-1]*2-tempmat_skyRMS[,yend]
    
    if(dim(tempmat_sky)[1]>1){
      
      old = FALSE
      if (old) {
        print("OLD SPLINE")
        #expand out map here!! and then use akima::bilinear function
        
        bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
        bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
        
        if(type=='bilinear'){
          if(useakima){
            tempgrid=expand.grid(xseq, yseq)
            temp_bi_sky=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_sky),xo=bigridx, yo=bigridy)$z
            temp_bi_skyRMS=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_skyRMS),xo=bigridx, yo=bigridy)$z
          }else{
            temp_bi_sky=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_sky))
            temp_bi_skyRMS=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_skyRMS))
          }
        }else if(type=='bicubic'){
          temp_bi_sky=akima::bicubic(xseq, yseq, tempmat_sky, bigridx, bigridy)$z
          temp_bi_skyRMS=akima::bicubic(xseq, yseq, tempmat_skyRMS, bigridx, bigridy)$z
        }else{
          stop('type must be one of bilinear / bicubic !')
        }
        
        rm(bigridx)
        rm(bigridy)
        
        temp_bi_sky=matrix(temp_bi_sky, dim(image)[1], dim(image)[2])
        temp_bi_skyRMS=matrix(temp_bi_skyRMS, dim(image)[1], dim(image)[2])
      } else {
      print("NEW SPLINE")
      #expand out map here!! and then use akima::bilinear function
      
      #bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
      #bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
      
      if(type=='bilinear'){
        if(useakima){
          tempgrid=expand.grid(xseq, yseq)
          temp_bi_sky=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_sky),xo=bigridx, yo=bigridy)$z
          temp_bi_skyRMS=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_skyRMS),xo=bigridx, yo=bigridy)$z
        }else{
          temp_bi_sky=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_sky))
          temp_bi_skyRMS=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_skyRMS))
        }
      }else if(type=='bicubic'){
        xxx = dim(image)[1];
        yyy = dim(image)[2];
        temp_bi_sky = scratch[['scratchSKY']]
        interpolateAkimaGrid(xseq,yseq,tempmat_sky,xxx,yyy,temp_bi_sky);
        temp_bi_skyRMS = scratch[['scratchSKYRMS']]
        interpolateAkimaGrid(xseq,yseq,tempmat_skyRMS,xxx,yyy,temp_bi_skyRMS);
      }else{
        stop('type must be one of bilinear / bicubic !')
      }
      }
    }else{
      temp_bi_sky=matrix(tempmat_sky[1,1], dim(image)[1], dim(image)[2])
      temp_bi_skyRMS=matrix(tempmat_skyRMS[1,1], dim(image)[1], dim(image)[2])
    }
    
    if(!is.null(mask)){
      temp_bi_sky[mask>0]=NA
      temp_bi_skyRMS[mask>0]=NA
    }
    
    invisible(list(sky=temp_bi_sky, skyRMS=temp_bi_skyRMS))
  }
}
