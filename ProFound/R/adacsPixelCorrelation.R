adacs_PixelCorrelation=function(image=NULL, bobjects=NULL, bmask=NULL, sky=0, skyRMS=1, lag=c(1:9,1:9*10,1:9*100,1:9*1e3,1:9*1e4), fft=TRUE, plot=FALSE, ylim=c(-1,1), log='x', grid=TRUE, ...){
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  
  lag=lag[lag<xlen & lag<ylen]
  
  image=(image-sky)/skyRMS
  
  flag_neg=image<0
  flag_pos=image>0
  
  if(!is.null(bobjects)){
    image[bobjects$trues()]=NA
  }
  
  image[bmask$trues()]=NA
  
  image_neg=image
  image_neg[flag_neg==FALSE]=NA #image is left with negative pixels only
  image_pos=image
  image_pos[flag_pos==FALSE]=NA #image is left with positive pixels only
  
  
  corx=cory=corx_neg=cory_neg=corx_pos=cory_pos=relsdx=relsdy=cor_err=rep(0,length(lag))
  
  for(i in 1:length(lag)){
    laguse=lag[i]
    corx[i]=cor(as.numeric(image[1:(xlen-laguse),]), as.numeric(image[1:(xlen-laguse)+laguse,]), use="na.or.complete")
    cory[i]=cor(as.numeric(image[,1:(ylen-laguse)]), as.numeric(image[,1:(ylen-laguse)+laguse]), use="na.or.complete")
    corx_neg[i]=cor(as.numeric(image_neg[1:(xlen-laguse),]), as.numeric(image[1:(xlen-laguse)+laguse,]), use="na.or.complete")
    cory_neg[i]=cor(as.numeric(image_neg[,1:(ylen-laguse)]), as.numeric(image[,1:(ylen-laguse)+laguse]), use="na.or.complete")
    corx_pos[i]=cor(as.numeric(image_pos[1:(xlen-laguse),]), as.numeric(image[1:(xlen-laguse)+laguse,]), use="na.or.complete")
    cory_pos[i]=cor(as.numeric(image_pos[,1:(ylen-laguse)]), as.numeric(image[,1:(ylen-laguse)+laguse]), use="na.or.complete")
    relsdx[i]=sd(as.numeric(image[1:(xlen-laguse),])-as.numeric(image[1:(xlen-laguse)+laguse,]), na.rm=doRMNA)/sqrt(2)
    relsdy[i]=sd(as.numeric(image[,1:(ylen-laguse)])-as.numeric(image[,1:(ylen-laguse)+laguse]), na.rm=doRMNA)/sqrt(2)
  }
  
  output_cortab=data.frame(lag=lag, corx=corx, cory=cory, corx_neg=corx_neg, cory_neg=cory_neg, corx_pos=corx_pos, cory_pos=cory_pos, corx_diff=corx_pos-corx_neg, cory_diff=cory_pos-cory_neg, relsdx=relsdx, relsdy=relsdy)
  
  tempapprox=approxfun(c(0,output_cortab$lag), abs(c(1,output_cortab$corx*output_cortab$cory)*c(0,output_cortab$lag)*pi*2), yleft=0, yright=0)
  
  for(i in 1:length(lag)){
    cor_err[i]=integrate(tempapprox,0,lag[i], stop.on.error=FALSE)$value/(pi*lag[i]^2)
  }
  cor_err[lag==1]=1
  cor_err_func=approxfun(pi*lag^2,cor_err, yleft=1, yright=1)
  
  if(fft){
    xlenpad=xlen+xlen%%2
    ylenpad=ylen+ylen%%2
    
    if(!is.null(bobjects)){
      image[bobjects$trues()]=rnorm(length(which(bobjects$trues())))
    }
    
    if(!is.null(mask)){
      image[bmask$trues()]=rnorm(length(bmask$trues()))
    }
    
    centre=matrix(c(rep(c(-1,1),ylenpad/2), rep(c(1,-1),ylenpad/2)),xlenpad,ylenpad)[1:xlen,1:ylen]
    FFT_Mod=Mod(fft(z=image*centre))
    dimx=dim(FFT_Mod)[1]
    dimy=dim(FFT_Mod)[2]
    if(dimx %% 2 == 0){
      x=c(-(dimx/2+0.5):(dimx/2-0.5))
    }else{
      x=c(-floor(dimx/2):floor(dimx/2))
    }
    if(dimy %% 2 == 0){
      y=c(-(dimy/2+0.5):(dimy/2-0.5))
    }else{
      y=c(-floor(dimy/2):floor(dimy/2))
    }
    output_FFT=list(x=x, y=y, z=FFT_Mod)
  }else{
    output_FFT=NULL
    image=NULL
  }
  
  if(plot){
    magplot(output_cortab[,c('lag','corx')], xlim=c(1,max(lag)), xlab='Pixel Lag', ylab='Correlation', type='l', col='blue', ylim=ylim, log=log, grid=grid, ...)
    lines(output_cortab[,c('lag','cory')], col='red')
    lines(output_cortab[,c('lag','corx_diff')], col='blue', lty=2)
    lines(output_cortab[,c('lag','cory_diff')], col='red', lty=2)
    lines(output_cortab[,c('lag','relsdx')], col='blue', lty=3)
    lines(output_cortab[,c('lag','relsdy')], col='red', lty=3)
    legend('topright', legend=c('x-cor','y-cor','x-cor-diff','y-cor-diff','x-rel-sd','y-rel-sd'), col=c('blue','red'), lty=c(1,1,2,2,3,3), bg='white')
  }
  
  invisible(list(cortab=output_cortab, fft=output_FFT, image_sky=image, cor_err_func=cor_err_func))
}
