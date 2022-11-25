# This script fills up an OM object with variables (many of which are dummies) just to pass MSEtool checks
# for OM completeness
# Note that only OM@h, OM@Perr, and OM@cpars: Wt_age, M_ageArray, age95, Mat_age are used by SRA_scope

DHOM<-function(Name,h,Perr,nsim,ny,na,np,interval=2,a,b){
  
  OM<-new('OM')
  OM@Name <- Name
  OM@nsim<-nsim<-nsim
  OM@nyears<-ny
  OM@maxage<-na
  OM@proyears<-np
  OM@interval<-interval
  OM@M<-rep(0.2,2)
  OM@h=rep(h,2)
  OM@Perr<-rep(Perr,2)
  OM@Size_area_1<-rep(0.5,2)
  OM@Frac_area_1<-rep(0.5,2)
  OM@Prob_staying<-rep(0.5,2)
  OM@K<-rep(0.2,2)
  OM@Linf <-rep(30,2)
  OM@LenCV<-rep(0.05,2)
  OM@EffYears<-c(1,5,OM@nyears)
  OM@EffLower<-rep(1,3)
  OM@EffUpper<-rep(1.1,3)
  OM@L50<-rep(22,2)
  OM@L50_95<-rep(1.5,2)
  OM@L5<-rep(20.5,2)
  OM@LFS<-rep(25,2)
  OM@Vmaxlen<-rep(1,2)
  OM@SRrel<-1
  OM@Msd<-OM@Ksd<-OM@Linfsd<-OM@qcv<-OM@Esd<-rep(0,2)
  OM@pstar<-0.5
  OM@reps<-1
  OM@Spat_targ<-c(1,1)
  OM@qinc<-rep(0,2)
  OM@a<-a
  OM@b<-b
  OM@D<-c(0.2,2)
  OM@isRel<-F
  OM@CurrentYr<-2018
  OM@AC<-c(0,0)
  OM@R0<-rep(1E3,2)
  
  # Invent an OM with full observation error model for replacing
  temp<-new('OM',Albacore,Generic_Fleet,Imprecise_Unbiased,Perfect_Imp)
  OM<-Replace(OM,temp,Sub="Obs")
  OM<-Replace(OM,temp,Sub="Imp")
  
  OM
  
}


# Wrapper for the SRA fitting function that uses IN lists of data, OMs and misc settings

SRA_wrap<-function(IN){ # wrapper function allows alternative names for Input lists (IN)
  
  SRA_scope(OM=IN$OM, data=IN$data, ESS=IN$misc$ESS, max_F=IN$misc$max_F,
            selectivity=IN$misc$selectivity, s_selectivity=IN$misc$s_selectivity,
            condition="catch2", mean_fit = T, rescale=1, LWT=IN$misc$LWT, cores=1,
            resample=T,comp_like=IN$misc$comp_like, map_log_rec_dev = c(1:(IN$OM@nyears-2), rep(NA, 2)),
            vul_par = IN$misc$vul_par, map_vul_par = IN$misc$map_vul_par)
  
}

SRA_parallel<-function(i,OMfolders){
  
  In<-readRDS(file=paste0(OMfolders[i],'/In.rda'))
  Fit<-SRA_wrap(In)
  saveRDS(Fit,file=paste0(OMfolders[i],'/Fit.rda'))
  
}

# Wrapper for the SRA plotting function that uses SRA objects and IN lists to plot standardized reports

plot_wrap<-function(out,IN,filename="SRA_scope",title="Operating model conditioning report",dir=tempdir()){ # wrapper function for plotting SRA results for alternative input and output lists
  
  if(out@conv[1]){
    plot(out, compare=T, f_name=IN$misc$f_name, s_name=IN$misc$s_name,
         filename=filename, title=title, dir=dir)
  }else{
    warning("Did not converge!")
  }
  
}


FitInfo<-function(i,OMfolders,OMcode,simno=1){
  
  Fit<-readRDS(paste0(OMfolders[i],'/Fit.rda')) # OM SSB NAA CAA CAL conv Misc mean_fit data config
  
  if(!is.na(Fit@conv[1])){
    temp<-unlist(c(OM=i,Code=OMcode[i],Fit@Misc[[simno]][c('conv','nll','nll_Catch','nll_Ceq','nll_CAA','nll_Index','q')]))
    temp[4:length(temp)]<-round(as.numeric(temp[4:length(temp)]),2)
  }else{
    temp<-c(OM=i,Code=OMcode[i],rep(0,17))
  }
  
  temp
}

RefPoints<-function(i,OMfolders,OMcode,simno=1,p=NA){
  
  temp<-readRDS(paste0(OMfolders[i],'/MSEproj.rda'))
  if(is.na(p)){
    Refs<-apply(cbind(temp@Misc$MSYRefs$Refs,D=temp@OM$D),2,mean)
  }else{
    Refs<-apply(cbind(temp@Misc$MSYRefs$Refs,D=temp@OM$D),2,quantile,p=p)
  }
  Refs<-c(OM=i,Code=OMcode[i],Refs)
  Refs[3:length(Refs)]<-round(as.numeric(Refs[3:length(Refs)]),3)
  Refs
  
}

runProjMSEs<-function(i, OMfolders){
  
  Fit<-readRDS(paste0(OMfolders[i],'/Fit.rda')) # OM SSB NAA CAA CAL conv Misc mean_fit data config
  # Ditch non feasible samples
  keep<-unlist(lapply(Fit@Misc,function(x)all(!is.na(x$F))))
  OM<-Sub_cpars(Fit@OM,sims=keep)
  OM@EffLower<-rep(1,OM@nyears)
  OM@EffUpper<-rep(2,OM@nyears)
  MSEproj<-runMSE(OM=OM,MPs=c("curE","CurC"),silent=T)
  saveRDS(MSEproj,file=paste0(OMfolders[i],'/MSEproj.rda'))
  
}

SSBs<-function(i,OMfolders,OMcode,simno=1,p=NA){
  
  Fit<-readRDS(paste0(OMfolders[i],'/Fit.rda')) # OM SSB NAA CAA CAL conv Misc mean_fit data config
  #temp<-runMSE(OM=Fit@OM,MPs="curE",Hist=T,silent=T)
  if(is.na(p)){
    SSB<-apply(Fit@SSB,2,mean,na.rm=T)
  }else{
    SSB<-apply(Fit@SSB,2,quantile,p=p,na.rm=T)
  }
  SSB
  
}

B_BMSY_proj<-function(i,OMfolders,p=NA,MP="curE"){
  
  temp<-readRDS(paste0(OMfolders[i],'/MSEproj.rda'))
  MPno<-match(MP,temp@MPs)
  
  if(is.na(p)){
    B_BMSY<-apply(temp@B_BMSY[,MPno,],2,mean)
  }else{
    B_BMSY<-apply(temp@B_BMSY[,MPno,],2,quantile,p=p)
  }
  
  B_BMSY
  
}

plotcomp<-function(ts,Design,refline=NA,ylim=NA,yrproj=F){
  
  coly<-c("#ff000050","#0000ff50","#00ff0050")
  colyt<-c('red','blue','green')
  yrs<-1967:2018
  if(yrproj) yrs<-2018+(1:50)
  nplot<-ncol(Design)
  ncol<-2
  nrow=ceiling(nplot/ncol)
  par(mfrow=c(nrow,ncol),mai=c(0.4,0.6,0.3,0.1),omi=c(0.01,0.01,0.5,0.01))
  if(is.na(ylim[1]))ylim=c(0,max(ts))
  
  for(i in 1:nplot){
    
    levs<-unique(Design[,i])
    cols<-coly[match(Design[,i],levs)]
    matplot(yrs,t(ts),col=cols,lty=1,type='l',main=names(Design)[i],ylab="",ylim=ylim,yaxs='i')
    legend('top',legend=levs,text.col=colyt,bty='n')
    if(!is.na(refline[1]))abline(h=refline,col="#99999950",lty=2,lwd=2)
    
  }
  mtext("Year",1,outer=T,line=0.3)
  mtext(deparse(substitute(ts)),2,outer=T,line=0.3)
  mtext("Impact on estimates among axes of uncertainty for reference OM fits",3,outer=T,line=0.3)
  
}

sense_comp_ts<-function(ts,sense_codes,ylim=NA,refline=NA,yrproj=F){
  
  yrs<-1967:2018
  if(yrproj) yrs<-2018+(1:50)
  type<-unlist(lapply(sense_codes,function(X)strsplit(X,"_")[[1]][1]))
  ncomp<-nrow(ts)
  types<-unique(type)
  types<-types[types!="RefCase"]
  nplots<-length(types)
  ncol<-min(2,nplots)
  nrow=ceiling(nplots/ncol)
  cols<-c('black','red','green','blue','grey','orange','pink')
  par(mfrow=c(nrow,ncol),mai=c(0.2,0.2,0.1,0.1),omi=c(0.4,0.4,0.4,0.01))
  if(is.na(ylim[1]))ylim=c(0,max(ts,na.rm=T))
  for(i in 1:nplots){
    ind<-c(1,(1:ncomp)[types[i]==type])
    temp<-ts[ind,]
    matplot(yrs,t(temp),type='l',col=cols,lty=1,ylab="",ylim=ylim,yaxs='i')
    legend('topright',legend=sense_codes[ind],text.col=cols,bty='n')
    if(!is.na(refline[1]))abline(h=refline,col="#99999950",lty=2,lwd=2)
    
  }
  mtext("Year",1,outer=T,line=1)
  mtext(deparse(substitute(ts)),2,outer=T,line=0.9)
  mtext("Impact on estimates among axes of uncertainty for sensitivity analyses",3,outer=T,line=0.3)
  
}

stoch_quants<-function(i,OMfolders){
  
  temp<-readRDS(paste0(OMfolders[i],'/MSEproj.rda'))
  D=temp@OM$D
  SSB_SSBMSY=D/Refs$SSBMSY_SSB0
  SSB<-Refs$SSB0*D
  list(cbind(temp@Misc$MSYRefs$Refs,D=D,SSB=SSB,SSB_SSBMSY))
  
}

plot_stoch<-function(Stoch,Design,quants=c("MSY","SSB_SSBMSY","SSB","UMSY","RefY")){
  
  getquant<-function(Stoch,quantnam)lapply(Stoch,function(x)x[,match(quantnam,names(x))])
  
  ncol<-ncol(Design)
  nrow=length(quants)
  fillcols<-c('#ff000050','#0000ff50','#00ff0050')
  textcols<-c('red','blue','green')
  par(mfrow=c(nrow,ncol),mai=c(0.5,0.2,0.05,0.1),omi=c(0.05,0.3,0.5,0.01))
  
  for(i in 1:nrow){
    
    quanty<-matrix(unlist(getquant(Stoch,quantnam=quants[i])),ncol=nrow(Design))
    breaks=seq(min(quanty),max(quanty),length.out=20)
    for(j in 1:ncol){
      hist(as.vector(quanty),breaks=breaks,col='#99999940',border=NA,main="",xlab="",ylab="")
      if(i==1)mtext(names(Design)[j],3,line=0.7,cex=0.9,font=2)
      if(j==2)mtext(quants[i],1,line=2.5,cex=0.9)
      levs<-unique(Design[,j])
      legend('topright',legend=c("All",levs),text.col=c('grey',textcols),bty='n')
      
      for(k in 1:length(levs)){
        vec<-as.vector(quanty[,Design[,j]%in%levs[k]])
        hist(as.vector(vec),breaks=breaks,col=fillcols[k],border=NA,main="",xlab="",ylab="",add=T)
        
      }
    }
  }
  
  mtext("No. Simulations",2,line=0.7,outer=T,font=2,cex=0.9)
  mtext("Reference Set Factor",3,line=1.8,outer=T,font=2)
  
}


plot_RobGrid<-function(OMcode,Selected=c("1 A H -","3 B H +","2 A L -"),B_BMSY,B_BMSY_p,ylim=c(0,9)){
  
  cols<-c('black','blue','red','green')
  
  par(mfrow=c(2,1),mai=c(0.5,0.2,0.05,0.1),omi=c(0.4,0.6,0.01,0.01))
  
  coly<-rep('#99999940',nrow(B_BMSY))
  coly[match(Selected,OMcode)]<-cols[1:length(Selected)]
  yrs<-1967:2018
  
  matplot(yrs,t(B_BMSY),type='l',col=coly,lty=1,ylab="",ylim=ylim,yaxs='i',lwd=2)
  abline(h=c(0.5,1),col='grey',lty=2)
  legend('topright',legend=c(Selected),text.col=c(cols),bty='n',cex=0.9)
  yrs<-2018+(1:50)
  matplot(yrs,t(B_BMSY_p),type='l',col=coly,lty=1,ylab="",ylim=ylim,yaxs='i',lwd=2)
  #legend('topleft',legend="MSE projection - current exp. rate",bty='n',cex=0.9)
  abline(h=c(0.5,1),col='grey',lty=2)
  mtext("Spawning Stock Biomass relative to MSY levels (B_BMSY)",2,line=1.2,outer=T)
  mtext("Year",1,outer=T)
  
  
}

SSBsOM<-function(i,OMfolders,OMcode){
  
  temp<-readRDS(paste0(OMfolders[i],'/MSEproj.rda'))
  c(NA,apply(apply(temp@SSB_hist,c(1,3),sum),2,mean))
  
}