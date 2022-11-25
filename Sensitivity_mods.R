# =================================================
# === Sensitivity OM modifiers ====================
# =================================================


# Catchability (scale)
Sense_Q<-function(IN,lev=1,verbose=F){
  if(lev==1){
    if(verbose) message("Q level qe: Acoustic survey catchability estimated")
    IN$data$abs_I<-c(0,0)
    IN$data$I_sd[AS_ind,1]<-AS[,3]/100 # s
  }else if(lev==2){
    if(verbose) message("Q Level q1: Acoustic survey catchability fixed to 1")
    IN$data$abs_I<-c(1,0)
    IN$data$I_sd[AS_ind,1]<-AS[,3]/100  # 1x precision
  }else if(lev==3){
    if(verbose) message("Q Level q1: Acoustic survey catchability fixed to 1 and 4x precision in Acoustic Index")
    IN$data$abs_I<-c(1,0)
    IN$data$I_sd[AS_ind,1]<-AS[,3]/200  # 4x precision
  }else if(lev==4){
    if(verbose) message("Q Level q1: Acoustic survey catchability fixed to 1 and 16x precision in Acoustic Index")
    IN$data$abs_I<-c(1,0)
    IN$data$I_sd[AS_ind,1]<-AS[,3]/400  # 16x precision
  }else{
    if(verbose) message("Q Level q1: Acoustic survey catchability fixed to 1 and 36x precision in Acoustic Index")
    IN$data$abs_I<-c(1,0)
    IN$data$I_sd[AS_ind,1]<-AS[,3]/600  # 36x precision
  }
  IN
}

# Recent recruitment
Sense_Rec<-function(Fit, lev=1, verbose=F){
  
  if(lev==4){
    if(verbose)message("Sampling recruitment 2017-2018")
    CV<-apply(log(Fit@OM@cpars$Perr_y[,(1993-1968):(Fit@OM@nyears-2)]),1,sd)
    ind<-(Fit@OM@maxage-1)+ Fit@OM@nyears-(1:0)
    Fit@OM@cpars$Perr_y[,ind]<-trlnorm(2*Fit@OM@nsim,1,CV)
  }else if(lev==5){
    if(verbose)message("Sampling recruitment 2017-2018")
    CV<-apply(log(Fit@OM@cpars$Perr_y[,(1993-1968):(Fit@OM@nyears-2)]),1,sd)
    mu_r<-mean(Fit@OM@cpars$Perr_y[,(1993-1968):(Fit@OM@nyears-2)])
    ind<-(Fit@OM@maxage-1)+ Fit@OM@nyears-(1:0)
    Fit@OM@cpars$Perr_y[,ind]<-matrix(data=mu_r,nrow=nsim,ncol=2)
  }
  Fit
}

class(Sense_Q)<-class(Sense_HC)<-class(Sense_Rec)<-"Sense"
