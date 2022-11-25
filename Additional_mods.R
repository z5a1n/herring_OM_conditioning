# =================================================
# === Additional OM modifiers ====================
# =================================================

Add_M<-function(IN,lev=1,verbose=F){
  # M_ageArray[nsim,na,np+ny]
  maxage<-IN$OM@maxage
  if(lev==1){
    if(verbose)message("Natural mortality rate (M) - Level 1:  value of 0.35 used for all ages")
    IN$OM@cpars$M_ageArray[]<-0.35
  }else if(lev==2){
    if(verbose)message("Natural mortality rate (M) - Level 2:  value of 0.2 used for all ages")
    IN$OM@cpars$M_ageArray[]<-0.2
  }
  IN
}

Add_NLS<-function(IN,lev=1,verbose=F){
  if(lev==1){
    IN$misc$LWT$Index[2]<-1
  }else if(lev==2){
    IN$misc$LWT$Index[2]<-0
  }
  IN
}

Add_NLcheck<-function(IN,lev=1,verbose=F){
  if(lev==1){
    IN$misc$LWT$Index[2]<-1
    
  }else if(lev==2){
    IN$misc$LWT$Index[2]<-0
    todo<-!is.na(IN$data$Index[,2])
    IN$data$Index[todo,2]<-1
  }
  IN
}


proj_rec<-function(Perr_y,yind,maxage,nyears,proyears){
  
  getlag1ac<-function(x)stats::acf(x)$acf[2,1,1]
  
  CV<-apply(Perr_y[,yind],1,sd)
  AC<-apply(Perr_y[,yind],1,getlag1ac)
  mu<-apply(Perr_y[,yind],1,mean)
  
  for(s in 1:nsim){
    for(y in (maxage-1+nyears-2)+(1:(proyears+2))){ # -1 is 2019, -3 is 2017 so 1:proyears if 2018-proyears
      Perr_y[s,y] <- AC[s]*Perr_y[s,y-1] + ((1-AC[s]^2)^0.5)* rnorm(n=1, mean=mu[s], sd=CV[s]) - (1-AC[s])* (CV[s]^2)/2
    }
  }
  
  Perr_y
  
}


# Regime shift or not
Add_Regime<-function(Fit, lev=1, verbose=F){
  
  Perr_y<-log(Fit@OM@cpars$Perr_y) # the process error matrix [nsim, (maxage-1)+nyears+proyears]
  nyears<-Fit@OM@nyears
  proyears<-Fit@OM@proyears
  maxage<-Fit@OM@maxage
  nsim<-Fit@OM@nsim
  
  if(lev==2){ # 1994 and later
    
    yind<-(maxage-1)+((1994-1968+1):(nyears-2))
    if(verbose) {message("Future recruitment from 1994-2016")}
    Fit@OM@cpars$Perr_y<-exp(proj_rec(Perr_y,yind,maxage,nyears,proyears))
    
  }else if(lev==3){ # pre 1994
    
    yind<-(maxage-1)+(1:(1994-1968))
    if(verbose)message("Future recruitment from 1968-1993")
    Fit@OM@cpars$Perr_y<-temp<-exp(proj_rec(Perr_y,yind,maxage,nyears,proyears))
    
  }else if(lev==4){ # after 2010
    
    yind<-(maxage-1)+((2010-1968+1):(nyears-2))
    if(verbose)message("Future recruitment from 2010-2019")
    Fit@OM@cpars$Perr_y<-temp<-exp(proj_rec(Perr_y,yind,maxage,nyears,proyears))
    
  }else if(lev==5){ # 1994-2009
    
    yind<-(maxage-1)+((1994-1968+1):(2009-1968+1))
    if(verbose)message("Future recruitment from 1994-2009")
    Fit@OM@cpars$Perr_y<-temp<-exp(proj_rec(Perr_y,yind,maxage,nyears,proyears))
    
  }else if(lev==6){ # linear decrease from 1994-2016
    
    yind<-(maxage-1)+((1994-1968+1):(2016-1968+1))
    if(verbose)message("Future recruitment linear decrease")
    Fit@OM@cpars$Perr_y<-temp<-exp(proj_rec(Perr_y,yind,maxage,nyears,proyears))
  }
  
  Fit
  
}



Add_G<-function(IN,lev="A",verbose=F){
  
  Wt_age<-Len_age<-array(NA,c(nsim,na,ny+np)) # Length at age by year is needed to impute Maturity at age from Length at maturity vectors
  G_ind<-(1:ny)[Year %in% rownames(LAA)] # by observed LAA data
  
  yamat<-as.matrix(expand.grid(1:nsim,1:na,G_ind))
  refmat<-as.matrix(expand.grid(1:nsim,1:na,1:(length(G_ind))))
  yind<-refmat[,3]
  
  # by Observed LAA data
  Len_age[yamat]<-LAA[refmat[,3:2]]
  Wt_age[yamat]<-WAA[refmat[,3:2]]
  
  
  if(lev=="A"){
    if(verbose)message("Growth - Level A:  projected growth is average of last 3 years")
    
    nyrsm<-3
    early<-G_ind[1:nyrsm]
    late<-G_ind[length(G_ind)-((nyrsm-1):0)]
    
    Wt_age_e<-apply(Wt_age[,,early,drop=F],1:2,mean)
    Wt_age_l<-apply(Wt_age[,,late,drop=F],1:2,mean)
    
    Wt_age[,,1:(G_ind[1]-1)]<-rep(Wt_age_e,G_ind[1]-1)
    Wt_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Wt_age_l,ny+np -G_ind[length(G_ind)])
    
    Len_age_e<-apply(Len_age[,,early,drop=F],1:2,mean)
    Len_age_l<-apply(Len_age[,,late,drop=F],1:2,mean)
    
    Len_age[,,1:(G_ind[1]-1)]<-rep(Len_age_e,G_ind[1]-1)
    Len_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Len_age_l,ny+np -G_ind[length(G_ind)])
    
  }
  if(lev=="B"){# CASE 2 --- now impute Wt_age and Len_age for missing years (projected using mean of last three years and regression for Wt-at-age projections) ---
    if(verbose)message("Growth - Level B:  projected growth is based on linear regression")
    
    nyrsm<-3
    early<-G_ind[1:nyrsm]
    late<-G_ind[length(G_ind)-((nyrsm-1):0)]
    
    Wt_age_e<-apply(Wt_age[,,early,drop=F],1:2,mean)
    
    Wt_age[,,1:(G_ind[1]-1)]<-rep(Wt_age_e,G_ind[1]-1)
    
    Len_age_e<-apply(Len_age[,,early,drop=F],1:2,mean)
    Len_age_l<-apply(Len_age[,,late,drop=F],1:2,mean)
    
    Len_age[,,1:(G_ind[1]-1)]<-rep(Len_age_e,G_ind[1]-1)
    Len_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Len_age_l,ny+np -G_ind[length(G_ind)])
    
    for(i in 1:dim(Wt_age)[3]) # add projections based on linear trend in WAA from waa_coeff.csv
    {
      if(is.na(Wt_age[1,1,i]))
      {
        for(j in 1:dim(Wt_age)[1])
        {
          yr <- i+1968-1
          Wt_age[j,,i] <- 10^(coeff$INTERCEPT+coeff$SLOPE*yr)
        }
      }
    }
  }  
  if(lev=="C"){# linear increase in weight at age based on a percentage of mean of last 3 years
    if(verbose)message("Growth - Level C:  projected growth is linear increase 0.5%")
    
    growthinc <- 0.5
    nyrsm<-3
    early<-G_ind[1:nyrsm]
    late<-G_ind[length(G_ind)-((nyrsm-1):0)]
    
    Wt_age_e<-apply(Wt_age[,,early,drop=F],1:2,mean)
    Wt_age_l<-apply(Wt_age[,,late,drop=F],1:2,mean)
    
    Wt_age[,,1:(G_ind[1]-1)]<-rep(Wt_age_e,G_ind[1]-1)  #fill in early missing Weights
    
    firstmissyrID <- 2019-1968+1
    Wt_age[,,firstmissyrID]<-Wt_age_l  #first missing year is average of last 3 years
    Wt_age[,,firstmissyrID+1]<-Wt_age[,,firstmissyrID]*(1 + growthinc/100) #next year is an increase of growthinc percent
    
    Delta <- Wt_age[,,firstmissyrID+1] - Wt_age[,,firstmissyrID]
    
    for(i in (firstmissyrID+2) : (ny+np))
    {
      Wt_age[,,i] <- Wt_age[,,i-1] + Delta
    }
    
    Len_age_e<-apply(Len_age[,,early,drop=F],1:2,mean)
    Len_age_l<-apply(Len_age[,,late,drop=F],1:2,mean)
    
    Len_age[,,1:(G_ind[1]-1)]<-rep(Len_age_e,G_ind[1]-1)
    Len_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Len_age_l,ny+np -G_ind[length(G_ind)])
  }
  if(lev=="D"){# linear increase in weight at age based on a percentage of mean of last 3 years
    if(verbose)message("Growth - Level C:  projected growth is linear increase 1%")
    
    growthinc <- 1
    nyrsm<-3
    early<-G_ind[1:nyrsm]
    late<-G_ind[length(G_ind)-((nyrsm-1):0)]
    
    Wt_age_e<-apply(Wt_age[,,early,drop=F],1:2,mean)
    Wt_age_l<-apply(Wt_age[,,late,drop=F],1:2,mean)
    
    Wt_age[,,1:(G_ind[1]-1)]<-rep(Wt_age_e,G_ind[1]-1)  #fill in early missing Weights
    
    firstmissyrID <- 2019-1968+1
    Wt_age[,,firstmissyrID]<-Wt_age_l  #first missing year is average of last 3 years
    Wt_age[,,firstmissyrID+1]<-Wt_age[,,firstmissyrID]*(1 + growthinc/100) #next year is an increase of growthinc percent
    
    Delta <- Wt_age[,,firstmissyrID+1] - Wt_age[,,firstmissyrID]
    
    for(i in (firstmissyrID+2) : (ny+np))
    {
      Wt_age[,,i] <- Wt_age[,,i-1] + Delta
    }
    
    Len_age_e<-apply(Len_age[,,early,drop=F],1:2,mean)
    Len_age_l<-apply(Len_age[,,late,drop=F],1:2,mean)
    
    Len_age[,,1:(G_ind[1]-1)]<-rep(Len_age_e,G_ind[1]-1)
    Len_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Len_age_l,ny+np -G_ind[length(G_ind)])
  }
  if(lev=="E"){# step back in growth
    if(verbose)message("Growth - Level E: step back in growth")
    
    nyrsm<-3
    early<-G_ind[1:nyrsm]
    late<-G_ind[length(G_ind)-((nyrsm-1):0)]
    
    Wt_age_e<-apply(Wt_age[,,early,drop=F],1:2,mean)
    Wt_age_l<-apply(Wt_age[,,late,drop=F],1:2,mean)
    
    Wt_age[,,1:(G_ind[1]-1)]<-rep(Wt_age_e,G_ind[1]-1)  #fill in early missing Weights
    
    for (i in 1:np)
    {
      Wt_age[,,ny+i] <- Wt_age[,,ny-i]
    }
    
    Len_age_e<-apply(Len_age[,,early,drop=F],1:2,mean)
    Len_age_l<-apply(Len_age[,,late,drop=F],1:2,mean)
    
    Len_age[,,1:(G_ind[1]-1)]<-rep(Len_age_e,G_ind[1]-1)
    Len_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Len_age_l,ny+np -G_ind[length(G_ind)])
  }
  if(lev=="X1"){# linear decrease in weight at age based on a percentage of mean of last 3 years
    if(verbose)message("Growth - Level X1:  projected growth is linear decrease 0.5%")
    
    growthinc <- 0.5
    nyrsm<-3
    early<-G_ind[1:nyrsm]
    late<-G_ind[length(G_ind)-((nyrsm-1):0)]
    
    Wt_age_e<-apply(Wt_age[,,early,drop=F],1:2,mean)
    Wt_age_l<-apply(Wt_age[,,late,drop=F],1:2,mean)
    
    Wt_age[,,1:(G_ind[1]-1)]<-rep(Wt_age_e,G_ind[1]-1)  #fill in early missing Weights
    
    firstmissyrID <- min(which(is.na(Wt_age[1,1,])))
    Wt_age[,,firstmissyrID]<-Wt_age_l  #first missing year is average of last 3 years
    Wt_age[,,firstmissyrID+1]<-Wt_age[,,firstmissyrID]*(1 - growthinc/100) #next year is an increase of growthinc percent
    
    Delta <- Wt_age[,,firstmissyrID+1] - Wt_age[,,firstmissyrID]
    
    for(i in (firstmissyrID+2) : (ny+np))
    {
      Wt_age[,,i] <- Wt_age[,,i-1] + Delta
    }
    
    Len_age_e<-apply(Len_age[,,early,drop=F],1:2,mean)
    Len_age_l<-apply(Len_age[,,late,drop=F],1:2,mean)
    
    Len_age[,,1:(G_ind[1]-1)]<-rep(Len_age_e,G_ind[1]-1)
    Len_age[,,(G_ind[length(G_ind)]+1):(ny+np)]<-rep(Len_age_l,ny+np -G_ind[length(G_ind)])
  }
  IN$OM@cpars$Wt_age <- Wt_age
  IN$OM@cpars$Len_age <-Len_age
  IN
}



