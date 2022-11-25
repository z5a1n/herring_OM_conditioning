# Reference OM modifiers

# Natural mortality rate, Ref Factor 1
Ref_M<-function(IN,lev=1,verbose=F){
  # M_ageArray[nsim,na,np+ny]
  maxage<-IN$OM@maxage
  if(lev==1){
    if(verbose)message("Natural mortality rate (M) - Level 1:  value of 0.35 used for all ages")
    IN$OM@cpars$M_ageArray[]<-0.35
  }else if(lev==2){
    if(verbose)message("Natural mortality rate (M) - Level 2:  0.49 for ages 1 & 2, 0.26 for ages 3+")
    IN$OM@cpars$M_ageArray[,1:2,]<-0.49
    IN$OM@cpars$M_ageArray[,3:maxage,]<-0.26
  }else{
    if(verbose)message("Natural mortality rate (M) - Level 3:  0.72 for ages 1 & 2, 0.45 for ages 3+")
    IN$OM@cpars$M_ageArray[,1:2,]<-0.72
    IN$OM@cpars$M_ageArray[,3:maxage,]<-0.45
  }
  IN
}

# Future growth, Ref Factor 2
Ref_G<-function(IN,lev="A",verbose=F){

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
 
  IN$OM@cpars$Wt_age <- Wt_age
  IN$OM@cpars$Len_age <-Len_age
  IN
}

# Steepness (Resilience), Ref Factor 3
Ref_R<-function(IN,lev="H",verbose=F){
  if(lev=="H"){
    if(verbose)message("Resilience - Level H:  steepness=0.95")
    IN$OM@h<-c(0.95,0.95)
  }else{
    if(verbose)message("Resilience - Level L:  steepness=0.75")
    IN$OM@h<-c(0.7,0.7)
  }
  IN
}

# Magnitude of weir catches 
Ref_C<-function(IN,lev="-",verbose=F){ # Weir catche / on or off
  IN$data$Chist<-as.matrix(Cat)/1000
  if(lev=='-'){
    if(verbose) message("Catch - Level -:  no weir or historical equilibrium catches included")
    IN$data$Chist[,4]<- 0# IN$data$Chist[,4]/1000 # one 1-thousandth of the historical catches.
    IN$data$C_eq[4]<- 0# IN$data$C_eq[4]/1000 # one 1-thousandth of the equilibrium catches.
  }else{ 
    if(verbose) message("Catch - Level +:  all weir catches and historical equilibrium catches included")
  }
  IN
}


class(Ref_M)<-class(Ref_G)<-class(Ref_R)<-class(Ref_C)<-"Ref"

