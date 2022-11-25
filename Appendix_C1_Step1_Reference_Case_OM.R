
# ===========================================================================================================================
# === 1. Specifying data, parameters and conditioning settings for a reference case operating model ===============================
# ===========================================================================================================================
#
# Atlantic Herring MSE
#
# Tom Carruthers t.carruthers@oceans.ubc.ca
# Updated July 24, 2020 by Tim Barrett [recruitment regime shift in 1994 used for projections; equilibrium catch added; larval survey changed to index of abundance]
#
# This script:
# a) Loads the most recent data and growth parameters, and formats these for OM conditioning
# b) Creates a default operating model object for herring
# c) Specified default settings for conditioning a Reference Case model
#
# The product of this script: the Reference Case operating model from which all other OMs are based 'OMs/OM_RC.rda'

# === Prerequisites =========================================================================================

library(MSEtool) #Version 1.7.0 is available as of June 24, 2020
library(DLMtool)

# --- Directories ---

OMdir<-"C:/Users/BARRETTTJ/Desktop/MSE"
setwd(OMdir)
datdir<-paste0(getwd(),"/Data Inputs/")

# --- Source code ---

source(paste0(OMdir,"/Source.r"))

# --- For plotting data ---

cols<-rep(c("black","red","green","blue","grey","orange","pink"),2)
ltys<-rep(1:2,each=7)


# === Load and Format data ===================================================================================
#
# Data prepared by Tim Barrett tim.barrett@dfo-mpo.gc.ca
# Accessible from: https://drive.google.com/open?id=18ZRW34Ea_kEJxvPNhg3LIJJzOQPhMc83
# Units are THOUSAND TONNES

data<-new('list')

# --- Catch data -----

Cat<-readRDS(file=paste0(datdir,"Chist_1_200228.rds"))

nf<-ncol(Cat)                                                 # fleet dimension
ny<-nrow(Cat)                                                 # year dimension
yind<-1:ny
Year<-as.numeric(row.names(Cat))
YearRNG<-range(Year)
matplot(Year,as.matrix(Cat),type='l',col=cols,lty=ltys,ylab="Catch (t)")
legend('topright',legend=names(Cat),text.col=cols)

data$Chist<-as.matrix(Cat)/1000 #Catch in 1000s of t

# --- Historical equilibrium catch (30 years)

HCat<-readRDS(file=paste0(datdir,"Chist_1_Inital_cond_200303.rds"))
data$C_eq <- unlist(HCat)/1000  #Catch in 1000s of t

# --- CAA data ------

CAA<-readRDS(file=paste0(datdir,"CAA_1_200228.rds")) # age, fleet, year (needs to be year, age, fleet)
na<-dim(CAA)[1]                                               # age dimension
data$CAA<-array(NA,c(ny,na,nf))
CAA_yind<-match(dimnames(CAA)[[3]],Year)
CAA_ind<-as.matrix(expand.grid(1:na,1:nf,1:(dim(CAA)[3])))
CAA_p_ind<-CAA_ind[,c(3,1,2)]
CAA_p_ind[,1]<-CAA_yind[CAA_p_ind[,1]]
data$CAA[CAA_p_ind]<-CAA[CAA_ind]

# a check: y<-44; y1=match(y,CAA_yind); a<-3; f=2; data$CAA[y,a,f]==CAA[a,f,y1]; data$CAA[,,2] # Test of reshaping and transposition


# --- CAL data ------

CAL<-readRDS(file=paste0(datdir,"CAL_1_200402.rds")) # age, fleet, year (needs to be year, age, fleet)
nl<-dim(CAL)[1]                                               # length class dimension
data$CAL<-array(NA,c(ny,nl, nf))
CAL_yind<-match(dimnames(CAL)[[3]],Year)
CAL_ind<-as.matrix(expand.grid(1:nl,1:nf,1:(dim(CAL)[3])))
CAL_p_ind<-CAL_ind[,c(3,1,2)]
CAL_p_ind[,1]<-CAL_yind[CAL_p_ind[,1]]
data$CAL[CAL_p_ind]<-CAL[CAL_ind]

# a check: y<-44; y1=match(y,CAL_yind); l<-11; f=2; data$CAL[y,l,f]==CAL[l,f,y1]; data$CAL[,,2] # Test of reshaping and transposition

ns<-2   # number of surveys

# --- Index data -----

data$Index <- data$I_sd<-matrix(NA,nrow=ny,ncol=ns)  # year, n surveys
data$s_CAA <- array(NA,c(ny,na,ns))              # year, ages, n surveys
data$s_CAL <- array(NA,c(ny,nl,ns))               # year, n lengths, n surveys

# Acoustic Survey
AS<-as.matrix(readRDS(file=paste0(datdir,"A_INDEX_200228.rds"))) # 1999-2018
AS_ind<-(1:ny)[Year %in% (1999:2018)]
data$Index[AS_ind,1]<-AS[,2] / 1000 # thousand tonnes
data$I_sd[AS_ind,1]<-AS[,3]/100 # starting with a 20% CV on acoustic surveys

# Acoustic Survey composition
AS_CAA<-readRDS(paste0(datdir,"ADAIndex_200224.rds"))
data$s_CAA[AS_ind,,1]<-t(AS_CAA[,1,])

ACAL <- read.csv(paste0(datdir,"ACAL_200318.csv"))
rownames(ACAL)<-ACAL$Year
ACAL<-ACAL[,colnames(ACAL)!="Year"]
data$s_CAL[AS_ind,,1]<-as.matrix(ACAL)

# Larval Survey
LS<-readRDS(file=paste0(datdir,"L_INDEX_200228.rds")) # 1999-2018
LS_ind<-(1:ny)[Year %in% LS$YEAR]
data$Index[LS_ind,2]<-LS[,2]
data$I_sd[LS_ind,2]<-LS[,4]/100 # sd are really CVs

#                Acoustic    Larval
data$I_type <- c('est',      'SSB')   # 'est' for estimating selectivity or otherwise selectivity follows maturity (SSB) or is 1 for all ages (B)
data$abs_I  <- c(0,           0)      # q is estmated for the larval survey
data$I_units <-c(1,           0)      # Biomass for acoustic, abundance for larval


data$length_bin<-as.numeric(dimnames(CAL)[[1]])


# === Build OM with time varying parameters ========================================================================================================================

# --- Additional operating model dimensions ---

np<-50     # Number of projected years
nsim<-48    # Number of simulations

# --- Growth stuff ---

Gpars<-read.csv(paste0(datdir,"VBGrowth_191106.csv"))
LAA<-readRDS(paste0(datdir,"LAA_200218.rds")) # Observed length at age
WAA<-readRDS(paste0(datdir,"WAA_200218.rds")) # Observed weight at age

Wt_age<-Len_age<-array(NA,c(nsim,na,ny+np)) # Length at age by year is needed to impute Maturity at age from Length at maturity vectors
G_ind<-(1:ny)[Year %in% rownames(LAA)] # by observed LAA data

yamat<-as.matrix(expand.grid(1:nsim,1:na,G_ind))
refmat<-as.matrix(expand.grid(1:nsim,1:na,1:(length(G_ind))))
yind<-refmat[,3]

# by Observed LAA data
Len_age[yamat]<-LAA[refmat[,3:2]]
Wt_age[yamat]<-WAA[refmat[,3:2]]

# --- now impute Wt_age for missing years (historical and projected) ---

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

# needed for level B
coeff <- read.csv(file=paste0(datdir,"waa_coeff_200504.csv"))

# --- Natural mortality rate ----------------------

M_ageArray <- array(0.35,c(nsim,na,np+ny))

# --- Maturity ------------------------------------

Mat <- read.csv(paste0(datdir,"MAT_200320.csv"))

Mat_age <- array(NA, dim = c(nsim, na, ny+np))
Mat_ind <-(1:ny)[Year %in% Mat$YEAR]

for(y in 1:length(Mat_ind))
{
  yind<-Mat_ind[y]
  Mat_age[,,yind] <- do.call("rbind", replicate(nsim, 1/(1+exp(-(Mat$Ab0[y] + Mat$Ab1[y] * c(1:na)))), simplify = FALSE))
}

nyrsm<-3
early<-Mat_ind[1:nyrsm]
late<-Mat_ind[length(Mat_ind)-((nyrsm-1):0)]

Mat_age_e<-apply(Mat_age[,,early,drop=F],1:2,mean)
Mat_age_l<-apply(Mat_age[,,late,drop=F],1:2,mean)

Mat_age[,,1:(Mat_ind[1]-1)]<-rep(Mat_age_e,Mat_ind[1]-1)
Mat_age[,,(Mat_ind[length(Mat_ind)]+1):(ny+np)]<-rep(Mat_age_l,ny+np -Mat_ind[length(Mat_ind)])


# === Make OM =========================================================================================

# --- Default herring OM builder (from 'Make_default_OM.R' script) -------

OM <- DHOM(Name="OM_Fleet_1", nsim=nsim, ny=ny, na=na, np=50, h=0.95,
           Perr=1, interval=2, a=mean(Gpars$a)*1E-6, b=mean(Gpars$b))

# --- Add custom parameters to the OM object -----------------------------

OM@cpars <- list(Wt_age=Wt_age, Len_age=Len_age, Mat_age=Mat_age, M_ageArray=M_ageArray)


# === Define miscellaneous setting and labels for fitting and plotting ================================

misc <- list(ESS=c(1000,1000), max_F=10, selectivity=c('dome','dome','logistic','dome'),
          s_selectivity=c("dome","dome"), f_name=names(Cat), s_name<-c("Aco. Surv.","Larval"),
          LWT=list(Index=rep(1,ns), s_CAA = rep(0.05,ns), s_CAL = rep(0,ns), CAA = c(0.05,0.05,0.05,0), CAL=rep(0,4)),
          comp_like="multinomial") #c(1,1,0.5,1)))


# === Create an OM conditioning input list IN that can be modified by other functions =================

source('Reference_mods.R')          # Modifiers for input object to make reference OMs
In_RC<-Ref_M(Ref_G(Ref_R(Ref_C(list(data=data, OM=OM, misc=misc))))) # apply reference level 1 for all factors

saveRDS(In_RC,"Reference_Case_OM/In_RC.rda")


# === Fit and report ==================================================================================

source("Additional_mods.R")
if(!exists('nofit')){ # later scripts need the objects above but don't need to rerun this
  Fit<-SRA_wrap(In_RC)
  Fit <- Add_Regime(Fit,lev=2,verbose=T)    ##add regime shift for recruitment
  saveRDS(Fit,"Reference_Case_OM/Fit.rda")
  plot_wrap(Fit, In_RC, dir=paste0(getwd(),"/Reference_Case_OM/"), filename="/RC_Report",  title="Reference Case Operating Model Fitting Report")
}

# ====================================================================================================================
# ======= END OF SCRIPT ==============================================================================================
# ====================================================================================================================


