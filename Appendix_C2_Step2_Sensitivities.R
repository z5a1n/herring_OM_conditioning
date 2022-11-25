# ===========================================================================================================================
# === 2. Sensitivity analyses ===============================================================================================
# ===========================================================================================================================
#
# Atlantic Herring MSE
#
# March 22nd 2020
# Tom Carruthers t.carruthers@oceans.ubc.ca
#
# This script:
# a) Loads a set of input object modifiers for sensitivities, reference set OMs and robustness set OMs
# b) Fits operating models for all
# c) Analyzes the importance of the various single factor modifications from the reference case model
#
# The product of this script: a set of OM diagnostics used for specifying refernece and robustness set OMs

# === Prerequisites =========================================================================================

rm(list=ls())
library(MSEtool)
library(DLMtool)

# --- Directories ---

OMdir<-"C:/Users/BARRETTTJ/Desktop/MSE"

setwd(OMdir)
n_sense<-30
outdirs<-paste0(OMdir,"/Sensitivity_OMs/",1:n_sense)
for(i in 1:n_sense)if(!dir.exists(outdirs[i]))dir.create(outdirs[i])

# --- Source Code
nofit<-T
source('Appendix_C1_Step1_Reference_Case_OM.R') # Contains various data objects and the reference OM required for OMs and returns the reference case OM
source('Sensitivity_mods.R')        # Modifiers for input object to make sensitivity OMs
source('Reference_mods.R')          # Modifiers for input object to make reference OMs
source('Additional_mods.R')         # Modifiers for input object to make robustness OMs

# --- objects

In_RC<-readRDS(paste0(OMdir,"/Reference_Case_OM/In_RC.rda"))

# === 1 Reference Case ==============================================================

saveRDS(In_RC,paste0(outdirs[1],'/In.rda'))

# === 2-3 Natural Mortality =========================================================

In <- Ref_M(In_RC,lev=2,verbose=T) # low M, age varying
saveRDS(In,paste0(outdirs[2],'/In.rda'))

In <- Ref_M(In_RC,lev=3,verbose=T) # high M, age varying
saveRDS(In,paste0(outdirs[3],'/In.rda'))

# === 4 Growth ======================================================================

In <- Ref_G(In_RC,lev='B',verbose=T) # linear extrapolation of growth changes
saveRDS(In,paste0(outdirs[4],'/In.rda'))

# === 5-9 Resilience ======================================================================

In <- In_RC # steepness = 0.85
In$OM@h<-rep(0.85,2)
saveRDS(In,paste0(outdirs[5],'/In.rda'))

In <- In_RC # steepness = 0.75
In$OM@h<-rep(0.75,2)
saveRDS(In,paste0(outdirs[6],'/In.rda'))

In <- In_RC # steepness = 0.7
In$OM@h<-rep(0.7,2)
saveRDS(In,paste0(outdirs[7],'/In.rda'))

In <- In_RC # # steepness = 0.65
In$OM@h<-rep(0.65,2)
saveRDS(In,paste0(outdirs[8],'/In.rda'))

In <- In_RC # steepness = 0.55
In$OM@h<-rep(0.55,2)
saveRDS(In,paste0(outdirs[9],'/In.rda'))


# === 10 Weir catches ======================================================================

In <- In_RC
In$data$Chist[,4] <- In$data$Chist[,4]*1000 # Ref case is /1000
In$data$C_eq[4] <- In$data$C_eq[4]*1000     # Ref case is /1000
In$misc$LWT$CAA[4] <- In$misc$LWT$CAA[1]    # put CAA weight back on Weir
saveRDS(In,paste0(outdirs[10],'/In.rda'))

# === 11 - 14 Acoustic survey q ===========================================================

for(lev in 2:5){
  In <- Sense_Q(In_RC,lev=lev,verbose=T) # All weir catches included
  saveRDS(In,paste0(outdirs[9+lev],'/In.rda'))
}

# === 15 Fit to catch at length comp ===================================================

In <- In_RC # both CAL and CAA
In$misc$LWT$s_CAL<-0.05
In$misc$LWT$CAL<-rep(0.05,4)
saveRDS(In,paste0(outdirs[15],'/In.rda'))


# === 16 Two Purse seinse fleet ==============================================================

source(paste0(OMdir,"/Sensitivity_OMs/Two_Purse_Seine.R")) # return
saveRDS(In,paste0(outdirs[16],'/In.rda'))

# === 17-20 Composition downweighting ===========================================================
In_RC<-readRDS(paste0(OMdir,"/Reference_Case_OM/In_RC.rda"))

In <- In_RC #
In$misc$LWT$s_CAA<-0.01
In$misc$LWT$CAA<-rep(0.01,4)
saveRDS(In,paste0(outdirs[17],'/In.rda'))

In <- In_RC #
In$misc$LWT$s_CAA<-0.1
In$misc$LWT$CAA<-rep(0.1,4)
saveRDS(In,paste0(outdirs[18],'/In.rda'))

In <- In_RC #
In$misc$LWT$s_CAA<-0.2
In$misc$LWT$CAA<-rep(0.2,4)
saveRDS(In,paste0(outdirs[19],'/In.rda'))

In <- In_RC #
In$misc$LWT$s_CAA<-0.5
In$misc$LWT$CAA<-rep(0.5,4)
saveRDS(In,paste0(outdirs[20],'/In.rda'))


# === 21 Old Natural Mortality = 0.2 =========================================================

In <- Add_M(In_RC,lev=2,verbose=T) # low M, age varying
saveRDS(In,paste0(outdirs[21],'/In.rda'))

# === 22 Old Natural Mortality = 0.2 and steepness 0.75 =========================================================

In <- In_RC
In$OM@cpars$M_ageArray[]<-0.2
In$OM@h<-c(0.75,0.75)
saveRDS(In,paste0(outdirs[22],'/In.rda'))

# === 23 No larval survey  =======================================================

In <- Add_NLS(In_RC,lev=2,verbose=T) # zero weight on larval survey
saveRDS(In,paste0(outdirs[23],'/In.rda'))

# === 24 No larval survey check all 1s  =======================================================

In <- Add_NLcheck(In_RC,lev=2,verbose=T) # zero weight on larval survey
saveRDS(In,paste0(outdirs[24],'/In.rda'))

# === 25 No Regime shift in recruitment POST FIT ============================================

Fit0<-Fit<-SRA_wrap(In_RC)
saveRDS(Fit, paste0(outdirs[25],'/Fit.rda')) # 
saveRDS(In_RC, paste0(outdirs[25],'/In.rda')) # 

# === 26 Regime shift to 1994- recruitment POST FIT ============================================

Fit <- Add_Regime(Fit0,lev=3,verbose=T) # pre 1994 process error mean, sd and ac
saveRDS(Fit, paste0(outdirs[26],'/Fit.rda')) # 
saveRDS(In_RC, paste0(outdirs[26],'/In.rda')) # 

# === 27 Regime shift after 2010 POST FIT ============================================

Fit <- Add_Regime(Fit0,lev=4,verbose=T) # post 2010 process error mean, sd and ac
saveRDS(Fit, paste0(outdirs[27],'/Fit.rda')) # 
saveRDS(In_RC, paste0(outdirs[27],'/In.rda')) # 

# === 28 recruitment 1994-2009  ==============================================================
Fit <- Add_Regime(Fit0,lev=5,verbose=T) # post process error mean, sd and ac
saveRDS(Fit, paste0(outdirs[28],'/Fit.rda')) # 
saveRDS(In_RC, paste0(outdirs[28],'/In.rda')) # 

# === 29 increase in growth ==============================================================

In <- Add_G(In_RC,lev='E',verbose=T) # step change in growth
saveRDS(In,paste0(outdirs[29],'/In.rda'))

# === 30  ==============================================================

In <- In_RC # steepness = 0.5
In$OM@h<-rep(0.9,2)
saveRDS(In,paste0(outdirs[30],'/In.rda'))

# === 21 Time varying M - Predator consumption in fleet 5 and change to M ==============================================================

#source(paste0(OMdir,"/Sensitivity_OMs/TVM.R")) # return
#saveRDS(In,paste0(outdirs[21],'/In.rda'))

# === create names and codes ==================================================================

sense_names<-c("Reference Case","Low age-varying M","High age-varying M","Linear extrapolation of growth",
               "Steepness 0.85","Steepness 0.75","Steepness 0.7","Steepness 0.65", "Steepness 0.55",
               "All weir catch", #10
               "Acoustic q = 1","Acoustic q = 1, 4x precision","Acoustic q = 1, 16x precision",
               "Acoustic q = 1, 36x precision","Fit to both age and length composition", 
               "Two Purse seine fleets",  ###16
               "Age composition weight = 0.01","Age composition weight = 0.1",
               "Age composition weight = 0.2","Age composition weight = 0.5",
               "M 0.2","M 0.2 Low h", "No larval survey", "No larval survey_check", #24
               "Rec 1968-2016","Rec 1968-1993","Rec 2010-2016", "Rec 1994-2009", "Stepback in growth", 
               "Steepness 0.90")

sense_codes<-c("RefCase","M_LowMv","M_HighMv","G_ChangeGrowth",
               "h_Steep85","h_Steep75","h_Steep70","h_Steep65","h_Steep55",
               "C_WeirCat", #10
               "q_1","q_1x4","q_1x16","q_1x36",
               "Comp_AgeLength","PS_Two", #16
               "CompWt_01","CompWt_1","CompWt_2","CompWt_5",
               "M_0.2", "M_0.2Lh", "LS_none", "LS_none1", #24
               "R_1968-2016", "R_1968-1993","R_2010-2016","R_1994-2009", "G_stepback", 
               "h_Steep90")

write.csv(cbind(1:length(sense_codes),sense_codes),paste0(OMdir,"/Tables/Step2_Sensitivities/sense_codes.csv"))


# === Fit OMs =================================================================================

#setup()               # set up parallel processing for MSEtool
#sfExport('SRA_wrap')  # export SRA_wrap function to the cluster
#system.time({sfSapply(1:n_sense,SRA_parallel,OMfolders=outdirs)}) # takes about 2 minutes

#don't fit for recruitment shifts i 25-28
for(i in c(1:n_sense)[-(25:28)]) { # takes about 20 minutes
  In<-readRDS(paste0(outdirs[i],"/In.rda"))
  Fit<-SRA_wrap(In)
  Fit<-Add_Regime(Fit,lev=2,verbose=T)    ##add regime shift for recruitment
  saveRDS(Fit,paste0(outdirs[i],"/Fit.rda"))
}  
# === Build OM reports =========================================================================

for(i in 1:n_sense) { # takes about 20 minutes
  In<-readRDS(paste0(outdirs[i],"/In.rda"))
  Fit<-readRDS(paste0(outdirs[i],"/Fit.rda"))
  plot_wrap(Fit, In, dir=outdirs[i], filename="/OM_Report",  title=paste("Sensitivity Operating Model Fitting Report:",sense_codes[i]))
}

# === Do MSE projections with status quo effort and catch ======================================

#sfSapply(1:n_sense,fun=runProjMSEs,OMfolders=outdirs) # Run projection MSEs with current effort (curE MP) and current catch (CurC MP), takes about 2 minutes

for(i in 1:n_sense) { 
 runProjMSEs(i,outdirs)
}

# === Extract model estimates ==================================================================

#FitDat<-as.data.frame(t(sfSapply(c(1:15,17:n_sense),FitInfo,OMfolders=outdirs,sense_codes)),stringsAsFactors=F) # Fitting data
FitDat<-as.data.frame(t(sapply(c(1:15,17:n_sense),FitInfo,OMfolders=outdirs,sense_codes)),stringsAsFactors=F) # Fitting data
FitDat_2PS<-as.data.frame(t(sapply(16,FitInfo,OMfolders=outdirs,sense_codes)),stringsAsFactors=F) # Fitting data
#FitDat_TVM<-as.data.frame(t(sapply(16,FitInfo,OMfolders=outdirs,sense_codes)),stringsAsFactors=F) # Fitting data

Fit<-readRDS(paste0(outdirs[16],'/Fit.rda')) # OM SSB NAA CAA CAL conv Misc mean_fit data config
Fit@Misc[[1]]['nll_CAA']

#for(i in 1:n_sense)runProjMSEs(i,outdirs)
Refs<-as.data.frame(t(sapply(1:n_sense,RefPoints,OMfolders=outdirs,sense_codes)),stringsAsFactors=F)     # MSY reference points
SSB<-t(sapply(1:n_sense,SSBs,OMfolders=outdirs,sense_codes))                                             # Spawning Stock Biomass extraction
D<-SSB/as.numeric(Refs$SSB0)
B_BMSY<-SSB/as.numeric(Refs$SSBMSY)

#B_BMSY_p<-t(sfSapply(1:n_sense,B_BMSY_proj,OMfolders=outdirs))
#B_BMSY_p_C<-t(sfSapply(1:n_sense,B_BMSY_proj,OMfolders=outdirs,MP="CurC"))

B_BMSY_p<-t(sapply(1:n_sense,B_BMSY_proj,OMfolders=outdirs))
B_BMSY_p_C<-t(sapply(1:n_sense,B_BMSY_proj,OMfolders=outdirs,MP="CurC"))

write.csv(FitDat,"Tables/Step2_Sensitivities/Fitting_data.csv")
write.csv(FitDat_2PS,"Tables/Step2_Sensitivities/Fitting_data_2PS.csv")
write.csv(Refs,"Tables/Step2_Sensitivities/Ref_points.csv")


# === Vizualize sensitivities ==================================================================

jpeg("Images/Step2_Sensitivities/SSB_sense.jpg",res=300,height=8,width=7,units='in')
  sense_comp_ts(SSB,sense_codes)
dev.off()

jpeg("Images/Step2_Sensitivities/Depletion_sense.jpg",res=300,height=8,width=7,units='in')
  sense_comp_ts(D,sense_codes,ylim=c(0,1.5),refline=c(0.5,1))
dev.off()

jpeg("Images/Step2_Sensitivities/B_BMSY_sense.jpg",res=300,height=8,width=7,units='in')
  sense_comp_ts(B_BMSY,sense_codes,ylim=c(0,4.5),refline=c(0.5,1))
dev.off()

jpeg("Images/Step2_Sensitivities/B_BMSY_proj_E.jpg",res=300,height=8,width=7,units='in')
  sense_comp_ts(B_BMSY_p,sense_codes,ylim=c(0,4.5),refline=c(0.5,1),yrproj=T)
dev.off()

jpeg("Images/Step2_Sensitivities/B_BMSY_proj_C.jpg",res=300,height=8,width=7,units='in')
  sense_comp_ts(B_BMSY_p_C,sense_codes,ylim=c(0,4.5),refline=c(0.5,1),yrproj=T)
dev.off()


jpeg("Images/Step2_Sensitivities/LLHprofile_h.jpg",res=300,height=3.7,width=7,units='in')
 
 par(mfrow=c(1,2),mai=c(0.5,0.5,0.01,0.01),omi=c(0.3,0.3,0.01,0.01))
 plot(c(0.95,0.85,0.75,0.7,0.65,0.55,0.9),FitDat$nll[FitDat$OM%in%c(1,5:9,30)],xlab="Steepness",ylab="",type='l')
 legend('topright',c("Profile for steepness","Reference set values"),text.col=c('black','blue'),cex=0.9)
 abline(v=c(0.95,0.70),col='blue',lty=2)
 lim <- c(4050,25+max(as.numeric(FitDat$nll[FitDat$OM%in%c(1:3,5:9,30,21,22)])))
         
 plot(c(0.95,0.85,0.75,0.7,0.65,0.55,0.9),FitDat$nll[c(1,5:9,30)],ylim=lim,xlab="Steepness",ylab="Total Neg. Log. Like.",type='l')
 legend('topright',c("Profile for steepness","Profile for Nat. Mort","Reference set values"),text.col=c('black','blue','red'),cex=0.9)
 
 abline(v=c(0.95,0.7),col='blue',lty=2)
 lines(seq(0.55,0.95,0.45/4),FitDat$nll[FitDat$OM%in%c(1:3,21,22)],col='red')
 points(seq(0.55,0.95,0.45/4),FitDat$nll[FitDat$OM%in%c(1:3,21,22)],col='red',pch=19)
 text(seq(0.55,0.95,0.45/4)+c(0.03,0.03,0.03,0.03,-0.03),as.numeric(FitDat$nll[FitDat$OM%in%c(1:3,21,22)])+c(5,5,5,0,0),FitDat$Code[FitDat$OM%in%c(1:3,21,22)],col='red',cex=0.9)
 mtext("Total Neg. Log. Like.",2,outer=T)
 mtext("Steepness (Bev-Holt SR)",1,outer=T)
 
dev.off()
#dev.new()
