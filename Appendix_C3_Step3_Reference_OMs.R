# ===========================================================================================================================
# === 3. Specifying Reference Case OMs ======================================================================================
# ===========================================================================================================================
#
# Atlantic Herring MSE
#
# March 22nd 2020
# Tom Carruthers t.carruthers@oceans.ubc.ca
#
# This script:
# a) Loads reference set input modifiers (to make inputs for reference set OMs)
# b) Makes the various input files for the reference set OMs
# c) Fits the reference set OMs
# d) Extracts summary data about model fit and estimates and visualizes reference set OM quantities
#
# The product of this script: the Reference set operating models

# === Prerequisites =========================================================================================
rm(list=ls())
library(MSEtool)
library(DLMtool)

# --- Directories ---

OMdir<-"C:/Users/BARRETTTJ/Desktop/MSE"
setwd(OMdir)

# --- Source Code

nofit<-TRUE   
source('Appendix_C1_Step1_Reference_Case_OM.R') # Contains various data objects and the reference OM required for reference grid OMs
source('Reference_mods.R') # Contains various data objects and the reference OM required for reference grid OMs


# === Make Reference OM grid Input objects =========================================================================================================

# --- Design of the Reference Grid ------------------------------

Design<-expand.grid(1:3,c("A","B"),c("H","L"),c("-","+"),stringsAsFactors=F)
names(Design)<-c("Mort","Growth","Resil","Weir")
saveRDS(Design,paste0(OMdir,"/Reference_OMs/Design.rda"))
write.csv(Design,paste0(OMdir,"/Tables/Step3_Reference_Set/Design.csv"))


# --- Make OM input objects

OMcode<-apply(Design,1,FUN=function(x)paste(x,collapse=" "))
nOMs<-length(OMcode)
OMfolders<-paste0(OMdir,"/Reference_OMs/",1:nOMs)
verbose<-F

for(i in 1:nOMs){
  
  print(paste("----------",i,":",OMcode[i],"--------------------"))
  if(!file.exists(OMfolders[i]))dir.create(OMfolders[i])
  In<-In_RC # copy the reference case input object as a template
  In<-Ref_M(In,lev=Design[i,1],verbose=verbose)  # Factor 1, M, Natural mortality rate
  In<-Ref_G(In,lev=Design[i,2],verbose=verbose)  # Factor 2, G, Projected growth
  In<-Ref_R(In,lev=Design[i,3],verbose=verbose)  # Factor 3, R, Resilience (steepness of the Stock-Recruitment relationship)
  In<-Ref_C(In,lev=Design[i,4],verbose=verbose)  # Factor 4, C, Catch of the weir fishery (0 or 100%)
  saveRDS(In,file=paste0(OMfolders[i],'/In.rda'))

}


# === Fit Reference OMs =========================================================================================================================

#setup()
#sfExport('SRA_wrap')
#system.time({sfSapply(1:nOMs,SRA_parallel,OMfolders=OMfolders)}) # takes about 2 minutes
system.time({
for(i in 1:nOMs) { # takes about 20 minutes
  In<-readRDS(paste0(OMfolders[i],"/In.rda"))
  Fit<-SRA_wrap(In)
  Fit<-Add_Regime(Fit,lev=2,verbose=T)    ##add regime shift for recruitment
  saveRDS(Fit,paste0(OMfolders[i],"/Fit.rda"))
} 
})  
# === Build OM reports =========================================================================

for(i in 1:nOMs){ # takes about 20 minutes
  In<-readRDS(paste0(OMfolders[i],"/In.rda"))
  Fit<-readRDS(paste0(OMfolders[i],"/Fit.rda"))
  keep<-unlist(lapply(Fit@Misc,function(x)all(!is.na(x$F))))
  Fit@OM<-Sub_cpars(Fit@OM,sims=keep)
  Fit@OM@EffLower<-rep(1,OM@nyears)
  Fit@OM@EffUpper<-rep(2,OM@nyears)
  plot_wrap(Fit, In, dir=OMfolders[i], filename="/OM_Report",  title=paste("Reference Case Operating Model Fitting Report:",OMcode[i]))
}


# === Gather diagnostics =======================================================================================================================

FitDat<-as.data.frame(t(sapply(1:nOMs,FitInfo,OMfolders=OMfolders,OMcode)),stringsAsFactors=F)
#FitDat[,c(4:ncol(FitDat))] <- as.numeric(unlist(FitDat[,c(4:ncol(FitDat))]))

sapply(1:nOMs,runProjMSEs,OMfolders=OMfolders) # Run projection MSEs with current effort (curE MP) and current catch (CurC MP), takes about 2 minutes
Refs<-as.data.frame(t(sapply(1:nOMs,RefPoints,OMfolders=OMfolders,OMcode)),stringsAsFactors=F) 
#Refs[,c(3:ncol(Refs))] <- as.numeric(unlist(Refs[,c(3:ncol(Refs))]))

SSB<-t(sapply(1:nOMs,SSBs,OMfolders=OMfolders,OMcode))
D<-SSB/as.numeric(Refs$SSB0)
B_BMSY<-SSB/as.numeric(Refs$SSBMSY)


B_BMSY_p<-t(sapply(1:nOMs,B_BMSY_proj,OMfolders=OMfolders))
B_BMSY_p_C<-t(sapply(1:nOMs,B_BMSY_proj,OMfolders=OMfolders,MP="CurC"))

Refs[,c(3:ncol(Refs))] <- as.numeric(unlist(Refs[,c(3:ncol(Refs))]))

Stoch<-sapply(1:nOMs,stoch_quants,OMfolders=OMfolders)

write.csv(FitDat,"Tables/Step3_Reference_Set/Fitting_data.csv")
Refs$SSB_SSBMSY<-Refs$D/Refs$SSBMSY_SSB0
write.csv(Refs,"Tables/Step3_Reference_Set/Ref_points.csv")


# === Vizualize outcomes =====================================================================

# Model estimates

jpeg("Images/Step3_Reference_Set/SSB_refOM.jpg",res=300,height=5,width=7,units='in')
  plotcomp(SSB,Design)
dev.off()

jpeg("Images/Step3_Reference_Set/Depletion_refOM.jpg",res=300,height=5,width=7,units='in')
  plotcomp(D,Design,refline=1)
dev.off()

jpeg("Images/Step3_Reference_Set/B_BMSY_refOM.jpg",res=300,height=5,width=7,units='in')
  plotcomp(B_BMSY,Design,refline=c(0.5,1),ylim=c(0,4))
dev.off()

# Status quo projections

jpeg("Images/Step3_Reference_Set/B_BMSY_proj_E.jpg",res=300,height=5,width=7,units='in')
  plotcomp(B_BMSY_p,Design,ylim=c(0,5.5),refline=c(0.5,1),yrproj=T)
dev.off()

jpeg("Images/Step3_Reference_Set/B_BMSY_proj_C.jpg",res=300,height=5,width=7,units='in')
  plotcomp(B_BMSY_p_C,Design,ylim=c(0,5.5),refline=c(0.5,1),yrproj=T)
dev.off()

# Stochastic estimates

jpeg("Images/Step3_Reference_Set/Stoch.jpg",res=300,height=10,width=7,units='in')
  plot_stoch(Stoch,Design)
dev.off()

# Selected Robustness OMs


jpeg("Images/Step3_Reference_Set/RobGrid.jpg",res=300,height=5,width=7,units='in')
   plot_RobGrid(OMcode,Selected=c("1 A H -","3 B H +","2 A L -"),B_BMSY,B_BMSY_p,ylim=c(0,8))
dev.off()




# ====================================================================================================================
# ======= END OF SCRIPT ==============================================================================================
# ====================================================================================================================
