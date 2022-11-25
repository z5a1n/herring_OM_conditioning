
i=1 #   i=23

RDF <- data.frame(Year=1968:2018)

for(i in 1:n_sense)
{
In<-readRDS(paste0(outdirs[i],"/In.rda"))
Fit<-readRDS(paste0(outdirs[i],"/Fit.rda"))

maxage <- Fit@OM@maxage
nyears <- Fit@OM@nyears
RecDevs <- colSums(Fit@OM@cpars$Perr_y)/nrow(Fit@OM@cpars$Perr_y)
RecDevs <- RecDevs[(maxage):(maxage+nyears-1)]
Rec <- RecDevs * Fit@OM@cpars$R0[1]
SSB <- colSums(Fit@SSB[,1:nyears])/nrow(Fit@SSB)

assign(paste0("R",i),Rec)
assign(paste0("SSB",i),SSB)

RDF <- cbind(RDF,get(paste0("SSB",i)),get(paste0("R",i)))

temp<-readRDS(paste0(outdirs[i],'/MSEproj.rda'))
mean(temp@Misc$MSYRefs$Refs$SSB0)
Fit@OM@cpars$R0[1]

plot(x=get(paste0("SSB",i)),y=get(paste0("R",i)))
}

write.csv(RDF,"RDF.csv")




index <- data$Index[AS_ind,1]
length(index)

estSSB <- colSums(Fit@SSB[,(nyears-20+1):nyears])/nrow(Fit@SSB)

length(estSSB)

plot(x=estSSB,y=index)

