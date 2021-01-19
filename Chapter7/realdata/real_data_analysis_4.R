#this was generated using a bash script, it brings together the results from the parallel loop to estimate 95% credible intervals

load(paste('renvironment_2_',fileref,'.RData',sep=''))
library(plyr)
library(tidyr)
library(dplyr)
library(sf)
library(chron)
library(multiplex)
library(R2OpenBUGS)
library(igraph)
library(splines2)
library(rootSolve)
library(ggplot2)
library(stringr)
library(Deriv)
library(ggpubr)

#collate all the information for 95% CIs
avspeed<-read.csv(paste("avspeed_",fileref,"_",r,".csv",sep=""))
max.val.sim<-read.csv(paste("maxvalsim_",fileref,"_",r,".csv",sep=""))
max.time.sim<-read.csv( paste("maxtimesim_",fileref,"_",r,".csv",sep=""))
maxmin.sim<-read.csv( paste("maxminsim_",fileref,"_",r,".csv",sep=""))



for(r in 2:boots){
  avspeed<-rbind(avspeed,read.csv(paste("avspeed_",fileref,"_",r,".csv",sep="")))
  max.val.sim<-rbind(max.val.sim,read.csv(paste("maxvalsim_",fileref,"_",r,".csv",sep="")))
  max.time.sim<-rbind(max.time.sim,read.csv( paste("maxtimesim_",fileref,"_",r,".csv",sep="")))
  maxmin.sim<-rbind(maxmin.sim,read.csv( paste("maxminsim_",fileref,"_",r,".csv",sep="")))
}

#get estimate of 95% CIs for each group
input$max.val.975<-apply(max.val.sim, 2, quantile, p=0.975, na.rm=T)
input$max.time.975<-apply(max.time.sim, 2, quantile, p=0.975, na.rm=T)
input$max.val.025<-apply(max.val.sim, 2, quantile, p=0.025, na.rm=T)
input$max.time.025<-apply(max.time.sim, 2, quantile, p=0.025, na.rm=T)
input$maxmin.025<-apply(maxmin.sim, 2, quantile, p=0.025, na.rm=T)
input$maxmin.974<-apply(maxmin.sim, 2, quantile, p=0.975, na.rm=T)


#average speed CIs
avspeed.conf<-data.frame(time=seq(from=start, to=end, length.out=100),
                         qlow=apply(avspeed, 2, quantile, p=0.025, na.rm=T),
                         qhigh=apply(avspeed, 2, quantile, p=0.975, na.rm=T))





#plots
avspeed_mlm<-ggplot(data=avspeed.conf)+
  geom_line(aes(x=time, y=qlow), linetype="dashed", colour="grey30")+
  geom_line(aes(x=time, y=qhigh), linetype="dashed", colour="grey30")+
  stat_function(fun=averagespeed,xlim=c(start, end), lwd=1, colour="black")+
  theme_minimal()+
  xlab(label="Time")+
  ylab(label="Average journey speed (km/h)")+
  scale_x_continuous(limits=c(14,21), breaks=c(14,15,16,17,18,19,20,21), labels=paste(c(14:21),"00",sep=":"))

avspeed_noci<-ggplot(data=avspeed.conf)+
  stat_function(fun=averagespeed, lwd=1, xlim=c(start,end),colour="black")+
  theme_minimal()+
  xlab(label="Time")+
  ylab(label="Average journey speed (km/h)")+
  scale_x_continuous(limits=c(14,21), breaks=c(14,15,16,17,18,19,20,21), labels=paste(c(14:21),"00",sep=":"))

timedensity<-ggplot(data=journeytimes, aes(x=start.time))+
  geom_density()+
  theme_minimal()+
  xlab(label="Time")+
  ylab(label="Density")+
  scale_x_continuous(limits=c(14,21), breaks=c(14,15,16,17,18,19,20,21), labels=paste(c(14:21),"00",sep=":"))


width=119
height=195

png(paste("speedfunction_noci",fileref, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
ggarrange(avspeed_mlm, avspeed_noci, ncol=2, nrow=1)
dev.off()
png(paste("speedfunction_timedensity",fileref, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
ggarrange(avspeed_mlm+theme(axis.text.x=element_blank(), axis.title.x=element_blank()),timedensity , ncol=1, nrow=2)
dev.off()

png(paste("speedfunction",fileref, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
avspeed_mlm
dev.off()
png(paste("timedensity",fileref, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
timedensity
dev.off()

png(paste("maxpoints",fileref, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
ggplot(data=input, aes(colour=as.factor(id)))+
  geom_point(aes(x=max.time, y=max.val))+
  geom_errorbar(aes(x=max.time, ymin=max.val.025, ymax=max.val.975))+
  geom_errorbarh(aes(y=max.val, xmin=max.time.025, xmax=max.time.975))+
  theme_minimal()+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(14,21), breaks=c(14,15,16,17,18,19,20,21), labels=paste(c(14:21),"00",sep=":"))+
  xlab("Time of maximum journey length")+
  ylab("Maximum journey length (hours)")
dev.off()

basisplot<-data.frame(time=journeytimes$start.time, basis)
colnames(basisplot)<-c("time",paste("basis", c(1:nbasis), sep="."))
basisplot<-reshape(basisplot, direction="long", varying=paste("basis",c(1:nbasis),sep="."), timevar="basisno",idvar="measure")
png(paste("basisplot",fileref,".png", sep="_"),width=width, height=width*5/7, units="mm",res=600)
ggplot(data=basisplot, aes(x=time, y=basis, group=basisno, colour=as.factor(basisno), linetype=as.factor(basisno)), lwd=1)+
  geom_line()+
  theme_minimal()+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(14,21), breaks=c(14,15,16,17,18,19,20,21), labels=paste(c(14:21),"00",sep=":"))+
  xlab(label="Time")+
  ylab(label="Basis Function Value")
dev.off()



write.csv(input, paste("pattern_feature_estimates",fileref,".csv",sep=""))

#save fixed effects coefficients
sdcoefs<-bugs.out$sd
meancoefslengths<-sapply(meancoefs, length)
fixedmean<-unlist(meancoefs[meancoefslengths==1])
fixedsd<-unlist(sdcoefs[meancoefslengths==1])

write.csv(rbind(fixedmean, fixedsd),  paste("model_coefficients",fileref,".csv",sep=""), row.names=T)
