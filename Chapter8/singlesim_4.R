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


load(paste('renvironment_2_',fileref,'.RData',sep=''))

max.val.sim<-read.csv(paste("maxvalsim_",fileref,"_",1,".csv",sep=""))
max.time.sim<-read.csv( paste("maxtimesim_",fileref,"_",1,".csv",sep=""))
maxmin.sim<-read.csv( paste("maxminsim_",fileref,"_",1,".csv",sep=""))

for(r in 2:40){
  max.val.sim<-rbind(max.val.sim,read.csv(paste("maxvalsim_",fileref,"_",r,".csv",sep="")))
  max.time.sim<-rbind(max.time.sim,read.csv( paste("maxtimesim_",fileref,"_",r,".csv",sep="")))
  maxmin.sim<-rbind(maxmin.sim,read.csv( paste("maxminsim_",fileref,"_",r,".csv",sep="")))
}

#get estimate for each group
input$max.val.975<-apply(max.val.sim, 2, quantile, p=0.975, na.rm=T)
input$max.time.975<-apply(max.time.sim, 2, quantile, p=0.975, na.rm=T)
input$max.val.025<-apply(max.val.sim, 2, quantile, p=0.025, na.rm=T)
input$max.time.025<-apply(max.time.sim, 2, quantile, p=0.025, na.rm=T)
input$maxmin.025<-apply(maxmin.sim, 2, quantile, p=0.025, na.rm=T)
input$maxmin.975<-apply(maxmin.sim, 2, quantile, p=0.975, na.rm=T)
input$max.val.sd<-apply(max.val.sim, 2, sd,  na.rm=T)
input$max.time.sd<-apply(max.time.sim, 2, sd,  na.rm=T)
input$maxmin.sd<-apply(maxmin.sim, 2, sd, na.rm=T)

#real values

#collate information from real simulation
a<-simdat$section1.function
input.real<-unique(journeytimes[,c("id","start.station_id","end.station_id")])
input.real<-input.real[order(input.real$id),]
colnames(input.real)<-c("id","from","to")
input.real<-cbind(input.real,max.val=NA, max.time=NA, maxmin=NA)
b<-list(length=length(a))

for (i in 1:length(a)){
  b[[i]]<-Deriv(a[[i]])
  roots<-uniroot.all(b[[i]], c(start,end))
  input.real$max.val[i]<-max(a[[i]](roots))/60
  input.real$max.time[i]<-roots[which.max(a[[i]](roots))]
  min<-min(a[[i]](roots))/60
  if(min!=input.real$max.val[i]){
    input.real$maxmin[i]<-(input.real$max.val[i]-min)
  } else{
    input.real$maxmin[i]<-(input.real$max.val[i]-optimize(a[[i]], c(start,end))$objective )
  }
  
}





colnames(input.real)<-paste(colnames(input.real), "real", sep=".")
savedata<-merge(input, input.real, by.x=c("from","to"), by.y=c("from.real","to.real"))
savedata$max.timedif<-savedata$max.time-savedata$max.time.real
savedata$max.valdif<-savedata$max.val-savedata$max.val.real
savedata$maxmin.dif<-savedata$maxmin-savedata$maxmin.real

savedata$max.val.ci.width<-savedata$max.val.975-savedata$max.val.025
savedata$max.time.ci.width<-savedata$max.time.975-savedata$max.time.025
savedata$maxmin.ci.width<-savedata$maxmin.975-savedata$maxmin.025


savedata$max.time.in.ci<-savedata$max.time.975>=savedata$max.time.real & savedata$max.time.025<=savedata$max.time.real
savedata$max.val.in.ci<-savedata$max.val.975>=savedata$max.val.real & savedata$max.val.025<=savedata$max.val.real
savedata$maxmin.in.ci<-savedata$maxmin.975>=savedata$maxmin.real & savedata$maxmin.025<=savedata$maxmin.real


savedata$max.realval.as.z<-(savedata$max.val.real-savedata$max.val)/savedata$max.val.sd
savedata$max.realtime.as.z<-(savedata$max.time.real-savedata$max.time)/savedata$max.time.sd
savedata$real.maxmin.as.z<-(savedata$maxmin.real-savedata$maxmin)/savedata$maxmin.sd



savedata$fileref<-fileref

result<- savedata
