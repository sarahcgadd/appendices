#this file was generated using a bash script - this brings together the estimates from the 4 parallel chains

fileref<-1

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

#load bugs objects

load(file=paste("renvironment",fileref,".RData",sep=""))
for(i in 2:4){
  load(paste("bugsout_",fileref,"_",i,".RData",sep=""))
}

###########################################
## COLLECT AVERAGE COEFFICIENTS ###########
###########################################

meancoefs1<-bugs.out$mean
meancoefs2<-bugs.out2$mean
meancoefs3<-bugs.out3$mean
meancoefs4<-bugs.out4$mean

meancoefs<-list()
for(i in 1:length(meancoefs1)){
  meancoefs[[i]]<-rbind(meancoefs1[[i]],meancoefs2[[i]],meancoefs3[[i]],meancoefs4[[i]])
  meancoefs[[i]]<-apply(meancoefs[[i]],2,mean)
}
names(meancoefs)<-names(meancoefs1)

sampcoefs1<-bugs.out$sims.list
sampcoefs2<-bugs.out2$sims.list
sampcoefs3<-bugs.out3$sims.list
sampcoefs4<-bugs.out4$sims.list

sampcoefs<-list()
for(i in 1:length(sampcoefs1)){
  if(class(sampcoefs1[[i]])=="numeric"){
    sampcoefs[[i]]<-c(sampcoefs1[[i]],sampcoefs2[[i]],sampcoefs3[[i]],sampcoefs4[[i]])
  } else{
    sampcoefs[[i]]<-rbind(sampcoefs1[[i]],sampcoefs2[[i]],sampcoefs3[[i]],sampcoefs4[[i]])
    
  }
}
names(sampcoefs)<-names(sampcoefs1)

nfrom<-length(unique(journeytimes$start.station_id))
nto<-length(unique(journeytimes$end.station_id))

coefsfrom<-meancoefs[c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))]
coefsfrom<-data.frame(matrix(nrow=nfrom, data=unlist(coefsfrom), byrow=F))
colnames(coefsfrom)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
rownames(coefsfrom)<-c(1:nfrom)

coefsto<-meancoefs[c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))]
coefsto<-data.frame(matrix(nrow=nto, data=unlist(coefsto), byrow=F))
colnames(coefsto)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
rownames(coefsto)<-c(1:nto)

fixedcoefs<-meancoefs[c("alpha.1",paste("betabs",c(1:nbasis),".1", sep=""))]
fixedcoefs<-unlist(fixedcoefs)
names(fixedcoefs)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))




########################################
## PROCESS TO GET 95% CRED INTERVALS  ##
########################################

sample<-sample(c(1:((ni2-nb2)*4)), boots)
iters<-(ni2-nb2)*4
#extract estimates of coefficients from samples

random.from<-sampcoefs[c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))]
random.from<-unlist(random.from)
random.from<-array(random.from, dim=c(iters,nfrom,(nbasis+1)), dimnames=list(c(1:iters), c(1:nfrom),
                                                                             c("(Intercept)",paste("basis",c(1:nbasis), sep=""))))
random.from<-random.from[sample,,]

random.to<-sampcoefs[c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))]
random.to<-unlist(random.to)
random.to<-array(random.to, dim=c(iters,nto,(nbasis+1)), dimnames=list(c(1:iters), c(1:nto),
                                                                       c("(Intercept)",paste("basis",c(1:nbasis), sep=""))))
random.to<-random.to[sample,,]

fixedboot<-sampcoefs[c("alpha.1",paste("betabs",c(1:nbasis),".1", sep=""))]
fixedboot<-data.frame(matrix(ncol=nbasis+1,data=unlist(fixedboot), byrow=F))
fixedboot<-fixedboot[sample,]
colnames(fixedboot)<- c("(Intercept)",paste("basis",c(1:nbasis), sep=""))



#function to predict value of model for group and daynumber at x
pred.val<-function(x, togrp, from) {
  togrp<-as.character(togrp)
  from<-as.character(from)
  newdata<-as.data.frame(predict(basis, newx=x))
  colnames(newdata)<-paste(rep("basis", ncol(newdata)), c(1:ncol(newdata)), sep="")
  coefficients<-fixedcoefs
  coefficients[colnames(coefsto)]<-as.numeric(coefsfrom[from,])+as.numeric(coefsto[togrp,])+coefficients[colnames(coefsto)]
  result<-apply(coefficients[2:length(coefficients)]*t(newdata), 2, sum)
  return(result<-result+coefficients[1])
}
#function to predict value of differential (with respect to time) for group and daynumber at x                          
pred.dif<-function(x, togrp,from) {
  togrp<-as.character(togrp)
  from<-as.character(from)
  newdata<-as.data.frame(predict(dspline, newx=x))
  colnames(newdata)<-paste(rep("basis", ncol(newdata)), c(1:ncol(newdata)), sep="")
  coefficients<-fixedcoefs
  coefficients[colnames(coefsto)]<-as.numeric(coefsfrom[from,])+as.numeric(coefsto[togrp,])+coefficients[colnames(coefsto)]
  
  return(apply(coefficients[c(2:length(coefficients))]*t(newdata), 2, sum)) 
}

###########################
#function for average speed

averagespeed<-function(time){
  nlinks<-nrow(input)
  timecent<-time-centre
  output<-rep(0, length(time))
  for (i in 1:nrow(input)){
    distance<-input$distancekm[i]
    output<-output+distance/pred.val(x=timecent, from=input$from[i], togrp=input$to[i])
  }
  output<-output/(nlinks)
  return(output)
}


#add distance to journeytimes
distance<-function(x){
  rad <- pi/180
  lat1<-as.numeric(x["start.lat"])
  lat2<-as.numeric(x["end.lat"])
  long1<-as.numeric(x["start.long"])
  long2<-as.numeric(x["end.long"])
  a1 <- lat1*rad
  a2 <- long1*rad
  b1 <- lat2*rad
  b2 <- long2*rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1)*cos(b1)*(sin(dlon/2))^2
  c <- 2*atan2(sqrt(a), sqrt(1 - a))
  R <- 6378137
  return( R*c)
}
journeytimes$distancekm<-apply(journeytimes, 1, distance)/1000



input<-unique(journeytimes[,c("start.station_id","end.station_id","distancekm")])
input<-cbind(c(1:nrow(input)), input)
colnames(input)<-c("id","from","to","distancekm")
input$max.val<-NA
input$max.time<-NA
input$maxmin<-NA

#collect average pattern features
for(i in 1:nrow(input)){
  
  roots<-uniroot.all(pred.dif,  interval=c(start,end)-centre, togrp=input$to[i], from=input$from[i])
  if(length(roots)>0){
    input$max.val[i]<-max(pred.val(roots, from=input$from[i], togrp=input$to[i]))
    input$max.time[i]<-roots[which.max(pred.val(roots, from=input$from[i], togrp=input$to[i]))]
    min<-min(pred.val(roots, from=input$from[i], togrp=input$to[i]))
    if(min!=input$max.val[i]){
      input$maxmin[i]<-input$max.val[i]-min
    } else{
      input$maxmin[i]<-input$max.val[i]-optimize(pred.val, c(start,end)-centre, from=input$from[i], togrp=input$to[i])$objective 
    }
  }
}
input$max.time<-input$max.time+centre

#save the image for next file
save.image(file=paste("renvironment_2_",fileref,".RData",sep=""))
