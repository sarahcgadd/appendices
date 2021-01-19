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


(file=paste("renvironment",fileref,".RData",sep=""))
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

################################################################################

###########################################
## COLLECT AVERAGE COEFFICIENTS ###########
###########################################


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

iters<-(ni2-nb2)*4
#extract estimates of coefficients from samples

random.from<-sampcoefs[c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))]
random.from<-unlist(random.from)
random.from<-array(random.from, dim=c(iters,nfrom,(nbasis+1)), dimnames=list(c(1:iters), c(1:nfrom),
                                                                             c("(Intercept)",paste("basis",c(1:nbasis), sep=""))))

random.to<-sampcoefs[c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))]
random.to<-unlist(random.to)
random.to<-array(random.to, dim=c(iters,nto,(nbasis+1)), dimnames=list(c(1:iters), c(1:nto),
                                                                       c("(Intercept)",paste("basis",c(1:nbasis), sep=""))))

fixedboot<-sampcoefs[c("alpha.1",paste("betabs",c(1:nbasis),".1", sep=""))]
fixedboot<-data.frame(matrix(ncol=nbasis+1,data=unlist(fixedboot), byrow=F))

colnames(fixedboot)<- c("(Intercept)",paste("basis",c(1:nbasis), sep=""))

#predict value of model for group and daynumber at x
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
#predict value of differential (with respect to time) for group and daynumber at x                          
pred.dif<-function(x, togrp,from) {
  togrp<-as.character(togrp)
  from<-as.character(from)
  newdata<-as.data.frame(predict(dspline, newx=x))
  colnames(newdata)<-paste(rep("basis", ncol(newdata)), c(1:ncol(newdata)), sep="")
  coefficients<-fixedcoefs
  coefficients[colnames(coefsto)]<-as.numeric(coefsfrom[from,])+as.numeric(coefsto[togrp,])+coefficients[colnames(coefsto)]
  
  return(apply(coefficients[c(2:length(coefficients))]*t(newdata), 2, sum)) 
}

#predict value of second derivative (with respect to time) for group and daynumber at x                          
pred.dif2<-function(x, togrp,from) {
  togrp<-as.character(togrp)
  from<-as.character(from)
  newdata<-as.data.frame(predict(dspline2, newx=x))
  colnames(newdata)<-paste(rep("basis", ncol(newdata)), c(1:ncol(newdata)), sep="")
  coefficients<-fixedcoefs
  coefficients[colnames(coefsto)]<-as.numeric(coefsfrom[from,])+as.numeric(coefsto[togrp,])+coefficients[colnames(coefsto)]
  
  return(apply(coefficients[c(2:length(coefficients))]*t(newdata), 2, sum)) 
}






input<-unique(journeytimes[,c("id","start.station_id","end.station_id")])
input<-input[order(input$id),]
colnames(input)<-c("id","from","to")

input$max.val<-NA
input$max.time<-NA
input$maxmin<-NA


for(i in 1:nrow(input)){
  
  tryCatch(expr={roots<-uniroot.all(pred.dif,  interval=c(start,end)-centre, togrp=input$to[i], from=input$from[i])},
           error={roots<-NULL}
  )
  if(length(roots)>0){
    input$max.val[i]<-max(pred.val(roots, from=input$from[i], togrp=input$to[i]))
    input$max.time[i]<-roots[which.max(pred.val(roots, from=input$from[i], togrp=input$to[i]))]
    min<-min(pred.val(roots, from=input$from[i], togrp=input$to[i]))
    if(min!=input$max.val[i]){
      input$maxmin[i]<-input$max.val[i]-min
    } else{
      input$maxmin[i]<-input$max.val[i]-optimize(pred.val, c(start,end)-centre, from=input$from[i], togrp=input$to[i])$objective 
    }
  } else {
    input$max.val[i]<-NA
    input$max.time[i]<-NA
    input$maxmin[i]<-NA
  }
}

input$max.time<-input$max.time+centre





####################################################

save.image(file=paste("renvironment_2_",fileref,".RData",sep=""))
