#this file was generated using a bash script, 4000 were generated for each tau parameter
#it goes through the loop process of estimating 95% credible itnervals for pattern features, but with each loop in parallel (through a task array)

fileref<-1
r<-1

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

#set up matrix for estimates (if in a loop, this would have more rows)
max.val.sim<-matrix(nrow=1, ncol=nrow(input))
max.time.sim<-matrix(nrow=1, ncol=nrow(input))
maxmin.sim<-matrix(nrow=1, ncol=nrow(input))
avspeed<-matrix(nrow=0, ncol=100)


#get random coefficients for loop r
coefsfromboot<-data.frame(random.from[r,,])
colnames(coefsfromboot)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
coefstoboot<-data.frame(random.to[r,,])
colnames(coefstoboot)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
fixedcoefssim<-as.numeric(fixedboot[r,])
names(fixedcoefssim)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))

inputboot<-input[,c("from","to","distancekm") ]

#function to predict value for group and daynumber at x
pred.val.boot<-function(x, togrp, from) {
  togrp<-as.character(togrp)
  from<-as.character(from)
  newdata<-as.data.frame(predict(basis, newx=x))
  colnames(newdata)<-paste("basis", c(1:ncol(newdata)), sep="")
  coefficients<-fixedcoefssim
  coefficients[colnames(coefstoboot)]<-as.numeric(coefsfromboot[from,])+as.numeric(coefstoboot[togrp,])+coefficients[colnames(coefstoboot)]
  result<-apply(coefficients[2:length(coefficients)]*t(newdata), 2, sum)
  return(result<-result+coefficients[1])
}
#predict value of differential (with respect to time) for group and daynumber at x                          
pred.dif.boot<-function(x, togrp,from) {
  togrp<-as.character(togrp)
  from<-as.character(from)
  newdata<-as.data.frame(predict(dspline, newx=x))
  colnames(newdata)<-paste(rep("basis", ncol(newdata)), c(1:ncol(newdata)), sep="")
  coefficients<-fixedcoefssim
  coefficients[colnames(coefstoboot)]<-as.numeric(coefsfromboot[from,])+as.numeric(coefstoboot[togrp,])+coefficients[colnames(coefstoboot)]
  
  return(apply(coefficients[c(2:length(coefficients))]*t(newdata), 2, sum)) 
}

inputboot$max.val<-NA
inputboot$max.time<-NA
inputboot$maxmin<-NA

#loop over each od pair and save pattern feature estimates for this model parameter sample
for(i in 1:nrow(inputboot)){
  
  roots<-uniroot.all(pred.dif.boot,  interval=c(start,end)-centre, togrp=inputboot$to[i], from=inputboot$from[i])
  if(length(roots)>0){
    inputboot$max.val[i]<-max(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))
    inputboot$max.time[i]<-roots[which.max(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))]
    min<-min(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))
    if(min!=inputboot$max.val[i]){
      inputboot$maxmin[i]<-inputboot$max.val[i]-min
    }else{
      input$maxmin[i]<-input$max.val[i]-optimize(pred.val.boot, c(start,end)-centre, from=input$from[i], togrp=input$to[i])$objective 
    }
  }
  
}
inputboot$max.time<-inputboot$max.time+centre

#put your estimate of maximum time and value in a matrix
max.val.sim[1,]<-inputboot$max.val
max.time.sim[1,]<-inputboot$max.time
maxmin.sim[1,]<-inputboot$maxmin

#calculate average speed
averagespeed<-function(time){
  nlinks<-nrow(inputboot)
  timecent<-time-centre
  output<-rep(0, length(time))
  for (i in 1:nrow(inputboot)){
    distance<-inputboot$distancekm[i]
    output<-output+distance/pred.val.boot(x=timecent, from=inputboot$from[i], togrp=inputboot$to[i])
  }
  output<-output/(nlinks)
  return(output)
}

avspeed<-rbind(avspeed, averagespeed(time=seq(from=start, to=end, length.out=100)))

#save all the information to collate later  
write.csv(avspeed, paste("avspeed_",fileref,"_",r,".csv",sep=""),row.names=F)

write.csv(max.val.sim, paste("maxvalsim_",fileref,"_",r,".csv",sep=""),row.names=F)
write.csv(max.time.sim, paste("maxtimesim_",fileref,"_",r,".csv",sep=""),row.names=F)
write.csv(maxmin.sim, paste("maxminsim_",fileref,"_",r,".csv",sep=""),row.names=F)



