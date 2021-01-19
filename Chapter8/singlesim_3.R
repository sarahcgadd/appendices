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



##############################
max.val.sim<-matrix(nrow=100, ncol=nrow(input))
max.time.sim<-matrix(nrow=100, ncol=nrow(input))
maxmin.sim<-matrix(nrow=100, ncol=nrow(input))

for(r in 1:100){
bootnum<-(loopnum-1)*100+r
#get random coefficients for boot number
coefsfromboot<-data.frame(random.from[bootnum,,])
colnames(coefsfromboot)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
coefstoboot<-data.frame(random.to[bootnum,,])
colnames(coefstoboot)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
fixedcoefssim<-as.numeric(fixedboot[bootnum,])
names(fixedcoefssim)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))

inputboot<-input[,c("id","from","to") ]

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


for(i in 1:nrow(inputboot)){
  tryCatch(expr={roots<-uniroot.all(pred.dif.boot,  interval=c(start,end)-centre, togrp=inputboot$to[i], from=inputboot$from[i])},
           error={roots<-NULL}
  )
  if(length(roots)>0){
    inputboot$max.val[i]<-max(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))
    inputboot$max.time[i]<-roots[which.max(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))]
    min<-min(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))
    if(min!=inputboot$max.val[i]){
      inputboot$maxmin[i]<-inputboot$max.val[i]-min
    }else{
      input$maxmin[i]<-input$max.val[i]-optimize(pred.val.boot, c(start,end)-centre, from=input$from[i], togrp=input$to[i])$objective 
    }
  }else{
    inputboot$max.val[i]<-NA
    inputboot$max.time[i]<-NA
    inputboot$maxmin[i]<-NA
  }
  
}
inputboot$max.time<-inputboot$max.time+centre

#put your estimate of maximum time and value in a matrix
max.val.sim[r,]<-inputboot$max.val
max.time.sim[r,]<-inputboot$max.time
maxmin.sim[r,]<-inputboot$maxmin
}

write.csv(max.val.sim, paste("maxvalsim_",fileref,"_",loopnum,".csv",sep=""),row.names=F)
write.csv(max.time.sim, paste("maxtimesim_",fileref,"_",loopnum,".csv",sep=""),row.names=F)
write.csv(maxmin.sim, paste("maxminsim_",fileref,"_",loopnum,".csv",sep=""),row.names=F)
