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






#############################################################################################
# RUN SECOND ANALYSIS WITH PRIORS & inits BASED ON FIRST MODEL + TEMPORAL COR
################################################################################

bugsdata<-list(length=journeytimes$length,
               SID=journeytimes$start.station_id, EID=journeytimes$end.station_id,
               adj1=adj1, weights1=weights1, num1=num1,
               Nobs=nrow(journeytimes), Nstart=length(unique(journeytimes$start.station_id)),
               Nend=length(unique(journeytimes$end.station_id)),timedif=journeytimes$timedif, timesetup=journeytimes$timesetup)

a<-length(bugsdata)
for(i in 1:nbasis){
  bugsdata[[i+a]]<-as.numeric(basis[,i])
}
names(bugsdata)[(a+1):(a+nbasis)]<-paste("bs",c(1:nbasis), sep="")

initials<-bugs.out$mean

inits<-function(){l<-initials[1:(length(initials)-1)]
l<-lapply(l,as.numeric)
l<-lapply(l, function(x){if(length(x)==117){return(c(x[1:19],0,x[20:117]))}else{return(x)}})
for(d in names(l)[which(str_detect(names(l),"space")==T & str_detect(names(l),"tau")==F)]){
  l[[d]]<-rep(0, length(l[[d]]))
}
l$acparam = 0.5
return(l)
}

fixedmeans<-list()
fixedsd<-list()

coefs<-c(paste("alpha.",c(1),sep=""), paste(paste("betabs",c(1:nbasis),sep=""),".",c(1),sep=""))
for(i in coefs){
  eval(parse(text=paste("fixedmeans$",i,"<-bugs.out$mean[['",i,"']]",sep="")))
  eval(parse(text=paste("fixedsd$",i,"<-bugs.out$sd[['",i,"']]",sep="")))
  
}

sink(paste("model",fileref,"_",chainnum,".txt",sep=""))
cat("model{


       
	      
	      length[1]~dnorm(a[1],tau.a)
		 
        a[1]<- mu[1] + ")

for(i in 1:(nbasis-1)){
  cat(paste("mu",i,"[1] + ",sep=""))
}
cat(paste("mu",nbasis,"[1]\n",sep=""))

cat("
    mu[1] <- alpha.1 + alpha.2[SID[1]] + alpha.3[EID[1]]
    ")




for (i in 1:nbasis){
  for(j in c(".1*bs", ".2[SID[1]]*bs", ".3[EID[1]]*bs")){
    if(j==".1*bs"){
      cat(paste("mu",i,"[1] <- ",sep=""))
    }
    if(j==".3[EID[1]]*bs"){
      cat(paste("betabs", i,j,i,"[1]\n", sep=""))
    }else{
      cat(paste("betabs", i,j,i,"[1] + ", sep=""))
    }
  }
}


cat("
    for(i in 2:Nobs){
   
        arcor[i]<-acparam*timesetup[i]/(timedif[i]+1)
    length[i]~dnorm(a[i],tau.a) 
        a[i] <-arcor[i]*length[i-1]+  mu[i]-arcor[i]*mu[i] + ")

for(i in 1:(nbasis-1)){
  cat(paste("mu",i,"[i] + ",sep=""))
}
cat(paste("mu",nbasis,"[i]\n",sep=""))

cat("
    mu[i] <- alpha.1 + alpha.2[SID[i]] + alpha.3[EID[i]]
    ")




for (i in 1:nbasis){
  for(j in c(".1", ".2[SID[i]]", ".3[EID[i]]")){
    if(j==".1"){
      cat(paste("mu",i,"[i] <- ",sep=""))
    }
    if(j==".3[EID[i]]"){
      cat(paste("betabs",i,j,"*bs",i,"[i]-arcor[i]*betabs", i,j,"*bs",i,"[i-1]\n", sep=""))
    }else{
      cat(paste("betabs",i,j,"*bs",i,"[i]-arcor[i]*betabs", i,j,"*bs",i,"[i-1] + ", sep=""))
    }
  }
}

cat("}
   
    for (j in 1:Nstart){
    alpha.2[j]<-alpha.2.space[j]+alpha.2.unstr[j]
    alpha.2.unstr[j]~dnorm(0,tau.int.start.unstr)
    ")
for ( i in 1:nbasis){
  cat(paste("betabs",i,".2[j]<-betabs",i,".2.space[j]+betabs",i,".2.unstr[j]\n", sep=""))
  cat(paste("betabs",i,".2.unstr[j]~dnorm(0, tau.bs",i,".start.unstr)\n",sep=""))
}


cat("}
    for (k in 1:Nend){
    alpha.3[k]<-alpha.3.space[k]+alpha.3.unstr[k]
    alpha.3.unstr[k]~dnorm(0,tau.int.end.unstr)
    ")
for ( i in 1:nbasis){
  cat(paste("betabs",i,".3[k]<-betabs",i,".3.space[k]+betabs",i,".3.unstr[k]\n", sep=""))
  cat(paste("betabs",i,".3.unstr[k]~dnorm(0, tau.bs",i,".end.unstr)\n",sep=""))
}

cat("}
    
    #set up all random intercept and basis coefficients following car distribution given spatial (weights1) or network (weights3) weights
    alpha.2.space[1:Nstart]~car.normal(adj1[], weights1[], num1[], tau.int.start.space)
    
    ")
for (i in 1:nbasis){
  for (j in c("space")){
    
    if(j=="space"){w<-1}else{w<-3}
    cat(paste("betabs",i,".2.",j,"[1:Nstart]~car.normal(adj",w,"[], weights",w,"[], num",w,"[], tau.bs",i,".start.",j,")\n", sep=""))
  }
}

cat("alpha.3.space[1:Nend]~car.normal(adj1[], weights1[], num1[], tau.int.end.space)
    
    ")

for (i in 1:nbasis){
  for (j in c("space")){
    
    if(j=="space"){w<-1}else{w<-3}
    cat(paste("betabs",i,".3.",j,"[1:Nend]~car.normal(adj",w,"[], weights",w,"[], num",w,"[], tau.bs",i,".end.",j,")\n", sep=""))
  }
}


cat("
    ### priors on regression coefficients and variances
    
    tau.a ~ dgamma(0.5, 0.005)
    sigma2.a <- 1/tau.a # residual error variance
    

    ### residual intercept variance
    tau.int.start.space ~ dgamma(0.5, 0.005)
    sigma2.int.start.space <- 1/tau.int.start.space #between start station, with spatial cor
    tau.int.start.unstr ~ dgamma(0.5, 0.005)
    sigma2.int.start.unstr <- 1/tau.int.start.unstr
    
    tau.int.end.space ~ dgamma(0.5, 0.005)
    sigma2.int.end.space <- 1/tau.int.end.space #between end station, with spatial cor
    tau.int.end.unstr ~ dgamma(0.5, 0.005)
    sigma2.int.end.unstr <- 1/tau.int.end.unstr
    
    ### variance in basis 1 coefficients
    ")

for(i in 1:nbasis){
  for(j in c("start","end")){
    for (k in c("space","unstr")){
      cat(paste("tau.bs",i,".",j,".",k," ~ dgamma(0.5, 0.005)\nsigma2.bs",i,".",j,".",k," <- 1/tau.bs",i,".",j,".",k,"\n", sep=""))
    }
  }
}



cat( paste("alpha.1 ~ dnorm(",gsub("e","E",as.character(fixedmeans$alpha.1)),",",gsub("e","E",as.character(1/(fixedsd$alpha.1^2))),") # intercept
     ",sep=""))
for(i in 1:nbasis){
  cat(paste("betabs",i,".1 ~ dnorm(",gsub("e","E",as.character(fixedmeans[[paste("betabs",i,".1",sep="")]])),",",gsub("e","E",as.character(1/(fixedsd[[paste("betabs",i,".1",sep="")]]^2))),")\n", sep=""))
}


cat("acparam ~ dgamma(1,1)I(0,1)
    
    }")
sink()   

sigmas<-expand.grid(c("sigma2.bs"),c(1:nbasis),c(".start",".end"),c(".space",".unstr"))
sigmas<-apply(sigmas, 1, paste, collapse="")
betas<-expand.grid(c("betabs"),c(1:nbasis),c(".2",".3"),c(".space",".unstr"))
betas<-apply(betas, 1, paste, collapse="")
sigmas<-gsub(" ","",sigmas)
betas<-gsub(" ","",betas)

params=c("acparam","alpha.1","alpha.2","alpha.3","sigma2.a",
         paste("betabs",c(1:nbasis),".1",sep=""),
         paste("betabs",c(1:nbasis),".2",sep=""),
         paste("betabs",c(1:nbasis),".3",sep=""),
         paste("sigma2.int.",c("start","start","end","end"),".",c("space","unstr"),sep=""),
         sigmas,betas)
rand_vect_cont <- function(N, sd = 1) {
  vec <- rnorm(N/2, 0, sd)
  vec<-c(vec, vec*(-1))
  vec[order(sample(N,N))]
}


bugs.out <- bugs(data=bugsdata, inits=inits, 
                 parameters.to.save=params, model.file=paste("model",fileref,"_",chainnum,".txt",sep=""), 
                 n.chains=nc, n.iter=ni2, n.burnin=nb2, n.thin=nt, debug=F, DIC=TRUE, bugs.seed=chainnum)






if (chainnum==1){
  save.image(file=paste("renvironment",fileref,".RData",sep=""))
} else{
  eval(parse(text=paste("bugs.out",chainnum,"<-bugs.out",sep="")))
  com<-paste( "save(list=c('bugs.out",chainnum,"'), file=paste('bugsout_',fileref,'_',chainnum,'.RData',sep=''))" ,sep="")
  eval(parse(text=com))
}


file.remove(paste("model",fileref,"_",chainnum,".txt",sep=""))
