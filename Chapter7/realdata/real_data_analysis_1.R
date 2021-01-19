#this file was generated from a bash script - 4 were generated for each tau parameter to run parallel chains for the same model.

#these two parameters will change depending on the file reference (determines initial tau) and chain number
fileref<-1
chainnum<-1

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

knotconfig<-23
sds<-c(1,20,10,5,0.5,0.1)
taus<-1/(sds^2)
tau.i<-taus[fileref]

#set up chain parameters
nc <- 1 #number of chains
ni <- 4000 #number of iterations
nb <- 3000 #number for burnin
nt <- 1 #thinning parameter
ni2<-20000 #iterations for model 2
nb2<-19000 #burn in for model 2
boots<-4000 #n. estimates to use for pattern feature CIs
nt2<-1 #thinning parameter model 2

#read journey data (downloaded from TFL open data)
journeys<-read.csv("Nov09JnyExport.csv")

#select relevant information
journeys<-journeys%>%dplyr::select(c(daytype, SubSystem, StartStn, EndStation, EntTimeHHMM, EXTimeHHMM))%>%
  filter(SubSystem=="LUL" & StartStn!="Unstarted" & EndStation!="Unfinished" & StartStn!="Not Applicable" & EndStation!="Not Applicable")

journeys$StartStn<-as.character(journeys$StartStn)
journeys$EndStation<-as.character(journeys$EndStation)
journeys$daytype<-factor(journeys$daytype, levels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
apply(journeys,2,class)
head(journeys)

#read in location information (downloaded from TFL open data)
locs<-st_read("stations.kml")
locs$Name<-as.character(locs$Name)
locs$Name<-gsub("\n\t\t\t","",locs$Name)
locs$Name<-gsub("\t","",locs$Name)
colnames(locs)

#add heathrow terminal 5, london city airport and woodlane
locadd<-locs[1:3,]
colnames(locadd)
locadd$Name<-c("Heathrow Terminal 5", "London City Airport", "Wood Lane")
locadd$description<-NA
locadd$geometry[[1]]<-st_point(x=c(-0.4880,51.4723,0), dim="XYZ")
locadd$geometry[[2]]<-st_point(x=c(-0.0532,51.5032,0), dim="XYZ")
locadd$geometry[[3]]<-st_point(x=c( -0.2212824482,51.5054663115,0), dim="XYZ")
colnames(locs)
colnames(locadd)
locs<-rbind(locs, locadd)
locs<-locs[order(locs$Name),]
rownames(locs)<-c(1:nrow(locs))

#shepherd's bush hammersmith and city should probably be shepherd's bush market as this isn't listed otherwise
locs$geometry[which(locs$Name=="Shepherd's Bush Hammersmith & City Station") ]<-st_point(x=c(-0.222499, 51.503497986, 0), dim="XYZ")

#read in the data describing the tube map connections and zones
netdata<-read.csv("station_names.csv")
netdata$Name<-as.character(netdata$Name)
netdata<-netdata[order(netdata$Name),]
netdata$Neighbours<-as.character(netdata$Neighbours)
colnames(netdata)[1]<-"station_id"

#limit to zone 2 stations and tube stations
locs<-locs[locs$Name%in%netdata$Name[netdata$Zone<=2 & netdata$Neighbours!=""],]
netdata<-netdata[netdata$Zone<=2 & netdata$Neighbours!="",]
locs$id<-row.names(locs)
preservelocs<-locs

#get longitude and latitude as separate variables
long<-sapply(locs$geometry, FUN=function(x){return(x[1])})
lat<-sapply(locs$geometry, FUN=function(x){return(x[2])})
locs<-data.frame(fullname=locs$Name, Name=locs$Name, station_id=locs$id, lat=lat, long=long)


#matching of stations from station locations to journey data
locs$Name<-gsub(" Station", "", locs$Name)
locs$Name[locs$Name=="Shepherd's Bush Hammersmith & City"]<-"shepherd's Bush Market"
locs$Name[locs$Name=="Shepherd's Bush Central"]<-"shepherd's Bush"

locs$Name<-tolower(locs$Name)
locs$Name<-gsub('[[:punct:]]+', '', locs$Name)
locs$Name<-gsub(' ', '', locs$Name)

journeys$startname<-gsub(' ', '', journeys$StartStn)
journeys$startname<-tolower(journeys$startname)
journeys$startname<-gsub('[[:punct:]]+', '', journeys$startname)

journeys$startname[journeys$startname=="canarywharfdlr" | journeys$startname=="canarywharfe2"]<-"canarywharf"
journeys$startname[journeys$startname=="edgwareroadb"]<-"edgwareroadbakerloo"
journeys$startname[journeys$startname=="edgwareroadm"]<-"edgwareroadcircle"
journeys$startname[journeys$startname=="greatportlandst"]<-"greatportlandstreet"
journeys$startname[journeys$startname=="hammersmithd"]<-"hammersmithdp"
journeys$startname[journeys$startname=="hammersmithm"]<-"hammersmith"
journeys$startname[journeys$startname=="heathrowterm4"]<-"heathrowterminal4"
journeys$startname[journeys$startname=="heathrowterm5"]<-"heathrowterminal5"
journeys$startname[journeys$startname=="heathrowterms123"]<-"heathrowterminals123"
journeys$startname[journeys$startname=="highbury"]<-"highburyislington"
journeys$startname[journeys$startname=="highstreetkens"]<-"highstreetkensington"
journeys$startname[journeys$startname=="kingscrossm"]<-"kingscrossstpancras"
journeys$startname[journeys$startname=="kingscrosst"]<-"kingscrossstpancras"
journeys$startname[journeys$startname=="shadwelldlr"]<-"shadwell"
journeys$startname[journeys$startname=="shepherdsbushund"]<-"shepherdsbush"
journeys$startname[journeys$startname=="shepherdsbushmkt"]<-"shepherdsbushmarket"
journeys$startname[journeys$startname=="tottenhamcourtrd"]<-"tottenhamcourtroad"
journeys$startname[journeys$startname=="totteridge"]<-"totteridgewhetstone"
journeys$startname[journeys$startname=="waterloojle"]<-"waterloo"
journeys$startname[journeys$startname=="watfordmet"]<-"watford"
journeys$startname[journeys$startname=="shepherdsbushund"]<-"shepherdsbush"

journeys$endname<-gsub(' ', '', journeys$EndStation)
journeys$endname<-tolower(journeys$endname)
journeys$endname<-gsub('[[:punct:]]+', '', journeys$endname)

journeys$endname[journeys$endname=="canarywharfdlr" | journeys$endname=="canarywharfe2"]<-"canarywharf"
journeys$endname[journeys$endname=="edgwareroadb"]<-"edgwareroadbakerloo"
journeys$endname[journeys$endname=="edgwareroadm"]<-"edgwareroadcircle"
journeys$endname[journeys$endname=="greatportlandst"]<-"greatportlandstreet"
journeys$endname[journeys$endname=="hammersmithd"]<-"hammersmithdp"
journeys$endname[journeys$endname=="hammersmithm"]<-"hammersmith"
journeys$endname[journeys$endname=="heathrowterm4"]<-"heathrowterminal4"
journeys$endname[journeys$endname=="heathrowterm5"]<-"heathrowterminal5"
journeys$endname[journeys$endname=="heathrowterms123"]<-"heathrowterminals123"
journeys$endname[journeys$endname=="highbury"]<-"highburyislington"
journeys$endname[journeys$endname=="highstreetkens"]<-"highstreetkensington"
journeys$endname[journeys$endname=="kingscrossm"]<-"kingscrossstpancras"
journeys$endname[journeys$endname=="kingscrosst"]<-"kingscrossstpancras"
journeys$endname[journeys$endname=="shadwelldlr"]<-"shadwell"
journeys$endname[journeys$endname=="shepherdsbushund"]<-"shepherdsbush"
journeys$endname[journeys$endname=="shepherdsbushmkt"]<-"shepherdsbushmarket"
journeys$endname[journeys$endname=="tottenhamcourtrd"]<-"tottenhamcourtroad"
journeys$endname[journeys$endname=="totteridge"]<-"totteridgewhetstone"
journeys$endname[journeys$endname=="waterloojle"]<-"waterloo"
journeys$endname[journeys$endname=="watfordmet"]<-"watford"
journeys$endname[journeys$endname=="shepherdsbushund"]<-"shepherdsbush"

#join location and journey information
a<-ncol(journeys)
apply(journeys,2,class)
nrow(journeys)
journeys<-inner_join(journeys, locs, by=c("startname"="Name"))
colnames(journeys)[(a+1):ncol(journeys)]<-paste("start.", colnames(journeys)[(a+1):ncol(journeys)], sep="")
a<-ncol(journeys)
nrow(journeys)
journeys$endname<-as.character(journeys$endname)

journeys<-inner_join(journeys, locs, by=c("endname"="Name"))
colnames(journeys)[(a+1):ncol(journeys)]<-paste("end.", colnames(journeys)[(a+1):ncol(journeys)], sep="")

#format start times
journeys$start.time<-times(paste(as.numeric(substr(journeys$EntTimeHHMM, start=1, stop=2))-2, substr(journeys$EntTimeHHMM, start=4, stop=5), '00', sep=':'))
journeys$start.time<-journeys$start.time+times("02:00:00")

journeys$end.time<-times(paste(as.numeric(substr(journeys$EXTimeHHMM, start=1, stop=2))-2, substr(journeys$EXTimeHHMM, start=4, stop=5), '00', sep=':'))
journeys$end.time<-journeys$end.time+times("02:00:00")

#remove journeys with same stard and finish
journeys<-journeys[journeys$endname!=journeys$startname,]

#keep useful columns
journeytimes<-journeys[,c("daytype","start.time","end.time","start.fullname","end.fullname","start.station_id",
                          "end.station_id", "start.lat","start.long","end.lat","end.long")]

#record journey length
journeytimes$length<-journeytimes$end.time-journeytimes$start.time
journeytimes$journeyid<-paste(journeytimes$start.station_id, journeytimes$end.station_id, sep=":")

#convert journey times to hours
journeytimes$start.time<-as.numeric(journeytimes$start.time)*24
journeytimes$end.time<-as.numeric(journeytimes$end.time)*24
journeytimes$length<-as.numeric(journeytimes$length)*24

journeytimes$start.station_id<-as.numeric(as.character(journeytimes$start.station_id))
journeytimes$end.station_id<-as.numeric(as.character(journeytimes$end.station_id))

#select only Wednesday and journeys between 2pm and 10pm
journeytimes<-journeytimes[journeytimes$daytype=="Wed" & journeytimes$start.time>=14 & journeytimes$start.time<=21,]

#renumber the stations to include only the ones used - this was done on the basis it might have been causing some issues with winbugs in early tests
usedstations<-unique(c(journeytimes$start.station_id, journeytimes$end.station_id))
usedstations<-usedstations[order(usedstations)]
usedstations<-data.frame(oldid=usedstations, newid=c(1:length(usedstations)))
#locs data
preservelocs$id<-as.numeric(preservelocs$id)
usedlocs<-left_join(usedstations, preservelocs, by=c("oldid"="id"))
usedlocs<-usedlocs[,2:ncol(usedlocs)]
colnames(usedlocs)[1]<-"station_id"
#network data
netdata<-left_join(usedstations, netdata, by=c("oldid"="station_id"))
netdata<-netdata[,2:ncol(netdata)]
netdata$Neighbours<-as.character(netdata$Neighbours)
netneighbours<-strsplit(netdata$Neighbours, "_")
netneighbours<-lapply(netneighbours, as.numeric)
maxlength<-max(sapply(netneighbours, length))
a<-ncol(netdata)
for (i in 1:maxlength){
  neighbourlist<-data.frame(id=sapply(netneighbours, function(x){return(x[i])}))
  neighbourlist<-left_join(neighbourlist, usedstations, by=c("id"="oldid"))
  netdata<-cbind(netdata, neighbourlist$newid)
  colnames(netdata)[a+i]<-paste("networkconnection.", i, sep="")
}
netdata<-netdata[,-c(a)]
colnames(netdata)[1]<-"station_id"
#journey info
a<-ncol(journeytimes)
journeytimes<-left_join(journeytimes, usedstations, by=c("start.station_id"="oldid"))
colnames(journeytimes)[a+1]<-paste("start.",colnames(journeytimes)[a+1], sep="")
a<-ncol(journeytimes)
journeytimes<-left_join(journeytimes, usedstations, by=c("end.station_id"="oldid"))
colnames(journeytimes)[a+1]<-paste("end.",colnames(journeytimes)[a+1], sep="")
journeytimes$start.station_id<-journeytimes$start.newid
journeytimes$end.station_id<-journeytimes$end.newid


###########################################################
## SET UP DISTANCE MATRICES                              ##
###########################################################

maxdist=2000 #this is in metres

usedlocs<-st_sf(usedlocs, sf_column_name="geometry")
sfdist<-matrix(nrow=nrow(usedlocs), ncol=nrow(usedlocs))

for(i in 1:nrow(usedlocs)){
  for(j in 1:nrow(usedlocs)){
    rad <- pi/180
    lat1<-usedlocs$geometry[[i]][2]
    lat2<-usedlocs$geometry[[j]][2]
    long1<-usedlocs$geometry[[i]][1]
    long2<-usedlocs$geometry[[j]][1]
    a1 <- lat1*rad
    a2 <- long1*rad
    b1 <- lat2*rad
    b2 <- long2*rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1)*cos(b1)*(sin(dlon/2))^2
    c <- 2*atan2(sqrt(a), sqrt(1 - a))
    R <- 6378137
    sfdist[i,j] <- R*c
  }
}
adj1<-apply(sfdist, 2, function(x){which(x<=maxdist & x>0)})
num1<-sapply(adj1, length)
adj1<-unlist(adj1)
sumNumNeigh1<-length(adj1)
#weights
weights1<-apply(sfdist, 2, function(x){x[which(x<=maxdist & x>0)]})
weights1<-unlist(weights1)
weights1[which(weights1!=0)]<-maxdist-weights1[which(weights1!=0)]
#scale weights
weights1[which(weights1!=0)]<-scale(weights1[which(weights1!=0)], center=F, scale=T)
names(adj1)<-NULL
names(weights1)<-NULL
names(num1)<-NULL


#centre time and set up for time lags
centre<-mean(journeytimes$start.time)
journeytimes$start.timec<-journeytimes$start.time-centre
journeytimes<-journeytimes[order(journeytimes$journeyid, journeytimes$start.time),]
journeytimes<-journeytimes%>%group_by(journeyid)%>%
  mutate(timedif=(start.timec-dplyr::lag(start.timec,n=1,default=NA))*60)
journeytimes$timesetup<-0
journeytimes$timesetup[is.na(journeytimes$timedif)==F]<-1
journeytimes$timedif[is.na(journeytimes$timedif)==T]<-0

#choose basis knots
start<-min(journeytimes$start.time)
end<-max(journeytimes$start.time)
splineref<-list()
for(nknots in 1:20){
  splineref[[nknots]]<-seq(from=start, to=end, length.out=nknots+2)[2:(nknots+1)]
}
for(nknots in 1:20){
  i<-nknots+20
  percentiles<-seq(from=0, to=1, length.out=nknots+2)[2:(nknots+1)]
  splineref[[i]]<-quantile(journeytimes$start.time, p=percentiles)
}

knots<-splineref[[knotconfig]]




nbasis<-length(knots)+3

basis<-bSpline(journeytimes$start.timec,  knots=knots-centre, Boundary.knots = c(start,end)-centre)
dspline<-dbs(journeytimes$start.timec, knots=knots-centre, Boundary.knots = c(start,end)-centre)


##################################
#run first model with no temporal
####################################

#list of data
bugsdata<-list(length=journeytimes$length,
               SID=journeytimes$start.station_id, EID=journeytimes$end.station_id,
               adj1=adj1, weights1=weights1, num1=num1,
               Nobs=nrow(journeytimes), Nstart=length(unique(journeytimes$start.station_id)),
               Nend=length(unique(journeytimes$end.station_id)))

a<-length(bugsdata)
for(i in 1:nbasis){
  bugsdata[[i+a]]<-as.numeric(basis[,i])
}
names(bugsdata)[(a+1):(a+nbasis)]<-paste("bs",c(1:nbasis), sep="")

#gen model file
sink(paste("model",fileref,"_",chainnum,".txt",sep=""))
cat("model{
    for(i in 1:Nobs){
    length[i]~dnorm(a[i], tau.a)
    a[i] <- mu[i] + ")

for(i in 1:(nbasis-1)){
  cat(paste("mu",i,"[i] + ",sep=""))
}
cat(paste("mu",nbasis,"[i]\n",sep=""))

cat("
    mu[i] <- alpha.1 + alpha.2[SID[i]] + alpha.3[EID[i]]
    ")




for (i in 1:nbasis){
  for(j in c(".1*bs", ".2[SID[i]]*bs", ".3[EID[i]]*bs")){
    if(j==".1*bs"){
      cat(paste("mu",i,"[i] <- ",sep=""))
    }
    if(j==".3[EID[i]]*bs"){
      cat(paste("betabs", i,j,i,"[i]\n", sep=""))
    }else{
      cat(paste("betabs", i,j,i,"[i] + ", sep=""))
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
    
    tau.a ~ dgamma(1,1)
    sigma2.a <- 1/tau.a # residual error variance
    
    ### residual intercept variance
    tau.int.start.space ~ dgamma(1,1)
    sigma2.int.start.space <- 1/tau.int.start.space #between start station, with spatial cor
    tau.int.start.unstr ~ dgamma(1,1)
    sigma2.int.start.unstr <- 1/tau.int.start.unstr
    
    tau.int.end.space ~ dgamma(1,1)
    sigma2.int.end.space <- 1/tau.int.end.space #between end station, with spatial cor
    tau.int.end.unstr ~ dgamma(1,1)
    sigma2.int.end.unstr <- 1/tau.int.end.unstr
    
    ### variance in basis 1 coefficients
    ")

for(i in 1:nbasis){
  for(j in c("start","end")){
    for (k in c("space","unstr")){
      cat(paste("tau.bs",i,".",j,".",k," ~ dgamma(1,1)\nsigma2.bs",i,".",j,".",k," <- 1/tau.bs",i,".",j,".",k,"\n", sep=""))
    }
  }
}



cat( "alpha.1 ~ dflat() # intercept
     ")
for(i in 1:nbasis){
  cat(paste("betabs",i,".1 ~ dflat()\n", sep=""))
}


cat("}")
sink()


#set up initial values

initscom<-paste("inits<-function() {list(
            alpha.1 = rnorm(1),
            
            tau.a = ",tau.i,",
            
            tau.int.start.space = ",tau.i,",
            tau.int.start.unstr = ",tau.i,",
            tau.int.end.space = ",tau.i,",
            tau.int.end.unstr = ",tau.i,",
            ",sep="")
for(i in 1:nbasis){
  initscom<-c(initscom, paste("betabs",i,".1=rnorm(1),",sep=""))
}
for(j in c("start","end")){
  for(k in c("space","unstr")){
    if(j=="start"){nn<-2}else{nn<-3}
    
    initscom<-c(initscom, paste("alpha.",nn,".",k,"=rep(0,length(unique(journeys$",j,".station_id))),", sep=""))
    
  }
}

for(i in 1:nbasis){
  for(j in c("start","end")){
    if(j=="start"){nn<-2}else{nn<-3}
    for(k in c("space","unstr")){
      
      if(i==nbasis & j=="end" & k=="unstr"){
        initscom<-c(initscom, paste("betabs",i,".",nn,".",k,"=rep(0,length(unique(journeys$",j,".station_id))),", sep=""))
        initscom<-c(initscom,paste("tau.bs",i,".",j,".",k," = ",tau.i,")}", sep=""))
        
      }else{
        initscom<-c(initscom, paste("betabs",i,".",nn,".",k,"=rep(0,length(unique(journeys$",j,".station_id))),", sep=""))
        initscom<-c(initscom,paste("tau.bs",i,".",j,".",k," = ",tau.i,",", sep=""))
      }
    }
  }
}

initscom<-paste(initscom,collapse="")
initscom<-gsub("\n"," ",initscom)

eval(parse(text=paste(initscom)))
p<-proc.time()

#parameters to save
params<-names(inits())
bugs.out <- bugs(data=bugsdata, inits=inits,  
                 parameters.to.save=params, model.file=paste("model",fileref,"_",chainnum,".txt",sep=""), 
                 n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=F, DIC=TRUE)


#############################################################################################
# RUN SECOND ANALYSIS WITH PRIORS & inits BASED ON FIRST MODEL + TEMPORAL COR
################################################################################

#data
bugsdata<-list(length=journeytimes$length,
               SID=journeytimes$start.station_id, EID=journeytimes$end.station_id,
               adj1=adj1, weights1=weights1, num1=num1,
               Nobs=nrow(journeytimes), Nstart=length(unique(journeytimes$start.station_id)),
               Nend=length(unique(journeytimes$end.station_id)), timedif=journeytimes$timedif, timesetup=journeytimes$timesetup)

a<-length(bugsdata)
for(i in 1:nbasis){
  bugsdata[[i+a]]<-as.numeric(basis[,i])
}
names(bugsdata)[(a+1):(a+nbasis)]<-paste("bs",c(1:nbasis), sep="")

#initial values from last model results
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

#model file
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
    
    tau.a ~ dgamma(2,2)
    sigma2.a <- 1/tau.a # residual error variance
    

    ### residual intercept variance
    tau.int.start.space ~ dgamma(2,2)
    sigma2.int.start.space <- 1/tau.int.start.space #between start station, with spatial cor
    tau.int.start.unstr ~ dgamma(2,2)
    sigma2.int.start.unstr <- 1/tau.int.start.unstr
    
    tau.int.end.space ~ dgamma(2,2)
    sigma2.int.end.space <- 1/tau.int.end.space #between end station, with spatial cor
    tau.int.end.unstr ~ dgamma(2,2)
    sigma2.int.end.unstr <- 1/tau.int.end.unstr
    
    ### variance in basis 1 coefficients
    ")

for(i in 1:nbasis){
  for(j in c("start","end")){
    for (k in c("space","unstr")){
      cat(paste("tau.bs",i,".",j,".",k," ~ dgamma(2,2)\nsigma2.bs",i,".",j,".",k," <- 1/tau.bs",i,".",j,".",k,"\n", sep=""))
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

#parameters to keep
params=c("acparam","alpha.1","alpha.2","alpha.3","sigma2.a",
         paste("betabs",c(1:nbasis),".1",sep=""),
         paste("betabs",c(1:nbasis),".2",sep=""),
         paste("betabs",c(1:nbasis),".3",sep=""),
         paste("sigma2.int.",c("start","start","end","end"),".",c("space","unstr"),sep=""),
         sigmas,betas)

#run model
p<-proc.time()
bugs.out <- bugs(data=bugsdata, inits=inits,bugs.seed=chainnum, 
                 parameters.to.save=params, model.file=paste("model",fileref,"_",chainnum,".txt",sep=""), 
                 n.chains=nc, n.iter=ni2, n.burnin=nb2, n.thin=nt2, debug=F, DIC=TRUE)

#save some info about model running
elapsed<-proc.time()-p
iters<-bugs.out$sims.list
iters<-length(iters[[1]])
return<-data.frame(fileref=fileref, time=elapsed[3], DIC=bugs.out$DIC,meandeviance=bugs.out$mean[["deviance"]], pD=bugs.out$pD,savediters=iters)

write.csv(return, paste("test_knots",fileref,"_",chainnum,".csv",sep=""),row.names=F)

# save trace data

tracearray<-bugs.out$sims.array
params<-c(params,"deviance")

tracemat<-as.data.frame(tracemat)

write.csv(tracemat,paste("tracedata",fileref,"_",chainnum,".csv",sep=""), row.names=F)

#save the r environment with the data loaded etc, for future files.
#or save the bugs object to reload
if (chainnum==1){
  save.image(file=paste("renvironment",fileref,".RData",sep=""))
} else{
  eval(parse(text=paste("bugs.out",chainnum,"<-bugs.out",sep="")))
  com<-paste( "save(list=c('bugs.out",chainnum,"'), file=paste('bugsout_',fileref,'_',chainnum,'.RData',sep=''))" ,sep="")
  eval(parse(text=com))
}



