#this file loops over all the results from simulated_analysis_1.R
#files with different file references

library(plyr)
library(tidyr)
library(dplyr)
library(sf)
library(chron)
library(multiplex)
library(igraph)
library(splines2)
library(rootSolve)
library(ggplot2)
library(stringr)
library(Deriv)
library(R2OpenBUGS)
library(ggpubr)
library(geomnet)
library(blandr)

#get coord info
'%ni%'<-Negate('%in%')

#read in station location information
locs<-st_read("stations.kml")
locs$Name<-as.character(locs$Name)
#deal with odd characters in this (only seemed to be a problem on HPC)
locs$Name<-gsub("\n\t\t\t","",locs$Name)
locs$Name<-gsub("\t","",locs$Name)

#add heathrow terminal 5, london city airport and woodlane - these were missing, the long lat coordinates were looked up online
locadd<-locs[1:3,]
locadd$Name<-c("Heathrow Terminal 5", "London City Airport", "Wood Lane")
locadd$description<-NA #I found on a windows 10 PC I had to change "description" in this line to "Description" (capital D)
locadd$geometry[[1]]<-st_point(x=c(-0.4880,51.4723,0), dim="XYZ")
locadd$geometry[[2]]<-st_point(x=c(-0.0532,51.5032,0), dim="XYZ")
locadd$geometry[[3]]<-st_point(x=c( -0.2212824482,51.5054663115,0), dim="XYZ")
locs<-rbind(locs, locadd)
locs<-locs[order(locs$Name),]
rownames(locs)<-c(1:nrow(locs))

#shepherd's bush hammersmith and city should probably be shepherd's bush market as this isn't listed otherwise
locs$geometry[which(locs$Name=="Shepherd's Bush Hammersmith & City Station") ]<-st_point(x=c(-0.222499, 51.503497986, 0), dim="XYZ")

#load file detailing underground network connections between stations and zone - these data were collected from 2009 and 2010 tube maps.
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

#get longitude and latitude as sepearate variables
long<-sapply(locs$geometry, FUN=function(x){return(x[1])})
lat<-sapply(locs$geometry, FUN=function(x){return(x[2])})
locs<-data.frame(fullname=locs$Name, Name=locs$Name, station_id=locs$id, lat=lat, long=long)

#limit to stations in simulation region
locs<-locs[locs$fullname%in%c("King's Cross St. Pancras Station","Euston Station","Euston Square Station","Great Portland Street Station",
                              "Mornington Crescent Station","Camden Town Station","Angel Station","Caledonian Road Station","Warren Street Station",
                              "Tottenham Court Road Station","Russell Square Station","Goodge Street Station"),]
preservelocs<-preservelocs[preservelocs$Name%in%locs$Name,]
locs$station_id<-as.numeric(as.character(locs$station_id))
netdata<-netdata[netdata$station_id%in%locs$station_id,]


#renumber the stations to include only the ones used
usedstations<-unique(locs$station_id)
usedstations<-usedstations[order(usedstations)]
usedstations<-data.frame(oldid=usedstations, newid=c(1:length(usedstations)))
#locs data
preservelocs$id<-as.numeric(as.character(preservelocs$id))
usedlocs<-left_join(usedstations, preservelocs, by=c("oldid"="id"))
usedlocs<-usedlocs[,2:ncol(usedlocs)]
colnames(usedlocs)[1]<-"station_id"
#locs data
locs<-left_join(usedstations, locs, by=c("oldid"="station_id"))
locs<-locs[,2:ncol(locs)]
colnames(locs)[1]<-"station_id"
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
  netdata[which(netdata[,paste("networkconnection.", i, sep="")]%ni%netdata$newid),paste("networkconnection.", i, sep="")]<-NA
}
netdata<-netdata[,-c(a)]
colnames(netdata)[1]<-"station_id"

mapedges<-reshape(netdata, direction="long", varying=colnames(netdata)[4:ncol(netdata)], idvar="station_id")
mapedges<-mapedges[is.na(mapedges$networkconnection)==F,c("station_id","networkconnection")]
colnames(mapedges)<-c("from","to")
a<-ncol(mapedges)+1
mapedges<-merge(mapedges, locs, by.x="from",by.y="station_id")
colnames(mapedges)[a:ncol(mapedges)]<-paste("from.",colnames(mapedges)[a:ncol(mapedges)], sep="")
a<-ncol(mapedges)+1
mapedges<-merge(mapedges, locs, by.x="to", by.y="station_id")
colnames(mapedges)[a:ncol(mapedges)]<-paste("to.",colnames(mapedges)[a:ncol(mapedges)], sep="")

#generate connections

redges<-expand.grid(locs$station_id, locs$station_id)
colnames(redges)<-c("from","to")
redges<-redges[redges$from!=redges$to,]
redges<-redges[order(redges$from,redges$to),]
redges.coord<-merge(redges, locs[,c("station_id","lat","long")], by.x="from",by.y="station_id")
colnames(redges.coord)[3:4]<-c("from.y","from.x")
redges.coord<-merge(redges.coord, locs[,c("station_id","lat","long")], by.x="to",by.y="station_id")
colnames(redges.coord)[5:6]<-c("to.y","to.x")
redges.coord<-redges.coord[,c("from","from.x","from.y","to","to.x","to.y")]

#add fixed effects for from station
numneighbours<-data.frame(station_id=netdata$station_id,num.neighbours=apply(netdata[,4:ncol(netdata)],1,function(x){length(which(is.na(x)==F))})) 

width=119
height=195

png("smalltubemap.png", width=width, height=width*5/7, units="mm", res=600)

print(ggplot()+
  geom_segment(data=mapedges, aes(x=from.long, y=from.lat, xend=to.long, yend=to.lat), colour="grey50")+
  geom_point(data=locs, aes(x=long, y=lat), shape=21, size=3, fill="white")+
  geom_text(data=locs, aes(x=long, y=lat, label=gsub(" Station","",fullname)), size=2, 
            hjust=c(0.3, -0.15, -0.15, -0.15, 1.5, -0.15, 0.1, -0.15, -0.15, -0.15, -0.15, -0.15),
            vjust=c(2, 0, 0, 0.5, 0, 0, 2, 1, 0, 0, 0, 0.6))+
  theme_minimal()+
  coord_quickmap(xlim=c(-0.145,-0.105))+
  theme( axis.title=element_blank(), axis.text=element_blank()))
dev.off()

g<-graph_from_edgelist(as.matrix(mapedges[,c("to","from")]))

tubenetinfo<-data.frame(centrality_closeness=centr_clo(g, mode = "all")$centralization,
                        edge_density=edge_density(g),
                        av_shortest_path=mean_distance(g))
write.csv(tubenetinfo , "tubenetinfo_small.csv", row.names=F)

png("smalltubenet_degreedist.png", width=width, height=width*5/7, units="mm", res=600)

print(ggplot(data=data.frame(degree=degree(g, mode="in")), aes(x=degree))+
  geom_histogram(binwidth=1, fill="white", colour="black")+
  theme_minimal()+
  scale_x_continuous(breaks=c(0:10), labels=c(0:10))+
  xlab("Node degree")+
  ylab("Count")+
  theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
dev.off()

rankcor<-data.frame(id=1:12, maxmin=NA, maxmin.lower=NA, maxmin.upper=NA, maxtime=NA, maxtime.lower=NA, maxtime.upper=NA)
blandaltstats<-data.frame(id=1:12, delaybias=NA, delaybiasupper=NA, delaybiaslower=NA, 
                          delayupperLOA=NA, delayupperLOAupper=NA, delayupperLOAlower=NA, 
                          delaylowerLOA=NA, delaylowerLOAupper=NA, delaylowerLOAlower=NA,
                          timebias=NA, timebiasupper=NA, timebiaslower=NA, 
                          timeupperLOA=NA, timeupperLOAupper=NA, timeupperLOAlower=NA, 
                          timelowerLOA=NA, timelowerLOAupper=NA, timelowerLOAlower=NA)
for(i in 1:12){

  patternfeatures<-read.csv(paste("pattern_feature_estimates",i,".csv",sep=""))
  
  patternfeatures<-merge(patternfeatures, numneighbours, by.x="from",by.y="station_id")
  colnames(patternfeatures)[ncol(patternfeatures)]<-paste(colnames(patternfeatures)[ncol(patternfeatures)],".from",sep="")
  patternfeatures<-merge(patternfeatures, numneighbours, by.x="to",by.y="station_id")
  colnames(patternfeatures)[ncol(patternfeatures)]<-paste(colnames(patternfeatures)[ncol(patternfeatures)],".to",sep="")

  png(paste("degree_v_time_sim",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=num.neighbours.from, y=max.time_real, group=num.neighbours.from))+
    geom_boxplot()+
    xlab("Degree of origin station")+
    ylab("Simulated time of maximum delay")+
    theme_minimal()+
    scale_y_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep="")))+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
  dev.off()
  
  png(paste("degree_v_time_mod",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=num.neighbours.from, y=max.time, group=num.neighbours.from))+
    geom_boxplot()+
    xlab("Degree of origin station")+
    ylab("Modelled time of maximum delay")+
    theme_minimal()+
    scale_y_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep="")))+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
  dev.off()
  
  png(paste("degree_v_delay_sim",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  a<-ggplot(data=patternfeatures,aes(x=num.neighbours.from, y=maxmin_real, group=num.neighbours.from))+
    geom_boxplot()+
    xlab("Degree of origin station")+
    ylab("Simulated maximum delay (minutes)")+
    theme_minimal()+
    scale_y_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
  print(a)
  dev.off()
  
  png(paste("degree_v_delay_mod",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  b<-ggplot(data=patternfeatures,aes(x=num.neighbours.from, y=maxmin, group=num.neighbours.from))+
    geom_boxplot()+
    xlab("Degree of origin station")+
    ylab("Modelled maximum delay (minutes)")+
    theme_minimal()+
    scale_y_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
  print(b)
  
  dev.off()
  
  png(paste("degree_v_delay_all",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  a<-a+ggtitle("A")
  b<-b+ggtitle("B")
  print(ggarrange(a,b,ncol=2,nrow=1))
  
  dev.off()
  
  png(paste("time_v_time",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=max.time_real, y=max.time))+
    geom_point(shape=1)+
    geom_errorbar(aes(ymin=max.time.025, ymax=max.time.975))+
    xlab("Simulated time of maximum delay")+
    ylab("Modelled time of maximum delay")+
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
    scale_y_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
      theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  png(paste("delay_v_delay",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=maxmin_real, y=maxmin))+
    geom_point(shape=1)+
    geom_errorbar(aes(ymin=maxmin.025, ymax=maxmin.974))+
    xlab("Simulated maximum delay (minutes)")+
    ylab("Modelled maximum delay (minutes)")+
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
    scale_y_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
      theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  png(paste("time_v_time_noci",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=max.time_real, y=max.time))+
    geom_point(shape=1)+
    xlab("Simulated time of maximum delay")+
    ylab("Modelled time of maximum delay")+
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
    scale_y_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
      theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  png(paste("delay_v_delay_noci",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=maxmin_real, y=maxmin))+
    geom_point(shape=1)+
    xlab("Simulated maximum delay (minutes)")+
    ylab("Modelled maximum delay (minutes)")+
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
    scale_y_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
      theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  png(paste("difdelay_v_delay",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=maxmin_real, y=maxmin-maxmin_real))+
          geom_point(shape=1)+
          xlab("Simulated maximum delay (minutes)")+
          ylab("Modelled - simulated maximum delay (minutes)")+
          theme_minimal()+
          scale_x_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
          scale_y_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
          theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  png(paste("timedelay_v_time",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=max.time_real, y=max.time-max.time_real))+
          geom_point(shape=1)+
          xlab("Simulated time of maximum delay")+
          ylab("Modelled - simulated time of maximum delay (hours)")+
          theme_minimal()+
          scale_x_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
          scale_y_continuous(breaks=seq(-24,24,by=0.5), labels=seq(-24,24,by=0.5))+
          theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  
  
  png(paste("timedelay_v_dif",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=maxmin_real, y=max.time-max.time_real))+
          geom_point(shape=1)+
          xlab("Simulated maximum delay (minutes)")+
          ylab("Modelled - simulated time of maximum delay (hours)")+
          theme_minimal()+
          scale_x_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
          scale_y_continuous(breaks=seq(-24,24,by=0.5), labels=seq(-24,24,by=0.5))+
          theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  png(paste("difdelay_v_time",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=max.time_real, y=maxmin-maxmin_real))+
          geom_point(shape=1)+
          xlab("Simulated time of maximum delay")+
          ylab("Modelled - simulated maximum delay (minutes)")+
          theme_minimal()+
          scale_x_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
          scale_y_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
          theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  

  
  
  
  png(paste("blandaltman_delay1",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  plot<-blandr.draw(patternfeatures$maxmin, patternfeatures$maxmin_real, method2name="Simulated maximum delay (minutes)",
                    method1name="Modelled maximum delay (minutes)", ciShading=FALSE)
  plot<-plot+ xlab("Mean of modelled and simulated maximum delay (minutes)")+
    ylab("Modelled - simulated maximum delay (minutes)")+
    theme_minimal()+
    scale_x_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
    scale_y_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
    theme(plot.title = element_blank(), axis.title=element_text(size=8), axis.text=element_text(size=8))
  print(plot)
  dev.off()
  
 
  png(paste("blandaltman_time",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  plot<-blandr.draw(patternfeatures$max.time, patternfeatures$max.time_real, method2name="Simulated time of maximum delay (hours)",
                    method1name="Modelled time of maximum delay (hours)", ciShading=FALSE)
  plot<-plot+ xlab("Mean of modelled and simulated time of maximum delay")+
          ylab("Modelled - simulated time of maximum delay (hours)")+
          theme_minimal()+
          scale_x_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
          scale_y_continuous(breaks=seq(-24,24,by=0.5), labels=seq(-24,24,by=0.5))+
          theme(plot.title = element_blank(), axis.title=element_text(size=8), axis.text=element_text(size=8))
  print(plot)
  dev.off()
  
 
  timeba<-blandr.statistics(patternfeatures$max.time, patternfeatures$max.time_real,sig.level=0.95)
  delayba<-blandr.statistics(patternfeatures$maxmin, patternfeatures$maxmin_real,sig.level=0.95)
  
  blandaltstats$delaybias[i]<-delayba$bias 
  blandaltstats$delaybiasupper[i]<-delayba$biasUpperCI 
  blandaltstats$delaybiaslower[i]<- delayba$biasLowerCI
    
  blandaltstats$delayupperLOA[i]<- delayba$upperLOA
    blandaltstats$delayupperLOAupper[i]<- delayba$upperLOA_upperCI
    blandaltstats$delayupperLOAlower[i]<- delayba$upperLOA_lowerCI
    
  blandaltstats$delaylowerLOA[i]<- delayba$lowerLOA
    blandaltstats$delaylowerLOAupper[i]<- delayba$lowerLOA_upperCI
    blandaltstats$delaylowerLOAlower[i]<- delayba$lowerLOA_lowerCI
  
    blandaltstats$timebias[i]<-timeba$bias 
    blandaltstats$timebiasupper[i]<-timeba$biasUpperCI 
    blandaltstats$timebiaslower[i]<- timeba$biasLowerCI
    
    blandaltstats$timeupperLOA[i]<- timeba$upperLOA
    blandaltstats$timeupperLOAupper[i]<- timeba$upperLOA_upperCI
    blandaltstats$timeupperLOAlower[i]<- timeba$upperLOA_lowerCI
    
    blandaltstats$timelowerLOA[i]<- timeba$lowerLOA
    blandaltstats$timelowerLOAupper[i]<- timeba$lowerLOA_upperCI
    blandaltstats$timelowerLOAlower[i]<- timeba$lowerLOA_lowerCI  

  
  rankcor$maxmin[i]<- cor(patternfeatures$maxmin, patternfeatures$maxmin_real, method="spearman")
  rankcor$maxtime[i]<- cor(patternfeatures$max.time, patternfeatures$max.time_real, method="spearman")
  ciwidth<-1.96/sqrt(length(patternfeatures)-3)
  rankcor$maxmin.lower[i]<-tanh(atanh(rankcor$maxmin[i])-ciwidth)
  rankcor$maxmin.upper[i]<-tanh(atanh(rankcor$maxmin[i])+ciwidth)
  rankcor$maxtime.lower[i]<-tanh(atanh(rankcor$maxtime[i])-ciwidth)
  rankcor$maxtime.upper[i]<-tanh(atanh(rankcor$maxtime[i])+ciwidth)
  
  
  balines<-data.frame(type=c("Mean Bias","95% Limit of agreement","95% Limit of agreement", rep("CI",6)),
                      val=c(delayba$bias,delayba$upperLOA,delayba$lowerLOA, delayba$biasUpperCI,delayba$biasLowerCI,delayba$upperLOA_upperCI,delayba$upperLOA_lowerCI,delayba$lowerLOA_upperCI,delayba$lowerLOA_lowerCI ))
  balines$type<-factor(balines$type, levels=c("Mean Bias","95% Limit of agreement", "CI"))
  png(paste("ba_difdelay_v_delay",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot()+
          geom_hline(data=balines[1:3,], aes(yintercept=val, linetype=type))+
          geom_errorbar(aes(ymin=balines$val[c(5,7,9)], ymax=balines$val[c(4,6,8)], x=49/60, width=1/60))+
          geom_point(data=patternfeatures,aes(x=maxmin_real, y=maxmin-maxmin_real),shape=1)+
          xlab("Simulated maximum delay (minutes)")+
          ylab("Modelled - simulated maximum delay (minutes)")+
          theme_minimal()+
          scale_linetype_manual(name="",values=c("solid","dashed"))+
          scale_x_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
          scale_y_continuous(breaks=seq(-0.1,24,by=0.1), labels=c(-6,seq(0,24,by=0.1)*60))+
          theme(legend.position="none",plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  balines<-data.frame(type=c("Mean Bias","95% Limit of agreement","95% Limit of agreement", rep("CI",6)),
                      val=c(timeba$bias,timeba$upperLOA,timeba$lowerLOA, timeba$biasUpperCI,timeba$biasLowerCI,timeba$upperLOA_upperCI,timeba$upperLOA_lowerCI,timeba$lowerLOA_upperCI,timeba$lowerLOA_lowerCI ))
  balines$type<-factor(balines$type, levels=c("Mean Bias","95% Limit of agreement", "CI"))
  
  png(paste("ba_timedelay_v_time",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  
  print(ggplot()+
          geom_hline(data=balines[1:3,], aes(yintercept=val, linetype=type))+
          geom_errorbar(aes(ymin=balines$val[c(5,7,9)], ymax=balines$val[c(4,6,8)], x=18, width=1/10))+
          geom_point(data=patternfeatures,aes(x=max.time_real, y=max.time-max.time_real),shape=1)+
          xlab("Simulated time of maximum delay")+
          ylab("Modelled - simulated time of maximum delay (hours)")+
          theme_minimal()+
          scale_linetype_manual(name="", values=c("solid","dashed"))+
          scale_x_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep=""))+
          scale_y_continuous(breaks=seq(-24,24,by=0.5), labels=seq(-24,24,by=0.5))+
          theme(legend.position="none",plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8)))
  dev.off()
  
  #summary of pattern features
  summary<-data.frame(type=c("Simulated","Modelled","Modelled - Simulated","95% CI Width"), 
                      Mean.maxmin=c(mean(patternfeatures$maxmin_real), mean(patternfeatures$maxmin), mean(patternfeatures$maxmin-patternfeatures$maxmin_real), mean(patternfeatures$maxmin.974-patternfeatures$maxmin.025)), 
                      SD.maxmin=c(sd(patternfeatures$maxmin_real), sd(patternfeatures$maxmin), sd(patternfeatures$maxmin-patternfeatures$maxmin_real), sd(patternfeatures$maxmin.974-patternfeatures$maxmin.025)), 
                      Mean.max.time=c(mean(patternfeatures$max.time_real), mean(patternfeatures$max.time), mean(patternfeatures$max.time-patternfeatures$max.time_real), mean(patternfeatures$max.time.975-patternfeatures$max.time.025)), 
                      SD.max.time=c(sd(patternfeatures$max.time_real), sd(patternfeatures$max.time), sd(patternfeatures$max.time-patternfeatures$max.time_real), sd(patternfeatures$max.time.975-patternfeatures$max.time.025)))
  
  write.csv(summary, paste("summaryfeatures",i,".csv",sep=""), row.names=F)
  
  tracedata<-read.csv(paste("tracedata",i,".csv",sep=""))
  
  fixedtrace<-tracedata[,which(str_count(colnames(tracedata), fixed("."))<=1 & str_count(colnames(tracedata), "sigma")==0 | colnames(tracedata)=="sigma2.a")]
 
   sigmastartunstrtrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")>0 & str_count(colnames(tracedata), "start")>0 | colnames(tracedata)%in%c("iter","chain"))]
  sigmaendunstrtrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")>0 & str_count(colnames(tracedata), "end")>0 | colnames(tracedata)%in%c("iter","chain"))]
  
  sigmastartspacetrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")==0 & str_count(colnames(tracedata), "start")>0 | colnames(tracedata)%in%c("iter","chain"))]
  sigmaendspacetrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")==0 & str_count(colnames(tracedata), "end")>0 | colnames(tracedata)%in%c("iter","chain"))]
  
  plotlist<-list()
  for(j in 1:(ncol(fixedtrace)-2)){
    data<-fixedtrace[,c(j, ncol(fixedtrace)-1, ncol(fixedtrace))]
    colnames(data)[1]<-"ydat"
    plotlist[[j]]<-ggplot(data=data, aes(x=iter, y=ydat, group=chain, colour=as.factor(chain)))+
      geom_line()+
      theme_minimal()+
      xlab("Sample")+
      ylab(colnames(fixedtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_fixed",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  
  plotlist<-list()
  for(j in 1:(ncol(sigmastartunstrtrace)-2)){
    data<-sigmastartunstrtrace[,c(j, ncol(sigmastartunstrtrace)-1, ncol(sigmastartunstrtrace))]
    colnames(data)[1]<-"ydat"
    plotlist[[j]]<-ggplot(data=data, aes(x=iter, y=ydat, group=chain, colour=as.factor(chain)))+
      geom_line()+
      theme_minimal()+
      xlab("Sample")+
      ylab(colnames(sigmastartunstrtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_random_unstr_start",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  
  plotlist<-list()
  for(j in 1:(ncol(sigmaendunstrtrace)-2)){
    data<-sigmaendunstrtrace[,c(j, ncol(sigmaendunstrtrace)-1, ncol(sigmaendunstrtrace))]
    colnames(data)[1]<-"ydat"
    plotlist[[j]]<-ggplot(data=data, aes(x=iter, y=ydat, group=chain, colour=as.factor(chain)))+
      geom_line()+
      theme_minimal()+
      xlab("Sample")+
      ylab(colnames(sigmaendunstrtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_random_unstr_end",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  
  plotlist<-list()
  for(j in 1:(ncol(sigmastartspacetrace)-2)){
    data<-sigmastartspacetrace[,c(j, ncol(sigmastartspacetrace)-1, ncol(sigmastartspacetrace))]
    colnames(data)[1]<-"ydat"
    plotlist[[j]]<-ggplot(data=data, aes(x=iter, y=ydat, group=chain, colour=as.factor(chain)))+
      geom_line()+
      theme_minimal()+
      xlab("Sample")+
      ylab(colnames(sigmastartspacetrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  
  png(paste("trace_random_space_start",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  
  plotlist<-list()
  for(j in 1:(ncol(sigmaendspacetrace)-2)){
    data<-sigmaendspacetrace[,c(j, ncol(sigmaendspacetrace)-1, ncol(sigmaendspacetrace))]
    colnames(data)[1]<-"ydat"
    plotlist[[j]]<-ggplot(data=data, aes(x=iter, y=ydat, group=chain, colour=as.factor(chain)))+
      geom_line()+
      theme_minimal()+
      xlab("Sample")+
      ylab(colnames(sigmaendspacetrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_random_space_end",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  
  }


write.csv(rankcor, "rankcorrelations.csv", row.names=F)
write.csv(blandaltstats, "blandaltman.csv", row.names=F)
