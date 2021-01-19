#this was run on a PC to generate plots, summarise results etc

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

#read in data, etc, as with real data annalysis1
journeys<-read.csv("Nov09JnyExport.csv")

journeys<-journeys%>%dplyr::select(c(daytype, SubSystem, StartStn, EndStation, EntTimeHHMM, EXTimeHHMM))%>%
  filter(SubSystem=="LUL" & StartStn!="Unstarted" & EndStation!="Unfinished" & StartStn!="Not Applicable" & EndStation!="Not Applicable")

journeys$StartStn<-as.character(journeys$StartStn)
journeys$EndStation<-as.character(journeys$EndStation)
journeys$daytype<-factor(journeys$daytype, levels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
apply(journeys,2,class)
head(journeys)
locs<-st_read("stations.kml")
locs$Name<-as.character(locs$Name)
locs$Name<-gsub("\n\t\t\t","",locs$Name)
locs$Name<-gsub("\t","",locs$Name)
colnames(locs)

#add heathrow term 5, london city airport and woodlane
locadd<-locs[1:3,]
colnames(locadd)
locadd$Name<-c("Heathrow Terminal 5", "London City Airport", "Wood Lane")
locadd$Description<-NA
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

#get longitude and latitude as seperate variables
long<-sapply(locs$geometry, FUN=function(x){return(x[1])})
lat<-sapply(locs$geometry, FUN=function(x){return(x[2])})
locs<-data.frame(fullname=locs$Name, Name=locs$Name, station_id=locs$id, lat=lat, long=long)



#matching stations
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

#join with journey data
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
head(journeys)
journeys$start.time<-times(paste(as.numeric(substr(journeys$EntTimeHHMM, start=1, stop=2))-2, substr(journeys$EntTimeHHMM, start=4, stop=5), '00', sep=':'))
journeys$start.time<-journeys$start.time+times("02:00:00")

journeys$end.time<-times(paste(as.numeric(substr(journeys$EXTimeHHMM, start=1, stop=2))-2, substr(journeys$EXTimeHHMM, start=4, stop=5), '00', sep=':'))
journeys$end.time<-journeys$end.time+times("02:00:00")

journeys<-journeys[journeys$endname!=journeys$startname,]


journeytimes<-journeys[,c("daytype","start.time","end.time","start.fullname","end.fullname","start.station_id",
                          "end.station_id", "start.lat","start.long","end.lat","end.long")]

journeytimes$length<-journeytimes$end.time-journeytimes$start.time
journeytimes$journeyid<-paste(journeytimes$start.station_id, journeytimes$end.station_id, sep=":")

journeytimes$start.time<-as.numeric(journeytimes$start.time)*24
journeytimes$end.time<-as.numeric(journeytimes$end.time)*24
journeytimes$length<-as.numeric(journeytimes$length)*24

journeytimes$start.station_id<-as.numeric(as.character(journeytimes$start.station_id))
journeytimes$end.station_id<-as.numeric(as.character(journeytimes$end.station_id))

#select only wednesday and journeys between 2pm and 10pm
journeytimes<-journeytimes[journeytimes$daytype=="Wed" & journeytimes$start.time>=14 & journeytimes$start.time<=21,]

#renumber the stations to include only the ones used
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

width=119
height=195

#generate graph of tube map to plot
long<-sapply(usedlocs$geometry, FUN=function(x){return(x[1])})
lat<-sapply(usedlocs$geometry, FUN=function(x){return(x[2])})
locs<-data.frame(fullname=usedlocs$Name, Name=usedlocs$Name, station_id=usedlocs$station_id, lat=lat, long=long)
mapedges<-reshape(netdata, direction="long", varying=colnames(netdata)[4:ncol(netdata)], idvar="station_id")
mapedges<-mapedges[is.na(mapedges$networkconnection)==F,c("station_id","networkconnection")]
colnames(mapedges)<-c("from","to")
a<-ncol(mapedges)+1
mapedges<-merge(mapedges, locs, by.x="from",by.y="station_id")
colnames(mapedges)[a:ncol(mapedges)]<-paste("from.",colnames(mapedges)[a:ncol(mapedges)], sep="")
a<-ncol(mapedges)+1
mapedges<-merge(mapedges, locs, by.x="to", by.y="station_id")
colnames(mapedges)[a:ncol(mapedges)]<-paste("to.",colnames(mapedges)[a:ncol(mapedges)], sep="")
g<-graph_from_edgelist(as.matrix(mapedges[,c("to","from")]))

#information about centrality based on tube map
numneighbours<-data.frame(station_id=netdata$station_id,
                          centrality=centr_clo(g, mode="all")$res,
                          num.neighbours=apply(netdata[,4:ncol(netdata)],1,function(x){length(which(is.na(x)==F))})) 


#to loop over other files if needed
for(i in 1){
  #read in pattern feature information
  patternfeatures<-read.csv(paste("pattern_feature_estimates",i,".csv",sep=""))
  patternfeatures<-merge(patternfeatures, numneighbours, by.x="from",by.y="station_id")
  colnames(patternfeatures)[(ncol(patternfeatures)-1):ncol(patternfeatures)]<-paste(colnames(patternfeatures)[(ncol(patternfeatures)-1):ncol(patternfeatures)],".from",sep="")
  patternfeatures<-merge(patternfeatures, numneighbours, by.x="to",by.y="station_id")
  colnames(patternfeatures)[(ncol(patternfeatures)-1):ncol(patternfeatures)]<-paste(colnames(patternfeatures)[(ncol(patternfeatures)-1):ncol(patternfeatures)],".to",sep="")
  
  
  #plots with centrality measures
  png(paste("degree_v_time_mod",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=num.neighbours.from, y=max.time, group=num.neighbours.from))+
          geom_boxplot()+
          xlab("Degree of origin station")+
          ylab("Modelled time of maximum delay")+
          theme_minimal()+
          scale_y_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep="")))+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
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
  
  png(paste("centr_v_time_mod",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  print(ggplot(data=patternfeatures,aes(x=centrality.from, y=max.time))+
          geom_point(shape=21)+
          xlab("Centrality (closeness) of origin station")+
          ylab("Modelled time of maximum delay")+
          theme_minimal()+
          scale_y_continuous(breaks=seq(0,24,by=2), labels=paste(seq(0,24,by=2),":00",sep="")))+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
  dev.off()
  
  
  
  png(paste("centr_v_delay_mod",i, ".png", sep="_"), width=width, height=width*5/7, units="mm", res=600)
  b<-ggplot(data=patternfeatures,aes(x=centrality.from, y=maxmin))+
    geom_point(shape=21)+
    xlab("Centrality (closeness) of origin station")+
    ylab("Modelled maximum delay (minutes)")+
    theme_minimal()+
    scale_y_continuous(breaks=seq(0,24,by=0.1), labels=seq(0,24,by=0.1)*60)+
    theme(plot.title = element_text(size=8), axis.title=element_text(size=8), axis.text=element_text(size=8))
  print(b)
  
  dev.off()
  
  #read in trace data
  tracedata<-read.csv(paste("tracedata",i,"_1.csv",sep=""))
  
  for(j in 2:4){
    num<-nrow(tracedata)
    tracedata<-rbind(tracedata,read.csv(paste("tracedata",i,"_",j,".csv",sep="")))
    tracedata$chain[num:nrow(tracedata)]<-j
  }
  
  #plots of trace data
  fixedtrace<-tracedata[,which(str_count(colnames(tracedata), fixed("."))<=1 & str_count(colnames(tracedata), "sigma")==0 | colnames(tracedata)=="sigma2.a")]
  
  sigmastartunstrtrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")>0 & str_count(colnames(tracedata), "start")>0 | colnames(tracedata)%in%c("iter","chain"))]
  sigmaendunstrtrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")>0 & str_count(colnames(tracedata), "end")>0 | colnames(tracedata)%in%c("iter","chain"))]
  
  sigmastartspacetrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")==0 & str_count(colnames(tracedata), "start")>0 | colnames(tracedata)%in%c("iter","chain"))]
  sigmaendspacetrace<-tracedata[ which(str_count(colnames(tracedata), "sigma")>0 & str_count(colnames(tracedata), "unstr")==0 & str_count(colnames(tracedata), "end")>0 | colnames(tracedata)%in%c("iter","chain"))]
  
  plotlist<-list()
  plotlist2<-list()
  for(j in 1:(ncol(fixedtrace)-2)){
    data<-fixedtrace[,c(j, ncol(fixedtrace)-1, ncol(fixedtrace))]
    colnames(data)[1]<-"ydat"
    plotlist[[j]]<-ggplot(data=data, aes(x=iter, y=ydat, group=chain, colour=as.factor(chain)))+
      geom_line()+
      theme_minimal()+
      xlab("Sample")+
      ylab(colnames(fixedtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
    plotlist2[[j]]<-ggplot(data=data, aes(x=ydat))+
      geom_density()+
      theme_minimal()+
      xlab(colnames(fixedtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
    
  }
  png(paste("trace_fixed",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  png(paste("density_fixed",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist2))
  dev.off()
  
  plotlist2<-list()
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
    plotlist2[[j]]<-ggplot(data=data, aes(x=ydat))+
      geom_density()+
      theme_minimal()+
      xlab(colnames(sigmastartunstrtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_random_unstr_start",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  png(paste("density_random_unstr_start",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist2))
  dev.off()
  
  plotlist2<-list()
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
    plotlist2[[j]]<-ggplot(data=data, aes(x=ydat))+
      geom_density()+
      theme_minimal()+
      xlab(colnames(sigmaendunstrtrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_random_unstr_end",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  png(paste("density_random_unstr_end",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist2))
  dev.off()
  
  plotlist2<-list()
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
    plotlist2[[j]]<-ggplot(data=data, aes(x=ydat))+
      geom_density()+
      theme_minimal()+
      xlab(colnames(sigmastartspacetrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  
  png(paste("trace_random_space_start",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  png(paste("density_random_space_start",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist2))
  dev.off()
  
  plotlist2<-list()
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
    plotlist2[[j]]<-ggplot(data=data, aes(x=ydat))+
      geom_density()+
      theme_minimal()+
      xlab(colnames(sigmaendspacetrace)[j])+
      theme(legend.position="none", axis.title=element_text(size=8), axis.text=element_text(size=8))
  }
  png(paste("trace_random_space_end",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist))
  dev.off()
  png(paste("density_random_space_end",i, ".png", sep="_"), width=3*width, height=3*width*5/7, units="mm", res=600)
  
  print(ggarrange(plotlist=plotlist2))
  dev.off()
  
  plotlist2<-list()
  
  #summary of pattern features
  summary<-data.frame(type=c("Modelled","95% CI Width"), 
                      Mean.maxmin=c( mean(patternfeatures$maxmin, na.rm=T),  mean(patternfeatures$maxmin.974-patternfeatures$maxmin.025, na.rm=T)), 
                      SD.maxmin=c( sd(patternfeatures$maxmin, na.rm=T), sd(patternfeatures$maxmin.974-patternfeatures$maxmin.025, na.rm=T)), 
                      Mean.max.time=c( mean(patternfeatures$max.time, na.rm=T), mean(patternfeatures$max.time.975-patternfeatures$max.time.025, na.rm=T)), 
                      SD.max.time=c( sd(patternfeatures$max.time, na.rm=T),  sd(patternfeatures$max.time.975-patternfeatures$max.time.025, na.rm=T)))
  
  write.csv(summary, paste("summaryfeatures",i,".csv",sep=""), row.names=F)
  
  
}
