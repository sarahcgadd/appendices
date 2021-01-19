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




knotconfig<-19

tau.i<-1
#CHAIN PARAMETERS
#temporary chain parameters for testing if models run
nc <- 1 #number of chains
ni <- 4000 #number of iterations
nb <- 3000 #number for burnin
nt <- 1 #thinning parameter
ncores<-1 #n cores to use
ni2<-20000
nb2<-19000
boots<-4000
#read data

'%ni%'<-Negate('%in%')
simulate_trajectory <- function(n.obs=1 , max.time=20 , record.times=NULL, record.times.vary=1,
                                part.1.average=character() , part.2.average=character() , 
                                part.1.slope=c("rnorm", "0", "1", "1") , part.2.slope=c("rnorm", "0", "1", "1") , 
                                time.step.change=c("rnorm", "0", "1") , value.step.change=character() , 
                                sd.error=c("rlnorm","0.5", "0.5"), between.cor=0, within.cor=0,
                                dayreps=1, dayreps.cor=0, no.day.var=F, fixed.spacegroups=NULL,
                                fixed.effect.part.1.slope=NULL,
                                spacegroups=NULL, spacegroup.size=NULL, 
                                origin.dest=FALSE, area.x=NULL, area.y=NULL, sd.ratio=0.5, 
                                sd.ratio.day=0.5, plot=TRUE, offset.amount=0){ 
  #load libraries required
  library(ggplot2)
  library(gridExtra)
  library(Matrix)
  library(matrixcalc)
  library(MASS)
  
  #load covariance function
  Covar <- function(n=2,SD=data.frame(1,1),Cor=data.frame(0.5)) {
    # cov(i,j) = cor(i,j)*sd(i)*sd(j)
    Cov <- matrix(nrow=n,ncol=n)
    for (i in 1:n) { for (j in 1:n) { Cov[i,j] <- Cor[i,j]*SD[i]*SD[j] }} 
    Cov <- as.matrix(forceSymmetric(Cov))
    if (!is.positive.definite(Cov)) {
      print("Warning: covariance matrix made Positive Definite")
      Cov <- as.matrix(nearPD(Cov)$mat) }
    return(Cov) 
  }
  
  if (length(part.1.slope)==3){
    part.1.slope[4]<-"1"
  }  
  if (length(part.2.slope)==3){
    part.2.slope[4]<-"1"
  } 
  
  #warning messages for if incorrect number of things are specified -- only two of part.2.average, part.1.average and value.step.change can be specified
  if (length(part.2.average)>0 && length(part.1.average)>0 && length(value.step.change)>0) {
    print("Warning: Only 3 distributions of part.1.average, part.2.average and value.step.change can be used. part.2.average input ignored.")
  }
  if(length(part.2.average)==0 && (length(part.1.average)==0 | length(value.step.change)==0) | (length(part.1.average)==0 && length(value.step.change)==0)){
    stop("Too few distributions specified. Two distributions out of part.1.average, part.2.average and value.step.change must be specified.")
  }
  if (abs(between.cor)>1){
    stop("between.cor incorrectly specified. between.cor must lie between -1 and 1")
  }
  if (abs(within.cor)>1){
    stop("within.cor incorrectly specified. within.cor must lie between -1 and 1")
  }
  if(is.null(fixed.spacegroups)==F){
    if(class(fixed.spacegroups)!="matrix" & class(fixed.spacegroups)!="data.frame"){
      stop("If specified, fixed.spacegroups should be a *matrix or data frame* with 3 or 6 columns (group1, x coordinate 1, y coordinate 1, group2, x coordinate 2, y coordinate 2) and n.obs rows.")
    }
    if(ncol(fixed.spacegroups)!=3 & ncol(fixed.spacegroups)!=6){
      stop("If specified, fixed.spacegroups should be a matrix or data frame with *3 or 6 columns* (group1, x coordinate 1, y coordinate 1, group2, x coordinate 2, y coordinate 2) and n.obs rows.")
    }
    if(nrow(fixed.spacegroups)!=n.obs){
      stop("If specified, fixed.spacegroups should be a matrix or data frame with 3 or 6 columns (group1, x coordinate 1, y coordinate 1, group2, x coordinate 2, y coordinate 2) and *n.obs rows*.")
    }
    if(ncol(fixed.spacegroups)==3){
      origin.dest=FALSE
    }
  }
  #need some error messages for things in the spacegroups2 dept
  #set up some values and things
  maxt<-max.time 
  obs<-n.obs 
  if (dayreps>1){
    obs<-n.obs*dayreps
  }
  slope.1.mult<-part.1.slope[4]
  slope.2.mult<-part.2.slope[4]
  
  
  #generate coordinate matrix depending on the number of groups
  if(is.null(fixed.spacegroups)==T){
    if (length(spacegroups)==1){
      spacegroup1<-round(runif(min=0.51, max=spacegroups+0.49, n=obs/dayreps))
      x.coord1<-numeric(length=obs/dayreps)
      y.coord1<-numeric(length=obs/dayreps)
      spacegroup1.mean.x<-runif(spacegroups, 0, area.x)
      spacegroup1.mean.y<-runif(spacegroups, 0, area.y)
      for(f in 1:spacegroups){
        x.coord1[which(spacegroup1==f)]<-rnorm(length(spacegroup1[which(spacegroup1==f)]), spacegroup1.mean.x[f], spacegroup.size)
        y.coord1[which(spacegroup1==f)]<-rnorm(length(spacegroup1[which(spacegroup1==f)]), spacegroup1.mean.y[f], spacegroup.size)
      }
      x.coord1<-rep(x.coord1,dayreps)
      y.coord1<-rep(y.coord1,dayreps)
      
      
      #matrix of relative distances 1
      distance1<-matrix(ncol=obs, nrow=obs)
      for (id1 in 1:obs){
        for(id2 in 1:obs){
          distance1[id1,id2]<-((x.coord1[id1]-x.coord1[id2])^2+(y.coord1[id1]-y.coord1[id2])^2)^0.5
        }
      }
      spacegroup.info<-cbind(group.no=c(1:spacegroups), spacegroup1.mean.x, spacegroup1.mean.y, coord.sd=rep(spacegroup.size, spacegroups))
      individual.coords1<-cbind(id=c(1:n.obs), spacegroup1=spacegroup1[1:n.obs],x.coord1=x.coord1[1:n.obs],y.coord1=y.coord1[1:n.obs])
      groupings<-unique(spacegroup1)
      groupref<-rep(spacegroup1, dayreps)
    } else {
      distance1<-matrix(ncol=obs, nrow=obs, rep(1, obs^2))
    }
    
    #repeat for destination coords
    if (origin.dest==TRUE){
      sample<-sample(c(1:(obs/dayreps)), obs/dayreps, replace=T)
      x.coord2<-x.coord1[sample]
      y.coord2<-y.coord1[sample]
      spacegroup2<-spacegroup1[sample]
      x.coord2<-rep(x.coord2,dayreps)
      y.coord2<-rep(y.coord2,dayreps)
      
      #matrix of relative distances 2
      distance2<-matrix(ncol=obs, nrow=obs)
      for (id1 in 1:obs){
        for(id2 in 1:obs){
          distance2[id1,id2]<-((x.coord2[id1]-x.coord2[id2])^2+(y.coord2[id1]-y.coord2[id2])^2)^0.5
        }
      }
      individual.coords2<-cbind(id=c(1:n.obs),spacegroup2=spacegroup2[1:n.obs],x.coord2=x.coord2[1:n.obs],y.coord2=y.coord2[1:n.obs])
      routegroup<-as.numeric(paste(spacegroup1, spacegroup2,sep=""))
      groupings<-unique(routegroup)
      groupref<-rep(routegroup, dayreps)
      odinfo<-rep(paste(spacegroup1, spacegroup2,sep=":"),dayreps)
      
    } else {
      distance2<-matrix(ncol=obs, nrow=obs, rep(1, obs^2))
    }
  } 
  if(is.null(fixed.spacegroups)==F){
    
    spacegroup1<-as.numeric(as.character(fixed.spacegroups[,1]))
    x.coord1<-as.numeric(fixed.spacegroups[,2])
    y.coord1<-as.numeric(fixed.spacegroups[,3])
    x.coord1<-rep(x.coord1,dayreps)
    y.coord1<-rep(y.coord1,dayreps)
    
    
    distance1<-matrix(ncol=obs, nrow=obs)
    for (id1 in 1:obs){
      for(id2 in 1:obs){
        distance1[id1,id2]<-((x.coord1[id1]-x.coord1[id2])^2+(y.coord1[id1]-y.coord1[id2])^2)^0.5
      }
    }
    spacegroup1.mean.x<-numeric()
    spacegroup1.mean.y<-numeric()
    spacegroup1.sd.x<-numeric()
    spacegroup1.sd.y<-numeric()
    for(f in unique(spacegroup1)){
      spacegroup1.mean.x<-c(spacegroup1.mean.x, mean(x.coord1[spacegroup1==f]))
      spacegroup2.mean.y<-c(spacegroup1.mean.y, mean(y.coord1[spacegroup1==f]))
      spacegroup1.sd.x<-c(spacegroup1.sd.x, sd(x.coord1[spacegroup1==f]))
      spacegroup2.sd.y<-c(spacegroup1.sd.y, sd(y.coord1[spacegroup1==f]))
    }
    spacegroup.info<-cbind(group.no=unique(spacegroup1), spacegroup1.mean.x, spacegroup1.mean.y, spacegroup1.sd.x, spacegroup1.sd.y)
    individual.coords1<-cbind(id=c(1:n.obs), spacegroup1=spacegroup1[1:n.obs],x.coord1=x.coord1[1:n.obs],y.coord1=y.coord1[1:n.obs])
    groupings<-unique(spacegroup1)
    groupref<-rep(spacegroup1, dayreps)
    
    if(ncol(fixed.spacegroups)==6){
      spacegroup2<-as.numeric(as.character(fixed.spacegroups[,4]))
      x.coord2<-as.numeric(fixed.spacegroups[,5])
      y.coord2<-as.numeric(fixed.spacegroups[,6])
      x.coord2<-rep(x.coord2,dayreps)
      y.coord2<-rep(y.coord2,dayreps)
      distance2<-matrix(ncol=obs, nrow=obs)
      for (id1 in 1:obs){
        for(id2 in 1:obs){
          distance2[id1,id2]<-((x.coord2[id1]-x.coord2[id2])^2+(y.coord2[id1]-y.coord2[id2])^2)^0.5
        }
      }
      individual.coords2<-cbind(id=c(1:n.obs),spacegroup2=spacegroup2[1:n.obs],x.coord2=x.coord2[1:n.obs],y.coord2=y.coord2[1:n.obs])
      if(all.equal(spacegroup1, spacegroup2)!=TRUE){
        routegroup<-as.numeric(paste(spacegroup1, spacegroup2,sep=""))
        groupref<-rep(routegroup, dayreps)
        groupings<-unique(routegroup)
        odinfo<-rep(paste(spacegroup1, spacegroup2,sep=":"),dayreps)
      } else {
        routegroup<-spacegroup1
      }
      
      
    } else {
      distance2<-matrix(ncol=obs, nrow=obs, rep(1, obs^2))
    }
  }
  #matrix of day differences
  daynumbers<-rep(c(1:dayreps), n.obs)
  daynumbers<-daynumbers[order(daynumbers)]
  days<-matrix(ncol=obs,nrow=obs)
  for (id1 in 1:obs){
    for(id2 in 1:obs){
      days[id1,id2]<-abs(daynumbers[id1]-daynumbers[id2])
    }
  }
  
  #simulate add on for day reps (mean 0)
  if (dayreps>1){
    slope1.day<-rnorm(dayreps, 0, as.numeric(part.1.slope[3])*sd.ratio.day)
    slope2.day<-rnorm(dayreps, 0, as.numeric(part.2.slope[3])*sd.ratio.day)
    jumptime.day<-rnorm(dayreps, 0, as.numeric(time.step.change[3])*sd.ratio.day)
    if(length(value.step.change)>0){
      jumpval.day<-rnorm(dayreps, 0, as.numeric(value.step.change[3])*sd.ratio.day)
    }
    if(length(part.1.average)>0){
      average1.day<-rnorm(dayreps, 0, as.numeric(part.1.average[3])*sd.ratio.day)
    }
    if(length(part.2.average)>0){
      average2.day<-rnorm(dayreps, 0, as.numeric(part.2.average[3])*sd.ratio.day)
    }
  } else{
    slope1.day<-0
    slope2.day<-0
    jumptime.day<-0
    jumpval.day<-0
    average1.day<-0
    average2.day<-0
  }
  
  daygroups<-rep(c(1:dayreps), length(groupings))
  daygroups2<-daygroups[order(daygroups)]
  spacedaygroups<-cbind(spatial.groups=rep(groupings, dayreps), daygroups2)
  groupref<-as.numeric(paste(groupref, daynumbers, sep=""))
  
  
  #then add this on for each thing
  #set group level mean parameter values 
  if(exists("groupings")==TRUE){
    eval(parse(text=paste("slope1.groups.num<-" , part.1.slope[1] , "(", length(groupings), "," , part.1.slope[2] , "," , part.1.slope[3] , ")" , sep="")))
    slope1.groups.num<-rep(slope1.groups.num, dayreps)+rep(slope1.day,length(groupings))[order(daygroups)]
    slope1.groups<-paste(slope1.groups.num, "*", slope.1.mult, "*t", sep="")
    eval(parse(text=paste("slope2.groups.num<-" , part.2.slope[1] , "(", length(groupings), "," , part.2.slope[2] , "," , part.2.slope[3] , ")" , sep="")))
    slope2.groups.num<-rep(slope2.groups.num, dayreps)+rep(slope2.day,length(groupings))[order(daygroups)]
    slope2.groups<-paste(slope2.groups.num, "*", slope.2.mult, "*t", sep="") 
    slope1.funct.groups<-list(length=length(slope1.groups)*dayreps)
    slope2.funct.groups<-list(length=length(slope2.groups)*dayreps)
    for (i in 1:length(slope2.groups)){
      eval(parse(text=paste("slope2.funct.groups[[i]]<-as.function(alist(t=, ", slope2.groups[i],"))", sep="")))
      eval(parse(text=paste("slope1.funct.groups[[i]]<-as.function(alist(t=, ", slope1.groups[i],"))", sep="")))
    }
    
    eval(parse(text=paste("jumptime.groups<-" , time.step.change[1] , "(", length(groupings), "," , time.step.change[2] , "," , time.step.change[3] , ")" , sep="")))
    jumptime.groups<-rep(jumptime.groups, dayreps)+rep(jumptime.day,length(groupings))[order(daygroups)]
    
    if (length(part.2.average)==0 | (length(part.2.average)>0 && length(part.1.average)>0 && length(value.step.change)>0)) {
      eval(parse(text=paste("jumpval.groups<-" , value.step.change[1] , "(", length(groupings), "," , value.step.change[2] , "," , value.step.change[3] , ")" , sep="")))
      jumpval.groups<-rep(jumpval.groups, dayreps)+rep(jumpval.day,length(groupings))[order(daygroups)]
      eval(parse(text=paste("average1.groups<-" , part.1.average[1] , "(", length(groupings), "," , part.1.average[2] , "," , part.1.average[3] , ")" , sep=""))) 
      average1.groups<-rep(average1.groups, dayreps)+rep(average1.day,length(groupings))[order(daygroups)]
      intercept.2.groups<-numeric(length=length(groupings)*dayreps)
      true.start.groups<-numeric(length=length(groupings)*dayreps)
      average2.groups<-numeric(length=length(groupings)*dayreps)
      for (i in 1:(length(groupings)*dayreps)){
        true.start.groups[i]<-(average1.groups[i]*jumptime.groups[i]-integrate(slope1.funct.groups[[i]], lower=0, upper=jumptime.groups[i] )$value)/jumptime.groups[i]
        intercept.2.groups[i]<-slope1.funct.groups[[i]](jumptime.groups[i])+true.start.groups[i]+jumpval.groups[i]-slope2.funct.groups[[i]](jumptime.groups[i])
        average2.groups[i]<-(integrate(slope2.funct.groups[[i]], upper=maxt, lower=jumptime.groups[i] )$value+intercept.2.groups[i]*(maxt-jumptime.groups[i]))/(maxt-jumptime.groups[i])
      }
      
    }
    if (length(part.1.average)==0) {
      eval(parse(text=paste("jumpval.groups<-" , value.step.change[1] , "(", length(groupings), "," , value.step.change[2] , "," , value.step.change[3] , ")" , sep="")))
      jumpval.groups<-rep(jumpval.groups, dayreps)+rep(jumpval.day,length(groupings))[order(daygroups)]
      eval(parse(text=paste("average2.groups<-" , part.2.average[1] , "(", length(groupings), "," , part.2.average[2] , "," , part.2.average[3] , ")" , sep=""))) 
      average2.groups<-rep(average2.groups, dayreps)+rep(average2.day,length(groupings))[order(daygroups)]
      intercept.2.groups<-numeric(length=length(groupings)*dayreps)
      true.start.groups<-numeric(length=length(groupings)*dayreps)
      average1.groups<-numeric(length=length(groupings)*dayreps)
      for (i in 1:(length(groupings)*dayreps)){
        intercept.2.groups[i]<-(average2.groups[i]*(maxt-jumptime.groups[i])-integrate(slope2.funct.groups[[i]], lower=jumptime.groups[i], upper=maxt )$value)/(maxt-jumptime.groups[i])
        true.start.groups[i]<-slope2.funct.groups[[i]](jumptime.groups[i])+intercept.2.groups[i]-jumpval.groups[i]-slope1.funct.groups[[i]](jumptime.groups[i])
        average1.groups[i]<-(integrate(slope1.funct.groups[[i]], upper=jumptime.groups[i], lower=0 )$value+true.start.groups[i]*jumptime.groups[i])/jumptime.groups[i]
      }
    }
    if (length(value.step.change)==0){
      eval(parse(text=paste("average1.groups<-" , part.1.average[1] , "(", length(groupings), "," , part.1.average[2] , "," , part.1.average[3] , ")" , sep=""))) 
      average1.groups<-rep(average1.groups, dayreps)+rep(average1.day,length(groupings))[order(daygroups)]
      eval(parse(text=paste("average2.groups<-" , part.2.average[1] , "(", length(groupings), "," , part.2.average[2] , "," , part.2.average[3] , ")" , sep=""))) 
      average2.groups<-rep(average2.groups, dayreps)+rep(average2.day,length(groupings))[order(daygroups)]
      intercept.2.groups<-numeric(length=length(groupings)*dayreps)
      true.start.groups<-numeric(length=length(groupings)*dayreps)
      jumpval.groups<-numeric(length=length(groupings)*dayreps)
      for (i in 1:(length(groupings)*dayreps)){
        true.start.groups[i]<-(average1.groups[i]*jumptime.groups[i]-integrate(slope1.funct.groups[[i]], lower=0, upper=jumptime.groups[i] )$value)/jumptime.groups[i]
        intercept.2.groups[i]<-(average2.groups[i]*(maxt-jumptime.groups[i])-integrate(slope2.funct.groups[[i]], lower=jumptime.groups[i], upper=maxt )$value)/(maxt-jumptime.groups[i])
        jumpval.groups[i]<-slope2.funct.groups[[i]](jumptime.groups[i])+intercept.2.groups[i]-slope1.funct.groups[[i]](jumptime.groups[i])-true.start.groups[i]
      }
    }
    #record group values
    group.info<-data.frame(cbind(spacedaygroups, average1.groups, average2.groups, slope1.groups, slope2.groups, jumptime.groups, jumpval.groups))
    
    groupings<-as.numeric(paste(spacedaygroups[,1], spacedaygroups[,2], sep=""))
    if(offset.amount>0){
      if(length(groupings)>1){
        offset.sd<-rnorm(length(groupings), 0, offset.amount)
        group.info$slope1.groups<-mapply(gsub, replacement=paste("(t+",offset.sd,")",sep=""), x=group.info$slope1.groups, MoreArgs=list(pattern="t"))
        group.info$slope2.groups<-mapply(gsub, replacement=paste("(t+",offset.sd,")",sep=""), x=group.info$slope2.groups, MoreArgs=list(pattern="t"))
        group.info$offset.sd<-offset.sd
      } else{
        offset.sd<-0
      }
    } else {
      if(length(groupings)>1){
        offset.sd<-rep(0, length(groupings))
        group.info$offset.sd<-offset.sd
      } else {
        offset.sd<-0
      }
      
    }
    #set true parameter values for individuals
    #set up vectors
    slope1<-character(length=obs)
    slope2<-character(length=obs)
    slope1.funct<-vector("list", obs)
    slope2.funct<-vector("list", obs)
    jumptime<-numeric(length=obs)
    errorsd<-numeric(length=obs)
    jumpval<-numeric(length=obs)
    average1<-numeric(length=obs)
    true.start<-numeric(length=obs)
    intercept.2<-numeric(length=obs)
    average2<-numeric(length=obs)
    for(q in groupings){
      number<-length(which(groupref==q))
      eval(parse(text=paste("slope1[which(groupref==q)]<-as.character(" , part.1.slope[1] , "(", number, "," , slope1.groups.num[which(groupings==q)] , "," , as.numeric(part.1.slope[3])*sd.ratio , "))" , sep="")))
      
      slope1[which(groupref==q)]<-paste(slope1[which(groupref==q)], "*", slope.1.mult, "*t", sep="")
      
      eval(parse(text=paste("slope2[which(groupref==q)]<-as.character(" , part.2.slope[1] , "(", number, "," , slope2.groups.num[which(groupings==q)] , "," , as.numeric(part.2.slope[3])*sd.ratio , "))" , sep="")))
      
      slope2[which(groupref==q)]<-paste(slope2[which(groupref==q)], "*", slope.2.mult, "*t", sep="") 
      if(is.null(fixed.effect.part.1.slope)==F){
        slope1[which(groupref==q)]<-paste(fixed.effect.part.1.slope[which(groupref==q)],"*(",slope1[which(groupref==q)],")",sep="")
      }
      
      offset1<-paste("(t+", rnorm(n=number, offset.sd[groupings==q], offset.amount/2), ")", sep="")
      offset2<-paste("(t+", rnorm(n=number, offset.sd[groupings==q], offset.amount/2), ")", sep="")
      slope1[which(groupref==q)]<-mapply(gsub, replacement=offset1, x=slope1[which(groupref==q)], MoreArgs=list(pattern="t"))
      slope2[which(groupref==q)]<-mapply(gsub, replacement=offset2, x=slope2[which(groupref==q)], MoreArgs=list(pattern="t"))
      
      for (i in c(which(groupref==q))){
        eval(parse(text=paste("slope2.funct[[i]]<-as.function(alist(t=, ", slope2[i],"))", sep="")))
        eval(parse(text=paste("slope1.funct[[i]]<-as.function(alist(t=, ", slope1[i],"))", sep="")))
      }
      
      eval(parse(text=paste("jumptime[which(groupref==q)]<-" , time.step.change[1] , "(", number, "," , jumptime.groups[which(groupings==q)] , "," , as.numeric(time.step.change[3])*sd.ratio , ")" , sep="")))
      eval(parse(text=paste("errorsd[which(groupref==q)]<-abs(" , sd.error[1] , "(", number, "," , sd.error[2] , "," , sd.error[3] , "))" , sep="")))
      
      if (length(part.2.average)==0 | (length(part.2.average)>0 && length(part.1.average)>0 && length(value.step.change)>0)) {
        eval(parse(text=paste("jumpval[which(groupref==q)]<-" , value.step.change[1] , "(", number, "," , jumpval.groups[which(groupings==q)] , "," , as.numeric(value.step.change[3])*sd.ratio , ")" , sep="")))
        eval(parse(text=paste("average1[which(groupref==q)]<-" , part.1.average[1] , "(", number, "," , average1.groups[which(groupings==q)] , "," , as.numeric(part.1.average[3])*sd.ratio , ")" , sep=""))) 
        for (i in which(groupref==q)){
          true.start[i]<-(average1[i]*jumptime[i]-integrate(slope1.funct[[i]], lower=0, upper=jumptime[i] )$value)/jumptime[i]
          intercept.2[i]<-slope1.funct[[i]](jumptime[i])+true.start[i]+jumpval[i]-slope2.funct[[i]](jumptime[i])
          average2[i]<-(integrate(slope2.funct[[i]], upper=maxt, lower=jumptime[i] )$value+intercept.2[i]*(maxt-jumptime[i]))/(maxt-jumptime[i])
        }
        
      }
      if (length(part.1.average)==0) {
        eval(parse(text=paste("jumpval[which(groupref==q)]<-" , value.step.change[1] , "(", number, "," , jumpval.groups[which(groupings==q)] , "," , as.numeric(value.step.change[3])*sd.ratio, ")" , sep="")))
        eval(parse(text=paste("average2[which(groupref==q)]<-" , part.2.average[1] , "(", number, "," , average2.groups[which(groupings==q)] , "," , as.numeric(part.2.average[3])*sd.ratio , ")" , sep=""))) 
        for (i in which(groupref==q)){
          intercept.2[i]<-(average2[i]*(maxt-jumptime[i])-integrate(slope2.funct[[i]], lower=jumptime[i], upper=maxt )$value)/(maxt-jumptime[i])
          true.start[i]<-slope2.funct[[i]](jumptime[i])+intercept.2[i]-jumpval[i]-slope1.funct[[i]](jumptime[i])
          average1[i]<-(integrate(slope1.funct[[i]], upper=jumptime[i], lower=0 )$value+true.start[i]*jumptime[i])/jumptime[i]
        }
      }
      if (length(value.step.change)==0){
        eval(parse(text=paste("average1[which(groupref==q)]<-" , part.1.average[1] , "(", number, "," , average1.groups[which(groupings==q)] , "," , as.numeric(part.1.average[3])*sd.ratio , ")" , sep=""))) 
        eval(parse(text=paste("average2[which(groupref==q)]<-" , part.2.average[1] , "(", number, "," , average2.groups[which(groupings==q)] , "," , as.numeric(part.2.average[3])*sd.ratio , ")" , sep=""))) 
        for (i in which(groupref==q)){
          true.start[i]<-(average1[i]*jumptime[i]-integrate(slope1.funct[[i]], lower=0, upper=jumptime[i] )$value)/jumptime[i]
          intercept.2[i]<-(average2[i]*(maxt-jumptime[i])-integrate(slope2.funct[[i]], lower=jumptime[i], upper=maxt )$value)/(maxt-jumptime[i])
          jumpval[i]<-slope2.funct[[i]](jumptime[i])+intercept.2[i]-slope1.funct[[i]](jumptime[i])-true.start[i]
        }
      }
    }
    if(no.day.var==T){
      slope1<-rep(slope1[1:n.obs], dayreps)
      slope2<-rep(slope2[1:n.obs], dayreps)
      offset.sd<-rep(offset.sd[1:n.obs], dayreps)
      for(i in 2:dayreps){
        for(j in 1:n.obs){
          slope1.funct[[n.obs*(i-1)+j]]<-slope1.funct[[j]]
          slope2.funct[[n.obs*(i-1)+j]]<-slope2.funct[[j]]
          
        }
      }
      jumptime<-rep(jumptime[1:n.obs], dayreps)
      errorsd<-rep(errorsd[1:n.obs], dayreps)
      jumpval<-rep(jumpval[1:n.obs], dayreps)
      average1<-rep(average1[1:n.obs], dayreps)
      true.start<-rep(true.start[1:n.obs], dayreps)
      intercept.2<-rep(intercept.2[1:n.obs], dayreps)
      average2<-rep(average2[1:n.obs], dayreps)
    }
  } 
  if(exists("groupings")==FALSE){ #if no spatial groupings specified
    
    #set true parameter values for individuals depending on what has been specified
    eval(parse(text=paste("slope1<-" , part.1.slope[1] , "(", obs, "," , part.1.slope[2] , "," , part.1.slope[3] , ")" , sep="")))
    slope1<-paste(slope1, "*", slope.1.mult, "*t", sep="")
    eval(parse(text=paste("slope2<-" , part.2.slope[1] , "(", obs, "," , part.2.slope[2] , "," , part.2.slope[3] , ")" , sep="")))
    slope2<-paste(slope2, "*", slope.2.mult, "*t", sep="") 
    #add fixed effect if present
    if(is.null(fixed.effect.part.1.slope)==F){
      slope1<-paste(fixed.effect.part.1.slope,"*(",slope1,")",sep="")
    }
    
    slope1.funct<-vector("list", length(slope1))
    slope2.funct<-vector("list", length(slope2))
    
    offset1<-paste("(t+", rnorm(n=obs, 0, offset.amount), ")", sep="")
    offset2<-paste("(t+", rnorm(n=obs, 0, offset.amount), ")", sep="")
    slope1<-mapply(gsub, replacement=offset1, x=slope1, MoreArgs=list(pattern="t"))
    slope2<-mapply(gsub, replacement=offset2, x=slope2, MoreArgs=list(pattern="t"))
    
    
    
    for (i in 1:length(slope2)){
      eval(parse(text=paste("slope2.funct[[i]]<-as.function(alist(t=, ", slope2[i],"))", sep="")))
      eval(parse(text=paste("slope1.funct[[i]]<-as.function(alist(t=, ", slope1[i],"))", sep="")))
    }
    
    eval(parse(text=paste("jumptime<-" , time.step.change[1] , "(", obs, "," , time.step.change[2] , "," , time.step.change[3] , ")" , sep="")))
    eval(parse(text=paste("errorsd<-" , sd.error[1] , "(", obs, "," , sd.error[2] , "," , sd.error[3] , ")" , sep="")))
    errorsd<-abs(errorsd)
    
    
    #set true parameter values for non-specified parameter
    if (length(part.2.average)==0 | (length(part.2.average)>0 && length(part.1.average)>0 && length(value.step.change)>0)) {
      eval(parse(text=paste("jumpval<-" , value.step.change[1] , "(", obs, "," , value.step.change[2] , "," , value.step.change[3] , ")" , sep="")))
      eval(parse(text=paste("average1<-" , part.1.average[1] , "(", obs, "," , part.1.average[2] , "," , part.1.average[3] , ")" , sep=""))) 
      intercept.2<-numeric(length=obs)
      true.start<-numeric(length=obs)
      average2<-numeric(length=obs)
      for (i in 1:obs){
        true.start[i]<-(average1[i]*jumptime[i]-integrate(slope1.funct[[i]], lower=0, upper=jumptime[i] )$value)/jumptime[i]
        intercept.2[i]<-slope1.funct[[i]](jumptime[i])+true.start[i]+jumpval[i]-slope2.funct[[i]](jumptime[i])
        average2[i]<-(integrate(slope2.funct[[i]], upper=maxt, lower=jumptime[i] )$value+intercept.2[i]*(maxt-jumptime[i]))/(maxt-jumptime[i])
      }
      
    }
    if (length(part.1.average)==0) {
      eval(parse(text=paste("jumpval<-" , value.step.change[1] , "(", obs, "," , value.step.change[2] , "," , value.step.change[3] , ")" , sep="")))
      eval(parse(text=paste("average2<-" , part.2.average[1] , "(", obs, "," , part.2.average[2] , "," , part.2.average[3] , ")" , sep=""))) 
      intercept.2<-numeric(length=obs)
      true.start<-numeric(length=obs)
      average1<-numeric(length=obs)
      for (i in 1:obs){
        intercept.2[i]<-(average2*(maxt-jumptime[i])-integrate(slope2.funct[[i]], lower=jumptime[i], upper=maxt )$value)/(maxt-jumptime[i])
        true.start[i]<-slope2.funct[[i]](jumptime[i])+intercept.2[i]-jumpval[i]-slope1.funct[[i]](jumptime[i])
        average1[i]<-(integrate(slope1.funct[[i]], upper=jumptime[i], lower=0 )$value+true.start[i]*jumptime[i])/jumptime[i]
      }
    }
    if (length(value.step.change)==0){
      eval(parse(text=paste("average1<-" , part.1.average[1] , "(", obs, "," , part.1.average[2] , "," , part.1.average[3] , ")" , sep=""))) 
      eval(parse(text=paste("average2<-" , part.2.average[1] , "(", obs, "," , part.2.average[2] , "," , part.2.average[3] , ")" , sep=""))) 
      intercept.2<-numeric(length=obs)
      true.start<-numeric(length=obs)
      jumpval<-numeric(length=obs)
      for (i in 1:obs){
        true.start[i]<-(average1[i]*jumptime[i]-integrate(slope1.funct[[i]], lower=0, upper=jumptime[i] )$value)/jumptime[i]
        intercept.2[i]<-(average2*(maxt-jumptime[i])-integrate(slope2.funct[[i]], lower=jumptime[i], upper=maxt )$value)/(maxt-jumptime[i])
        jumpval[i]<-slope2.funct[[i]](jumptime[i])+intercept.2[i]-slope1.funct[[i]](jumptime[i])-true.start[i]
      }
    }
  }
  #set up more values
  
  
  
  
  #sample measurement times
  
  #if uniform sample uniformly for each column and order
  #if waves, sample each wave based on distribution and order -- resample if outside max time
  #if no sampling just have the measurement times
  if (length(record.times)==1) {
    times<-matrix(ncol=obs, nrow=record.times)
    for (i in 1:obs){
      times[,i]<-runif(record.times, 0, maxt)
      times[,i]<-times[order(times[,i]), i]
    }
  } 
  if (length(record.times)>1 & length(record.times)>0) {
    times<-matrix(ncol=obs, nrow=length(record.times))
    for (i in 1:obs){
      for (j in 1:length(record.times)){
        while(is.na(times[j,i])==T | times[j,i]<0 | times [j,i]>maxt){
          times[j,i]<-rnorm(1, record.times[j], record.times.vary)
        }
      }
      times[,i]<-times[order(times[,i]), i]
    }
  }
  if (length(record.times)==0){
    times<-matrix(data=rep(c(1:maxt), obs), ncol=obs, nrow=maxt)
  }
  
  
  #vector of error values for each measure for each individual together
  sd<-numeric(length=(nrow(times))*obs)
  for (x in 1:(nrow(times))) {
    for (y in 1:obs) { 
      sd[(nrow(times))*(y-1)+x]<-errorsd[y]
    }
  }
  
  #matrix of phases
  phases<-matrix(ncol=obs, nrow=nrow(times))
  for (i in 1:obs){
    for (j in 1:nrow(phases)){
      if (times[j,i]<=jumptime[i]){
        phases[j,i]<-1
      } else {
        phases[j,i]<-2
      }
    }
  }
  
  #matrix of true values
  value.true<-matrix(ncol=obs, nrow=nrow(times))
  for (i in 1:obs){
    for (j in 1:nrow(value.true)){
      if (phases[j,i]==1){
        value.true[j,i]<-true.start[i]+slope1.funct[[i]](times[j,i])#function 1 of time[j,i]
      } else {
        value.true[j,i]<-intercept.2[i]+slope2.funct[[i]](times[j,i])
      }
    }
  }
  
  
  #correlation and then covariance matrix using time and space differences
  covmat<-matrix(ncol=nrow(times)*obs, nrow=nrow(times)*obs)
  for (id1 in 1:obs){
    for (time1 in 1:nrow(times)){
      for(id2 in 1:obs){
        for(time2 in 1:nrow(times)){
          time.dif<-abs(times[time1,id1]-times[time2,id2])
          if (id1==id2){
            covmat[(nrow(times)*(id1-1))+time1, (nrow(times)*(id2-1))+time2]<-within.cor^time.dif
          } else {
            covmat[(nrow(times)*(id1-1))+time1, (nrow(times)*(id2-1))+time2]<-within.cor^time.dif*between.cor^(distance1[id1, id2]+distance2[id1, id2])*dayreps.cor^days[id1,id2]
          }
        }
      }
    }
  }
  
  covmat<-Covar(n=obs*(nrow(times)), sd<-sd, Cor=covmat)
  
  #sim data and add times to the data in wide format
  value<-matrix(as.numeric(mvrnorm(n=1, Sigma=data.frame(covmat), mu<-as.numeric(value.true))), nrow=(nrow(times)))
  value<-data.frame(t(value))
  times<-data.frame(t(times))
  for (i in 1:ncol(times)){
    colnames(times)[i]<-c(paste("time", i, sep="."))
    colnames(value)[i]<-c(paste("value", i, sep="."))
  }
  
  
  
  #generate summary stats for each individual
  summary<-data.frame()
  prejumpt<-numeric(length=obs)
  postjumpt<-numeric(length=obs)
  prepostchange<-numeric(length=obs)
  for (i in 1:obs) {
    
    if (jumptime[i]<min(times[i,])){
      
      prejumpt[i]<-NA
      postjumpt[i]<-0
      mean1<-NA
      realmean1<-NA
      jump<-NA
      prepostchange[i]<-NA
      postjumpt[i]<-min(times[i,which(times[i,]>jumptime[i])])
      mean2<-mean(as.numeric(value[i, which(times[i,]>=postjumpt[i])], na.rm=TRUE))
      realmean2<-(integrate(slope2.funct[[i]], lower=postjumpt[i], upper=maxt )$value+intercept.2[i]*(maxt-postjumpt[i]))/(maxt-postjumpt[i])
    } else{
      
      prejumpt[i]<-max(times[i,which(times[i,]<=jumptime[i])])
      mean1<-mean(as.numeric(value[i, which(times[i,]<=prejumpt[i])], na.rm=TRUE))
      realmean1<-(integrate(slope1.funct[[i]], lower=0, upper=prejumpt[i] )$value+true.start[i]*prejumpt[i])/prejumpt[i]
      if (jumptime[i]>=max(times[i,])){
        
        postjumpt[i]<-NA
        mean2<-NA
        realmean2<-NA
        jump<-NA
        prepostchange<-NA
      }
      if (jumptime[i]<max(times[i,])){
        
        postjumpt[i]<-min(times[i,which(times[i,]>jumptime[i])])
        prepostchange[i]<-slope2.funct[[i]](postjumpt[i])+intercept.2[i]-slope1.funct[[i]](prejumpt[i])-true.start[i]
        mean2<-mean(as.numeric(value[i, which(times[i,]>prejumpt[i])], na.rm=TRUE))
        realmean2<-(integrate(slope2.funct[[i]], lower=postjumpt[i], upper=maxt )$value+intercept.2[i]*(maxt-postjumpt[i]))/(maxt-postjumpt[i])
        jump<-value[i, which(times[i,]==postjumpt[i])]-value[i, which(times[i,]==prejumpt[i])]
        
      }
    }
    
    summary<-rbind(summary, c(i, 
                              average1[i], 
                              realmean1, 
                              mean1,
                              average2[i],  
                              realmean2, 
                              mean2,
                              jumptime[i], 
                              jumpval[i], 
                              prepostchange[i], 
                              jump,
                              errorsd[i]))
  }
  
  colnames(summary)<-c( "id",  
                        "average1.input", "average1.expected", "average1.output",
                        "average2.input", "average2.expected", "average2.output",
                        "jumptime.input", "jump.input", "jump.expected", "jump.output", "error.sd")
  #reshape data to long format
  value$id<-c(1:nrow(value))
  times$id<-c(1:nrow(times))
  value<-reshape(value, varying=colnames(value[,1:(ncol(value)-1)]), idvar<-"meas.no", direction="long", timevar<-"id", sep=".", v.names=c("value"))
  times<-reshape(times, varying=colnames(times[,1:(ncol(times)-1)]), idvar<-"meas.no", direction="long", timevar<-"id", sep=".", v.names=c("time"))
  
  value<-cbind(value, times$time)
  if (dayreps>1){
    value<-cbind(value, daynumbers)
    value$id<-rep(c(1:n.obs), dayreps)
    value$idday<-as.factor(paste(value$id, value$daynumbers, sep=":"))
  }
  if(exists("groupings")==T){
    value<-cbind(value, groupref)
  }
  if(exists("odinfo")==T){
    value<-cbind(value, odinfo)
  }
  colnames(value)[4]<-c("time")
  value<-value[order(value$id, value$meas.no),]
  #plots and stuff
  
  if (obs>10){
    if(dayreps==1){
      if(length(groupings)>1 & length(groupings)<20){
        p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(groupref)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Spatial group")+
          ggtitle("Simulated data patterns")
      } else{
        p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(id)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(legend.position="none", plot.title=element_text(size=10))+
          ggtitle("Simulated data patterns")
      }
    } else{
      if(length(groupings)>1 & length(groupings)<20){
        p0<-ggplot(data=value, aes(x=time, group=idday, colour=as.factor(groupref)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme( plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Spatial group")+
          ggtitle("Simulated data patterns")  
      } else{
        p0<-ggplot(data=value, aes(x=time, group=idday, colour=as.factor(id)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(legend.position="none", plot.title=element_text(size=10))+
          ggtitle("Simulated data patterns")  
      }
    }
  } else {
    if(dayreps==1){
      if(length(groupings)>1){
        p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(groupref)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Spatial group")+
          ggtitle("Simulated data patterns")
      } else{
        p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(id)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Observation ID")+
          ggtitle("Simulated data patterns")
      }
    } else{
      if(length(groupings)>1){
        p0<-ggplot(data=value, aes(x=time, group=idday, colour=as.factor(groupref)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Spatial group")+
          ggtitle("Simulated data patterns")  
      } else{
        p0<-ggplot(data=value, aes(x=time, group=idday, colour=as.factor(id)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Observation ID")+
          ggtitle("Simulated data patterns")  
      }
    }
  }
  
  
  summary$average1.dist<-summary$average1.output-summary$average1.expected
  summary$average1.dist.2<-summary$average1.output-summary$average1.input
  
  summary$average2.dist<-summary$average2.output-summary$average2.expected
  summary$average2.dist.2<-summary$average2.output-summary$average2.input
  
  summary$jump.dist<-summary$jump.output-summary$jump.expected
  summary$jump.dist.2<-summary$jump.output-summary$jump.input
  
  p3<-ggplot(summary, aes(average1.dist.2))+
    geom_density()+
    xlab("Distances from average 1 input")+
    theme(axis.title=element_text(size=10))
  
  p6<-ggplot(summary, aes(average2.dist))+
    geom_density()+
    xlab("Distances from average 2 input")+
    theme(axis.title=element_text(size=10))
  
  p7<-ggplot(summary, aes(jump.dist))+
    geom_density()+
    xlab("Distances from expected step change")+
    theme(axis.title=element_text(size=10))
  
  p8<-ggplot(summary, aes(jump.dist.2))+
    geom_density()+
    xlab("Distances from step change input")+
    theme(axis.title=element_text(size=10))
  
  summary.plot<-grid.arrange(grobs=list(p3,p6,p7,p8,p0),  layout_matrix=rbind(c(1,1,2,2), c(3,3,4,4), c(5,5,5,5), c(5,5,5,5)))
  if(plot==TRUE){
    summary.plot
  }
  #return list of data, summary and plot
  
  slope1<-paste(true.start, "+", slope1, sep="")
  slope2<-paste(intercept.2, "+", slope2, sep="") 
  slope1.funct<-list(length=length(slope1))
  slope2.funct<-list(length=length(slope2))
  for (i in 1:length(slope2)){
    eval(parse(text=paste("slope2.funct[[i]]<-as.function(alist(t=, ", slope2[i],"))", sep="")))
    eval(parse(text=paste("slope1.funct[[i]]<-as.function(alist(t=, ", slope1[i],"))", sep="")))
  }
  if (length(spacegroups==1) | is.null(fixed.spacegroups)==F){
    if (origin.dest==TRUE){
      ret.list<-list(data=value, 
                     simulation.summary=summary[,c(1:12)],
                     plot=summary.plot, 
                     section1.function=slope1.funct, 
                     section2.function=slope2.funct,
                     coordinates.info=list(spacegroup.info, individual.coords1, individual.coords2),
                     group.info=group.info,
                     distmat=list(distance1,distance2))
    }else{
      ret.list<-list(data=value, 
                     simulation.summary=summary[,c(1:12)],
                     plot=summary.plot, 
                     section1.function=slope1.funct, 
                     section2.function=slope2.funct,
                     coordinates.info=list(spacegroup.info, individual.coords1),
                     group.info=group.info, distmat=distance1)
    }
  } else {
    ret.list<-list(data=value, 
                   simulation.summary=summary[,c(1:12)],
                   plot=summary.plot, 
                   section1.function=slope1.funct, 
                   section2.function=slope2.funct)
  }
  return(ret.list)
  
}




locs<-st_read("stations.kml")
locs$Name<-as.character(locs$Name)
locs$Name<-gsub("\n\t\t\t","",locs$Name)
locs$Name<-gsub("\t","",locs$Name)

#add heathrow term 5, london city airport and woodlane
locadd<-locs[1:3,]
locadd$Name<-c("Heathrow Terminal 5", "London City Airport", "Wood Lane")
locadd$description<-NA
locadd$geometry[[1]]<-st_point(x=c(-0.4880,51.4723,0), dim="XYZ")
locadd$geometry[[2]]<-st_point(x=c(-0.0532,51.5032,0), dim="XYZ")
locadd$geometry[[3]]<-st_point(x=c( -0.2212824482,51.5054663115,0), dim="XYZ")
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

#limit to certain area
# preservelocs<-locs[locs$lat<=locs$lat[locs$fullname=="King's Cross St. Pancras Station"] &
#              locs$lat>=locs$lat[locs$fullname=="Waterloo Station"] &
#              locs$long>=locs$long[locs$fullname=="Marble Arch Station"] &
#              locs$long<=locs$long[locs$fullname=="Aldgate Station"],]
# 
# locs<-locs[locs$lat<=locs$lat[locs$fullname=="King's Cross St. Pancras Station"] &
#           locs$lat>=locs$lat[locs$fullname=="Waterloo Station"] &
#           locs$long>=locs$long[locs$fullname=="Marble Arch Station"] &
#           locs$long<=locs$long[locs$fullname=="Aldgate Station"],]
locs<-locs[locs$fullname%in%c("King's Cross St. Pancras Station","Euston Station","Euston Square Station","Great Portland Street Station",
                              "Mornington Crescent Station","Camden Town Station","Angel Station","Caledonian Road Station","Warren Street Station",
                              "Tottenham Court Road Station","Russell Square Station","Goodge Street Station"),]
preservelocs<-preservelocs[preservelocs$Name%in%locs$Name,]
locs$station_id<-as.numeric(as.character(locs$station_id))
netdata<-netdata[netdata$station_id%in%locs$station_id,]


#renumber the stations to include only the ones used - I think winbugs will have issues otherwise
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

#####################################################
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
fixeff<-numeric()
for(i in 1:nrow(redges.coord)){
  fixeff<-c(fixeff, numneighbours$num.neighbours[numneighbours$station_id==redges.coord$from[i]])
}
fixeff<-fixeff/2

######################################################
#simulate data
simdat<-simulate_trajectory(n.obs=nrow(redges.coord), 
                            max.time=24, 
                            record.times=50,
                            part.1.average=c("rnorm","25", "3") , 
                            part.2.average=NULL,
                            part.1.slope=c("rnorm", "4", "0.1", "((((cos(t/(24/(6*pi)))+1)*(cos((t-3)/(24/(2*pi)))*-0.8+0.7)+(cos((t-3)/(24/(2*pi)))*-0.8+1)))/t)") , 
                            part.2.slope=c("rnorm", "0", "1", "1") , 
                            time.step.change=c("rnorm", "24", "0") , 
                            value.step.change=c("rnorm", "0","0") , 
                            sd.error=c("rlnorm", "-1", "0.2"), 
                            between.cor=0.3, 
                            within.cor=0.3, 
                            dayreps=1,
                            fixed.spacegroups=redges.coord,
                            fixed.effect.part.1.slope=fixeff,
                            origin.dest=T, 
                            sd.ratio=0.9,
                            sd.ratio.day=0.5, 
                            plot=F,
                            offset.amount=1)

nfrom<-length(unique(redges.coord$from))
nto<-length(unique(redges.coord$to))

#extract data
journeytimes<-simdat$data
od<-matrix(unlist(strsplit(as.character(journeytimes$odinfo),":")),ncol=2,byrow=T)
journeytimes$from<-as.numeric(od[,1])
journeytimes$to<-as.numeric(od[,2])

journeytimes<-merge(journeytimes, redges.coord, by=c("from","to"))
journeytimes$id<-as.numeric(journeytimes$id)
journeytimes<-journeytimes[order(journeytimes$id, journeytimes$meas.no),]




#journeytimes<-journeytimes[journeytimes$time>=14 & journeytimes$time<=21,]

journeytimes<-journeytimes[,c("id","from","to","from.x","from.y","to.x","to.y","time","value")]
colnames(journeytimes)<-c("id","start.station_id","end.station_id","start.long","start.lat","end.long","end.lat","start.time","length")
journeytimes$length<-journeytimes$length/60
journeytimes<-journeytimes[order(journeytimes$start.station_id, journeytimes$end.station_id, journeytimes$start.time),]

################################# SET UP FOR WINBUGS ANALYSIS







###########################################################
## SET UP DISTANCE MATRICES                              ##
###########################################################

################### winbugs neighbour vectors start stations distance
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


################## get the tube map and work with the network neighbours
#reshape the net data and turn into graph (igraph)
netdatalong<-pivot_longer(netdata,cols= -c(station_id, Name, Zone), values_to="end")
netdatalong<-netdatalong[is.na(netdatalong$end)==F, c("station_id","Name","Zone","end")]
colnames(netdatalong)[1]<-"start"
netdatalong<-netdatalong[c("start","end")]
locations<-data.frame(id=usedlocs$station_id, x=sapply(usedlocs$geometry, function(x){return(x[1])}),
                      y=sapply(usedlocs$geometry, function(x){return(x[2])}))
tubenetwork<-graph.data.frame(netdatalong, vertices=locations,directed=F)


# plot(tubenetwork, layout=layout.norm(as.matrix(locations[,2:3])))
#graph has come out with two links between each connected node, so going to generate adjacency matrix and then turn all the 2s to 1s, turn into a graph agan
a<-as_adjacency_matrix(tubenetwork, sparse=F)
a[which(a==2)]<-1
tubenetwork2<-graph_from_adjacency_matrix(a, mode="undirected")
# plot(tubenetwork2, layout=layout.norm(as.matrix(locations[,2:3])))

network_distances<-distances(tubenetwork2) #this creates a distance matrix which we can use as a weight matrix. 
#we want to keep only information on distances below a certain value, lets say <=3 for now.
maxdist<-3
#for each row, want a list of cols which are <=maxdist
adj3<-apply(network_distances, 2, function(x){which(x<=maxdist & x>0)})
#num vector is num neighbours (with distance<maxdist) for each one
num3<-sapply(adj3, length)
adj3<-unlist(adj3)
sumNumNeigh3<-length(adj3)
#weights
weights3<-apply(network_distances, 2, function(x){x[which(x<=maxdist & x>0)]})
weights3<-unlist(weights3)
weights3[which(weights3!=0)]<-1/weights3[which(weights3!=0)]

names(adj3)<-NULL
names(weights3)<-NULL
names(num3)<-NULL

######################################################
## CHOICE OF BASIS KNOTS                            ##
######################################################

#centre time
journeytimes$journeyid<-paste(journeytimes$start.station_id, journeytimes$end.station_id,sep=":")
centre<-mean(journeytimes$start.time)
journeytimes$start.timec<-journeytimes$start.time-centre
journeytimes<-journeytimes[order(journeytimes$start.station_id, journeytimes$end.station_id,journeytimes$start.time),]
journeytimes<-journeytimes%>%group_by(journeyid)%>%
  mutate(timedif=(start.timec-dplyr::lag(start.timec,n=1,default=NA))*60) #in minutes
journeytimes$timesetup<-0
journeytimes$timesetup[is.na(journeytimes$timedif)==F]<-1
journeytimes$timedif[is.na(journeytimes$timedif)==T]<-0

lengthsd<-sd(journeytimes$length)
lengthmean<-mean(journeytimes$length)
journeytimes$scalelength<-(journeytimes$length-lengthmean)/lengthsd

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

sink(paste("model",fileref,".txt",sep=""))
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



cat( "alpha.1 ~ dflat() # intercept
     ")
for(i in 1:nbasis){
  cat(paste("betabs",i,".1 ~ dflat()\n", sep=""))
}


cat("}")
sink()


rand_vect_cont <- function(N, sd = 1) {
  vec <- rnorm(N/2, 0, sd)
  vec<-c(vec, vec*(-1))
  vec[order(sample(N,N))]
}

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
    
    initscom<-c(initscom, paste("alpha.",nn,".",k,"=rep(0,length(unique(journeytimes$",j,".station_id))),", sep=""))
    
  }
}

for(i in 1:nbasis){
  for(j in c("start","end")){
    if(j=="start"){nn<-2}else{nn<-3}
    for(k in c("space","unstr")){
      
      if(i==nbasis & j=="end" & k=="unstr"){
        initscom<-c(initscom, paste("betabs",i,".",nn,".",k,"=rep(0,length(unique(journeytimes$",j,".station_id))),", sep=""))
        initscom<-c(initscom,paste("tau.bs",i,".",j,".",k," = ",tau.i,")}", sep=""))
        
      }else{
        initscom<-c(initscom, paste("betabs",i,".",nn,".",k,"=rep(0,length(unique(journeytimes$",j,".station_id))),", sep=""))
        initscom<-c(initscom,paste("tau.bs",i,".",j,".",k," = ",tau.i,",", sep=""))
      }
    }
  }
}

initscom<-paste(initscom,collapse="")
initscom<-gsub("\n"," ",initscom)

eval(parse(text=paste(initscom)))
p<-proc.time()
#, 



params<-names(inits())
bugs.out <- bugs(data=bugsdata, inits=inits, 
                 parameters.to.save=params, model.file=paste("model",fileref,".txt",sep=""), 
                 n.chains=4, n.iter=ni, n.burnin=nb, n.thin=nt, debug=F, DIC=TRUE)

save.image(file=paste("renvironment0_",fileref,".RData",sep=""))
