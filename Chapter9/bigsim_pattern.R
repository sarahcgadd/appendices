nc <- 4 #chains
ni <- 4000 #iterations
nb <- 3000 #discard as burn in
nt <- 1 #thinning parameter

bugsdir<-"C:\\Program Files\\WinBugs\\winbugs143_unrestricted\\winbugs14_full_patched\\WinBUGS14"

#result<-foreach(i=1:3) %do% {
  knotty<-list(c(4,10,12,14,20),
               c(4.5,6,7.5,12,16.5,18,19.5),
               c(2,4,6.5,9.5,12,14.5,17.5,20,22),
               c(2,4,6.5,9.5,12,14.5,17.5,20,22),
               c(4,10,12,14,20),
               c(2.5,4,5.5,10.5,16.5,18,19.5),
               c(2,4,6.5,9.5,12,14.5,18,19.5),
               c(2,4,6.5,9.5,12,14.5,18,19.5)
               ) #chose these parameters based on the functions specified below - I know where I'm expecting the peaks to be
  #this is the equaivalent of prior knowledge about the process going on e.g. you know that peak congestion is at 4pm or 
  #you know children grow most at 12 years of age or something
  
  
  patterns<-list()
  patterns[[1]]<-"((-1)*cos(t/(12/pi)))*2/t"
  patterns[[2]]<-"((-1)*cos(t/(6/pi))-0.1*t)/t"
  patterns[[3]]<-"((-1)*cos(t/(4/pi))-0.1*t)/t"
  patterns[[4]]<-"((-1)*cos(t/(4/pi))-0.01*(t-12)^2)*1.15/t"
  patterns[[5]]<- "(1/(((t-12)^2)+1))*4/t"
  patterns[[6]]<-"(0.5/(((t-18)^2)+1) + 1/(((t-4)^2)+1))*4/t"
  patterns[[7]]<-"(0.5/(((t-12)^2)+1) + 1/(((t-4)^2)+1) + 0.3/(((t-18)^2)+1))*4/t"
  patterns[[8]]<-"(1/(((t-12)^2)+1) + 0.7/(((t-4)^2)+1) + 0.7/(((t-18)^2)+1))*4/t"
  
  horoffs<-c(0, 2, 4)
  
  params<-expand.grid(pattern=c(1:8), horoff=c(1:3))
  knots<-knotty[[params$pattern[paramref]]]
  nbasis<-length(knots)+3
  
  
  library(plyr)
  library(dplyr)
  #load libraries
  library(igraph)
  library(lubridate)
  library(rootSolve)
  library(Deriv)
  library(sf)
  library(splines2)
  library(R2OpenBUGS)
#  options(mc.cores = 4)
 # options(mc.cores = 1)
 # nlinks=40
  
  
  
  
  simulate_trajectory <- function(n.obs=1 , max.time=20 , record.times=NULL, record.times.vary=1,
                                  part.1.average=character() , part.2.average=character() , 
                                  part.1.slope=c("rnorm", "0", "1", "1") , part.2.slope=c("rnorm", "0", "1", "1") , 
                                  time.step.change=c("rnorm", "0", "1") , value.step.change=character() , 
                                  sd.error=c("rlnorm","0.5", "0.5"), between.cor=0, within.cor=0,
                                  dayreps=1, dayreps.cor=0, no.day.var=F, fixed.spacegroups=NULL,
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
  
  coord<-as.data.frame(cbind(from=c(1:nstations), 
                                   x.from=runif(nstations, 0, 20),
                                   y.from=runif(nstations, 0, 20)))
  journeys<-expand.grid(from=c(1:nstations), to=c(1:nstations))
  journeys<-journeys[journeys$from!=journeys$to,]
  journeys<-merge(journeys, coord, by="from")
  colnames(coord)<-c("to","x.to","y.to")
  journeys<-merge(journeys, coord, by="to")
  journeys<-journeys[order(journeys$from, journeys$to),c("from","x.from","y.from","to","x.to","y.to")]
  
  pattern<-patterns[[params$pattern[paramref]]]

  simdat<-simulate_trajectory(n.obs=nrow(journeys), 
                              max.time=24, 
                              record.times=20,
                              part.1.average=c("rnorm","12", "3") , 
                              part.2.average=NULL,
                              part.1.slope=c("rnorm", "1", "0.1", pattern) , 
                              part.2.slope=c("rnorm", "0", "1", "1") , 
                              time.step.change=c("rnorm", "24", "0") , 
                              value.step.change=c("rnorm", "0","0") , 
                              sd.error=c("rlnorm", "-1", "0.5"), 
                              between.cor=0.3, 
                              within.cor=0.3, 
                              dayreps=1,
                              fixed.spacegroups=journeys,
                              origin.dest=T, 
                              sd.ratio=0.9,
                              sd.ratio.day=0.5, 
                              plot=F,
                              offset.amount=horoffs[params$horoff[paramref]])
  
  
  
  colnames(coord)<-c("id", "y", "x")
  coord<-coord[order(coord$id),]
  
  journeys$id<-c(1:nrow(journeys))
  
  redges<-journeys[order(journeys$id),c("id","from","to")]
  
  redges.coord<-journeys
  colnames(redges.coord)<-c("from","from.x","from.y","to","to.x","to.y","id")
  redges.coord<-redges.coord[,c("id","from","to","from.x","from.y","to.x","to.y")]
  nfrom<-length(unique(redges.coord$from))
  nto<-length(unique(redges.coord$to))
  
  #extract data
  data<-simdat$data
  
  data<-merge(data, redges, by="id", all.x=F)
  data<-data[,c("id","from","to","meas.no","time","value")]
  data<-data[order(data$id, data$meas.no),]
  
  
  centre<-mean(data$time)
  data$timec<-data$time-centre
  
  
  basis<-bSpline(data$timec,  knots=knots-centre, Boundary.knots = c(0,24)-centre)
  dspline<-dbs(data$timec, knots=knots-centre, Boundary.knots = c(0,24)-centre)
  dspline2<-dbs(data$timec, derivs=2L, knots=knots-centre, Boundary.knots = c(0,24)-centre)
  
  #generate relative distance matrix
  maxdist<-7 #cutoff - arbitrary, but remembering locations are simulated with coordinates random from 0-20 in both axes
  weights1<-list()
  adj1<-list()
  for(i in 1:nrow(coord)){
    weights1[[i]]<-apply(coord[-c(i),], 1, function(x){
      if(((x["x"]-coord$x[i])^2+(x["y"]-coord$y[i])^2)^0.5<=maxdist){
        return(((x["x"]-coord$x[i])^2+(x["y"]-coord$y[i])^2)^0.5)
      }else{
        return(NA)
      }
    })
    adj1[[i]]<-apply(coord[-c(i),], 1, function(x){
      if(((x["x"]-coord$x[i])^2+(x["y"]-coord$y[i])^2)^0.5<=maxdist){
        return(x["id"])
      }else{
        return(NA)
      }
    })
  }
  
  adj1<-lapply(adj1, function(x){return(x[which(is.na(x)==F)])})
  weights1<-lapply(weights1, function(x){return(x[which(is.na(x)==F)])})
  num1<-sapply(adj1, FUN=length)
  adj1<-unlist(adj1)
  weights1<-unlist(weights1)
  
  ########
  bugsdata<-list(value=data$value,
                 SID=data$from, EID=data$to,
                 adj1=adj1, weights1=weights1, num1=num1,
                 Nobs=nrow(data), Nstart=length(unique(data$from)),
                 Nend=length(unique(data$to)))
  
  a<-length(bugsdata)
  for(i in 1:nbasis){
    bugsdata[[i+a]]<-as.numeric(basis[,i])
  }
  names(bugsdata)[(a+1):(a+nbasis)]<-paste("bs",c(1:nbasis), sep="")
  
  sink("model1.txt")
  cat("model{
    for(i in 1:Nobs){
    value[i]~dnorm(a[i], tau.a)
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
      
        w<-1
        cat(paste("betabs",i,".2.",j,"[1:Nstart]~car.normal(adj",w,"[], weights",w,"[], num",w,"[], tau.bs",i,".start.",j,")\n", sep=""))
      }
  }
  
  cat("alpha.3.space[1:Nend]~car.normal(adj1[], weights1[], num1[], tau.int.end.space)
    ")
  
  for (i in 1:nbasis){
    for (j in c("space")){
      
        w<-1
        cat(paste("betabs",i,".3.",j,"[1:Nend]~car.normal(adj",w,"[], weights",w,"[], num",w,"[], tau.bs",i,".end.",j,")\n", sep=""))
      }
  }
  
  
  cat("
    ### priors on regression coefficients and variances
    
    tau.a ~ dgamma(0.1,0.1)
    sigma2.a <- 1/tau.a # residual error variance
    
    ### residual intercept variance
    tau.int.start.space ~ dgamma(0.1,0.1)
    sigma2.int.start.space <- 1/tau.int.start.space #between start station, with spatial cor
    tau.int.start.unstr ~ dgamma(0.1,0.1)
    sigma2.int.start.unstr <- 1/tau.int.start.unstr
    
    tau.int.end.space ~ dgamma(0.1,0.1)
    sigma2.int.end.space <- 1/tau.int.end.space #between end station, with spatial cor
    tau.int.end.unstr ~ dgamma(0.1,0.1)
    sigma2.int.end.unstr <- 1/tau.int.end.unstr
    
    ### variance in basis 1 coefficients
    ")
  
  for(i in 1:nbasis){
    for(j in c("start","end")){
      for (k in c("space","unstr")){
        cat(paste("tau.bs",i,".",j,".",k," ~ dgamma(0.1,0.1)\nsigma2.bs",i,".",j,".",k," <- 1/tau.bs",i,".",j,".",k,"\n", sep=""))
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

params=c("alpha.1","alpha.2","alpha.3",
         paste("betabs",c(1:nbasis),".1",sep=""),
         paste("betabs",c(1:nbasis),".2",sep=""),
         paste("betabs",c(1:nbasis),".3",sep=""))

initscom<-c("inits<-function() {list(
            alpha.1 = rnorm(1),
            
            tau.a = 1,
            
            tau.int.start.space = 1,
            tau.int.start.unstr = 1,
            tau.int.end.space = 1,
            tau.int.end.unstr = 1,")
for(i in 1:nbasis){
  initscom<-c(initscom, paste("betabs",i,".1=rnorm(1),",sep=""))
}
for(i in 1:nbasis){
  for(j in c("start","end")){
    for(k in c("space","unstr")){
      if(i==nbasis & j=="end" & k=="unstr"){
        initscom<-c(initscom,paste("tau.bs",i,".",j,".",k," = 1)}", sep=""))
        
      }else{
        initscom<-c(initscom,paste("tau.bs",i,".",j,".",k," = 1,", sep=""))
      }
    }
  }
}

initscom<-paste(initscom,collapse="")
initscom<-gsub("\n"," ",initscom)

eval(parse(text=paste(initscom)))




  errorcount<-NA
  warningcount<-NA
  tryCatch(expr={bugs.out <- pbugs(data=bugsdata, inits=inits, 
                                   parameters.to.save=params, model.file="model1.txt", 
                                   n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=F, DIC=TRUE,
                                   OpenBUGS.pgm = bugsdir, clusters=cores,
                                   working.directory=getwd())}, 
  error=function(e){errorcount<-paste(e)},
  warning=function(w){warningcount<-paste(w)})
  
  #get coefficients
  
  meancoefs<-bugs.out$mean
  sampcoefs<-bugs.out$sims.list
  indexes<-bugs.out$indexes.short
  indexes<-lapply(indexes, unlist)
  names(indexes)<-names(meancoefs)
  missing<-which(num1==0)
  
  #collect coefficients
  coefsfrom<-meancoefs[c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))]
  #order the coefficients according to their index
  for(i in c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))){
    coefsfrom[[i]]<-coefsfrom[[i]][order(indexes[[i]])]
  }
  #insert any missing ones (these are ones with no neighbours, they are not given random parameters)
  coefsfrom<-lapply(coefsfrom, function(x){
    x2<-x
    for(i in missing){
      if(i<=length(x2)){
      x2<-c(x2[1:(i-1)], 0, x2[i:length(x2)])
      }else{
        x2<-c(x2[1:(i-1)], 0)
      }
    }
    return(x2)
  })
  coefsfrom<-data.frame(matrix(nrow=nfrom, data=unlist(coefsfrom), byrow=F))
  colnames(coefsfrom)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
  rownames(coefsfrom)<-c(1:nfrom)
  
  coefsto<-meancoefs[c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))]
  for(i in c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))){
    coefsto[[i]]<-coefsto[[i]][order(indexes[[i]])]
  }
  coefsto<-lapply(coefsto, function(x){
    x2<-x
    for(i in missing){
      if(i<=length(x2)){
        x2<-c(x2[1:(i-1)], 0, x2[i:length(x2)])
      }else{
        x2<-c(x2[1:(i-1)], 0)
      }
    }
    return(x2)
  })
  coefsto<-data.frame(matrix(nrow=nto, data=unlist(coefsto), byrow=F))
  colnames(coefsto)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
  rownames(coefsto)<-c(1:nto)
  
  fixedcoefs<-meancoefs[c("alpha.1",paste("betabs",c(1:nbasis),".1", sep=""))]
  fixedcoefs<-unlist(fixedcoefs)
  names(fixedcoefs)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
  
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
  
  
  
  
  input<-redges
  input$max.val<-NA
  input$max.time<-NA
  input$rate.to.max<-NA
  input$min.val<-NA
  input$min.time<-NA
  input$max.rate<-NA
  input$max.rate.time<-NA
  input$min.rate<-NA
  input$min.rate.time<-NA
  
  
  for(i in 1:nrow(input)){
    
    tryCatch(expr={roots<-uniroot.all(pred.dif,  interval=c(0,24)-centre, togrp=input$to[i], from=input$from[i])},
             error={roots<-NULL}
    )
    tryCatch(expr={roots2<-uniroot.all(pred.dif2,  interval=c(0,24)-centre, togrp=input$to[i], from=input$from[i])},
             error={roots2<-NULL}
    )
    
    if(length(roots)>0){
      input$max.val[i]<-max(pred.val(roots, from=input$from[i], togrp=input$to[i]))
      input$max.time[i]<-roots[which.max(pred.val(roots, from=input$from[i], togrp=input$to[i]))]
      input$rate.to.max[i]<-(input$max.val[i]-pred.val(0, from=input$from[i], togrp=input$to[i]))/input$max.time[i]
      input$min.val[i]<-min(pred.val(roots, from=input$from[i], togrp=input$to[i]))
      input$min.time[i]<-roots[which.min(pred.val(roots, from=input$from[i], togrp=input$to[i]))]
    } else{
      input$max.val[i]<-NA
      input$max.time[i]<-NA
      input$rate.to.max[i]<-NA
      input$min.val[i]<-NA
      input$min.time[i]<-NA
    }
    if(length(roots2)>0){
      input$max.rate<-max(pred.dif(roots2, from=input$from[i], togrp=input$to[i]))
      input$max.rate.time[i]<-roots2[which.max(pred.dif(roots2, from=input$from[i], togrp=input$to[i]))]
      input$min.rate<-min(pred.dif(roots2, from=input$from[i], togrp=input$to[i]))
      input$min.rate.time[i]<-roots2[which.min(pred.dif(roots2, from=input$from[i], togrp=input$to[i]))]
    } else{
      input$max.rate<-NA
      input$max.rate.time[i]<-NA
      input$min.rate<-NA
      input$min.rate.time[i]<-NA
    }
  }
  input$max.time<-input$max.time+centre
  input$min.time<-input$min.time+centre
  input$max.rate.time<-input$max.rate.time+centre
  input$min.rate.time<-input$min.rate.time+centre
  
  sample<-sample(c(1:((ni-nb)/nc)), boots)
  
  #random.from wants to be 1000,15,7 (boots,groups,pars)
  
  random.from<-sampcoefs[c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))]
  #order the coefficients according to their index
  for(i in c("alpha.2",paste("betabs",c(1:nbasis),".2", sep=""))){
    random.from[[i]]<-random.from[[i]][,order(indexes[[i]])]
  }
  random.from<-lapply(random.from, function(x){
    x2<-x
    for(i in missing){
      if(i<=ncol(x2)){
        x2<-cbind(x2[,1:(i-1)], rep(0, nrow(x2)), x2[,i:ncol(x2)])
      }else{
        x2<-cbind(x2[,1:(i-1)], rep(0, nrow(x2)))
      }
    }
    return(x2)
  })
  dim1<-length(random.from[[1]])
  random.from<-unlist(random.from)
  random.from<-array(random.from, dim=c(dim1,nfrom,(nbasis+1)), dimnames=list(c(1:dim1), c(1:nfrom),
                                                                                       c("(Intercept)",paste("basis",c(1:nbasis), sep=""))))
  random.from<-random.from[sample,,]
  
   
  random.to<-sampcoefs[c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))]
  for(i in c("alpha.3",paste("betabs",c(1:nbasis),".3", sep=""))){
    random.to[[i]]<-random.to[[i]][,order(indexes[[i]])]
  }
  random.to<-lapply(random.to, function(x){
    x2<-x
    for(i in missing){
      if(i<=ncol(x2)){
        x2<-cbind(x2[,1:(i-1)], rep(0, nrow(x2)), x2[,i:ncol(x2)])
      }else{
        x2<-cbind(x2[,1:(i-1)], rep(0, nrow(x2)))
      }
    }
    return(x2)
  })
  random.to<-unlist(random.to)
  random.to<-array(random.to, dim=c(dim1,nto,(nbasis+1)), dimnames=list(c(1:dim1), c(1:nto),
                                                                               c("(Intercept)",paste("basis",c(1:nbasis), sep=""))))
  random.to<-random.to[sample,,]
  
  fixedboot<-sampcoefs[c("alpha.1",paste("betabs",c(1:nbasis),".1", sep=""))]
  fixedboot<-data.frame(matrix(ncol=length(fixedboot),data=unlist(fixedboot), byrow=F))
  fixedboot<-fixedboot[sample,]
  colnames(fixedboot)<- c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
  
  
  
  ##############################
  max.val.sim<-matrix(nrow=boots, ncol=nrow(redges))
  max.time.sim<-matrix(nrow=boots, ncol=nrow(redges))
  min.val.sim<-matrix(nrow=boots, ncol=nrow(redges))
  min.time.sim<-matrix(nrow=boots, ncol=nrow(redges))
  max.rate.sim<-matrix(nrow=boots, ncol=nrow(redges))
  max.rate.time.sim<-matrix(nrow=boots, ncol=nrow(redges))
  min.rate.sim<-matrix(nrow=boots, ncol=nrow(redges))
  min.rate.time.sim<-matrix(nrow=boots, ncol=nrow(redges))
  w.dens<-matrix(nrow=0, ncol=100)
  
  ##bootstrappy process
  
  
  for(r in 1:boots){
    if ((100*r/boots)%%10==0){
      cat(paste("...",100*r/boots, "%", sep=""))
    }
    
    
    #get random coefficients for boot number
    coefsfromboot<-data.frame(random.from[r,,])
    colnames(coefsfromboot)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
    coefstoboot<-data.frame(random.to[r,,])
    colnames(coefstoboot)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
    fixedcoefssim<-as.numeric(fixedboot[r,])
    names(fixedcoefssim)<-c("(Intercept)",paste("basis",c(1:nbasis), sep=""))
    
    inputboot<-redges
    
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
    pred.dif.boot2<-function(x, togrp,from) {
      togrp<-as.character(togrp)
      from<-as.character(from)
      newdata<-as.data.frame(predict(dspline2, newx=x))
      colnames(newdata)<-paste(rep("basis", ncol(newdata)), c(1:ncol(newdata)), sep="")
      coefficients<-fixedcoefssim
      coefficients[colnames(coefstoboot)]<-as.numeric(coefsfromboot[from,])+as.numeric(coefstoboot[togrp,])+coefficients[colnames(coefstoboot)]
      
      return(apply(coefficients[c(2:length(coefficients))]*t(newdata), 2, sum)) 
    }
    
    inputboot$max.val<-NA
    inputboot$max.time<-NA
    inputboot$rate.to.max<-NA
    inputboot$min.val<-NA
    inputboot$min.time<-NA
    inputboot$max.rate<-NA
    inputboot$max.rate.time<-NA
    inputboot$min.rate<-NA
    inputboot$min.rate.time<-NA
    
    
    
    for(i in 1:nrow(inputboot)){
      tryCatch(expr={roots<-uniroot.all(pred.dif.boot,  interval=c(0,24)-centre, togrp=inputboot$to[i], from=inputboot$from[i])},
               error={roots<-NULL}
      )
      tryCatch(expr={roots2<-uniroot.all(pred.dif.boot2,  interval=c(0,24)-centre, togrp=inputboot$to[i], from=inputboot$from[i])},
               error={roots2<-NULL}
      )
      
      if(length(roots)>0){
        inputboot$max.val[i]<-max(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))
        inputboot$max.time[i]<-roots[which.max(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))]
        inputboot$rate.to.max[i]<-(inputboot$max.val[i]-pred.val.boot(0, from=inputboot$from[i], togrp=inputboot$to[i]))/inputboot$max.time[i]
        inputboot$min.val[i]<-min(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))
        inputboot$min.time[i]<-roots[which.min(pred.val.boot(roots, from=inputboot$from[i], togrp=inputboot$to[i]))]
      } else{
        inputboot$max.val[i]<-NA
        inputboot$max.time[i]<-NA
        inputboot$rate.to.max[i]<-NA
        inputboot$min.val[i]<-NA
        inputboot$min.time[i]<-NA
      }
      if(length(roots2)>0){
        inputboot$max.rate[i]<-max(pred.dif.boot(roots2, from=inputboot$from[i], togrp=inputboot$to[i]))
        inputboot$max.rate.time[i]<-roots2[which.max(pred.dif.boot(roots2, from=inputboot$from[i], togrp=inputboot$to[i]))]
        inputboot$min.rate[i]<-min(pred.dif.boot(roots2, from=inputboot$from[i], togrp=inputboot$to[i]))
        inputboot$min.rate.time[i]<-roots2[which.min(pred.dif.boot(roots2, from=inputboot$from[i], togrp=inputboot$to[i]))]
      } else{
        inputboot$max.rate[i]<-NA
        inputboot$max.rate.time[i]<-NA
        inputboot$min.rate[i]<-NA
        inputboot$min.rate.time[i]<-NA
      }
      #plodadd<-data.frame(y=seq(from=0, to=24, length.out=200),
      #                    x=pred.val.boot(seq(from=0-centre, to=24-centre, length.out=200), from=inputboot$from[i], togrp=inputboot$to[i]),
      #                    grp=rep(r+1, 200),
      #                    col=rep(1,200))
      #eval(parse(text=paste("colnames(plodadd)<-colnames(plod",i,")", sep="")))
      #eval(parse(text=paste("plod",i,"<-rbind(plodadd,plod",i,")", sep="")))
    }
    inputboot$max.time<-inputboot$max.time+centre
    inputboot$min.time<-inputboot$min.time+centre
    inputboot$max.rate.time<-inputboot$max.rate.time+centre
    inputboot$min.rate.time<-inputboot$min.rate.time+centre
    
    #put your estimate of maximum time and value in a matrix
    max.val.sim[r,]<-inputboot$max.val
    max.time.sim[r,]<-inputboot$max.time
    min.val.sim[r,]<-inputboot$min.val
    min.time.sim[r,]<-inputboot$min.time
    max.rate.sim[r,]<-inputboot$max.rate
    max.rate.time.sim[r,]<-inputboot$max.rate.time
    min.rate.sim[r,]<-inputboot$min.rate
    min.rate.time.sim[r,]<-inputboot$min.rate.time
    
    #calculate weighted density
    wdensity<-function(time){
      nnodes<-length(unique(c(inputboot$from, inputboot$to)))
      timecent<-time-centre
      output<-rep(0, length(time))
      for (i in 1:nrow(inputboot)){
        output<-output+pred.val.boot(x=timecent, from=inputboot$from[i], togrp=inputboot$to[i])
      }
      output<-output/(nnodes*(nnodes-1))
      return(output)
    }
    
    w.dens<-rbind(w.dens, wdensity(time=seq(from=0, to=24, length.out=100)))
    
  }
  
  
  
  #get estimate for each group
  input$max.val.975<-apply(max.val.sim, 2, quantile, p=0.975, na.rm=T)
  input$max.time.975<-apply(max.time.sim, 2, quantile, p=0.975, na.rm=T)
  input$max.val.025<-apply(max.val.sim, 2, quantile, p=0.025, na.rm=T)
  input$max.time.025<-apply(max.time.sim, 2, quantile, p=0.025, na.rm=T)
  input$max.val.sd<-apply(max.val.sim, 2, sd,  na.rm=T)
  input$max.time.sd<-apply(max.time.sim, 2, sd,  na.rm=T)
  input$min.val.975<-apply(min.val.sim, 2, quantile, p=0.975, na.rm=T)
  input$min.time.975<-apply(min.time.sim, 2, quantile, p=0.975, na.rm=T)
  input$min.val.025<-apply(min.val.sim, 2, quantile, p=0.025, na.rm=T)
  input$min.time.025<-apply(min.time.sim, 2, quantile, p=0.025, na.rm=T)
  input$min.val.sd<-apply(min.val.sim, 2, sd,  na.rm=T)
  input$min.time.sd<-apply(min.time.sim, 2, sd,  na.rm=T)
  
  input$max.rate.975<-apply(max.rate.sim, 2, quantile, p=0.975, na.rm=T)
  input$max.rate.time.975<-apply(max.rate.time.sim, 2, quantile, p=0.975, na.rm=T)
  input$max.rate.025<-apply(max.rate.sim, 2, quantile, p=0.025, na.rm=T)
  input$max.rate.time.025<-apply(max.rate.time.sim, 2, quantile, p=0.025, na.rm=T)
  input$max.rate.sd<-apply(max.rate.sim, 2, sd,  na.rm=T)
  input$max.rate.time.sd<-apply(max.rate.time.sim, 2, sd,  na.rm=T)
  input$min.rate.975<-apply(min.rate.sim, 2, quantile, p=0.975, na.rm=T)
  input$min.rate.time.975<-apply(min.rate.time.sim, 2, quantile, p=0.975, na.rm=T)
  input$min.rate.025<-apply(min.rate.sim, 2, quantile, p=0.025, na.rm=T)
  input$min.rate.time.025<-apply(min.rate.time.sim, 2, quantile, p=0.025, na.rm=T)
  input$min.rate.sd<-apply(min.rate.sim, 2, sd,  na.rm=T)
  input$min.rate.time.sd<-apply(min.rate.time.sim, 2, sd,  na.rm=T)
  
  ids<-unique(data$id)[order(unique(data$id))]
  
  #real values
  
  a<-simdat$section1.function
  a<-a[input$id[order(input$id)]]
  input.real<-cbind(redges.coord[order(redges.coord$id),c("id","to","from")],max.val=NA, max.time=NA, rate.to.max=NA, min.val=NA, min.time=NA)
  b<-list(length=length(a))
  c<-list(length=length(a))
  for (i in 1:length(a)){
    b[[i]]<-Deriv(a[[i]])
    c[[i]]<-Deriv(b[[i]])
    
    tryCatch(expr={roots<-uniroot.all(b[[i]], c(0, 24))},
             error={roots<-NULL})
    tryCatch(expr={roots2<-uniroot.all(c[[i]], c(0, 24))},
             error={roots2<-NULL})
    
    if(length(roots)>0){
    input.real$max.val[i]<-max(a[[i]](roots))
    input.real$max.time[i]<-roots[which.max(a[[i]](roots))]
    input.real$rate.to.max[i]<-(input.real$max.val[i]-a[[i]](0))/(input.real$max.time[i]-(0))
    input.real$min.val[i]<-min(a[[i]](roots))
    input.real$min.time[i]<-roots[which.min(a[[i]](roots))]
    } else{
      input.real$max.val[i]<-NA
      input.real$max.time[i]<-NA
      input.real$rate.to.max[i]<-NA
      input.real$min.val[i]<-NA
      input.real$min.time[i]<-NA
    }
    if(length(roots2)>0){
    input.real$max.rate[i]<-max(b[[i]](roots2))
    input.real$max.rate.time[i]<-roots2[which.max(b[[i]](roots2))]
    input.real$min.rate[i]<-min(b[[i]](roots2))
    input.real$min.rate.time[i]<-roots2[which.min(b[[i]](roots2))]
    } else {
      input.real$max.rate[i]<-NA
      input.real$max.rate.time[i]<-NA
      input.real$min.rate[i]<-NA
      input.real$min.rate.time[i]<-NA
    }
  }
  
  colnames(input.real)<-paste(colnames(input.real), "real", sep=".")
  savedata<-merge(input, input.real, by.x="id", by.y="id.real")
  savedata$max.timedif<-savedata$max.time-savedata$max.time.real
  savedata$max.valdif<-savedata$max.val-savedata$max.val.real
  savedata$max.val.ci.width<-savedata$max.val.975-savedata$max.val.025
  savedata$max.time.ci.width<-savedata$max.time.975-savedata$max.time.025
  savedata$max.time.in.ci<-savedata$max.time.975>=savedata$max.time.real & savedata$max.time.025<=savedata$max.time.real
  savedata$max.val.in.ci<-savedata$max.val.975>=savedata$max.val.real & savedata$max.val.025<=savedata$max.val.real
  savedata$max.realval.as.z<-(savedata$max.val.real-savedata$max.val)/savedata$max.val.sd
  savedata$max.realtime.as.z<-(savedata$max.time.real-savedata$max.time)/savedata$max.time.sd
  savedata$min.timedif<-savedata$min.time-savedata$min.time.real
  savedata$min.valdif<-savedata$min.val-savedata$min.val.real
  savedata$min.val.ci.width<-savedata$min.val.975-savedata$min.val.025
  savedata$min.time.ci.width<-savedata$min.time.975-savedata$min.time.025
  savedata$min.time.in.ci<-savedata$min.time.975>=savedata$min.time.real & savedata$min.time.025<=savedata$min.time.real
  savedata$min.val.in.ci<-savedata$min.val.975>=savedata$min.val.real & savedata$min.val.025<=savedata$min.val.real
  savedata$min.realval.as.z<-(savedata$min.val.real-savedata$min.val)/savedata$min.val.sd
  savedata$min.realtime.as.z<-(savedata$min.time.real-savedata$min.time)/savedata$min.time.sd
  
  savedata$max.rate.timedif<-savedata$max.rate.time-savedata$max.rate.time.real
  savedata$max.ratedif<-savedata$max.rate-savedata$max.rate.real
  savedata$max.rate.ci.width<-savedata$max.rate.975-savedata$max.rate.025
  savedata$max.rate.time.ci.width<-savedata$max.rate.time.975-savedata$max.rate.time.025
  savedata$max.rate.time.in.ci<-savedata$max.rate.time.975>=savedata$max.rate.time.real & savedata$max.rate.time.025<=savedata$max.rate.time.real
  savedata$max.rate.in.ci<-savedata$max.rate.975>=savedata$max.rate.real & savedata$max.rate.025<=savedata$max.rate.real
  savedata$max.realrate.as.z<-(savedata$max.rate.real-savedata$max.rate)/savedata$max.rate.sd
  savedata$max.realrate.time.as.z<-(savedata$max.rate.time.real-savedata$max.time)/savedata$max.rate.time.sd
  savedata$min.rate.timedif<-savedata$min.rate.time-savedata$min.rate.time.real
  savedata$min.ratedif<-savedata$min.rate-savedata$min.rate.real
  savedata$min.rate.ci.width<-savedata$min.rate.975-savedata$min.rate.025
  savedata$min.rate.time.ci.width<-savedata$min.rate.time.975-savedata$min.rate.time.025
  savedata$min.rate.time.in.ci<-savedata$min.rate.time.975>=savedata$min.rate.time.real & savedata$min.rate.time.025<=savedata$min.rate.time.real
  savedata$min.rate.in.ci<-savedata$min.rate.975>=savedata$min.rate.real & savedata$min.rate.025<=savedata$min.rate.real
  savedata$min.realrate.as.z<-(savedata$min.rate.real-savedata$min.rate)/savedata$min.rate.sd
  savedata$min.realrate.time.as.z<-(savedata$min.rate.time.real-savedata$min.rate.time)/savedata$min.rate.time.sd
  savedata$modelerror<-errorcount
  savedata$modelwarning<-warningcount
  savedata$fileref<-fileref
  
 result<- savedata




