simulate_trajectory <- function(n.obs=1 , 
                                max.time=20 , 
                                record.times=NULL,
                                record.times.vary=1,
                                part.1.average=character() , 
                                part.2.average=character() , 
                                part.1.slope=c("rnorm", "0", "1", "1") , 
                                part.2.slope=c("rnorm", "0", "1", "1") , 
                                time.step.change=c("rnorm", "0", "1") , 
                                value.step.change=character() , 
                                sd.error=c("rlnorm","0.5", "0.5"),
                                between.cor=0, 
                                within.cor=0,
                                dayreps=1, 
                                dayreps.cor=0, 
                                fixed.spacegroups=NULL,
                                fixed.effect.part.1.slope=NULL,
                                spacegroups=NULL, 
                                spacegroup.size=NULL, 
                                origin.dest=FALSE, 
                                area.x=NULL, 
                                area.y=NULL, 
                                sd.ratio=0.5, 
                                sd.ratio.day=0.5, 
                                print.plot=TRUE, 
                                offset.amount=0){ 
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
    if(class(fixed.spacegroups)[1]!="matrix" & class(fixed.spacegroups)[1]!="data.frame"){
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
  if(is.null(fixed.effect.part.1.slope)==F & length(fixed.effect.part.1.slope)!=n.obs){
    stop("If specified, fixed.effect.part.1.slope should be a numeric vector of length n.obs")
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
      sep.group.ref<-data.frame(spacegroup=rep(spacegroup1, dayreps))
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
      sep.group.ref<-data.frame(spacegroup=rep(paste(spacegroup1, spacegroup2, sep=":"), dayreps))
      
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
    sep.group.ref<-data.frame(spacegroup=rep(spacegroup1, dayreps))
    
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
        sep.group.ref<-data.frame(spacegroup=rep(paste(spacegroup1, spacegroup2, sep=":"), dayreps))
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
  sep.group.ref$day<-daynumbers
  
  
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
        if(sd.ratio.day!=0){
          offset.sd<-rnorm(length(groupings), 0, offset.amount)
        } else{
          offset.sd<-rep(rnorm(length(groupings)/dayreps, 0, offset.amount), dayreps)
        }
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
      if(time.step.change[1]=="rnorm"){
        jumptime.sd<-as.numeric(time.step.change[3])*sd.ratio
      }else{
        jumptime.sd<-1/sd.ratio
      }
      eval(parse(text=paste("jumptime[which(groupref==q)]<-rnorm(", number, "," , jumptime.groups[which(groupings==q)] , "," , jumptime.sd , ")" , sep="")))
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
  if (dayreps>1){
    value$id<-rep(c(1:n.obs), dayreps)
    times$id<-rep(c(1:n.obs), dayreps)
  }
  if(exists("groupings")==T){
    value<-cbind(value, sep.group.ref)
    value<-reshape(value, 
                    varying=colnames(value[,1:(ncol(value)-1-ncol(sep.group.ref))]), 
                    idvar<-"meas.no", 
                    direction="long", 
                    timevar<-"added", 
                    sep=".", 
                    v.names=c("value"))
    value<-value[,-which(colnames(value)=="added")]
    times<-cbind(times, sep.group.ref)
    times<-reshape(times, 
                   varying=colnames(times[,1:(ncol(times)-1-ncol(sep.group.ref))]), 
                   idvar<-"meas.no", 
                   direction="long", 
                   timevar<-"added", 
                   sep=".", 
                   v.names=c("time"))
    times<-times[,-which(colnames(times)=="added")]
    
  } else{
  value<-reshape(value, varying=colnames(value[,1:(ncol(value)-1)]), idvar<-"meas.no", direction="long", timevar<-"id", sep=".", v.names=c("value"))
  times<-reshape(times, varying=colnames(times[,1:(ncol(times)-1)]), idvar<-"meas.no", direction="long", timevar<-"id", sep=".", v.names=c("time"))
  }
  
  value$time<-times$time
  
  value<-value[order(value$id, value$meas.no),]
  #plots and stuff
  
  if (obs>10){
    if(dayreps==1){
      if(length(groupings)>1){
        if(length(groupings)<=20){
          p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(spacegroup)))+
            geom_line(aes(y=value))+
            geom_point(aes(y=value))+
            theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
            scale_colour_discrete(name="Spatial group")+
            ggtitle("Simulated data patterns")
        }else{
          p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(spacegroup)))+
            geom_line(aes(y=value))+
            geom_point(aes(y=value))+
            theme(plot.title=element_text(size=10), legend.position="none")+
            ggtitle("Simulated data patterns") 
        }
      } else{
        p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(id)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(legend.position="none", plot.title=element_text(size=10))+
          ggtitle("Simulated data patterns")
      }
    } else{
      if(length(groupings)>1){
        if(length(groupings)<=20){
          p0<-ggplot(data=value, aes(x=time, group=interaction(id, day), colour=as.factor(spacegroup)))+
            geom_line(aes(y=value))+
            geom_point(aes(y=value))+
            theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
            scale_colour_discrete(name="Spatial group")+
            ggtitle("Simulated data patterns")
        }else{
          p0<-ggplot(data=value, aes(x=time, group=interaction(id,day), colour=as.factor(spacegroup)))+
            geom_line(aes(y=value))+
            geom_point(aes(y=value))+
            theme(plot.title=element_text(size=10), legend.position="none")+
            ggtitle("Simulated data patterns") 
        }
      } else{
        p0<-ggplot(data=value, aes(x=time, group=interaction(id,day), colour=as.factor(id)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(legend.position="none", plot.title=element_text(size=10))+
          ggtitle("Simulated data patterns")  
      }
    }
  } else {
    if(dayreps==1){
      if(length(groupings)>1){
        p0<-ggplot(data=value, aes(x=time, group=id, colour=as.factor(spatial.groups)))+
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
        p0<-ggplot(data=value, aes(x=time, group=interaction(id,day), colour=as.factor(spatial.groups)))+
          geom_line(aes(y=value))+
          geom_point(aes(y=value))+
          theme(plot.title=element_text(size=10), legend.title=element_text(size=10))+
          scale_colour_discrete(name="Spatial group")+
          ggtitle("Simulated data patterns")  
      } else{
        p0<-ggplot(data=value, aes(x=time, group=interaction(id,day), colour=as.factor(id)))+
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
  if(print.plot==TRUE){
    plot(summary.plot)
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
                     group.info=group.info, 
                     distmat=distance1)
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
