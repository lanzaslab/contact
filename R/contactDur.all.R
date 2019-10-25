#' Identify Inter-animal Contacts
#'
#' This function uses the output from dist.all to determine when and for how 
#'     long tracked individuals are in "contact" with one another. Individuals 
#'     are said to be in a "contact" event if they are observed within a given 
#'     distance (<= dist.threshold) at a given timestep. Contacts are broken 
#'     when individuals are observed outside the specified distance threshold 
#'     from one another for > sec.threshold seconds. Sec.threshold dictates the
#'     maximum amount of time between concurrent observations during which 
#'     potential "contact" events remain unbroken. For example, if 
#'     sec.threshold == 10, only "contacts" occurring within 10secs of one 
#'     another will be regarded as a single "contact" event of duration sum(h).
#'     If in this case, a time difference between contacts was 11 seconds, the 
#'     function will report two separate contact events.
#'     
#' The output of this function is a data frame containing a time-ordered 
#'     contact edge set detailing inter-animal contacts.
#' @param x Output from the dist.all function. Can be either a data frame or 
#'     non-data-frame list.
#' @param dist.threshold Numeric. Radial distance (in meters) within which 
#'     "contact" can be said to occur. Defaults to 1. Note: If you are 
#'     defining conttacts as occurring when polygons intersect, set 
#'     dist.threshold to 0.
#' @param sec.threshold Numeric. Dictates the maximum amount of time between 
#'     concurrent observations during which potential "contact" events remain 
#'     unbroken. Defaults to 10. 
#' @param blocking Logical. If TRUE, contacts will be evaluated for temporal 
#'     blocks spanning blockLength blockUnit (e.g., 6 hours) within the data 
#'     set. Defaults to FALSE.
#' @param blockUnit Numerical. Describes the number blockUnits within each 
#'     temporal block. Defaults to 1.
#' @param blockLength Character string taking the values, "secs," "mins," 
#'     "hours," "days," or "weeks." Describes the temporal unit associated with
#'     each block. Defaults to "hours."
#' @param equidistant.time Logical. If TRUE, location fixes in individuals' 
#'     movement paths are temporally equidistant (e.g., all fix intervals are 
#'     30 seconds). Defaults to FALSE. Note: This is a time-saving argument. 
#'     A sub-function here calculates the time difference (dt) between each 
#'     location fix. If all fix intervals in an individuals' path are 
#'     identical, it saves a lot of time.
#' @param parallel Logical. If TRUE, sub-functions within the contactDur.all 
#'     wrapper will be parallelized. Note that this can significantly speed up 
#'     processing of relatively small data sets, but may cause R to crash due 
#'     to lack of available memory when attempting to process large datasets. 
#'     Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'     parallel processes. Defaults to the maximum number of cores available
#'     (i.e., parallel::detectCores()).
#' @param reportParameters Logical. If TRUE, function argument values will be 
#'     appended to output data frame(s). Defaults to TRUE.
#' @keywords data-processing contact
#' @return Returns a data frame (or list of data frames if \code{x} is a 
#'    list of data frames) with the following columns:
#'    
#'    \item{dyadMember1}{The unique ID of an individual observed in contact 
#'    with a specified second individual.}
#'    \item{dyadMember2}{The unique ID of an individual observed in contact 
#'    with \code{dyadMember1}.}
#'    \item{dyadID}{The unique dyad ID used to identify the pair
#'    of individuals \code{dyadMember1} and \code{dyadMember2}.}    
#'    \item{contactDuration}{The number of sequential timepoints in \code{x} 
#'    that \code{dyadMember1} and \code{dyadMember2} were observed to be in 
#'    contact with one another.}
#'    \item{contactStartTime}{The timepoint in \code{x} at which contact
#'    between \code{dyadMember1} and \code{dyadMember2} begins.}
#'    \item{contactEndTime}{The timepoint in \code{x} at which contact
#'    between \code{dyadMember1} and \code{dyadMember2} ends.}
#'    
#'    If blocking == TRUE, the following columns are appended to the output
#'    data frame described above:
#'    
#'    \item{block}{Integer ID describing unique blocks of time during which 
#'    contacts occur.}
#'    \item{block.start}{The timepoint in \code{x} at which the \code{block}
#'    begins.}
#'    \item{block.end}{The timepoint in \code{x} at which the \code{block}
#'    ends.}
#'    \item{numBlocks}{Integer describing the total number of time blocks 
#'    observed within \code{x} at which the \code{block}}
#'     
#'    Finally, if reportParameters == TRUE function arguments 
#'    \code{distThreshold}, \code{secThreshold}, \code{equidistant.time},
#'    and \code{blockLength} (if applicable) will be appended to the 
#'    output data frame.
#' @export
#' @examples
#' 
#' data(calves)
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, time =
#'     calves$time) #create a dataframe with dateTime identifiers for location foxes
#'     
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag,
#'     dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x,
#'     point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE,
#'     extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE,
#'     na.rm = TRUE, smooth.type = 1) #smooth locations to 5-min fix intervals.
#' 
#' calves.dist<-dist2All_df(x = calves.agg, parallel = FALSE, dataType = "Point",
#'     lonlat = FALSE) #calculate distance between all individuals at each timepoint
#'     
#' calves.contact.block<-contactDur.all(x = calves.dist, dist.threshold=1,
#'     sec.threshold=10, blocking = TRUE, blockUnit = "hours", blockLength = 1,
#'     equidistant.time = FALSE, parallel = FALSE, reportParameters = TRUE)
#'     
#' calves.contact.NOblock<-contactDur.all(x = calves.dist, dist.threshold=1,
#'     sec.threshold=10, blocking = TRUE, blockUnit = "hours", blockLength = 1,
#'     equidistant.time = FALSE, parallel = FALSE, reportParameters = TRUE)

contactDur.all<-function(x,dist.threshold=1,sec.threshold=10, blocking = FALSE, blockUnit = "hours", blockLength = 1, equidistant.time = FALSE, parallel = FALSE, nCores = parallel::detectCores(), reportParameters = TRUE){ 
  
  timeDifference = function(x){
    t1 = unname(unlist(x[1]))
    t2 = unname(unlist(x[2]))
    dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
    return(dt)
  }
  datetime.append1 = function(x){
    
    timevec = x$dateTime
    daySecondList = lubridate::hour(timevec) * 3600 + lubridate::minute(timevec) * 60 + lubridate::second(timevec) #This calculates a day-second
    lub.dates = lubridate::date(x$dateTime)
    dateseq = unique(lub.dates)
    dayIDList = NULL
    dayIDseq = seq(1,(length(dateseq)),1)
    for(b in dayIDseq){
      ID = rep(b,length(which(lub.dates == dateseq[which(dayIDseq == b)])))
      dayIDList = c(dayIDList, ID)
    } 
    x$totalSecond = ((dayIDList - min(dayIDList))*86400) + daySecondList #This calculates the total second (the cumulative second across the span of the study's timeframe)
    
    return(x)
  }
  mat.breaker <-function(x, distthreshold, timebreakVec, dateTimeFrame){
    breakVec <- unname(distthreshold[,match(paste("dist.to.indiv_",as.character(x[2]),sep =""), colnames(distthreshold))])
    
    timebreakVec <- timebreakVec[timebreakVec != 1] #added/fixed 1/8, needed because if timebreakVec == 1, it means that the dt between an individuals' first point was >secThreshold seconds after the previous individuals' last point.
    
    if(length(which(breakVec == 1)) >0){ #If there are no 1s, then there's no reason to preoceed with calculations below
      timeVec1 <-dateTimeFrame[,1]
      
      #The for-loop below adjusts the values in distthreshold if dt indicates that enough time has passed between individuals to break a contact, but two observations describing sequential contacts exist across the span of this dt, the timebreakVec loop inserts a "0" between the contacts.
      
      if(length(timebreakVec) > 0){
        addtoi = 0 #everytime the vectors are adjusted, all observations below i will move down by one.
        for(i in timebreakVec){
          if(breakVec[(i + addtoi)] == 0){next}else{
            preBreak.dist <- breakVec[1:((i + addtoi) -1)]
            postBreak.dist <- breakVec[(i + addtoi):length(breakVec)]
            breakVec <-c(preBreak.dist,0,postBreak.dist)
            preBreak.time <- as.character(droplevels(timeVec1[1:((i + addtoi) -1)])) #fixed 1/8. if levels are not dropped here, there further calculations will be erroneous.						
            postBreak.time <- as.character(droplevels(timeVec1[(i + addtoi):length(timeVec1)])) #fixed 1/8
            timeVec2 <-c(preBreak.time,as.character(timeVec1[((i + addtoi) -1)]),postBreak.time) #fixed 1/8						
            timeVec1 <-as.factor(timeVec2) #fixed 1/8 ; had to make this a factor so that later "droplevels" commands will not trigger an error.
            addtoi = addtoi + 1
          }
        }
      }
      
      breakVal <- rle(breakVec)
      values <- unlist(breakVal[2])
      repTimes <- unlist(breakVal[1])
      
      finish<-unname(cumsum(unlist(breakVal[1])))
      start <-unname(((finish-repTimes)+1))
      contact.start <- start[which(values == 1)]
      contact.finish <- finish[which(values == 1)]
      durlengths = unname(repTimes[which(values == 1)])
      
      member1 = unlist(rep(x[1],length(contact.start)))
      member2 = unlist(rep(x[2],length(contact.start)))
      dyad = paste(member1,"-",member2,sep="")
      
      times.start <- timeVec1[contact.start]
      times.finish <- timeVec1[contact.finish]
      
      durationTab = data.frame("dyadMember1" = member1,"dyadMember2" = member2,"dyadID" = dyad, "contactDuration" = durlengths, "contactStartTime" = times.start, "contactEndTime" = times.finish)
      
    }else{ #If there were no recorded contacts
      
      durationTab <- data.frame(matrix(ncol = 6, nrow = 0))
      colnames(durationTab) <- c("dyadMember1","dyadMember2","dyadID", "contactDuration", "contactStartTime", "contactEndTime")
    }
    
    return(durationTab)
  }
  contactMatrix.maker <- function(x,idVec1, dist.all.reduced){
    spec.dist = dist.all.reduced[which(dist.all.reduced$id == as.character(unname(unlist(x[1])))),]
    dt = spec.dist$dt
    timebreakVec <-which(dt > unname(unlist(x[3]))) #needed here. fixed 1/8
    distmat = data.matrix(spec.dist[,c(match(paste("dist.to.indiv_",as.character(idVec1),sep =""), names(spec.dist)))])
    distmat.noNA <-ifelse(is.na(distmat) == TRUE,1000000000,distmat)
    distthreshold<-ifelse(distmat.noNA<=as.numeric(unname(unlist(x[2]))),1,0) 
    idVec.redac = idVec1[-c(1:which(idVec1 == unname(unlist(x[1]))))]
    idVecFrame = data.frame(unlist(rep(unname(unlist(x[1])),length(idVec.redac))),idVec.redac, unlist(rep(unname(unlist(x[2])),length(idVec.redac))),unlist(rep(unname(unlist(x[3])),length(idVec.redac))))  
    timeVec <- unname(spec.dist[,match("dateTime", colnames(spec.dist))])
    dateTimeFrame <- data.frame(timeVec)
    idDurations <- data.frame(data.table::rbindlist(apply(idVecFrame, 1, mat.breaker,distthreshold, timebreakVec, dateTimeFrame)))
    return(idDurations)
  }
  durFinder.noblock<-function(x,parallel, dist.threshold, sec.threshold, equidistant.time, nCores){
    dist.all<-x[order(x$id, x$dateTime),]
    idVec1 = unique(dist.all$id)
    dist.all.reduced <-dist.all[-which(dist.all$id == idVec1[length(idVec1)]),] #there's no need to process contacts associated with the last id values, because if they contacted any other individuals, the contacts would already be processed earlier on.
    
    if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
      dist.all.reduced$dt = 0
    }else{ #if equidistant.time == FALSE
      
      if(nrow(dist.all.reduced) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all.reduced is only created from blocks with observed contacts, there will never be a case where dist.all.reduced < 1)
        dist.all.reduced$dt = 0
      }else{ #if there's more than one row in dist.all.reduced
        timesFrame = data.frame(dist.all.reduced$dateTime[1:(nrow(dist.all.reduced) - 1)], dist.all.reduced$dateTime[2:nrow(dist.all.reduced)])
        if (parallel == TRUE){
          cl<-parallel::makeCluster(nCores)
          timedif<-parallel::parApply(cl, timesFrame, 1, timeDifference)
          dist.all.reduced$dt = c(0, timedif) #timedif represents the time it takes to move from location i-1 to location i	
          parallel::stopCluster(cl)
        }else{
          timedif = apply(timesFrame, 1, timeDifference) 
          dist.all.reduced$dt = c(0, timedif) 
        }
      }
    }
    comboFrame = data.frame(unique(dist.all.reduced$id),dist.threshold,sec.threshold)
    
    if (parallel == TRUE){
      cl<-parallel::makeCluster(nCores)
      duration<-parallel::parApply(cl, comboFrame, 1, contactMatrix.maker,idVec1, dist.all.reduced)
      parallel::stopCluster(cl)
    }else{
      duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all.reduced)	
    }
    durationTable = data.frame(data.table::rbindlist(duration))
    return(durationTable)
  }
  durFinder.block.List<-function(x,dist.threshold, sec.threshold, equidistant.time){
    dist.all<-x[order(x$id, x$dateTime),]
    idVec1 = unique(dist.all$id)
    dist.all.reduced <-dist.all[-which(dist.all$id == idVec1[length(idVec1)]),] #there's no need to process contacts associated with the last id values, because if they contacted any other individuals, the contacts would already be processed earlier on.
    
    if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
      dist.all.reduced$dt = 0
    }else{ #if equidistant.time == FALSE
      if(nrow(dist.all.reduced) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all.reduced is only created from blocks with observed contacts, there will never be a case where dist.all.reduced < 1)
        dist.all.reduced$dt = 0
      }else{ #if there's more than one row in dist.all.reduced
        timesFrame = data.frame(dist.all.reduced$dateTime[1:(nrow(dist.all.reduced) - 1)], dist.all.reduced$dateTime[2:nrow(dist.all.reduced)])
        timedif = apply(timesFrame, 1, timeDifference)
        dist.all.reduced$dt = c(0, timedif) 
      }
    }
    comboFrame = data.frame(unique(dist.all.reduced$id),dist.threshold,sec.threshold)
    duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all.reduced)	
    durationTable = data.frame(data.table::rbindlist(duration))
    if(nrow(durationTable) > 0){ #if there was at least one contact duration, block information is appended to the data frame. If there are no observations, durationTable becomes NULL
      durationTable$block<- unique(dist.all$block)
      durationTable$block.start<- unique(dist.all$block.start)
      durationTable$block.end<- unique(dist.all$block.end)
      durationTable$numBlocks<- unique(dist.all$numBlocks)
    }else{
      durationTable <- NULL
    }
    return(durationTable)
  }
  
  list.breaker2<-function(x,y,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores){ 
    input<- data.frame(y[unname(unlist(x[1]))])
    durationTable<-duration.generator1(input,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores)
    return(durationTable)    
  }
  duration.generator1<-function(x,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores){
    if(blocking == TRUE){
      
      if(length(x$block) == 0 ){ #If there's no "block" column in dist.all output, then we need to define blocks here. Note here that, if individuals wanted to append block information to the dist.all file (for whatever reason), they may do so using the timeblock.append function. If users did this, there's no need to waste time remaking blocks here. 
        
        if(blockUnit == "Secs" || blockUnit == "SECS" || blockUnit == "secs"){
          blockLength1 <- blockLength
        }
        if(blockUnit == "Mins" || blockUnit == "MINS" || blockUnit == "mins"){
          blockLength1 <- blockLength*60 #num seconds in a minute
        }
        if(blockUnit == "Hours" || blockUnit == "HOURS" || blockUnit == "hours"){
          blockLength1 <- blockLength*60*60 #num seconds in an hour
        }
        if(blockUnit == "Days" || blockUnit == "DAYS" || blockUnit == "days"){
          blockLength1 <- blockLength*60*60*24 #num seconds in a day
        }
        if(blockUnit == "Weeks" || blockUnit == "WEEKS" || blockUnit == "weeks"){
          blockLength1 <- blockLength*60*60*24*7 #num seconds in a week
        }
        
        x<-x[order(x$dateTime),] #Just in case the data wasn't already ordered in this way.
        x<-datetime.append1(x) #adds the total second column to the dataframe
        studySecond <- (x$totalSecond -min(x$totalSecond)) + 1
        numblocks <- ceiling((max(studySecond) - 1)/blockLength1)
        block <-rep(0,length(studySecond))
        for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
          block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
        }
        if(length(which(block == 0)) > 0){ #identifies the last block
          block[which(block == 0)] = numblocks
        }
        blockVec<-unique(block)
        dateTimeVec2<-x$dateTime
        minBlockTimeSeq <- rep(0, length(block)) #Added 2/4/2019. This vector will identify the minimum timepoint in each block.
        maxBlockTimeSeq <- rep(0, length(block)) #Added 2/4/2019. This vector will identify the maximum timepoint in each block.
        for(f in 1:length(blockVec)){
          minBlockTime<-as.character(dateTimeVec2[min(which(block == blockVec[f]))])
          minBlockTimeSeq[which(block == blockVec[f])] <- minBlockTime
          maxBlockTime<-as.character(dateTimeVec2[max(which(block == blockVec[f]))])
          maxBlockTimeSeq[which(block == blockVec[f])] <- maxBlockTime
        }
        x$block <- block
        x$block.start <- minBlockTimeSeq
        x$block.end <- maxBlockTimeSeq
        x$numBlocks <- max(blockVec) #the contactTest function will require thus information (i.e. the number of blocks in the dataset)

        blockList<-list()
        
        for(j in 1:length(blockVec)){
          blockL1<- list(x[which(block == blockVec[j]),])
          blockList<- c(blockList, blockL1)
        }
        
      }else{ #if importBlocks == TRUE and length(x$block > 0)
        block = x$block
        blockVec<-unique(block)
        blockList<-list()
        for(j in 1:length(blockVec)){
          blockL1<- list(x[which(block == blockVec[j]),])
          blockList<- c(blockList, blockL1)
        }  
      }
      duration.block = lapply(blockList, durFinder.block.List, dist.threshold, sec.threshold, equidistant.time)	
      durationTable <- data.frame(data.table::rbindlist(duration.block))
      if(nrow(durationTable) > 0){ #Here we change components of the durationTable to the appropriate data type. This only occurs if there was at least one observation, however (otherwise an error will be produced). 
        durationTable[,match("block", names(durationTable))]<- as.factor(durationTable[,match("block", names(durationTable))])
        durationTable[,match("block.start", names(durationTable))]<- as.factor(durationTable[,match("block.start", names(durationTable))])
        durationTable[,match("block.end", names(durationTable))]<- as.factor(durationTable[,match("block.end", names(durationTable))])
        durationTable[,match("numBlocks", names(durationTable))]<- as.factor(durationTable[,match("numBlocks", names(durationTable))])       
      }
    }else{ #If blocking == FALSE
      durationTable <- durFinder.noblock(x,parallel, dist.threshold, sec.threshold, equidistant.time, nCores)
    }
    if(nrow(durationTable) > 0){ 
      #Here we change the rest of the components of the durationTable to the appropriate data type.
      durationTable[,match("dyadMember1", names(durationTable))]<- as.factor(durationTable[,match("dyadMember1", names(durationTable))])
      durationTable[,match("dyadMember2", names(durationTable))]<- as.factor(durationTable[,match("dyadMember2", names(durationTable))])
      durationTable[,match("dyadID", names(durationTable))]<- as.factor(durationTable[,match("dyadID", names(durationTable))])
      durationTable[,match("contactDuration", names(durationTable))]<- as.integer(durationTable[,match("contactDuration", names(durationTable))])
      if(reportParameters == TRUE){
        durationTable$distThreshold = dist.threshold
        durationTable$secThreshold = sec.threshold
        durationTable$equidistant.time = equidistant.time
        if(blocking == TRUE){
          durationTable$blockLength <- paste(blockLength,blockUnit, sep = " ")
        }
      }
    }
    return(durationTable)
  }
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.dur <- apply(breakFrame, 1, list.breaker2,y = x,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.
    return(list.dur)
    
  }else{ #if(is.list(x) == FALSE)
    
    frame.dur<-duration.generator1(x,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores)
    return(frame.dur)
  } 
}
