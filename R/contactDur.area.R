#' Identify Environmental Contacts
#'
#' This function uses the output from distToArea to determine when tracked individuals are in "contact" with fixed locations. Individuals are said to be in a "contact" event (h) if they are observed within a given distance (<= dist.threshold) at a given timestep(i). Sec.threshold dictates the maximum amount of time a single, potential "contact" event should exist. For example, if sec.threshold=10, only "contacts" occurring within 10secs of one another will be regarded as a single "contact" event of duration sum(h). If in this case, a time difference between contacts was 11 seconds, the function will report two separate contact events.
#' The output of this function is a data frame containing a time-ordered contact edge set detailing animal-environment contacts.
#' @param x Output from the distToArea function (either df or sf variant). Can be either a data frame or non-data-frame list.
#' @param dist.threshold Numeric. Radial distance (in meters) within which "contact" can be said to occur. Defaults to 1. Note: If you are defining conttacts as occurring when polygons intersect, set dist.threshold to 0.
#' @param sec.threshold Numeric. Dictates the maximum amount of time between concurrent observations during which potential "contact" events remain unbroken. Defaults to 10. 
#' @param blocking Logical. If TRUE, contacts will be evaluated for temporal blocks spanning blockLength blockUnit (e.g., 6 hours) within the data set. Defaults to FALSE.
#' @param blockUnit Numerical. Describes the number blockUnits within each temporal block. Defaults to 1.
#' @param blockLength Character string taking the values, "secs," "mins," "hours," "days," or "weeks." Describes the temporal unit associated with each block. Defaults to "hours."
#' @param equidistant.time Logical. If TRUE, location fixes in individuals' movement paths are temporally equidistant (e.g., all fix intervals are 30seconds). Defaults to FALSE. Note: This is a time-saving argument. A sub-function here calculates the time difference (dt) between each location fix. If all fix intervals are identical, it saves a lot of time.
#' @param parallel Logical. If TRUE, sub-functions within the contactDur.all wrapper will be parallelized. Note that this can significantly speed up processing of relatively small data sets, but may cause R to crash due to lack of available memory when attempting to process large datasets. Defaults to TRUE.
#' @param reportParameters Logical. If TRUE, function argument values will be appended to output data frame(s). Defaults to TRUE.
#' @keywords data-processing contact
#' @export
#' @examples
#' #load the calves data set
#' data(calves)
#' 
#' #pre-process the data
#' calves.dateTime<-datetime.append(calves, date = calves$date, time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, point.y = calves.dateTime$y, secondAgg = 10, extrapolate.left = FALSE, extrapolate.right = FALSE, resolutionLevel = "Full", parallel = TRUE, na.rm = FALSE, smooth.type = 1) #smooth locations to 10-second fix intervals. Note that na.rm was set to "FALSE" because randomizing this data set according to Spiegel et al.'s method (see below) requires equidistant time points.
#'
#' #delineate the water trough polygon (showing where the water trough in the calves' feedlot pen is)
#' water_trough.x<- c(61.43315, 61.89377, 62.37518, 61.82622) #water x coordinates
#' water_trough.y<- c(62.44815 62.73341 61.93864 61.67411) #water y coordintates
#' water_poly<-data.frame(point1.x = water_trough.x[1], point1.y = water_trough.y[1], point2.x = water_trough.x[2], point2.y = water_trough.y[2], point3.x = water_trough.x[3], point3.y = water_trough.y[3], point4.x = water_trough.x[4], point4.y = water_trough.y[4])
#' 
#' #generate empirical time-ordered network edges.
#' water.dist<-dist2Area_df(x = calves.agg, y = water_poly, parallel = TRUE, x.id = calves.agg$id, y.id = "water", dateTime = calves.agg$dateTime, point.x = calves.agg$x, point.y = calves.agg$y, dataType = "Point", lonlat = FALSE) #calculate distance between all individuals and the water polygon at each timepoint.
#' water_contacts <- contactDur.area(water.dist, dist.threshold=1,sec.threshold=10, blocking = FALSE, blockUnit = "mins", blockLength = 10, equidistant.time = FALSE, parallel = TRUE, reportParameters = TRUE)
#' 
#' More examples coming later

contactDur.area<-function(x,dist.threshold=1,sec.threshold=10, blocking = TRUE, blockUnit = "mins", blockLength = 10, equidistant.time = FALSE, parallel = TRUE, reportParameters = TRUE){ 
  timeDifference = function(x){
    t1 = unlist(unname(x[1]))
    t2 = unlist(unname(x[2]))
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
    breakVec <- unname(distthreshold[,match(unname(unlist(x[2])), colnames(distthreshold))]) #fixed 1/8
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
      member1 = unlist(rep(unname(unlist(x[1])),length(contact.start)))
      member2 = unlist(rep(substring(unname(unlist(x[2])),9),length(contact.start))) #specific area id in "dist.to." colnames start at character 9
      dyad = paste(member1,"-",member2,sep="")
      times.start <- timeVec1[contact.start]
      times.finish <- timeVec1[contact.finish]
      durationTab = data.frame("indiv.id" = member1,"area.id" = member2,"contact.id" = dyad, "contactDuration" = durlengths, "contactStartTime" = times.start, "contactEndTime" = times.finish)
      
    }else{ #If there were no recorded contacts
      durationTab <- data.frame(matrix(ncol = 6, nrow = 0))
      colnames(durationTab) <- c("indiv.id","area.id","contact.id", "contactDuration", "contactStartTime", "contactEndTime")
    }
    return(durationTab)
  }
  contactMatrix.maker <- function(x,idVec1, dist.all){
    spec.dist = dist.all[which(dist.all$id == as.character(unlist(unname(x[1])))),]
    dt = spec.dist$dt
    timebreakVec <-which(dt > unname(unlist(x[3]))) #needed here. fixed 1/8
    distColumns<-colnames(dist.all)
    areaIDVec<- distColumns[grep("dist.to", distColumns)]
    distmat = data.matrix(spec.dist[,c(match(areaIDVec, names(spec.dist)))])
    distmat.noNA <-ifelse(is.na(distmat) == TRUE,1000000000,distmat)
    distthreshold<-ifelse(distmat.noNA<=as.numeric(unlist(unname(x[2]))),1,0) 
    idVecFrame = data.frame(unlist(rep(unlist(unname(x[1])),length(areaIDVec))),areaIDVec, unlist(rep(unlist(unname(x[2])),length(areaIDVec))),unlist(rep(unlist(unname(x[3])),length(areaIDVec))))  
    timeVec <- unname(spec.dist[,match("dateTime", colnames(spec.dist))])
    dateTimeFrame <- data.frame(timeVec)
    idDurations <- data.frame(data.table::rbindlist(apply(idVecFrame, 1, mat.breaker,distthreshold, timebreakVec, dateTimeFrame)))
    return(idDurations)
  }
  durFinder.noblock<-function(x,parallel, dist.threshold, sec.threshold, equidistant.time){
    dist.all<-x[order(x$id, x$dateTime),]
    idVec1 = unique(dist.all$id)
    if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
      dist.all.reduced$dt = 0
    }else{ #if equidistant.time == FALSE
      if(nrow(dist.all) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all is only created from blocks with observed contacts, there will never be a case where dist.all < 1)
        dist.all$dt = 0
      }else{ #if there's more than one row in dist.all
        timesFrame = data.frame(dist.all$dateTime[1:(nrow(dist.all) - 1)], dist.all$dateTime[2:nrow(dist.all)])
        
        if (parallel == TRUE){
          cl<-parallel::makeCluster(parallel::detectCores())
          timedif<-parallel::parApply(cl, timesFrame, 1, timeDifference); dist.all$dt = c(0, timedif) #timedif represents the time it takes to move from location i-1 to location i	
          parallel::stopCluster(cl)
        }else{
          timedif = apply(timesFrame, 1, timeDifference)
          dist.all$dt = c(0, timedif) 
        }
      }
    }
    comboFrame = data.frame(unique(dist.all$id),dist.threshold,sec.threshold)
    
    if (parallel == TRUE){
      cl<-parallel::makeCluster(parallel::detectCores())
      duration<-parallel::parApply(cl, comboFrame, 1, contactMatrix.maker,idVec1, dist.all)
      parallel::stopCluster(cl)
    }else{
      duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all)
    }
    durationTable = data.frame(data.table::rbindlist(duration))
    return(durationTable)
  }
  durFinder.block.List<-function(x,dist.threshold, sec.threshold, equidistant.time){
    dist.all<-x[order(x$id, x$dateTime),]
    idVec1 = unique(dist.all$id)
    if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
      dist.all.reduced$dt = 0
    }else{ #if equidistant.time == FALSE
      if(nrow(dist.all) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all is only created from blocks with observed contacts, there will never be a case where dist.all < 1)
        dist.all$dt = 0
      }else{ #if there's more than one row in dist.all
        timesFrame = data.frame(dist.all$dateTime[1:(nrow(dist.all) - 1)], dist.all$dateTime[2:nrow(dist.all)])
        timedif = apply(timesFrame, 1, timeDifference)
        dist.all$dt = c(0, timedif)
      }
    }
    comboFrame = data.frame(unique(dist.all$id),dist.threshold,sec.threshold)
    duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all)	
    durationTable = data.frame(data.table::rbindlist(duration))
    durationTable$block<- unique(dist.all$block)
    durationTable$block.start<- unique(dist.all$block.start)
    durationTable$block.end<- unique(dist.all$block.end)
    durationTable$numBlocks<- unique(dist.all$numBlocks)
    return(durationTable)
  }	
  list.breaker5 <-function(x,y,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters){
    
    input<- data.frame(y[unname(unlist(x[1]))])
    durationTable<-duration.generator2(input,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters)
    return(durationTable)     
  }
  duration.generator2 <-function(x,dist.threshold, sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters){
    
    if(blocking == TRUE){
      
      if(length(x$block) == 0){ #If there's no "block" column in dist.all output (i.e., blocking was "FALSE" when dist.all was run), then we need to define blocks here.
        
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
          minBlockTime<-dateTimeVec2[min(which(block == blockVec[f]))]
          minBlockTimeSeq[which(block == blockVec[f])] <- minBlockTime
          maxBlockTime<-dateTimeVec2[max(which(block == blockVec[f]))]
          maxBlockTimeSeq[which(block == blockVec[f])] <- maxBlockTime
        }
        x$block <- block
        x$block.start <- minBlockTimeSeq
        x$block.end <- maxBlockTimeSeq
        x$numBlocks <- max(blockVec) #the contactTest function will require thus information (i.e. the number of blocks in the dataset)
        
      }else{ #if length(x$block > 0)
        block = x$block
        blockVec<-unique(block)
      }
      blockList<-list()
      
      for(j in 1:length(blockVec)){
        blockL1<- list(x[which(block == blockVec[j]),])
        blockList<- c(blockList, blockL1)
      }
      
      originTab <-x
      duration.block = lapply(blockList, durFinder.block.List, dist.threshold, sec.threshold, equidistant.time)	
      durationTable <- data.frame(data.table::rbindlist(duration.block))
      
      #Here we change components of the durationTable to the appropriate data type.
      durationTable[,match("block", names(durationTable))]<- as.factor(durationTable[,match("block", names(durationTable))])
      durationTable[,match("block.start", names(durationTable))]<- as.factor(durationTable[,match("block.start", names(durationTable))])
      durationTable[,match("block.end", names(durationTable))]<- as.factor(durationTable[,match("block.end", names(durationTable))])
      durationTable[,match("numBlocks", names(durationTable))]<- as.factor(durationTable[,match("numBlocks", names(durationTable))])
      
    }else{ #If blocking == FALSE
      durationTable <- durFinder.noblock(x,parallel, dist.threshold, sec.threshold, equidistant.time)
    }
    
    #Here we change the rest of the components of the durationTable to the appropriate data type.
    durationTable[,match("indiv.id", names(durationTable))]<- as.factor(durationTable[,match("indiv.id", names(durationTable))])
    durationTable[,match("area.id", names(durationTable))]<- as.factor(durationTable[,match("area.id", names(durationTable))])
    durationTable[,match("contact.id", names(durationTable))]<- as.factor(durationTable[,match("contact.id", names(durationTable))])
    durationTable[,match("contactDuration", names(durationTable))]<- as.integer(durationTable[,match("contactDuration", names(durationTable))])
    
    if(nrow(durationTable) >0){ #added 1/11 to ensure that even if no contacts existed in the dataset, no error will be returned if reportParameters == TRUE
      if(reportParameters == TRUE){
        durationTable$distThreshold = dist.threshold
        durationTable$secThreshold = sec.threshold
        durationTable$equidistant.time = equidistant.time
        if(blocking == TRUE){
          if((length(blockLength) > 0) & (length(blockUnit) > 0)){
            durationTable$blockLength <- paste(blockLength,blockUnit, sep = " ")
          }
        }
      }
    }
    return(durationTable)
  }
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.dur <- apply(breakFrame, 1, list.breaker5,y = x,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.
    return(list.dur)
  }else{ #if(is.list(x) == FALSE)
    frame.dur<-duration.generator2(x,dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters)
    return(frame.dur)
  } 
} 