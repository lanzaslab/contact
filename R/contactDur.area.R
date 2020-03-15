#' Identify Environmental Contacts
#'
#' This function uses the output from dist2Area to determine when tracked 
#'    individuals are in "contact" with fixed locations. Individuals are said 
#'    to be in a "contact" event (h) if they are observed within a given 
#'    distance (<= dist.threshold) at a given timestep(i). Sec.threshold 
#'    dictates the maximum amount of time a single, potential "contact" event 
#'    should exist. For example, if sec.threshold=10, only "contacts" occurring
#'    within 10secs of one another will be regarded as a single "contact" event
#'    of duration sum(h). If in this case, a time difference between contacts 
#'    was 11 seconds, the function will report two separate contact events.
#'    
#' The output of this function is a data frame containing a time-ordered 
#'    contact edge set detailing animal-environment contacts.
#' @param x Output from the dist2Area function (either df or sf variant). Can 
#'    be either a data frame or non-data-frame list.
#' @param dist.threshold Numeric. Radial distance (in meters) within which 
#'    "contact" can be said to occur. Defaults to 1. Note: If you are defining 
#'    conttacts as occurring when polygons intersect, set dist.threshold to 0.
#' @param sec.threshold Numeric. Dictates the maximum amount of time between 
#'    concurrent observations during which potential "contact" events remain 
#'    unbroken. Defaults to 10. 
#' @param blocking Logical. If TRUE, contacts will be evaluated for temporal 
#'    blocks spanning blockLength blockUnit (e.g., 6 hours) within the data 
#'    set. Defaults to FALSE.
#' @param blockLength Integer. Describes the number blockUnits within each 
#'    temporal block. Defaults to 1.
#' @param blockUnit Character string taking the values, "secs," "mins," 
#'    "hours," "days," or "weeks." Describes the temporal unit associated with 
#'    each block. Defaults to "hours."
#' @param equidistant.time Logical. If TRUE, location fixes in individuals' 
#'    movement paths are temporally equidistant (e.g., all fix intervals are 30
#'    seconds). Defaults to FALSE. Note: This is a time-saving argument. A 
#'    sub-function here calculates the time difference (dt) between each 
#'    location fix. If all fix intervals are identical, it saves a lot of time.
#' @param parallel Logical. If TRUE, sub-functions within the contactDur.all 
#'    wrapper will be parallelized. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
#' @param reportParameters Logical. If TRUE, function argument values will be 
#'    appended to output data frame(s). Defaults to TRUE.
#' @keywords data-processing contact
#' @return Returns a data frame (or list of data frames if \code{x} is a 
#'    list of data frames) with the following columns:
#'    
#'    \item{indiv.id}{The unique ID of an individual observed in contact 
#'    with a specified fixed point/polygon.}
#'    \item{area.id}{The unique ID of a fixed point/polygon observed in contact
#'    with \code{indiv.id}.}
#'    \item{contact.id}{The unique ID used to identify contacts between the 
#'    \code{indiv.id} and \code{contact.id} pair.}    
#'    \item{contactDuration}{The number of sequential timepoints in \code{x} 
#'    that \code{indiv.id} and \code{area.id} were observed to be in 
#'    contact.}
#'    \item{contactStartTime}{The timepoint in \code{x} at which contact
#'    between \code{indiv.id} and \code{area.id} begins.}
#'    \item{contactEndTime}{The timepoint in \code{x} at which contact
#'    between \code{indiv.id} and \code{area.id} ends.}
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
#' @import foreach
#' @export
#' @examples
#' data(calves)
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'   time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#'
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth to 5-min fix intervals.
#'
#' water<- data.frame(x = c(61.43315, 61.89377, 62.37518, 61.82622),
#'                   y = c(62.44815, 62.73341, 61.93864, 61.67411)) 
#'
#' water_poly<-data.frame(matrix(ncol = 8, nrow = 1)) #(ncol = number of vertices)*2 #arrange data
#' colnum = 0
#' for(h in 1:nrow(water)){
#'  water_poly[1,colnum + h] <- water$x[h] #pull the x location for each vertex
#'  water_poly[1, (colnum + 1 + h)] <- water$y[h] #pull the y location for each vertex
#'  colnum <- colnum + 1
#' }
#'
#' water_dist<-dist2Area_df(x = calves.agg, y = water_poly, 
#'   x.id = calves.agg$id, y.id = "water", dateTime = "dateTime", point.x = calves.agg$x, 
#'   point.y = calves.agg$y, poly.xy = NULL, parallel = FALSE, dataType = "Point", 
#'   lonlat = FALSE, numVertices = NULL) #find distances to the water trough 
#'
#' water_contacts <- contactDur.area(water_dist, dist.threshold=1,
#'   sec.threshold=10, blocking = FALSE, blockUnit = "mins", blockLength = 10,
#'   equidistant.time = FALSE, parallel = FALSE, reportParameters = TRUE)
#'   

contactDur.area<-function(x,dist.threshold=1,sec.threshold=10, blocking = FALSE, blockLength = 1, blockUnit = "hours", equidistant.time = FALSE, parallel = FALSE, nCores = (parallel::detectCores()/2), reportParameters = TRUE){ 

  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  l <- NULL
  k <- NULL
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    
    listBreak_dur.generator<-function(x, dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores){ #this function just calls what happens when x is not a list o data frames
    
      durFinder.noblock<-function(parallel, dist.threshold, sec.threshold, equidistant.time, nCores, environmentTag){
        
        eval(expr = { #evaluate these function steps in the master-function environment to avoid cloning x (i.e., dist2... ouput), which can be quite large
          
          #write nested sub functions
          
          timeDifference = function(x){
            t1 = unlist(unname(x[1]))
            t2 = unlist(unname(x[2]))
            dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
            return(dt)
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
              durationTab = data.frame("indiv.id" = member1,"area.id" = member2,"contact.id" = dyad, "contactDuration" = durlengths, "contactStartTime" = times.start, "contactEndTime" = times.finish, stringsAsFactors = TRUE)
              
            }else{ #If there were no recorded contacts
              durationTab <- data.frame(matrix(ncol = 6, nrow = 0), stringsAsFactors = TRUE)
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
            colnames(distthreshold)<-areaIDVec #relate the distance threshold columns to each area id
            idVecFrame = data.frame(unlist(rep(unlist(unname(x[1])),length(areaIDVec))),areaIDVec, unlist(rep(unlist(unname(x[2])),length(areaIDVec))),unlist(rep(unlist(unname(x[3])),length(areaIDVec))), stringsAsFactors = TRUE)  
            timeVec <- unname(spec.dist[,match("dateTime", colnames(spec.dist))])
            dateTimeFrame <- data.frame(timeVec, stringsAsFactors = TRUE)
            idDurations <- data.frame(data.table::rbindlist(apply(idVecFrame, 1, mat.breaker,distthreshold, timebreakVec, dateTimeFrame)), stringsAsFactors = TRUE)
            return(idDurations)
          }
          
          dist.all<-x[order(x$id, x$dateTime),]
          idVec1 = unique(dist.all$id)
          if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
            dist.all.reduced$dt = 0
          }else{ #if equidistant.time == FALSE
            if(nrow(dist.all) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all is only created from blocks with observed contacts, there will never be a case where dist.all < 1)
              dist.all$dt = 0
            }else{ #if there's more than one row in dist.all
              timesFrame = data.frame(dist.all$dateTime[1:(nrow(dist.all) - 1)], dist.all$dateTime[2:nrow(dist.all)], stringsAsFactors = TRUE)
              
              if (parallel == TRUE){
                cl<-parallel::makeCluster(nCores)
                on.exit(parallel::stopCluster(cl))
                timedif<-parallel::parApply(cl, timesFrame, 1, timeDifference); dist.all$dt = c(0, timedif) #timedif represents the time it takes to move from location i-1 to location i	
              }else{
                timedif = apply(timesFrame, 1, timeDifference)
                dist.all$dt = c(0, timedif) 
              }
            }
          }
          comboFrame = data.frame(unique(dist.all$id),dist.threshold,sec.threshold, stringsAsFactors = TRUE)
          
          if (parallel == TRUE){
            duration<-parallel::parApply(cl, comboFrame, 1, contactMatrix.maker,idVec1, dist.all) #note that cl is defined above
          }else{
            duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all)
          }
          durationTable = data.frame(data.table::rbindlist(duration), stringsAsFactors = TRUE)
          
        }, envir = environmentTag)
        
        return(durationTable)
      }
      durFinder.block.List<-function(x,dist.threshold, sec.threshold, equidistant.time){
        
        timeDifference = function(x){
          t1 = unlist(unname(x[1]))
          t2 = unlist(unname(x[2]))
          dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
          return(dt)
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
            durationTab = data.frame("indiv.id" = member1,"area.id" = member2,"contact.id" = dyad, "contactDuration" = durlengths, "contactStartTime" = times.start, "contactEndTime" = times.finish, stringsAsFactors = TRUE)
            
          }else{ #If there were no recorded contacts
            durationTab <- data.frame(matrix(ncol = 6, nrow = 0), stringsAsFactors = TRUE)
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
          colnames(distthreshold)<-areaIDVec #relate the distance threshold columns to each area id
          idVecFrame = data.frame(unlist(rep(unlist(unname(x[1])),length(areaIDVec))),areaIDVec, unlist(rep(unlist(unname(x[2])),length(areaIDVec))),unlist(rep(unlist(unname(x[3])),length(areaIDVec))), stringsAsFactors = TRUE)  
          timeVec <- unname(spec.dist[,match("dateTime", colnames(spec.dist))])
          dateTimeFrame <- data.frame(timeVec, stringsAsFactors = TRUE)
          idDurations <- data.frame(data.table::rbindlist(apply(idVecFrame, 1, mat.breaker,distthreshold, timebreakVec, dateTimeFrame)), stringsAsFactors = TRUE)
          return(idDurations)
        }
        
        dist.all<-x[order(x$id, x$dateTime),]
        idVec1 = unique(dist.all$id)
        if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
          dist.all.reduced$dt = 0
        }else{ #if equidistant.time == FALSE
          if(nrow(dist.all) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all is only created from blocks with observed contacts, there will never be a case where dist.all < 1)
            dist.all$dt = 0
          }else{ #if there's more than one row in dist.all
            timesFrame = data.frame(dist.all$dateTime[1:(nrow(dist.all) - 1)], dist.all$dateTime[2:nrow(dist.all)], stringsAsFactors = TRUE)
            timedif = apply(timesFrame, 1, timeDifference)
            dist.all$dt = c(0, timedif)
          }
        }
        comboFrame = data.frame(unique(dist.all$id),dist.threshold,sec.threshold, stringsAsFactors = TRUE)
        duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all)	
        durationTable = data.frame(data.table::rbindlist(duration, fill = TRUE), stringsAsFactors = TRUE)
        
        if(nrow(durationTable) > 0){ #if there was at least one contact duration, block information is appended to the data frame. If there are no observations, durationTable becomes NULL
          durationTable$block<- unique(dist.all$block)
          durationTable$block.start<- unique(dist.all$block.start)
          durationTable$block.end<- unique(dist.all$block.end)
          durationTable$numBlocks<- unique(dist.all$numBlocks)
        }else{ #if nrow(durationTable) == 0
          durationTable <- NULL
        }
        
        return(durationTable)
      }	
      
      thisEnvironment<-environment() #tag the environment of the parent function so that sub-functions may work within it. This is done so that we don't have to clone the list of data frames in the various sub functions over and over again.
      
      if(blocking == TRUE){
        
        if(length(x$block) == 0){ #If there's no "block" column in dist.all output (i.e., blocking was "FALSE" when dist.all was run), then we need to define blocks here.
          
          #the code below comes from the contact::timeBlock.append function. We just refrain from calling it here to prevent the inputs being cloned within the function, and thus save memory.
          
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
          
          daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
          lub.dates = lubridate::date(x$dateTime)
          x<-x[order(lub.dates, daySecondList),] #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data
          
          #for some odd reason, difftime will output mostly zeroes (incorrectly) if there are > 1 correct 0 at the beginning. We use a crude fix here to address this. Basically, we create the zeroes first and combine it with other values afterwards
          totSecond <- rep(0, length(which(x$dateTime == x$dateTime[1])))
          totSecond2<-as.integer(difftime(x$dateTime[(length(totSecond) +1): nrow(x)] ,x$dateTime[1] , units = c("secs")))
          studySecond <- as.integer((c(totSecond, totSecond2) -min(c(totSecond, totSecond2))) + 1)
          
          numblocks <- ceiling((max(studySecond) - 1)/blockLength1)
          block <-rep(0,length(studySecond))
          for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
            block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
          }
          if(length(which(block == 0)) > 0){ #identifies the last block
            block[which(block == 0)] = numblocks
          }
          
          block.start<-as.character(as.POSIXct(x$dateTime[1]) + ((block - 1)*blockLength1)) #identify the timepoint where each block starts (down to the second resolution)
          block.end<-as.character(as.POSIXct(x$dateTime[1]) + ((block - 1)*blockLength1) + (blockLength1 -1)) #identify the timepoint where each block ends (down to the second resolution)
          
          x$block <- block
          x$block.start <- block.start
          x$block.end <- block.end
          x$numBlocks <- max(block) #the contactTest function will require this information (i.e. the number of blocks in the dataset)
          
          rm(list = c("daySecondList", "lub.dates", "totSecond", "totSecond2", "studySecond", "block", "numblocks", "block.start", "block.end")) #remove these objects because they are no longer needed.
          
          blockList<-list()
          blockVec <- unique(x$block)
          
          for(j in 1:length(blockVec)){
            blockL1<- list(x[which(x$block == blockVec[j]),])
            blockList<- c(blockList, blockL1)
          }
          
        }else{ #if length(x$block > 0)
          block = x$block
          blockVec<-unique(block)
          
          blockList<-list()
          
          for(j in 1:length(blockVec)){
            blockL1<- list(x[which(block == blockVec[j]),])
            blockList<- c(blockList, blockL1)
          }
          
        }
        
        rm(x) #remove x to free up local memory
        
        if(parallel == TRUE){#if parallel is true, we run blocks in parallel
          
          cl <- parallel::makeCluster(nCores)
          doParallel::registerDoParallel(cl)
          on.exit(parallel::stopCluster(cl))
          duration.block <-foreach(l = seq(from = 1, to = length(blockList), by = 1)) %dopar% {durFinder.block.List(x = blockList[[l]], dist.threshold, sec.threshold, equidistant.time)} #run the durFinder.block.list function for each list entry
          
        }else{ #if parallel == FALSE
          
          duration.block <-foreach(l = seq(from = 1, to = length(blockList), by = 1)) %do% {durFinder.block.List(x = blockList[[l]], dist.threshold, sec.threshold, equidistant.time)} #run the durFinder.block.list function for each list entry
          
        }
        
        durationTable <- data.frame(data.table::rbindlist(duration.block), stringsAsFactors = TRUE)
        
        if(nrow(durationTable) > 0){ #Here we change components of the durationTable to the appropriate data type. This only occurs if there was at least one observation, however (otherwise an error will be produced). 
          durationTable[,match("block", names(durationTable))]<- as.factor(durationTable[,match("block", names(durationTable))])
          durationTable[,match("block.start", names(durationTable))]<- as.factor(durationTable[,match("block.start", names(durationTable))])
          durationTable[,match("block.end", names(durationTable))]<- as.factor(durationTable[,match("block.end", names(durationTable))])
          durationTable[,match("numBlocks", names(durationTable))]<- as.factor(durationTable[,match("numBlocks", names(durationTable))])
        }
        
      }else{ #If blocking == FALSE
        durationTable <- durFinder.noblock(parallel, dist.threshold, sec.threshold, equidistant.time, nCores, environmentTag = thisEnvironment) #note that if blocking == FALSE, we can still improve speed a bit by parallelizing a later sub function. Further note that x isn't carried over because the function works in the master-function environment
      }
      
      if(nrow(durationTable) > 0){ 
        #Here we change the rest of the components of the durationTable to the appropriate data type.
        durationTable[,match("indiv.id", names(durationTable))]<- as.factor(durationTable[,match("indiv.id", names(durationTable))])
        durationTable[,match("area.id", names(durationTable))]<- as.factor(durationTable[,match("area.id", names(durationTable))])
        durationTable[,match("contact.id", names(durationTable))]<- as.factor(durationTable[,match("contact.id", names(durationTable))])
        durationTable[,match("contactDuration", names(durationTable))]<- as.integer(durationTable[,match("contactDuration", names(durationTable))])
        
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
      }else{ #if nrow(durationTable) == 0
        
        warning("no contacts observed at this dist.threshold value") #if no contacts were found, ensure that users are made aware.
        
      }
      return(durationTable)
      
    }
      
    list.dur<- foreach::foreach(k = 1:length(x), .packages = 'foreach') %do% listBreak_dur.generator(x[[k]],dist.threshold,sec.threshold, blocking, blockUnit, blockLength, equidistant.time, parallel, reportParameters, nCores) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops
    return(list.dur)
  }else{ #if(is.data.frame(x) == TRUE)
    
    durFinder.noblock<-function(parallel, dist.threshold, sec.threshold, equidistant.time, nCores, environmentTag){
      
      eval(expr = { #evaluate these function steps in the master-function environment to avoid cloning x (i.e., dist2... ouput), which can be quite large
        
      #write nested sub functions
      
      timeDifference = function(x){
        t1 = unlist(unname(x[1]))
        t2 = unlist(unname(x[2]))
        dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
        return(dt)
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
          durationTab = data.frame("indiv.id" = member1,"area.id" = member2,"contact.id" = dyad, "contactDuration" = durlengths, "contactStartTime" = times.start, "contactEndTime" = times.finish, stringsAsFactors = TRUE)
          
        }else{ #If there were no recorded contacts
          durationTab <- data.frame(matrix(ncol = 6, nrow = 0), stringsAsFactors = TRUE)
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
        colnames(distthreshold)<-areaIDVec #relate the distance threshold columns to each area id
        idVecFrame = data.frame(unlist(rep(unlist(unname(x[1])),length(areaIDVec))),areaIDVec, unlist(rep(unlist(unname(x[2])),length(areaIDVec))),unlist(rep(unlist(unname(x[3])),length(areaIDVec))), stringsAsFactors = TRUE)  
        timeVec <- unname(spec.dist[,match("dateTime", colnames(spec.dist))])
        dateTimeFrame <- data.frame(timeVec, stringsAsFactors = TRUE)
        idDurations <- data.frame(data.table::rbindlist(apply(idVecFrame, 1, mat.breaker,distthreshold, timebreakVec, dateTimeFrame)), stringsAsFactors = TRUE)
        return(idDurations)
      }
      
      dist.all<-x[order(x$id, x$dateTime),]
      idVec1 = unique(dist.all$id)
      if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
        dist.all.reduced$dt = 0
      }else{ #if equidistant.time == FALSE
        if(nrow(dist.all) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all is only created from blocks with observed contacts, there will never be a case where dist.all < 1)
          dist.all$dt = 0
        }else{ #if there's more than one row in dist.all
          timesFrame = data.frame(dist.all$dateTime[1:(nrow(dist.all) - 1)], dist.all$dateTime[2:nrow(dist.all)], stringsAsFactors = TRUE)
          
          if (parallel == TRUE){
            cl<-parallel::makeCluster(nCores)
            on.exit(parallel::stopCluster(cl))
            timedif<-parallel::parApply(cl, timesFrame, 1, timeDifference); dist.all$dt = c(0, timedif) #timedif represents the time it takes to move from location i-1 to location i	
          }else{
            timedif = apply(timesFrame, 1, timeDifference)
            dist.all$dt = c(0, timedif) 
          }
        }
      }
      comboFrame = data.frame(unique(dist.all$id),dist.threshold,sec.threshold, stringsAsFactors = TRUE)
      
      if (parallel == TRUE){
        duration<-parallel::parApply(cl, comboFrame, 1, contactMatrix.maker,idVec1, dist.all) #note that cl is defined above
      }else{
        duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all)
      }
      durationTable = data.frame(data.table::rbindlist(duration), stringsAsFactors = TRUE)
      
      }, envir = environmentTag)
      
      return(durationTable)
    }
    durFinder.block.List<-function(x,dist.threshold, sec.threshold, equidistant.time){
      
      timeDifference = function(x){
        t1 = unlist(unname(x[1]))
        t2 = unlist(unname(x[2]))
        dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
        return(dt)
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
          durationTab = data.frame("indiv.id" = member1,"area.id" = member2,"contact.id" = dyad, "contactDuration" = durlengths, "contactStartTime" = times.start, "contactEndTime" = times.finish, stringsAsFactors = TRUE)
          
        }else{ #If there were no recorded contacts
          durationTab <- data.frame(matrix(ncol = 6, nrow = 0), stringsAsFactors = TRUE)
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
        colnames(distthreshold)<-areaIDVec #relate the distance threshold columns to each area id
        idVecFrame = data.frame(unlist(rep(unlist(unname(x[1])),length(areaIDVec))),areaIDVec, unlist(rep(unlist(unname(x[2])),length(areaIDVec))),unlist(rep(unlist(unname(x[3])),length(areaIDVec))), stringsAsFactors = TRUE)  
        timeVec <- unname(spec.dist[,match("dateTime", colnames(spec.dist))])
        dateTimeFrame <- data.frame(timeVec, stringsAsFactors = TRUE)
        idDurations <- data.frame(data.table::rbindlist(apply(idVecFrame, 1, mat.breaker,distthreshold, timebreakVec, dateTimeFrame)), stringsAsFactors = TRUE)
        return(idDurations)
      }
      
      dist.all<-x[order(x$id, x$dateTime),]
      idVec1 = unique(dist.all$id)
      if(equidistant.time == TRUE){ #added 02/04/2019 to make the dt calculations a toggleable parameter that users may turn off if all data points in their data set are temporally equidistant. This saves a large amount of time (approx. 3.5 mins/day)
        dist.all.reduced$dt = 0
      }else{ #if equidistant.time == FALSE
        if(nrow(dist.all) ==1){ #if there's only one row, there cannot be any time difference (note, because dist.all is only created from blocks with observed contacts, there will never be a case where dist.all < 1)
          dist.all$dt = 0
        }else{ #if there's more than one row in dist.all
          timesFrame = data.frame(dist.all$dateTime[1:(nrow(dist.all) - 1)], dist.all$dateTime[2:nrow(dist.all)], stringsAsFactors = TRUE)
          timedif = apply(timesFrame, 1, timeDifference)
          dist.all$dt = c(0, timedif)
        }
      }
      comboFrame = data.frame(unique(dist.all$id),dist.threshold,sec.threshold, stringsAsFactors = TRUE)
      duration = apply(comboFrame, 1, contactMatrix.maker,idVec1, dist.all)	
      durationTable = data.frame(data.table::rbindlist(duration, fill = TRUE), stringsAsFactors = TRUE)
      
      if(nrow(durationTable) > 0){ #if there was at least one contact duration, block information is appended to the data frame. If there are no observations, durationTable becomes NULL
        durationTable$block<- unique(dist.all$block)
        durationTable$block.start<- unique(dist.all$block.start)
        durationTable$block.end<- unique(dist.all$block.end)
        durationTable$numBlocks<- unique(dist.all$numBlocks)
      }else{ #if nrow(durationTable) == 0
        durationTable <- NULL
      }
  
      return(durationTable)
    }	
    
    thisEnvironment<-environment() #tag the environment of the parent function so that sub-functions may work within it. This is done so that we don't have to clone the list of data frames in the various sub functions over and over again.
    
    if(blocking == TRUE){
      
      if(length(x$block) == 0){ #If there's no "block" column in dist.all output (i.e., blocking was "FALSE" when dist.all was run), then we need to define blocks here.
        
        #the code below comes from the contact::timeBlock.append function. We just refrain from calling it here to prevent the inputs being cloned within the function, and thus save memory.
        
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
        
        daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
        lub.dates = lubridate::date(x$dateTime)
        x<-x[order(lub.dates, daySecondList),] #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data
        
        #for some odd reason, difftime will output mostly zeroes (incorrectly) if there are > 1 correct 0 at the beginning. We use a crude fix here to address this. Basically, we create the zeroes first and combine it with other values afterwards
        totSecond <- rep(0, length(which(x$dateTime == x$dateTime[1])))
        totSecond2<-as.integer(difftime(x$dateTime[(length(totSecond) +1): nrow(x)] ,x$dateTime[1] , units = c("secs")))
        studySecond <- as.integer((c(totSecond, totSecond2) -min(c(totSecond, totSecond2))) + 1)
        
        numblocks <- ceiling((max(studySecond) - 1)/blockLength1)
        block <-rep(0,length(studySecond))
        for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
          block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
        }
        if(length(which(block == 0)) > 0){ #identifies the last block
          block[which(block == 0)] = numblocks
        }
        
        block.start<-as.character(as.POSIXct(x$dateTime[1]) + ((block - 1)*blockLength1)) #identify the timepoint where each block starts (down to the second resolution)
        block.end<-as.character(as.POSIXct(x$dateTime[1]) + ((block - 1)*blockLength1) + (blockLength1 -1)) #identify the timepoint where each block ends (down to the second resolution)
        
        x$block <- block
        x$block.start <- block.start
        x$block.end <- block.end
        x$numBlocks <- max(block) #the contactTest function will require this information (i.e. the number of blocks in the dataset)
        
        rm(list = c("daySecondList", "lub.dates", "totSecond", "totSecond2", "studySecond", "block", "numblocks", "block.start", "block.end")) #remove these objects because they are no longer needed.
        
        blockList<-list()
        blockVec <- unique(x$block)
        
        for(j in 1:length(blockVec)){
          blockL1<- list(x[which(x$block == blockVec[j]),])
          blockList<- c(blockList, blockL1)
        }
        
      }else{ #if length(x$block > 0)
        block = x$block
        blockVec<-unique(block)
        
        blockList<-list()
        
        for(j in 1:length(blockVec)){
          blockL1<- list(x[which(block == blockVec[j]),])
          blockList<- c(blockList, blockL1)
        }
        
      }
      
      rm(x) #remove x to free up local memory
      
      if(parallel == TRUE){#if parallel is true, we run blocks in parallel
        
        cl <- parallel::makeCluster(nCores)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))
        duration.block <-foreach(l = seq(from = 1, to = length(blockList), by = 1)) %dopar% {durFinder.block.List(x = blockList[[l]], dist.threshold, sec.threshold, equidistant.time)} #run the durFinder.block.list function for each list entry
        
      }else{ #if parallel == FALSE
        
        duration.block <-foreach(l = seq(from = 1, to = length(blockList), by = 1)) %do% {durFinder.block.List(x = blockList[[l]], dist.threshold, sec.threshold, equidistant.time)} #run the durFinder.block.list function for each list entry
        
      }
      
      durationTable <- data.frame(data.table::rbindlist(duration.block), stringsAsFactors = TRUE)
      
      if(nrow(durationTable) > 0){ #Here we change components of the durationTable to the appropriate data type. This only occurs if there was at least one observation, however (otherwise an error will be produced). 
        durationTable[,match("block", names(durationTable))]<- as.factor(durationTable[,match("block", names(durationTable))])
        durationTable[,match("block.start", names(durationTable))]<- as.factor(durationTable[,match("block.start", names(durationTable))])
        durationTable[,match("block.end", names(durationTable))]<- as.factor(durationTable[,match("block.end", names(durationTable))])
        durationTable[,match("numBlocks", names(durationTable))]<- as.factor(durationTable[,match("numBlocks", names(durationTable))])
      }
      
    }else{ #If blocking == FALSE
      durationTable <- durFinder.noblock(parallel, dist.threshold, sec.threshold, equidistant.time, nCores, environmentTag = thisEnvironment) #note that if blocking == FALSE, we can still improve speed a bit by parallelizing a later sub function. Further note that x isn't carried over because the function works in the master-function environment
    }
    
    if(nrow(durationTable) > 0){ 
      #Here we change the rest of the components of the durationTable to the appropriate data type.
      durationTable[,match("indiv.id", names(durationTable))]<- as.factor(durationTable[,match("indiv.id", names(durationTable))])
      durationTable[,match("area.id", names(durationTable))]<- as.factor(durationTable[,match("area.id", names(durationTable))])
      durationTable[,match("contact.id", names(durationTable))]<- as.factor(durationTable[,match("contact.id", names(durationTable))])
      durationTable[,match("contactDuration", names(durationTable))]<- as.integer(durationTable[,match("contactDuration", names(durationTable))])
    
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
    }else{ #if nrow(durationTable) == 0
      
      warning("no contacts observed at this dist.threshold value") #if no contacts were found, ensure that users are made aware.
      
    }
    return(durationTable)
  } 
} 
