#' Smooth Point-Locations Over Time
#'
#' Aggregate location data by secondAgg seconds over the course of each day 
#'    represented in the dataset. The function smooths xy data forwards 
#'    (smooth.type == 1) or backwards (smooth.type == 2) according to a 
#'    data-point-averaging loess smoothing methodology. As part of the 
#'    smoothing process, tempAggregate fills in any missing values (either due 
#'    to a lack of data transmission or faulty prior interpolation). We 
#'    recognize that this procedure is not sensitive to individual presence at 
#'    given timesteps (e.g., some individuals may be missing on certain days, 
#'    hours, etc., and therefore may produce inaccurate location aggregates if 
#'    days/hours exist where individuals are not present in the dataset (e.g., 
#'    they were purposefully removed, or moved outside of the monitoring 
#'    area)). To increase accuracy, package users may specify a resolutionLevel
#'    ("full" or "reduced") to process individuals' locations at different 
#'    resolutions. If resolution == "reduced", if no locations of individuals 
#'    exist over any secondAgg time block, NAs will be produced for the time 
#'    blocks of interest. 
#'    
#' This function is based on real-time-location-data-smoothing methods 
#'    presented by Dawson et al. 2019. 
#' @param x Data frame or list of data frames containing real-time-location 
#'    data.  
#' @param id Vector of length nrow(data.frame(x)) or singular character data, 
#'    detailing the relevant colname in x, that denotes what unique ids for 
#'    tracked individuals will be used. If argument == NULL, the function 
#'    assumes a column with the colname "id" exists in x. Defaults to NULL.
#' @param point.x Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what planar-x or 
#'    longitude coordinate information will be used. If argument == NULL, the 
#'    function assumes a column with the colname "x" exists in x. Defaults to 
#'    NULL.
#' @param point.y Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what planar-y or 
#'    lattitude coordinate information will be used. If argument == NULL, the 
#'    function assumes a column with the colname "y" exists in x. Defaults to 
#'    NULL.
#' @param dateTime Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what dateTime 
#'    information will be used. If argument == NULL, the function assumes a 
#'    column with the colname "dateTime" exists in x. Defaults to NULL.
#' @param secondAgg Numerical. The number of seconds over which 
#'    tracked-individuals' location will be averaged. Defaults to 10.
#' @param extrapolate.left Logical. If TRUE, individuals position at time 
#'    points prior to their first location fix will revert to their first 
#'    recorded location. If FALSE, NAs will be placed at these time points in 
#'    individuals' movement paths. Defaults to FALSE.
#' @param extrapolate.right Logical. If TRUE, individuals position at time 
#'    points following their last location fix will revert to their final 
#'    recorded location. If FALSE, NAs will be placed at these time points in 
#'    individuals' movement paths. Defaults to FALSE.
#' @param resolutionLevel Character string taking the value of "full" or 
#'    "reduced."  If "full," if no known locations of individuals exist over 
#'    any secondAgg time block, xy-coordinates revert to the last-known values 
#'    for that individual. If "reduced," if no known locations of individuals 
#'    exist over any secondAgg time block, NAs will be produced for the time 
#'    blocks of interest. Defaults to "full."
#' @param parallel Logical. If TRUE, sub-functions within the tempAggregate 
#'    wrapper will be parallelized. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
#' @param na.rm Logical. If TRUE, all unknown locations (i.e., xy-coordinate 
#'    pairs reported as NAs) will be removed from the output. Defaults to TRUE.
#'    Note that if na.rm == FALSE, all aggregated location fixes will be 
#'    temporally equidistant.
#' @param smooth.type Numerical, taking the values 1 or 2. Indicates the type 
#'    of smooting used to average individuals' xy-coordinates. If 
#'    smooth.type == 1, data are smoothed forwards. If smooth.type == 2, data 
#'    are smoothed backwards. Defaults to 1.
#' @references Dawson, D.E., Farthing, T.S., Sanderson, M.W., and Lanzas, C. 
#'    2019. Transmission on empirical dynamic contact networks is influenced by
#'    data processing decisions. Epidemics 26:32-42. 
#'    https://doi.org/10.1016/j.epidem.2018.08.003/
#' @keywords data-processing smoothing location point
#' @return Returns a data frame (or list of data frames if \code{x} is a 
#'    list of data frames) with the following columns:
#'    
#'    \item{id}{The unique ID of tracked individuals.}
#'    \item{x}{Smoothed x coordinates.}
#'    \item{y}{Smoothed y coordinates.}
#'    \item{dateTime}{Timepoint at which smoothed points were observed.}
#'
#' @import foreach        
#' @export
#' @examples
#' data("calves")
#' head(calves) #observe that fix intervals occur ever 4-5 seconds.
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'    time = calves$time) #add dateTime identifiers for location fixes.
#'    
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth to 5-min fix intervals.
#'    

tempAggregate <- function(x = NULL, id = NULL, point.x = NULL, point.y = NULL, dateTime = NULL, secondAgg = 10, extrapolate.left = FALSE, extrapolate.right = FALSE, resolutionLevel = "full", parallel = FALSE, nCores = (parallel::detectCores()/2), na.rm = TRUE, smooth.type = 1) { #removed totalSecond = NULL argument on 01102019

  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  l <- NULL
  
  #write sub-function
  Agg.generator<-function(x, id, point.x, point.y, dateTime, secondAgg, extrapolate.left, extrapolate.right, resolutionLevel, parallel, na.rm, smooth.type, nCores){
    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
      originTab = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime, stringsAsFactors = TRUE)
    }

    if(length(x) > 0){ #for some reason using an "else" statement would always result in an originTab table with 0 records...
      if(length(id) > 0){
        if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than id being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
          x$id <- x[,match(id, names(x))]
        }else{ #if length(id) > 1
          x$id = id
        }
      }
      idVec1 <- x$id

      if(length(point.x) > 0){
        if(length(point.x) == 1 && is.na(match(point.x[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$x <- x[,match(point.x, names(x))]
        }else{ #if length(point.x) > 1
          x$x = point.x
        }
      }
      if(length(point.y) > 0){
        if(length(point.y) == 1 && is.na(match(point.y[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$y <- x[,match(point.y, names(x))]
        }else{ #if length(point.y) > 1
          x$y = point.y
        }
      }
      xyFrame1<- data.frame(x = x$x, y = x$y, stringsAsFactors = TRUE)

      if(length(dateTime) > 0){

        if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$dateTime <- x[,match(dateTime, names(x))]
        }else{ #if length(dateTime) > 1
          x$dateTime = dateTime
        }
      }
      dateTimeVec<-x$dateTime
      bindlist1<-list(idVec1, xyFrame1, dateTimeVec)
      originTab <- do.call("cbind", bindlist1)
      names(originTab)[c(1,ncol(originTab))]<-c("id", "dateTime")
    }
    originTab$date = lubridate::date(x$dateTime)

    ###The following lines creating the totalSecond column in originTab (derived from datetime.append) were added to version 01102019 to remove the need for including the dayID and  totalSecond columns initially, and to fix there error where the total seconds at breakpoints (n) exceeds the number of seconds in the dataset (coord.tmp) (Note: the specific error produced was "Error in coord.tmp[n, 1] : subscript out of bounds")
    originTab<-originTab[order(originTab$dateTime),] #Just in case the data wasn't already ordered in this way.
    timevec <- originTab$dateTime

    daySecondVec = lubridate::hour(timevec) * 3600 + lubridate::minute(timevec) * 60 + lubridate::second(timevec) #This calculates a day-second
    dates = unique(originTab$date)
    dayIDVec = NULL
    dayIDseq = seq(1,(length(dates)),1)
    for(b in dayIDseq){
      dayID = rep(b,length(which(originTab$date == dates[b])))
      dayIDVec = c(dayIDVec, dayID)
    }

    originTab$totalSecond <- ((dayIDVec - min(dayIDVec))*86400) + daySecondVec + 1 #This calculates the total second (the cumulative second across the span of the study's timeframe) # We add the "+ 1" because lubridate::second goes from 0-59, but for our purposes, we need 1-60 b/c these these seconds will relate to rownumbers (n) later on, and cannot take a 0 value.
    leftExtrap = extrapolate.left
    rightExtrap = extrapolate.right
    
    if(resolutionLevel == "r" || resolutionLevel == "re" || resolutionLevel == "red" || resolutionLevel == "redu" || resolutionLevel == "reduc" || resolutionLevel == "reduce" || resolutionLevel == "REDUCED" || resolutionLevel == "R" || resolutionLevel == "RE" || resolutionLevel == "RED" || resolutionLevel == "REDU" || resolutionLevel == "REDUC" || resolutionLevel == "REDUCE"){ #if users input anything into the resolutionLevel argument that even resembles "reduced," the function will interpret that as what they were intending. Any other value will result in the function keeping resolutionLevel as "full."
      resolutionLevel = "reduced"
    }

    Processing_and_Aggregation_Procedures<-function(brk.point, originTab, extrapolate.left = leftExtrap, extrapolate.right = rightExtrap, resolutionLevel, smooth.type, secondAgg){

      Nseconds.aggregate=function(coord, secondAgg) {
        m=secondAgg
        tmp1<-coord[,1]
        tmp2<-coord[,2]
        tmp1.mx<-matrix(tmp1,nrow=m) # aggregate to m seconds
        tmp2.mx<-matrix(tmp2,nrow=m)
        tmp1.mmin<-apply(tmp1.mx,2,mean,na.rm = TRUE)
        tmp2.mmin<-apply(tmp2.mx,2,mean,na.rm = TRUE)
        tmp.mmin<-cbind(tmp1.mmin, tmp2.mmin)
        return(tmp.mmin)
      }

      coord.tmp<-matrix(rep(0,length(unique(lubridate::date(originTab$dateTime)))*24*60*60*2),ncol=2)

      data.one=originTab[c(brk.point[1]:brk.point[2]),]
      data.len=nrow(data.one)
      xcol<-match("x", names(data.one))
      ycol<-match("y", names(data.one))

      tmp=as.numeric(data.one[1,c(xcol,ycol)])

      n=data.one[1,match("totalSecond", names(data.one))] # get current position at 1st data point (previous point)

      if(extrapolate.left == TRUE){ #This section extrapolates the left tails of individuals' paths by filling in the missing time between the first recorded second and the first second of the day. It assumes that the individual starts the first second of the day at the same place it spent the first recorded second.
        coord.tmp[1:n,1]=tmp[1]
        coord.tmp[1:n,2]=tmp[2]

      }else{ #if extrapolate.left == FALSE

        if(n > 1){ #Note that n cannot be smaller than "1," as "totalSecond" calculated above cannot be smaller than 1.
          coord.tmp[1:(n-1),1]=NA
          coord.tmp[1:(n-1),2]=NA
        }
        coord.tmp[n,1]=tmp[1] #This is the timepoint where individuals first appear in the dataset
        coord.tmp[n,2]=tmp[2]
      }

      #Interpolate missing values: The point of the code below is to fill in any missing values (either due to a lack of data transmition or faulty prior interpolation). This code works by setting n initially to the first value of totalSecond's column of data.one. Using a loop, it sets the first value of n.current to the ith row(starting at 2) of "totalSecond". A temporal object "tmp" pulls the xy coordinates from the the ith row of data.one. Then, it sets the xy coordinates of coord.tmp, starting at n + 1 position to the n.current position to values of tmp. This works because n moves up by 1 second, while n.current(and data.len) moves at whatever interval spacing the data is in in data.one (in seconds). It inches along, filling in the slots small or large until it reaches the last data point.

      if(data.len > 1){
        for (i in 2:data.len) { #this grabs the second position 1:length of the subset; thats because it moves from the first second to last on in set looking for holes
          n.current=data.one[i,match("totalSecond", names(data.one))] #this grabs the next second from the subset
          if(smooth.type == 1){ #forwards smoothing
            if((n.current - n) > secondAgg & resolutionLevel == "reduced"){ #So if there are many secondAgg or greater seconds with no data, and resolutionLevel is set to "reduced," the function will set these coordinates to NA
              tmp1 = c(NA,NA)
            }else{ #if there  are fewer than secondAgg with no data, or resolutionLevel is "full," then the points will be smoothed normally.
              tmp1=as.numeric(data.one[(i-1),xcol:ycol]) #this pulls the value from the preceding second to tmp1
            }
            tmp2=as.numeric(data.one[i,xcol:ycol]) #this pulls the value from the current second to tmp2; #Because tmp2 represents a known data point, it would not be changed to NA
            if((n.current - n) > 1){ #This is needed here, because if the difference between n.current and n == 1, there is nothing to smooth (smoothing forward in these cases would overwrite actually-observed xy values)
              coord.tmp[(n+1):(n.current - 1),1]=tmp1[1] #this takes the value of tmp, and assigns it to all the slots beween n and n.current -1
              coord.tmp[(n+1):(n.current - 1),2]=tmp1[2]
            }
            coord.tmp[n.current,1]=tmp2[1] #this takes the value of tmp, and assigns it to all the n.current slot
            coord.tmp[n.current,2]=tmp2[2]
          }
          if(smooth.type == 2){ #backwards smoothing
            tmp=as.numeric(data.one[i,xcol:ycol]) #this pulls the value from that second to tmp
            if((n.current - n) > secondAgg & resolutionLevel == "reduced"){ #So if there are many secondAgg or greater seconds with no data, and resolutionLevel is set to "reduced," the function will set these coordinates to NA
              coord.tmp[(n+1):(n.current -1),1]=NA #this fills all the slots beween n and n.current with NAs
              coord.tmp[(n+1):(n.current -1),2]=NA
              coord.tmp[n.current,1]=tmp[1] #this takes the value of tmp, and assigns it to the n.current slot (because n.current will always be a known point.)
              coord.tmp[n.current,2]=tmp[2]
            }else{ #if there  are fewer than secondAgg with no data, or resolutionLevel is "full," then the points will be smoothed normally.
              coord.tmp[(n+1):n.current,1]=tmp[1] #this takes the value of tmp, and assigns it to all the slots beween n and n.current
              coord.tmp[(n+1):n.current,2]=tmp[2]
              }

          }
          n=n.current # assign n with most current position
        }
      }

      if(nrow(coord.tmp) > max(data.one$totalSecond)){ #This step will only proceed if there are empty spaces at the end of coord.tmp (i.e., the last rownumber does not exceed the maximum potential number of second in the dataset)
        if(extrapolate.right == TRUE){#Extrapolate Right Tail: this simply fills in the rest of the values using the location of the last datavalue
          tmp=as.numeric(data.one[length(data.one[,xcol]), c(xcol,ycol)])
          coord.tmp[(n+1):nrow(coord.tmp),1]=tmp[1]
          coord.tmp[(n+1):nrow(coord.tmp),2]=tmp[2]
          
          #Added 05/21/2019. For some reason, there was a bug that prevented the last row of coord.tmp from taking new values using the code above.... I'm not quite sure why that is, but to fix it, I've added the code below -tsf.
          coord.tmp[nrow(coord.tmp),1] <- tmp[1]
          coord.tmp[nrow(coord.tmp),2] <- tmp[2]

        }else{ #if extrapolate.right == FALSE

          coord.tmp[(n+1):nrow(coord.tmp),1]=NA #This reports NAs following the last instance of observed movement
          coord.tmp[(n+1):nrow(coord.tmp),2]=NA

          #Added 05/21/2019. For some reason, there was a bug that prevented the last row of coord.tmp from taking new values using the code above.... I'm not quite sure why that is, but to fix it, I've added the code below -tsf.
          coord.tmp[nrow(coord.tmp),1] <- NA
          coord.tmp[nrow(coord.tmp),2] <- NA
          
        }
      }

      loc.aggregate = Nseconds.aggregate(coord.tmp, brk.point[4])
      return(loc.aggregate)
    }

    originTab<-originTab[order(originTab$id, originTab$dateTime),] #This sorts the dataset by individuals' IDs and timestep
    rownames(originTab) <-seq(1,nrow(originTab),1) #This is necessary for the brk.points to be accurate.

    locmatrix <- NULL
    start.brk <- NULL
    locTable <- NULL
    dateSeq <- dates #just rename for the same of convenience
    indivSeq <- unique(originTab$id)

    for(i in indivSeq){ #This loop determines the rows in the dataset where each individual's path begins.
      start.brk = c(start.brk, min(which(originTab$id == i)))
    }

    if(length(start.brk) > 1){end.brk=c(start.brk[2:length(indivSeq)]-1, nrow(originTab))}
    if(length(start.brk) == 1){end.brk=nrow(originTab)}
    brk.point<-cbind(start.brk,end.brk,seq(1, length(end.brk),1), secondAgg) #these two vectors offset beginning and endpoints for obsercations from different individuals.

    if (parallel == TRUE){
      cl<-parallel::makeCluster(nCores)
      on.exit(parallel::stopCluster(cl))
      locmatrix<-parallel::parApply(cl, brk.point, 1, Processing_and_Aggregation_Procedures, originTab, leftExtrap, rightExtrap, resolutionLevel, smooth.type, secondAgg) ##Note, this procedure produces a table that in which the first half is the x axis, and the second is the y axis. You have to merge (see code below)
    }else{locmatrix<-apply(brk.point, 1, Processing_and_Aggregation_Procedures, originTab, leftExtrap, rightExtrap, resolutionLevel, smooth.type, secondAgg)}

    #The code below breaks the table into the x and y coordinates that we'd expect here.
    for (m in 1:length(indivSeq)){
      xvec<-locmatrix[1:(length(locmatrix[,1])/2),m]
      yvec<-locmatrix[((length(locmatrix[,1])/2)+1):length(locmatrix[,1]),m]
      id1 = rep(indivSeq[m],length(xvec))

      datesAgg = NULL

      for(n in 1:length(dateSeq)){
        date = data.frame(date = rep(dateSeq[n],(nrow(locmatrix)/2)/length(dateSeq)), stringsAsFactors = TRUE)
        datesAgg = data.frame(data.table::rbindlist(list(datesAgg,date)), stringsAsFactors = TRUE)

      }
      tempTable = data.frame(id = id1, x = xvec, y = yvec, date = datesAgg, stringsAsFactors = TRUE)
      tempTable1 = NULL

      for(o in 1:length(dateSeq)){
        datesub = subset(tempTable, date == unique(tempTable$date)[o])
        datesub$year = lubridate::year(datesub$date)
        datesub$month = lubridate::month(datesub$date)
        datesub$day = lubridate::day(datesub$date)
        rownames(datesub) = seq(1,nrow(datesub),1)
        datesub$UID = as.numeric(rownames(datesub))

        H = floor((datesub$UID - 1)/(60*(60/secondAgg)))
        M = floor((datesub$UID - 1)/(60/secondAgg)) - (H*60)
        S = round((datesub$UID/(60/secondAgg) - floor(datesub$UID/(60/secondAgg)))*60)
        S = c(0,S[1:(length(S) - 1)]) #I can't seem to get a zero-second starting point without this step
        Hour = ifelse(H < 10, paste(0,H, sep = ""), H)
        Minute = ifelse(M < 10, paste(0,M, sep = ""), M)
        Second = ifelse(S < 10, paste(0,S, sep = ""), S)

        datesub$time = paste(Hour,":", Minute, ":", Second, sep ="")
        date.time = paste(datesub[,match("date",names(datesub))], datesub[,match("time",names(datesub))], sep = " ")
        timevec = lubridate::ymd_hms(date.time)
        datesub$dateTime = timevec
        datesub$hour = H
        datesub$minute = M
        datesub$second = S

        tempTable1 = data.frame(data.table::rbindlist(list(tempTable1,datesub)), stringsAsFactors = TRUE)
      }
      locTable = data.frame(data.table::rbindlist(list(locTable,tempTable1)), stringsAsFactors = TRUE)
    }

    locTable <- locTable[,-match("UID",names(datesub))] #removes the UID column to save space

    id.finder<-function(id, idSeq.subset, data.sub){
      if(isTRUE(id%in%idSeq.subset) == FALSE){ #This looks to see if individuals present over the course of the dataset were present in a given time frame (i.e., had at least one recorded location).
        removeVec <- which(data.sub$id == id)
        return(removeVec)
      }
    }

    if(resolutionLevel == "reduced"){ #if resolutionLevel is "reduced," there may be NaN values in the dataset (if secondAgg timeblock consisting only of NAs was averaged in the Nseconds.Aggregate function). For later processing, we want them to be NAs
      NaNVec<-which(is.nan(locTable$x) == TRUE)
      if(length(NaNVec)> 0){
        locTable$x[NaNVec] <-NA
      }
      
    }
    
    if(na.rm == TRUE){
      navec <- which(is.na(locTable$x) == TRUE)
      if (length(navec) >= 1){
        locTable <- locTable[-navec,] #removes all observations related to specific individuals individual where an NA location was reported (i.e., where the function determined they were absent from the data at a specific resolution).
      }
    }

    locTable<-locTable[order(locTable$id, locTable$dateTime),]
    rownames(locTable)<-seq(1,nrow(locTable),1)

    locTable <- locTable[,-c(match("date",names(locTable)), match("year",names(locTable)), match("month",names(locTable)), match("day",names(locTable)), match("time",names(locTable)), match("hour",names(locTable)), match("minute",names(locTable)), match("second",names(locTable)))] #removes all the strictly unnecessary columns to reduce file size. All these columns can be later derrived again from the dateTime column.

    if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #if id is representative of a column name, this maintains that column name in the output file
      names(locTable)[match("id",names(locTable))] <- id
    }
    if(length(point.x) == 1 && is.na(match(point.x[1], names(x))) == FALSE){ #if point.x is representative of a column name, this maintains that column name in the output file
      names(locTable)[match("x",names(locTable))] <- point.x
    }
    if(length(point.y) == 1 && is.na(match(point.y[1], names(x))) == FALSE){ #if point.y is representative of a column name, this maintains that column name in the output file
      names(locTable)[match("y",names(locTable))] <- point.y
    }
    if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #if dateTime is representative of a column name, this maintains that column name in the output file
      names(locTable)[match("dateTime",names(locTable))] <- dateTime
    }

    return(locTable)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R treats dataframes as lists.
    list.Agg <- foreach::foreach(l = seq(from = 1, to = length(x), by = 1)) %do% Agg.generator(x[[l]], id, point.x, point.y, dateTime, secondAgg, extrapolate.left, extrapolate.right, resolutionLevel, parallel, na.rm, smooth.type, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here.
    return(list.Agg)

  }else{ #if x is a dataFrame
    frame.Agg<- Agg.generator(x, id, point.x, point.y, dateTime, secondAgg, extrapolate.left, extrapolate.right, resolutionLevel, parallel, na.rm, smooth.type, nCores)
    return(frame.Agg)
  }

}
