#' Identify and Remove Data Points Based on Observed Movement Speed
#'
#' mps (a.k.a. Meters-per-Second Filter) identifies and removes timepoints when
#'    tracked individuals were observed moving faster than a set distance 
#'    threshold (representing either the great-circle distance between two 
#'    points a planar distance metric, depending on whether or not lonlat == 
#'    TRUE or FALSE, respectively) per second. (i.e., if it is 
#'    impossible/highly unlikely that individuals moved faster than a given 
#'    speed (mps), we can assume that any instances when they were observed 
#'    doing so were the result of erroneous reporting, and should be removed). 
#'    When running the mps filter, users have the option of setting 
#'    lonlat == TRUE (by default lonlat == FALSE). lonlat is a logical 
#'    argument that tells the function to calculate the distance between 
#'    points on the WGS ellipsoid (if lonlat == TRUE), or on a plane 
#'    (lonlat == FALSE) (see raster::pointDistance). If lonlat == TRUE, 
#'    coordinates should be in degrees. Otherwise, coordinates should represent
#'    planar ('Euclidean') space (e.g. units of meters).
#'
#' If users are not actually interested in filtering datasets, but rather 
#'    determining what observations should be filtered, they may set 
#'    filterOutput == FALSE. By doing so, this function will append up an "mps"
#'    column to the dataset, which reports the avg distance per second 
#'    individuals moved to get from observation i-1 to observation i.
#' @param x List or data frame containing real-time location data that will be 
#'    filtered.
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
#' @param mpsThreshold Numerical. Distance (in meters) representing the maximum
#'    distance individuals can realistically travel over a single second.
#' @param lonlat Logical. If TRUE, point.x and point.y contain geographic 
#'    coordinates (i.e., longitude and lattitude). If FALSE, point.x and 
#'    point.y contain planar coordinates. Defaults to FALSE.
#' @param parallel Logical. If TRUE, sub-functions within the mps wrapper will 
#'    be parallelized. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
#' @param filterOutput Logical. If TRUE, output will be a data frame or list of
#'    data frames (depending on whether or not x is a data frame or not) 
#'    containing only points that adhere to the mpsThreshold rule. If FALSE, no
#'    observartions are removed and an "mps" column is appended to x,which 
#'    reports the avg distance per second individuals moved to get from 
#'    observation i-1 to observation i. Defaults to TRUE.
#' @keywords filter
#' @return If filterOutput == TRUE, returns \code{x} less observations 
#'    representing impossible/unlikely movements.
#'    
#'    If filterOutput == FALSE, returns \code{x} appended with an 
#'    "mps" column which reports the avg distance per second 
#'    individuals moved to get from observation i-1 to observation i.
#' @import foreach
#' @export
#' @examples
#' data(calves) #load calves data
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date,
#'    time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#' 
#' calves_filter1 <- mps(calves.dateTime, id = calves.dateTime$calftag,
#'    point.x = calves.dateTime$x, point.y = calves.dateTime$y, 
#'    dateTime = calves.dateTime$dateTime, mpsThreshold = 10, lonlat = FALSE, parallel = FALSE, 
#'    filterOutput = TRUE) 
#'
mps <- function(x, id = NULL, point.x = NULL, point.y = NULL, dateTime = NULL, mpsThreshold = 10, lonlat = FALSE, parallel = FALSE, nCores = (parallel::detectCores()/2), filterOutput = TRUE){

  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  l <- NULL
  
  #write the sub-functions
  filter3.func<-function(x, id, point.x, point.y, dateTime, mpsThreshold, lonlat, parallel, filterOutput, nCores){
    idVec <- NULL
    xVec <- NULL
    yVec <- NULL
    dateTimeVec <- NULL

    if(length(id) > 0){
      if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #added 1/21 Rather than id being a vector of length(nrow(x)), it may be more convenient to designate the colname for intended "id" values
        idVec <- x[,match(id, names(x))]
      }else{ #if length(id) > 1
        idVec <- id
      }
    }

    if(length(point.x) > 0){
      if(length(point.x) == 1 && is.na(match(point.x[1], names(x))) == FALSE){
        xVec <- x[,match(point.x, names(x))]
      }else{ #if length(point.x) > 1
        xVec <- point.x
      }
    }
    if(length(point.y) > 0){
      if(length(point.y) == 1 && is.na(match(point.y[1], names(x))) == FALSE){
        yVec <- x[,match(point.y, names(x))]
      }else{ #if length(point.y) > 1
        yVec <- point.y
      }
    }

    if(length(dateTime) > 0){

      if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){
        dateTimeVec <- x[,match(dateTime, names(x))]
      }else{ #if length(dateTime) > 1
        dateTimeVec <- dateTime
      }
    }

    if(length(xVec) > 0){ #if xVec is still NULL, then the function assumes an "x" column exists in x and assigns that column value to xVec1. In this case, if no "x" column exists, an error will be returned (Note: The error will look like this: Error in ans[!xVec1 & ok] <- rep(no, length.out = length(ans))[!xVec1 &  : replacement has length zero In addition: Warning message:In rep(no, length.out = length(ans)) : 'x' is NULL so the result will be NULL).
      xVec1 <- xVec
    }else{
      xVec1 <- x$x
    }

    if(length(yVec) > 0){
      yVec1 <- yVec
    }else{
      yVec1 <- x$y
    }

    if(length(idVec) > 0){
      idVec1 <- idVec
    }else{
      idVec1 <- x$id
    }

    if(length(dateTimeVec) > 0){
      dateTimeVec1 <- dateTimeVec
    }else{
      dateTimeVec1 <- x$dateTime
    }
    x$idVec1 <- idVec1
    x$dateTimeVec1 <- dateTimeVec1
    x$xVec1 <- xVec1
    x$yVec1 <- yVec1

    if(length(x$totalSecond) > 0){ #if there is a totalSecond column (indicating the data has likely been processed using dateTime.append), then the totSecVec will be the totalSecond-column values, otherwise we need to calculate the totalSecond for each dateTime (using the same process as in dateTime.append)
      x<-x[order(x$idVec1, x$dateTimeVec1),] #this sorts the data so that future processes will work
      totSecVec<-x$totalSecond
    }else{
      x<-x[order(x$dateTimeVec1),] #We need to ensure that the data is ordered in the same way as timeVec will be, so that all observations line up as they should.
      timeVec<-x$dateTimeVec1
      daySecondVec <- lubridate::hour(timeVec) * 3600 + lubridate::minute(timeVec) * 60 + lubridate::second(timeVec) #This calculates a day-second
      lub.dates = lubridate::date(timeVec)
      dateseq = unique(lub.dates)
      dayIDVec = NULL
      dayIDseq = seq(1,length(dateseq),1)
      for(b in dayIDseq){
        ID = rep(b,length(which(lub.dates == dateseq[which(dayIDseq == b)])))
        dayIDVec = c(dayIDVec, ID)
      } #This part of the function takes awhile (especially for large datasets)
      totSecVec <- ((dayIDVec - 1)*86400) + daySecondVec #This calculates the total second (the cumulative second across the span of the study's timeframe)
      id_dateTime.order<-order(x$idVec1, x$dateTimeVec1) #Identify the necessary order so that the data so that future processes will work
      x<-x[id_dateTime.order,] #this sorts the data so that future processes will work
      totSecVec<-totSecVec[id_dateTime.order] #We need to make sure these observations have been ordered in the same way as  by idVec1 and dateTimeVec1
    }
    euc=function(x, dist.measurement) {
      point1 = x.cor=unlist(c(x[1],x[2]))
      point2 = x.cor=unlist(c(x[3],x[4]))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = dist.measurement)
      return(euc.dis)
    }

    indivDist <- function(x, y, dist.measurement){
      xytab = y[which(y$idVec1 == x[1]),]
      if(nrow(xytab) > 1){
        distCoordinates = data.frame(xytab$xVec1[1:(nrow(xytab) - 1)], xytab$yVec1[1:(nrow(xytab) - 1)], xytab$xVec1[2:nrow(xytab)], xytab$yVec1[2:nrow(xytab)], stringsAsFactors = TRUE)
        dist = apply(distCoordinates,1,euc, dist.measurement)
        dist1 = c(NA, dist)
      }else{#if nrow(xytab) == 0 or 1.
        dist1 <- NA
      }
      return(dist1)
    }
    dist.measurement = lonlat
    rownames(x) <-seq(1,nrow(x),1) #This renames the rownames because they may no longer be sequential following the confinementFilter
    indivSeqFrame=data.frame(unique(x$idVec1), stringsAsFactors = TRUE) #The list of individual IDs.
    if(parallel == TRUE){
      cl<-parallel::makeCluster(nCores)
      on.exit(parallel::stopCluster(cl))
      distance = parallel::parApply(cl,indivSeqFrame,1,indivDist, y = x, dist.measurement) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point travelled to reach it from the subsequent point.
    }else{
      distance = apply(indivSeqFrame,1,indivDist, y = x, dist.measurement) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point travelled to reach it from the subsequent point.
    }
    distTotal = unlist(distance)
    if(length(totSecVec) > 1){
      totalSecond.prevStep <- c(NA,totSecVec[1:(length(totSecVec) - 1)]) #There's no need to subset this by id (to make sure there's an NA at the beginning of each individual's path) because NAs are already present in the appropriate places in distTotal.
    }else{ #if length(totSecVec) <= 1
      totalSecond.prevStep <- NA
    }
    secondDifference <- totSecVec - totalSecond.prevStep
    mps <- distTotal/secondDifference
    nanVec <- which(is.nan(mps) == TRUE) #NaN values in mps indicate that distance/seconds = 0/0. In otherwords, observation i was a duplicate of observation i-1. This is only applicable if duplicates were not removed in the previous step (i.e., duplicateFilter or filterOutput == FALSE). (Note: If a duplicate timepoint existed, but individuals were in a different location from the previous observation, an NaN value would not be produced. Instead, distance/seconds = (>0)/0 would return an Inf value)
    if(length(nanVec) > 0){ #This loop replaces the NaN values with the value at observation i-1. Due to the nature of this subfunction, it is impossible to determine a correct mps value if Inf is returned.
      nanVecPrev <- nanVec - 1
      mps[nanVec] <- mps[nanVecPrev]
    }
    remove = which(mps > mpsThreshold)
    if(filterOutput == TRUE){
      if(length(remove) > 0){
        x<-x[-remove,] #removes records where the distance traveled was greater than the mps threshold
      }
    }else{
      x$mps <- mps
    }
    if(nrow(x) > 0){
      rownames(x) <-seq(1,nrow(x),1)
    }
    x <- x[,-match("xVec1",names(x))]
    x <- x[,-match("yVec1",names(x))]
    x <- x[,-match("idVec1",names(x))]
    x <- x[,-match("dateTimeVec1",names(x))]
    return(x)
  }
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    list.mps <- foreach::foreach(l = seq(from = 1, to = length(x), by = 1)) %do% filter3.func(x[[l]], id, point.x, point.y, dateTime, mpsThreshold, lonlat, parallel, filterOutput, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here.
    return(list.mps)
  }else{ #if x is a dataFrame
    frame.mps<- filter3.func(x, id, point.x, point.y, dateTime, mpsThreshold, lonlat, parallel, filterOutput, nCores)
    return(frame.mps)
  }
}
