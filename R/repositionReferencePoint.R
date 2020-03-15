#' Move Data Point a Specified Distance
#'
#' Translates locations of a single rfid tag/gps transmitter to a different 
#'    location a fixed distance away, given a known angular offset (in 
#'    degrees), while maintaining orientations associated with observed
#'    movements (see vignette or Farthing et al. in Review (note: when this 
#'    manuscript is officially published, we will update this 
#'    citation/reference information)) For example, calves in our study (see 
#'    calves2018) were equiped with RFID tags on their left ear. With this 
#'    function, we can move this reference point somewhere else on the body of 
#'    each individual. This might be done for a number of reasons, but is very 
#'    useful for use in the referencePoint2Polygon function later on (for 
#'    delineating polygons representing entire individuals). Currently, this 
#'    function only supports input data with coordinates representing planar 
#'    ('Euclidean') space (e.g. units of meters).
#'
#' In this function, if the distance individuals moved was less than/equal to 
#'    the noted immobThreshold, individuals are said to be immobile ("immob"),
#'    and their position will not change relative to their previous one. (i.e.,
#'    you assume that any observed movement less than immobThreshold was due to
#'    errors or miniscule bodily movements (e.g., head shaking) that are not 
#'    indicative of actual movement.)
#'
#' If distance == NULL, then function will require information (dist, dx, dy) 
#'    from 2 points on an individual's path to work properly. Because of this, 
#'    when no gyroscopic data are provided, at least the first point in each 
#'    individual's path will be removed (the function will report NAs for 
#'    adjusted locations). Also note that if the distance between an 
#'    individual's first point in their path and the second one is 0, the 
#'    function will also report NAs for the second point's adjusted 
#'    coordinates. The first non-NA values will only be reported for the 
#'    instance where dist > 0.
#'
#' Note that populating the direction argument with gyroscopic accelerometer 
#'    data (or data collected using similar devices) collected concurrently 
#'    with point-locations allows us to overcome a couple of assumptions 
#'    associated with using point-locations alone.
#'    
#'    First, unless the direction argument is specifically given (i.e., 
#'    direction != NULL), new point-locations in output are subject to the 
#'    assumption that dt values are sufficiently small to capture individuals'
#'    orientations (i.e., individuals do not face unknown directions 
#'    inbetween observed relocations). If input was previously processed using 
#'    tempAggregate with resolutionLevel == "reduced," dt > secondAgg indicates
#'    that tracked individuals were missing in the original dataset for a 
#'    period of time. In this case, the assumption that individuals are facing 
#'    a given direction because they moved from the previous timepoint may not 
#'    be accurate. Consider removing these rows (rows following one with 
#'    dt > secondAgg; remember that dt indicates the time between recording xy 
#'    coordinates in row i to row i + 1) from your data set.
#'
#'    Second, unless the direction argument is specifically given (i.e., 
#'    direction != NULL), this function assumes tracked individuals are always 
#'    forward-facing. This is because by observing only a single point on each 
#'    individual, we cannot ascertain the true positioning of individuals' 
#'    bodies. For example, even if we know a point-location moved x distance in
#'    a 90-degree direction, from this information alone we cannot determine 
#'    what direction said individual was facing at the time (e.g., this could
#'    be an example of forward, bawckward, or sideward movement). However, 
#'    gyroscopic data (or data collected using similar devices) can tell us 
#'    absolute movement directions, as opposed to relative ones.
#'    
#' @param x Data frame or list of data frames containing real-time-location 
#'    point data.  
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
#' @param direction Numerical vector of length nrow(data.frame(x)) or singular 
#'    character data detailing the relevant colname in x, that denotes what 
#'    movement-direction information will be used. Observations in this vector
#'    represent the direction (in degrees) that tracked individuals moved to 
#'    reach their position at each time point, NOT the direction that they will
#'    move to reach their subsequent position (i.e., values represent known
#'    orientations at each time point). Note that for the purposes of this 
#'    function, observations of 0, 90, 180, and 270 degrees indicates that an 
#'    individual moved straight Eastward, Northward, Westward, and Southward,
#'    respectively. If NULL, direction will be calculated using observed 
#'    point-locations. Defaults to NULL.
#' @param repositionAngle Numerical. Describes the angle (in degrees) 
#'    between empirical point-locations and the desired vertex location as 
#'    represented in a planar model (see vignette or Farthing et al. in Review 
#'    (note: when this manuscript is officially published, we will update this 
#'    citation/reference information)). Essentially, this is the direction you
#'    want new points to be from orginal points. Note that for the purposes of 
#'    this function, observations of 0, 90, 180, and 270 degrees indicates that
#'    an individual moved straight Eastward, Northward, Westward, and 
#'    Southward, respectively. Defaults to 0.
#' @param repositionDist Numerical. Describes the distance from the empirical 
#'    point-locations to desired vertex locations in planar units (e.g., 
#'    meters) inherent to the real-time-location input. Defaults to 1.
#' @param immobThreshold Numerical. Describes what we call, the immobility 
#'    threshold, which is a movement distance (in planar units) within which we
#'    assume individuals’ physical locations and orientations remain unchanged.
#'    This immobility threshold allows us to discount observed movements so 
#'    miniscule that the majority of animals’ physical-space usage is likely 
#'    unaffected (e.g., head shaking). Defaults to 0.
#' @param parallel Logical. If TRUE, sub-functions within the 
#'    repositionReferencePoint wrapper will be parallelized. Note that this can
#'    significantly speed up processing of relatively small data sets, but may 
#'    cause R to crash due to lack of available memory when attempting to 
#'    process large datasets. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
#' @param modelOrientation Numerical. Describes the relative orientation (in 
#'   degrees) of a planar model (see vignette or Farthing et al. in Press 
#'   (note: when this manuscript is officially published, we will update this 
#'   citation/reference information)) describing vertex locations relative to
#'   tracking-device point-locations. Defaults to 90.
#' @keywords data-processing location point planar
#' @return Output is a data frame with the following columns:
#'       
#'    \item{id}{Unique ID of tracked individuals.}
#'    \item{x.original}{Original x coordinates from input.}
#'    \item{y.original}{Original y coordinates from input.}
#'    \item{distance.original}{Original planar distance (m) between 
#'    point-location i to point-location i + 1.}
#'    \item{dx.original}{Original difference between point-location 
#'    x-coordinate i to x-coordinate i + 1.}
#'    \item{dy.original}{Original difference between point-location 
#'    y-coordinate i to y-coordinate i + 1.}
#'    \item{x.adjusted}{Translated x coordinates.}
#'    \item{y.adjusted}{Translated y coordinates.}
#'    \item{dist.adjusted}{Planar distance (m) between translated
#'    point-location i to translated point-location i + 1.}
#'    \item{dx.adjusted}{Difference between translated point-location 
#'    x-coordinate i to translated x-coordinate i + 1.}
#'    \item{dy.adjusted}{Difference between translated point-location 
#'    y-coordinate i to translated y-coordinate i + 1.}
#'    \item{movementDirection}{Describes the angle of movement (in degrees) 
#'    required to translate point-locations to be congruent with planar-model
#'    adjustments.}
#'    \item{repositionAngle}{Describes the value \code{repositionAngle} of the
#'    argument.}
#'    \item{repositionDist}{Describes the value \code{repositionDist} of the
#'    argument.}
#'    \item{immob}{If "0", distance between observed movements is 
#'    < \code{immobThreshold}.}
#'    \item{immobThreshold}{Returns the value from the \code{immobThreshold}
#'    argument.}
#'    \item{dateTime}{Timepoint at which polygons were observed.}
#'    \item{dt}{The the time between reported xy coordinates in row i to row 
#'    i + 1 in each individuals' movement path.}
#' @references Farthing, T.S., Dawson, D.E., Sanderson, M.W., and Lanzas, 
#'    C. in Press. Accounting for space and uncertainty in real-time-location-
#'    system-derived contact networks. Ecology and Evolution.
#' @import foreach
#' @export
#' @examples
#' \donttest{
#' data("calves")
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'    time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#' 
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth to 5-min fix intervals.
#' 
#' leftShoulder.point<-repositionReferencePoint(x = calves.agg, 
#'    id = calves.agg$id, dateTime = calves.agg$dateTime, 
#'    point.x = calves.agg$x, point.y = calves.agg$y, direction = NULL, 
#'    repositionAngle = 180, repositionDist = 0.0835, immobThreshold = 0, parallel = FALSE, 
#'    modelOrientation = 90)
#'    }

repositionReferencePoint <- function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, direction = NULL, repositionAngle = 0, repositionDist = 1, immobThreshold = 0, parallel = FALSE, nCores = (parallel::detectCores()/2), modelOrientation = 90){ 

  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  i <- NULL
  j <- NULL
  k <- NULL
  l <- NULL
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    
    listBreak_reposition.generator <-function(x, id, dateTime, point.x, point.y, direction, repositionAngle, repositionDist, immobThreshold, parallel, nCores, modelOrientation){ #this function is exactly the same as what happens when the x input to the master function is a single data frame.
      
      #define sub-functions
      reposition.generator<-function(x, repositionAngle, repositionDist, immobThreshold, modelOrientation){
        euc=function(x) {
          point1 = x.cor=unlist(c(x[1],x[2]))
          point2 = x.cor=unlist(c(x[3],x[4]))
          euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
          return(euc.dis)
        }
        translate<-function(x, y, eta, etaStar, theta, distV){
          #(2*etaStar) - (theta + eta) #This is the angle of the new vertex relative to a clockwise direction from a point etaStar degrees from empirical points given in x
          vertex.angle<- (360 - (etaStar - (theta + eta))) #This is the relative angle of the new vertex from the empirical points given in x
          #simplify this value now
          vertex.angle360s <- floor(vertex.angle/360) # determine how many times individuals turn in complete circles
          vertex.angle.simplified <- vertex.angle - (vertex.angle360s*360) #this reduces direction values to between 0 and 359 
          vertex.angle.simplified.rad<- vertex.angle.simplified / (180/pi) #convert back to radians
          
          newX<- distV*cos(vertex.angle.simplified.rad) + x #identify the x coordinate that lies distV units in the vertex.angle.simplified.rad direction from the empirical x
          newY<- distV*sin((vertex.angle.simplified.rad)) + y #identify the y coordinate that lies distV units in the vertex.angle.simplified.rad direction from the empirical y
          
          output<-data.frame(x.adjusted = newX, y.adjusted = newY, stringsAsFactors = TRUE)
          return(output)
        }
        immobAdjustment.point = function(x, locMatrix, immobVec){
          dur = 0
          standDurTest <- FALSE
          while (!standDurTest) { #This means when standDurTest == FALSE #This part of the function determines how long an individual has been immob.
            standDurTest <- immobVec[(x[1] - (1 + dur))] == 0
            standDurTest <- ifelse(is.na(standDurTest) == TRUE, TRUE, standDurTest) #This is for when x[1] -1 relates to the first observation for any individual (because those observations are NA)
            dur = (dur + 1)
          }
          standFrame = data.frame(replaceRow = x[1], replacePoint.x = locMatrix$x.adjusted[(x[1] - dur)], replacePoint.y = locMatrix$y.adjusted[(x[1] - dur)], stringsAsFactors = TRUE) #if individuals are immob, their location at a given timestep will be the same as at the previous timestep.
          return(standFrame)
        }
        timeDifference = function(x){
          t1 = x[1]
          t2 = x[2]
          dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
          return(dt)
        }
        
        distCoordinates = data.frame(x$x[1:(nrow(x) - 1)], x$y[1:(nrow(x) - 1)], x$x[2:nrow(x)], x$y[2:nrow(x)], stringsAsFactors = TRUE)
        
        dist = apply(distCoordinates,1,euc) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
        dist = c(dist, NA) #to make dist the same length as nrow(x)
        
        dx = c((distCoordinates[,3] - distCoordinates[,1]),NA) #calculate differences in x (x2 - x1) and put an NA at the end to ensure it is of equal length to the input data
        dy = c((distCoordinates[,4] - distCoordinates[,2]),NA) #calculate differences in y (y2 - y1) and put an NA at the end to ensure it is of equal length to the input data
        
        timesFrame = data.frame(x$dateTime[1:(nrow(x) - 1)], x$dateTime[2:nrow(x)], stringsAsFactors = TRUE)
        
        dt = apply(timesFrame, 1, timeDifference)
        dt = c(dt, NA) #to make dt the same length as nrow(x)
        
        idVec = unique(x$id) ; idSeqVec = x$id; idVecForDtRemoval = idVec[-length(idVec)] #There's no need to remove the max point of the final id point b/c there's already an NA there.
        for(a in idVecForDtRemoval){ #This loop removes the dt value at points describing the last point in each individual's path
          dx[max(which(idSeqVec == a))] = NA
          dy[max(which(idSeqVec == a))] = NA
          dist[max(which(idSeqVec == a))] = NA
          dt[max(which(idSeqVec == a))] = NA
        }
        x$dx = dx
        x$dy = dy
        x$dist = dist
        x$dt = dt
        
        dataShift = data.frame(x = x$x, y = x$y, dx = c(NA,x$dx[1:(nrow(x) - 1)]), dy = c(NA,x$dy[1:(nrow(x) - 1)]), dist = c(NA,x$dist[1:(nrow(x) - 1)]), stringsAsFactors = TRUE) #This is necessary because of the way we calculated/listed these values above. These columns (dist, dx, and dt) refer to the changes a point must make to reach the subsequent point. However, later on in this function, we are not interested in future point alterations. Rather, we need to know how tracked individuals moved during the preceding time step to reach their current point (to determine directionality of movement) #Note that this code shifts values down, but remember that beacuse values are downshifted, the first observation for each id individual will be incorrect. Code below addresses this issue.
        for(b in idVec){ #This code replaces dist, dx, and dy values in the row when a new individual (id) is first observed with NAs to fix the problem noted above. #Note that this loop assumes that the data is sorted by id (i.e., a given id will not repeat in the dataset, after a new id has appeared)
          dataShift[min(which(idSeqVec == b)),] = NA
        }
        
        #eta was defined earlier and appended to x as part of the efficiency improvement implemented on 02/20. Here we vectorize it again here, and remove the appended column in x
        if(length(x$eta..) > 0){ #the other x inputs are required, but direction is not
          
          eta <- x$eta..
          x<- droplevels(x[,-match("eta..", colnames(x))])
          
        }else{ #if length(direction) == 0 in the master function input
          eta <- atan2(dataShift$dy,dataShift$dx)*(180/pi) + 360 #calculate the relative angle of movements given point (x1,y1) lies on the axes' origin. Multiplying by 180/pi converts radians to degrees and adding 360 ensures positive values.
        }
        
        #simplify eta (i.e., make all observations fall between 0 and 359)
        eta360s <- floor(eta/360) # determine how many times individuals turn in complete circles
        eta.simplified <- eta - (eta360s*360) #this reduces direction values to between 0 and 359
        
        #translate the empirical coordinates to create the new ones
        newCoordinates <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = repositionAngle, distV = repositionDist) #repositions the empirical point directly to the East (i.e., right) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the top-right corner of the polygon.
        
        RepositionMatrix <- data.frame(matrix(nrow = nrow(x), ncol = (18)), stringsAsFactors = TRUE) #One row for each observation, Columns listed below
        colnames(RepositionMatrix) <- c("id","x.original","y.original","dist.original", "dx.original", "dy.original","x.adjusted","y.adjusted", "dist.adjusted", "dx.adjusted", "dy.adjusted","movementDirection","repositionAngle","repositionDist","immob","immobThreshold", "dateTime","dt")
        RepositionMatrix$id = x$id
        RepositionMatrix$x.original = x$x
        RepositionMatrix$y.original = x$y
        RepositionMatrix$dist.original = x$dist
        RepositionMatrix$dx.original = x$dx
        RepositionMatrix$dy.original = x$dy
        RepositionMatrix$movementDirection = eta.simplified
        RepositionMatrix$repositionAngle = repositionAngle
        RepositionMatrix$repositionDist = repositionDist
        RepositionMatrix$immobThreshold = immobThreshold
        RepositionMatrix$dateTime = x$dateTime
        RepositionMatrix$dt = x$dt
        RepositionMatrix$x.adjusted = newCoordinates$x.adjusted
        RepositionMatrix$y.adjusted = newCoordinates$y.adjusted
        RepositionMatrix$immob = ifelse(dataShift$dist > immobThreshold, 0, 1) #if the distance individuals moved was less than / equal to the noted immobThreshold, individuals are said to be "immob," and their position will not change relative to their previous one. (i.e., you assume that any observed movement less than immobThreshold was due to errors or miniscule bodily movements (e.g., head shaking) that are not indicative of actual movement.)
        
        standlist = which(RepositionMatrix$immob == 1)
        
        if ((length(standlist) >= 1)){ #To save processing time and reduce chances of errors, this evaluation will not take place if there is only one observation.
          
          immobVec = RepositionMatrix$immob
          immob.list <- foreach::foreach(k = standlist) %do% immobAdjustment.point(k, locMatrix = RepositionMatrix, immobVec)
          immobFrame <- data.frame(data.table::rbindlist(immob.list), stringsAsFactors = TRUE)
          
          if(nrow(immobFrame) > 0){
            RepositionMatrix[immobFrame$replaceRow,match("x.adjusted", colnames(RepositionMatrix))] = immobFrame$replacePoint.x
            RepositionMatrix[immobFrame$replaceRow,match("y.adjusted", colnames(RepositionMatrix))] = immobFrame$replacePoint.y
          }
        }
        
        newDistCoordinates = data.frame(RepositionMatrix$x.adjusted[1:(nrow(RepositionMatrix) - 1)], RepositionMatrix$y.adjusted[1:(nrow(RepositionMatrix) - 1)], RepositionMatrix$x.adjusted[2:nrow(RepositionMatrix)], RepositionMatrix$y.adjusted[2:nrow(RepositionMatrix)], stringsAsFactors = TRUE)
        
        newDist = apply(newDistCoordinates,1,euc) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
        newDist = c(newDist, NA)
        
        RepositionMatrix$dist.adjusted = newDist
        
        newdx = c((newDistCoordinates[,3] - newDistCoordinates[,1]),NA)
        newdy = c((newDistCoordinates[,4] - newDistCoordinates[,2]),NA)
        RepositionMatrix$dx.adjusted = newdx
        RepositionMatrix$dy.adjusted = newdy
        
        return(RepositionMatrix)
      }
      
      day_listReposition <- function(x, data.list, repositionAngle, repositionDist, immobThreshold, modelOrientation){ #Because this function slows down when trying to process large data frames AND large list sets, we must concattenate both here. We did so to the former by breaking the data frame into hourly lists, and the latter by breaking these lists into daily subsets with this function.
        
        day_lists <- data.list[grep(unname(unlist(x[1])), names(data.list))] #pulls the hour lists within a given day
        names(day_lists)<-NULL #ensure that list names do not mess up column names
        list.reposition <- lapply(day_lists, reposition.generator, repositionAngle, repositionDist, immobThreshold, modelOrientation)
        reposition.bind <- data.frame(data.table::rbindlist(list.reposition, fill = TRUE), stringsAsFactors = TRUE) #bind these hours back together
        
        return(reposition.bind)
      }
      
      date_hourSub.func<-function(x, data, date_hour){ #added 02/20/2020 #This function will be used to break down data sets into hourly time blocks prior to further processing to increase speed. Admittedly, this is not a pretty fix for increasing efficiency of processing large data sets, but it's a working fix nonetheless. 
        date_hour1 <- droplevels(data[which(date_hour == unname(unlist(x[1]))),]) #subset data
        return(date_hour1)
      }
      
      if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors.
        x <- data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime, stringsAsFactors = TRUE)
      }
      
      if(length(x) > 0){ #for some reason using an "else" statement would always result in a table with 0 records...
        
        if(length(id) > 0){
          if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than id being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
            x$id <- x[,match(id, names(x))]
          }else{ #if length(id) > 1
            x$id = id
          }
        }
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
        if(length(dateTime) > 0){
          if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
            x$dateTime <- x[,match(dateTime, names(x))]
          }else{ #if length(dateTime) > 1
            x$dateTime = dateTime
          }
        }
      }
      
      if(length(direction) > 0){ #we create a disposable column in x to pull direction values into the generator function (this is required here because we break x into smaller lists below)
        if(length(direction) == 1 && is.na(match(direction[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$eta.. <- x[,match(direction, names(x))]
        }else{ #if length(direction) > 1
          x$eta.. <- direction
        }
      }
      
      x<-x[order(x$id, x$dateTime),] #data tables should already be ordered like this, but because this function assumes this ordering scheme, I wanted to ensure it was so.
      rownames(x) <- seq(1, nrow(x),1)
      
      data.dates<-lubridate::date(x$dateTime) #now we can start concattenating the data by subsetting it into smaller lists
      
      date_hour <- paste(data.dates, lubridate::hour(x$dateTime), sep = "_") #create a tag for each unique date_hour combination in the data set
      date_hour.vec <- unique(date_hour)
      date.vec <- unique(data.dates)
      
      if(length(date_hour.vec) == 1){ #the processing step requires a list of data frames. If there's only a single hour represented in x, we can just create the list using the "list" function
        data.list <- list(x)
      }else{
        data.list <- foreach::foreach(i = date_hour.vec) %do% date_hourSub.func(i, x, date_hour)
      }
      
      names(data.list)<-date_hour.vec #add names to list to pull for date lists below
      
      rm(list =  c("x", "data.dates", "date_hour.vec")) #remove the unneeded objects to free up memory
      
      if(parallel == TRUE){
        
        cl <- parallel::makeCluster(nCores)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))
        repositions<-foreach::foreach(j = date.vec, .packages = 'foreach') %dopar% day_listReposition(j, data.list, repositionAngle, repositionDist, immobThreshold, modelOrientation) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
        
      }else{ #if parallel == FALSE
        
        repositions<-foreach::foreach(j = date.vec, .packages = 'foreach') %do% day_listReposition(j, data.list, repositionAngle, repositionDist, immobThreshold, modelOrientation) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
        
      }
      
      frame.rep<- data.frame(data.table::rbindlist(repositions, fill = TRUE), stringsAsFactors = TRUE)
      
      return(frame.rep)
      
    }
    
    list.rep<- foreach::foreach(l = 1:length(x), .packages = 'foreach') %do% listBreak_reposition.generator(x[[l]], id, dateTime, point.x, point.y, direction, repositionAngle, repositionDist, immobThreshold, parallel, nCores, modelOrientation) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops
    
    return(list.rep)

  }else{ #if x is a dataFrame

    #define sub-functions
    reposition.generator<-function(x, repositionAngle, repositionDist, immobThreshold, modelOrientation){
      euc=function(x) {
        point1 = x.cor=unlist(c(x[1],x[2]))
        point2 = x.cor=unlist(c(x[3],x[4]))
        euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
        return(euc.dis)
      }
      translate<-function(x, y, eta, etaStar, theta, distV){
        #(2*etaStar) - (theta + eta) #This is the angle of the new vertex relative to a clockwise direction from a point etaStar degrees from empirical points given in x
        vertex.angle<- (360 - (etaStar - (theta + eta))) #This is the relative angle of the new vertex from the empirical points given in x
        #simplify this value now
        vertex.angle360s <- floor(vertex.angle/360) # determine how many times individuals turn in complete circles
        vertex.angle.simplified <- vertex.angle - (vertex.angle360s*360) #this reduces direction values to between 0 and 359 
        vertex.angle.simplified.rad<- vertex.angle.simplified / (180/pi) #convert back to radians
        
        newX<- distV*cos(vertex.angle.simplified.rad) + x #identify the x coordinate that lies distV units in the vertex.angle.simplified.rad direction from the empirical x
        newY<- distV*sin((vertex.angle.simplified.rad)) + y #identify the y coordinate that lies distV units in the vertex.angle.simplified.rad direction from the empirical y
        
        output<-data.frame(x.adjusted = newX, y.adjusted = newY, stringsAsFactors = TRUE)
        return(output)
      }
      immobAdjustment.point = function(x, locMatrix, immobVec){
        dur = 0
        standDurTest <- FALSE
        while (!standDurTest) { #This means when standDurTest == FALSE #This part of the function determines how long an individual has been immob.
          standDurTest <- immobVec[(x[1] - (1 + dur))] == 0
          standDurTest <- ifelse(is.na(standDurTest) == TRUE, TRUE, standDurTest) #This is for when x[1] -1 relates to the first observation for any individual (because those observations are NA)
          dur = (dur + 1)
        }
        standFrame = data.frame(replaceRow = x[1], replacePoint.x = locMatrix$x.adjusted[(x[1] - dur)], replacePoint.y = locMatrix$y.adjusted[(x[1] - dur)], stringsAsFactors = TRUE) #if individuals are immob, their location at a given timestep will be the same as at the previous timestep.
        return(standFrame)
      }
      timeDifference = function(x){
        t1 = x[1]
        t2 = x[2]
        dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
        return(dt)
      }
      
      distCoordinates = data.frame(x$x[1:(nrow(x) - 1)], x$y[1:(nrow(x) - 1)], x$x[2:nrow(x)], x$y[2:nrow(x)], stringsAsFactors = TRUE)

      dist = apply(distCoordinates,1,euc) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
      dist = c(dist, NA) #to make dist the same length as nrow(x)
      
      dx = c((distCoordinates[,3] - distCoordinates[,1]),NA) #calculate differences in x (x2 - x1) and put an NA at the end to ensure it is of equal length to the input data
      dy = c((distCoordinates[,4] - distCoordinates[,2]),NA) #calculate differences in y (y2 - y1) and put an NA at the end to ensure it is of equal length to the input data
      
      timesFrame = data.frame(x$dateTime[1:(nrow(x) - 1)], x$dateTime[2:nrow(x)], stringsAsFactors = TRUE)

      dt = apply(timesFrame, 1, timeDifference)
      dt = c(dt, NA) #to make dt the same length as nrow(x)

      idVec = unique(x$id) ; idSeqVec = x$id; idVecForDtRemoval = idVec[-length(idVec)] #There's no need to remove the max point of the final id point b/c there's already an NA there.
      for(a in idVecForDtRemoval){ #This loop removes the dt value at points describing the last point in each individual's path
        dx[max(which(idSeqVec == a))] = NA
        dy[max(which(idSeqVec == a))] = NA
        dist[max(which(idSeqVec == a))] = NA
        dt[max(which(idSeqVec == a))] = NA
      }
      x$dx = dx
      x$dy = dy
      x$dist = dist
      x$dt = dt
      
      dataShift = data.frame(x = x$x, y = x$y, dx = c(NA,x$dx[1:(nrow(x) - 1)]), dy = c(NA,x$dy[1:(nrow(x) - 1)]), dist = c(NA,x$dist[1:(nrow(x) - 1)]), stringsAsFactors = TRUE) #This is necessary because of the way we calculated/listed these values above. These columns (dist, dx, and dt) refer to the changes a point must make to reach the subsequent point. However, later on in this function, we are not interested in future point alterations. Rather, we need to know how tracked individuals moved during the preceding time step to reach their current point (to determine directionality of movement) #Note that this code shifts values down, but remember that beacuse values are downshifted, the first observation for each id individual will be incorrect. Code below addresses this issue.
      for(b in idVec){ #This code replaces dist, dx, and dy values in the row when a new individual (id) is first observed with NAs to fix the problem noted above. #Note that this loop assumes that the data is sorted by id (i.e., a given id will not repeat in the dataset, after a new id has appeared)
        dataShift[min(which(idSeqVec == b)),] = NA
      }
  
      #eta was defined earlier and appended to x as part of the efficiency improvement implemented on 02/20. Here we vectorize it again here, and remove the appended column in x
      if(length(x$eta..) > 0){ #the other x inputs are required, but direction is not
        
        eta <- x$eta..
        x<- droplevels(x[,-match("eta..", colnames(x))])
        
      }else{ #if length(direction) == 0 in the master function input
        eta <- atan2(dataShift$dy,dataShift$dx)*(180/pi) + 360 #calculate the relative angle of movements given point (x1,y1) lies on the axes' origin. Multiplying by 180/pi converts radians to degrees and adding 360 ensures positive values.
      }
      
      #simplify eta (i.e., make all observations fall between 0 and 359)
      eta360s <- floor(eta/360) # determine how many times individuals turn in complete circles
      eta.simplified <- eta - (eta360s*360) #this reduces direction values to between 0 and 359
      
      #translate the empirical coordinates to create the new ones
      newCoordinates <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = repositionAngle, distV = repositionDist) #repositions the empirical point directly to the East (i.e., right) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the top-right corner of the polygon.
      
      RepositionMatrix <- data.frame(matrix(nrow = nrow(x), ncol = (18)), stringsAsFactors = TRUE) #One row for each observation, Columns listed below
      colnames(RepositionMatrix) <- c("id","x.original","y.original","dist.original", "dx.original", "dy.original","x.adjusted","y.adjusted", "dist.adjusted", "dx.adjusted", "dy.adjusted","movementDirection","repositionAngle","repositionDist","immob","immobThreshold", "dateTime","dt")
      RepositionMatrix$id = x$id
      RepositionMatrix$x.original = x$x
      RepositionMatrix$y.original = x$y
      RepositionMatrix$dist.original = x$dist
      RepositionMatrix$dx.original = x$dx
      RepositionMatrix$dy.original = x$dy
      RepositionMatrix$movementDirection = eta.simplified
      RepositionMatrix$repositionAngle = repositionAngle
      RepositionMatrix$repositionDist = repositionDist
      RepositionMatrix$immobThreshold = immobThreshold
      RepositionMatrix$dateTime = x$dateTime
      RepositionMatrix$dt = x$dt
      RepositionMatrix$x.adjusted = newCoordinates$x.adjusted
      RepositionMatrix$y.adjusted = newCoordinates$y.adjusted
      RepositionMatrix$immob = ifelse(dataShift$dist > immobThreshold, 0, 1) #if the distance individuals moved was less than / equal to the noted immobThreshold, individuals are said to be "immob," and their position will not change relative to their previous one. (i.e., you assume that any observed movement less than immobThreshold was due to errors or miniscule bodily movements (e.g., head shaking) that are not indicative of actual movement.)
      
      standlist = which(RepositionMatrix$immob == 1)
      
      if ((length(standlist) >= 1)){ #To save processing time and reduce chances of errors, this evaluation will not take place if there is only one observation.
        
        immobVec = RepositionMatrix$immob
        immob.list <- foreach::foreach(k = standlist) %do% immobAdjustment.point(k, locMatrix = RepositionMatrix, immobVec)
        immobFrame <- data.frame(data.table::rbindlist(immob.list), stringsAsFactors = TRUE)
        
        if(nrow(immobFrame) > 0){
          RepositionMatrix[immobFrame$replaceRow,match("x.adjusted", colnames(RepositionMatrix))] = immobFrame$replacePoint.x
          RepositionMatrix[immobFrame$replaceRow,match("y.adjusted", colnames(RepositionMatrix))] = immobFrame$replacePoint.y
        }
      }
      
      newDistCoordinates = data.frame(RepositionMatrix$x.adjusted[1:(nrow(RepositionMatrix) - 1)], RepositionMatrix$y.adjusted[1:(nrow(RepositionMatrix) - 1)], RepositionMatrix$x.adjusted[2:nrow(RepositionMatrix)], RepositionMatrix$y.adjusted[2:nrow(RepositionMatrix)], stringsAsFactors = TRUE)
      
      newDist = apply(newDistCoordinates,1,euc) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
      newDist = c(newDist, NA)
      
      RepositionMatrix$dist.adjusted = newDist
      
      newdx = c((newDistCoordinates[,3] - newDistCoordinates[,1]),NA)
      newdy = c((newDistCoordinates[,4] - newDistCoordinates[,2]),NA)
      RepositionMatrix$dx.adjusted = newdx
      RepositionMatrix$dy.adjusted = newdy
      
      return(RepositionMatrix)
    }
    
    day_listReposition <- function(x, data.list, repositionAngle, repositionDist, immobThreshold, modelOrientation){ #Because this function slows down when trying to process large data frames AND large list sets, we must concattenate both here. We did so to the former by breaking the data frame into hourly lists, and the latter by breaking these lists into daily subsets with this function.
      
      day_lists <- data.list[grep(unname(unlist(x[1])), names(data.list))] #pulls the hour lists within a given day
      names(day_lists)<-NULL #ensure that list names do not mess up column names
      list.reposition <- lapply(day_lists, reposition.generator, repositionAngle, repositionDist, immobThreshold, modelOrientation)
      reposition.bind <- data.frame(data.table::rbindlist(list.reposition, fill = TRUE), stringsAsFactors = TRUE) #bind these hours back together
      
      return(reposition.bind)
    }
    
    date_hourSub.func<-function(x, data, date_hour){ #added 02/20/2020 #This function will be used to break down data sets into hourly time blocks prior to further processing to increase speed. Admittedly, this is not a pretty fix for increasing efficiency of processing large data sets, but it's a working fix nonetheless. 
      date_hour1 <- droplevels(data[which(date_hour == unname(unlist(x[1]))),]) #subset data
      return(date_hour1)
    }
    
    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors.
      x <- data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime, stringsAsFactors = TRUE)
    }
    
    if(length(x) > 0){ #for some reason using an "else" statement would always result in a table with 0 records...
      
      if(length(id) > 0){
        if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than id being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
          x$id <- x[,match(id, names(x))]
        }else{ #if length(id) > 1
          x$id = id
        }
      }
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
      if(length(dateTime) > 0){
        if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$dateTime <- x[,match(dateTime, names(x))]
        }else{ #if length(dateTime) > 1
          x$dateTime = dateTime
        }
      }
    }
    
    if(length(direction) > 0){ #we create a disposable column in x to pull direction values into the generator function (this is required here because we break x into smaller lists below)
      if(length(direction) == 1 && is.na(match(direction[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
        x$eta.. <- x[,match(direction, names(x))]
      }else{ #if length(direction) > 1
        x$eta.. <- direction
      }
    }
    
    x<-x[order(x$id, x$dateTime),] #data tables should already be ordered like this, but because this function assumes this ordering scheme, I wanted to ensure it was so.
    rownames(x) <- seq(1, nrow(x),1)
    
    data.dates<-lubridate::date(x$dateTime) #now we can start concattenating the data by subsetting it into smaller lists
    
    date_hour <- paste(data.dates, lubridate::hour(x$dateTime), sep = "_") #create a tag for each unique date_hour combination in the data set
    date_hour.vec <- unique(date_hour)
    date.vec <- unique(data.dates)
    
    if(length(date_hour.vec) == 1){ #the processing step requires a list of data frames. If there's only a single hour represented in x, we can just create the list using the "list" function
      data.list <- list(x)
    }else{
      data.list <- foreach::foreach(i = date_hour.vec) %do% date_hourSub.func(i, x, date_hour)
    }
    
    names(data.list)<-date_hour.vec #add names to list to pull for date lists below
    
    rm(list =  c("x", "data.dates", "date_hour.vec")) #remove the unneeded objects to free up memory
    
    if(parallel == TRUE){
      
      cl <- parallel::makeCluster(nCores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      repositions<-foreach::foreach(j = date.vec, .packages = 'foreach') %dopar% day_listReposition(j, data.list, repositionAngle, repositionDist, immobThreshold, modelOrientation) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
      
    }else{ #if parallel == FALSE
      
      repositions<-foreach::foreach(j = date.vec, .packages = 'foreach') %do% day_listReposition(j, data.list, repositionAngle, repositionDist, immobThreshold, modelOrientation) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
      
    }
    
    frame.rep<- data.frame(data.table::rbindlist(repositions, fill = TRUE), stringsAsFactors = TRUE)
    
    return(frame.rep)
  }

}
