#' Create a Rectangular Polygon Using Planar XY Coordinates
#'
#' This function creates a square/rectangular polygon from a single reference 
#'    point by translating its location multiple times using the same method 
#'    used in repositionReferencePoint. For example, even though calves in our 
#'    study see (data(calves2018)) were only equiped with RFID tags on their 
#'    left ear. With this function, we can create polygons that account for the
#'    total space used by each individual at each time step.This function is 
#'    different from similar point-to-polygon functions for two reasons:
#'    1.) It does not assume points lie within the center of the polygon. 
#'        Rather, the reference point must be a corner of the polygon 
#'        (Note: "UL" denotes that reference point lies on the upper-left 
#'        corner of the polygon, "UR" denotes that reference point lies on the 
#'        upper-right corner of the polygon,"DL" denotes that reference point 
#'        lies on the down-left corner of the polygon, "DR" denotes that 
#'        reference point lies on the down-left corner of the polygon). Note 
#'        that if you want the reference point to be at the center of the 
#'        polygon, you can first translate the reference point to a central 
#'        location on tracked individuals using repositionReferencePoint.
#'    2.) Polygon angles/directionality are based on observed movements of 
#'        tracked individuals.
#'
#' Currently, this function only supports input data with coordinates 
#'    representing planar ('Euclidean') space (e.g. units of meters).
#'
#' In the output, point1.x and point1.y represent the xy coordinates from the 
#'    input file. Point2-n coordinates move in a clockwise direction from 
#'    point1. For example: if point1 is located on the upper left ("UL") corner
#'    of the polygon, point2 would be on the upper right corner, point3 on the 
#'    bottom right, and point 4 on the bottom left.
#'
#' Because this function needs information (dist, dx, dy) from 2 points on an 
#'    individual's path to work, at least the first point in each individual's 
#'    path will be removed (the function will report NAs for adjusted 
#'    locations). Also note that if the distance between an individual's first 
#'    point in their path and the second one is 0, the function will also 
#'    report NAs for the second point's adjusted coordinates. The first non-NA 
#'    values will only be reported for the instance where dist > 0.
#'
#' Note that populating the direction argument with gyroscopic accelerometer 
#'    data (or data collected using similar devices) collected concurrently 
#'    with point-locations allows us to overcome a couple of assumptions 
#'    associtated with using point-locations alone.
#'    
#'    First, unless the direction argument is specifically given (i.e., 
#'    direction != NULL), vertex locations in output are subject to the 
#'    assumption that dt values are sufficiently small to capture individuals'
#'    orientations (i.e., individuals do not face unknown directions 
#'    inbetween observed relocations). If input was previously processed using 
#'    tempAggregate with resolutionLevel == "reduced," dt > secondAgg indicates
#'    that tracked individuals were missing in the original dataset for a 
#'    period of time. In this case, the assumption that individuals are facing 
#'    a given direction because they moved from the previous timepoint may not 
#'    be accurate. Consider removing these rows (rows following one with 
#'    dt > secondAgg; remember that dt indicates the time between reported xy 
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
#' @param StartLocation Character string taking the values "UL," "UR," "DL," or
#'    "DR" describing where the reference point (i.e., point corresponding to 
#'    xy-coordinates in the data set) lies on the rectangle that this function 
#'    will delineate. Defaults to "UL."
#' @param UpDownRepositionLen Numerical. Describes the height, in planar units 
#'    (e.g., meters) of the output polygon. Planar units are inherent to the 
#'    real-time-location input. Defaults to 1.
#' @param LeftRightRepositionLen Numerical. Describes the width, in planar 
#'    units (e.g., meters) of the output polygon. Planar units are inherent to 
#'    the real-time-location input. Defaults to 1.
#' @param CenterPoint Logical. If TRUE, in addition to the xy-coordinates for 
#'    each polygon vertex, xy-coordinates for centroid of each polygon will be 
#'    reported in the output. Defaults to FALSE.
#' @param MidPoints Logical. If TRUE, in addition to the xy-coordinates for 
#'    each polygon vertex, xy-coordinates for mid-point of each polygon edge 
#'    will be reported in the output. Defaults to FALSE.
#' @param immobThreshold Numerical. Describes what we call, the immobility 
#'    threshold, which is a movement distance (in planar units) within which 
#'    we assume individuals’ physical locations and orientations remain 
#'    unchanged. This immobility threshold allows us to discount observed 
#'    movements so miniscule that the majority of animals’ physical-space 
#'    usage is likely unaffected (e.g., head shaking). Defaults to 0.
#' @param parallel Logical. If TRUE, sub-functions within the 
#'    referencePoint2Polygon wrapper will be parallelized. Note that this can 
#'    significantly speed up processing of relatively small data sets, but may 
#'    cause R to crash due to lack of available memory when attempting to 
#'    process large datasets. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to the maximum number of cores available
#'    (i.e., parallel::detectCores()).
#' @param modelOrientation Numerical. Describes the relative orientation (in 
#'   degrees) of a planar model (see vignette or Farthing et al. in Review 
#'   (note: when this manuscript is officially published, we will update this 
#'   citation/reference information)) describing vertex locations relative to
#'   tracking-device point-locations. Defaults to 90.
#' @keywords data-processing polygon point location planar
#' @return Output is a data frame with the following columns:
#'       
#'    \item{id}{Unique ID of tracked individuals.}
#'    \item{cornerPoint...x}{Planar x coordinates of polygon-corner vertices.}
#'    \item{cornerPoint...y}{Planar y coordinates of polygon-corner vertices.}
#'    \item{startLocation}{Describes the location of input point-locations in 
#'    the vertex outputs. see \code{StartLocation} argument.}
#'    \item{upDownRepositionLength}{Describes the vertical movement of 
#'    point-locations on planar models. see \code{UpDownRepositionLen} 
#'    argument.}
#'    \item{leftRightRepositionLength}{Describes the horizontal movement of 
#'    point-locations on planar models. see \code{leftRightRepositionLen} 
#'    argument.}
#'    \item{immob}{If "0", distance between observed movements is 
#'    < \code{immobThreshold}.}
#'    \item{immobThreshold}{Returns the value from the \code{immobThreshold}
#'    argument.}
#'    \item{dateTime}{Timepoint at which polygons were observed.}
#'    \item{dt}{The the time between reported xy coordinates in row i to row 
#'    i + 1 in each individuals' movement path.}
#'    
#'    If MidPoints or CenterPoints == TRUE, additional columns will be appended
#'    to output data frame.
#'    
#' @references Farthing, T.S., Dawson, D.E., Sanderson, M.W., and Lanzas, 
#'    C. in Review. Accounting for space and uncertainty in real-time-location-
#'    system-derived contact networks. Ecology and Evolution.
#' @export
#' @examples
#' \donttest{
#' data("calves")
#' calves.dateTime<-datetime.append(calves, date = calves$date,
#'    time = calves$time) #add dateTime identifiers for location fixes.
#' 
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth to 5-min fix intervals.
#' 
#' calf_heads <- referencePoint2Polygon(x = calves.agg,
#'    id = calves.agg$id, dateTime = calves.agg$dateTime,
#'    point.x = calves.agg$x, point.y = calves.agg$y, direction = NULL,
#'    StartLocation = "DL", UpDownRepositionLen = 0.333, LeftRightRepositionLen = 0.333,
#'    CenterPoint = FALSE, MidPoints = FALSE, immobThreshold = 0.1, parallel = FALSE,
#'    modelOrientation = 90)
#'    }

referencePoint2Polygon <-function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, direction = NULL, StartLocation = "UL", UpDownRepositionLen = 1, LeftRightRepositionLen = 1, CenterPoint = FALSE, MidPoints = FALSE, immobThreshold = 0, parallel = FALSE, nCores = parallel::detectCores(), modelOrientation = 90){

  poly.generator<-function(x, id, dateTime, point.x, point.y, direction, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, immobThreshold, parallel, modelOrientation, nCores){
    euc=function(x) {
      point1 = x.cor=unlist(c(x[1],x[2]))
      point2 = x.cor=unlist(c(x[3],x[4]))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    timeDifference = function(x){
      t1 = x[1]
      t2 = x[2]
      dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
      return(dt)
    }
    immobAdjustment.polygon = function(x, locMatrix, immobVec){

      dur = 0
      standDurTest <- FALSE
      while (!standDurTest) { #This means when standDurTest == FALSE #This part of the function determines how long an individual has been immob.
        standDurTest <- immobVec[(x[1] - (1 + dur))] == 0
        standDurTest <- ifelse(is.na(standDurTest) == TRUE, TRUE, standDurTest) #This is for when x[1] -1 relates to the first observation for any individual (because those observations are NA)
        dur = (dur + 1)
      }

      standFrame = data.frame(replaceRow = x[1], replacePoint2.x = locMatrix$cornerPoint2.x[(x[1] - dur)], replacePoint2.y = locMatrix$cornerPoint2.y[(x[1] - dur)], replacePoint3.x = locMatrix$cornerPoint3.x[(x[1] - dur)], replacePoint3.y = locMatrix$cornerPoint3.y[(x[1] - dur)], replacePoint4.x = locMatrix$cornerPoint4.x[(x[1] - dur)], replacePoint4.y = locMatrix$cornerPoint4.y[(x[1] - dur)]) #if individuals are immob, their location at a given timestep will be the same as at the previous timestep. There's no need to replace point1 because it is the pre-translated (a.k.a. original) location of the individual.

      return(standFrame)
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
      
      output<-data.frame(x.adjusted = newX, y.adjusted = newY)
      return(output)
    }

    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors.
      x = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime)
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

    x<-x[order(x$id, x$dateTime),] #data tables should already be ordered like this, but because this function assumes this ordering scheme, I wanted to ensure it was so.
    rownames(x)<-seq(1,nrow(x),1)

    distCoordinates = data.frame(x1 = x$x[1:(nrow(x) - 1)], y1 = x$y[1:(nrow(x) - 1)], x2 = x$x[2:nrow(x)], y2 = x$y[2:nrow(x)])
    if(parallel == TRUE){ #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
      cl<-parallel::makeCluster(nCores)
      dist = parallel::parApply(cl, distCoordinates,1,euc)
      dist = c(dist, NA) #to make dist the same length as nrow(x)
      parallel::stopCluster(cl)
    }else{
      dist = apply(distCoordinates,1,euc)
      dist = c(dist, NA) #to make dist the same length as nrow(x)
    }
    
    dx = c((distCoordinates[,3] - distCoordinates[,1]),NA) #calculate differences in x (x2 - x1) and put an NA at the end to ensure it is of equal length to the input data
    dy = c((distCoordinates[,4] - distCoordinates[,2]),NA) #calculate differences in y (y2 - y1) and put an NA at the end to ensure it is of equal length to the input data
    
    timesFrame = data.frame(x$dateTime[1:(nrow(x) - 1)], x$dateTime[2:nrow(x)])
    if(parallel == TRUE){
      cl<-parallel::makeCluster(nCores)
      dt = parallel::parApply(cl, timesFrame, 1, timeDifference)
      dt = c(dt, NA)  #to make length(dt) == nrow(x)
      parallel::stopCluster(cl)
    }else{
      dt = apply(timesFrame, 1, timeDifference)
      dt = c(dt, NA)  #to make length(dt) == nrow(x)
    }

    idVec = unique(x$id)
    idSeqVec = x$id
    idVecForDtRemoval = idVec[-length(idVec)] #There's no need to remove the max point of the final id point b/c there's already an NA there.
    for(a in idVecForDtRemoval){ #This loop removes the dt value at points describing the last point in each individual's path
      dx[max(which(idSeqVec == a))] = NA
      dy[max(which(idSeqVec == a))] = NA
      dt[max(which(idSeqVec == a))] = NA
    }
    x$dx = dx
    x$dy = dy
    x$dt = dt
    x$dist<-dist

    dataShift = data.frame(x = x$x, y = x$y, dx = c(NA,x$dx[1:(nrow(x) - 1)]), dy = c(NA,x$dy[1:(nrow(x) - 1)]), dist = c(NA,x$dist[1:(nrow(x) - 1)])) #This is necessary because of the way we calculated/listed these values above. These columns (dist, dx, and dt) refer to the changes a point must make to reach the subsequent point. However, later on in this function, we are not interested in future point alterations. Rather, we need to know how tracked individuals moved during the preceding time step to reach their current point (to determine directionality of movement) #Note that this code shifts values down, but remember that beacuse values are downshifted, the first observation for each id individual will be incorrect. Code below addresses this issue.

    for(b in idVec){ #This code replaces dx, and dy values in the row when a new individual (id) is first observed with NAs to fix the problem noted above. #Note that this loop assumes that the data is sorted by id (i.e., a given id will not repeat in the dataset, after a new id has appeared)
      dataShift[min(which(idSeqVec == b)),] = NA
    }

    if(length(direction) > 0){ #the other x inputs are required, but direction is not
      if(length(direction) == 1 && is.na(match(direction[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
        eta <- x[,match(direction, names(x))]
      }else{ #if length(direction) > 1
        eta <- direction
      }
    }else{ #if length(direction) == 0
      eta <- atan2(dataShift$dy,dataShift$dx)*(180/pi) + 360 #calculate the relative angle of movements given point (x1,y1) lies on the axes' origin. Multiplying by 180/pi converts radians to degrees and adding 360 ensures positive values.
    }
    
    #simplify eta (i.e., make all observations fall between 0 and 359)
    eta360s <- floor(eta/360) # determine how many times individuals turn in complete circles
    eta.simplified <- eta - (eta360s*360) #this reduces direction values to between 0 and 359
    
    #repositionHypoteneuse<-sqrt((LeftRightRepositionLen^2) + (UpDownRepositionLen^2))
    
    if (StartLocation == "UL"){
      ##newCoordinates here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      newCoordinates1 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 0, distV = LeftRightRepositionLen) #repositions the empirical point directly to the East (i.e., right) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the top-right corner of the polygon.
      newCoordinates2 <-translate(x = newCoordinates1$x.adjusted, y = newCoordinates1$y.adjusted, eta = eta.simplified, etaStar = modelOrientation, theta = 270, distV = UpDownRepositionLen) #repositions the adjusted top-right-corner point directly to the South (i.e., down) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the Bottom-right corner of the polygon.
      newCoordinates3 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 270, distV = UpDownRepositionLen) #repositions the empirical point directly to the South (i.e., down) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the Bottom-left corner of the polygon.
    }

    if (StartLocation == "UR"){
      ##newCoordinates here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      newCoordinates1 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 270, distV = UpDownRepositionLen) #repositions the empirical point directly to the South (i.e., down) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the Bottom-right corner of the polygon.
      newCoordinates2 <-translate(x = newCoordinates1$x.adjusted, y = newCoordinates1$y.adjusted, eta = eta.simplified, etaStar = modelOrientation, theta = 180, distV = LeftRightRepositionLen) #repositions the adjusted bottom-right-corner point directly to the West (i.e., left) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the Bottom-left corner of the polygon.
      newCoordinates3 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 180, distV = LeftRightRepositionLen) #repositions the empirical point directly to the West (i.e., left) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the upper-left corner of the polygon.
      
    }

    if (StartLocation == "DR"){
      ##newCoordinates here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      newCoordinates1 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 180, distV = LeftRightRepositionLen) #repositions the empirical point directly to the West (i.e., left) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the bottom-left corner of the polygon.
      newCoordinates2 <-translate(x = newCoordinates1$x.adjusted, y = newCoordinates1$y.adjusted, eta = eta.simplified, etaStar = modelOrientation, theta = 90, distV = UpDownRepositionLen) #repositions the adjusted bottom-left-corner point directly to the North (i.e., up) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the upper-left corner of the polygon.
      newCoordinates3 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 90, distV = UpDownRepositionLen) #repositions the empirical point directly to the North (i.e., up) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the upper-right corner of the polygon.
    }
    
    if (StartLocation == "DL"){
      ##newCoordinates here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      newCoordinates1 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 90, distV = UpDownRepositionLen) #repositions the empirical point directly to the North (i.e., up) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the upper-left corner of the polygon.
      newCoordinates2 <-translate(x = newCoordinates1$x.adjusted, y = newCoordinates1$y.adjusted, eta = eta.simplified, etaStar = modelOrientation, theta = 0, distV = LeftRightRepositionLen) #repositions the adjusted bottom-right-corner point directly to the East (i.e., right) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the upper-right corner of the polygon.
      #newCoordinates2 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 45, distV = repositionHypoteneuse) #repositions the adjusted bottom-right-corner point directly to the East (i.e., right) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the upper-right corner of the polygon.
      newCoordinates3 <-translate(x = dataShift$x, y = dataShift$y, eta = eta.simplified, etaStar = modelOrientation, theta = 0, distV = LeftRightRepositionLen) #repositions the empirical point directly to the East (i.e., right) in a planar model oriented to 90-degrees (i.e., etaStar == 90). Forms the bottom-right corner of the polygon.
    }


    polygonMatrix = data.frame(matrix(nrow = nrow(x), ncol = 27)) #One row for each observation, Columns listed below
    colnames(polygonMatrix) = c("id","cornerPoint1.x","cornerPoint1.y","midPoint1.x","midPoint1.y","cornerPoint2.x","cornerPoint2.y","midPoint2.x","midPoint2.y","cornerPoint3.x","cornerPoint3.y","midPoint3.x","midPoint3.y","cornerPoint4.x","cornerPoint4.y","midPoint4.x","midPoint4.y","centroid.x","centroid.y","startLocation","movementDirection","upDownRepositionLength","leftRightRepositionLength","immob","immobThreshold", "dateTime","dt")

    polygonMatrix$id = x$id
    polygonMatrix$cornerPoint1.x = x$x
    polygonMatrix$cornerPoint1.y = x$y
    polygonMatrix$cornerPoint2.x = newCoordinates1$x.adjusted
    polygonMatrix$cornerPoint2.y = newCoordinates1$y.adjusted
    polygonMatrix$cornerPoint3.x = newCoordinates2$x.adjusted
    polygonMatrix$cornerPoint3.y = newCoordinates2$y.adjusted
    polygonMatrix$cornerPoint4.x = newCoordinates3$x.adjusted
    polygonMatrix$cornerPoint4.y = newCoordinates3$y.adjusted
    polygonMatrix$movementDirection = eta.simplified
    polygonMatrix$startLocation = StartLocation
    polygonMatrix$upDownRepositionLength = UpDownRepositionLen
    polygonMatrix$leftRightRepositionLength = LeftRightRepositionLen
    polygonMatrix$immobThreshold = immobThreshold
    polygonMatrix$dateTime = x$dateTime
    polygonMatrix$dt = x$dt

    polygonMatrix$immob = ifelse(dataShift$dist > immobThreshold, 0, 1) #if the distance individuals moved was less than / equal to the noted immobThreshold, individuals are said to be "immob," and their position will not change relative to their previous one. (i.e., you assume that any observed movement less than immobThreshold was due to errors or miniscule bodily movements (e.g., head shaking) that are not indicative of actual movement.)

    standlist = which(polygonMatrix$immob == 1)

    if ((length(standlist) >= 1)){ #To save processing time and reduce chances of errors, this evaluation will not take place if there is only one observation.

      immobVec = polygonMatrix$immob
      standlistFrame = data.frame(immob = standlist)

      if(parallel == TRUE){
        cl<-parallel::makeCluster(nCores)
        immobFrame = data.frame(data.table::rbindlist(parallel::parApply(cl, standlistFrame, 1, immobAdjustment.polygon, locMatrix = polygonMatrix, immobVec)))
        parallel::stopCluster(cl)
      }else{
        immobFrame = data.frame(data.table::rbindlist(apply(standlistFrame, 1, immobAdjustment.polygon, locMatrix = polygonMatrix, immobVec)))
      }
      if(nrow(immobFrame) > 0){
        polygonMatrix[immobFrame$replaceRow,match("cornerPoint2.x", names(polygonMatrix))] = immobFrame$replacePoint2.x
        polygonMatrix[immobFrame$replaceRow,match("cornerPoint2.y", names(polygonMatrix))] = immobFrame$replacePoint2.y
        polygonMatrix[immobFrame$replaceRow,match("cornerPoint3.x", names(polygonMatrix))] = immobFrame$replacePoint3.x
        polygonMatrix[immobFrame$replaceRow,match("cornerPoint3.y", names(polygonMatrix))] = immobFrame$replacePoint3.y
        polygonMatrix[immobFrame$replaceRow,match("cornerPoint4.x", names(polygonMatrix))] = immobFrame$replacePoint4.x
        polygonMatrix[immobFrame$replaceRow,match("cornerPoint4.y", names(polygonMatrix))] = immobFrame$replacePoint4.y
      }
    }

    if(CenterPoint == TRUE){
      polygonMatrix$centroid.x = ((polygonMatrix$cornerPoint1.x + polygonMatrix$cornerPoint3.x)/2)
      polygonMatrix$centroid.y = ((polygonMatrix$cornerPoint1.y + polygonMatrix$cornerPoint3.y)/2)
    }else{polygonMatrix <- polygonMatrix[,-c(match("centroid.x", names(polygonMatrix)),match("centroid.y", names(polygonMatrix)))]}

    if(MidPoints == TRUE){
      polygonMatrix$midPoint1.x = ((polygonMatrix$cornerPoint1.x + polygonMatrix$cornerPoint2.x)/2)
      polygonMatrix$midPoint1.y = ((polygonMatrix$cornerPoint1.y + polygonMatrix$cornerPoint2.y)/2)
      polygonMatrix$midPoint2.x = ((polygonMatrix$cornerPoint2.x + polygonMatrix$cornerPoint3.x)/2)
      polygonMatrix$midPoint2.y = ((polygonMatrix$cornerPoint2.y + polygonMatrix$cornerPoint3.y)/2)
      polygonMatrix$midPoint3.x = ((polygonMatrix$cornerPoint3.x + polygonMatrix$cornerPoint4.x)/2)
      polygonMatrix$midPoint3.y = ((polygonMatrix$cornerPoint3.y + polygonMatrix$cornerPoint4.y)/2)
      polygonMatrix$midPoint4.x = ((polygonMatrix$cornerPoint4.x + polygonMatrix$cornerPoint1.x)/2)
      polygonMatrix$midPoint4.y = ((polygonMatrix$cornerPoint4.y + polygonMatrix$cornerPoint1.y)/2)
    }else{polygonMatrix <- polygonMatrix[,-c(match("midPoint1.x", names(polygonMatrix)),match("midPoint1.y", names(polygonMatrix)),match("midPoint2.x", names(polygonMatrix)),match("midPoint2.y", names(polygonMatrix)),match("midPoint3.x", names(polygonMatrix)),match("midPoint3.y", names(polygonMatrix)),match("midPoint4.x", names(polygonMatrix)),match("midPoint4.y", names(polygonMatrix)))]}

    return(polygonMatrix)
  }

  list.breaker<-function(x,y,id, dateTime, point.x, point.y, direction, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, immobThreshold, parallel, modelOrientation, nCores){
    input<- data.frame(y[unname(unlist(x[1]))])
    newPolygons<-poly.generator(input, id, dateTime, point.x, point.y, direction, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, immobThreshold, parallel, modelOrientation, nCores)
    return(newPolygons)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.poly <- apply(breakFrame, 1, list.breaker, y = x,id, dateTime, point.x, point.y, direction, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, immobThreshold, parallel, modelOrientation, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.

    return(list.poly)

  }else{ #if x is a dataFrame
    frame.poly<- poly.generator(x,id, dateTime, point.x, point.y, direction, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, immobThreshold, parallel, modelOrientation, nCores)

    return(frame.poly)
  }
}
