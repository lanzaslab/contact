#' Calculate Distances Between All Individuals
#'
#' Calculates the distance between all tracked individuals at a given timestep.
#'    Users can choose whether to calculate distances based on a single point, 
#'    or polygons representative of individuals' locations. If individuals set 
#'    dataType == "Polygon", the distance matrix reported describes the 
#'    shortest distances between polygons' edges (Note that the 
#'    rgeos::gDistance function is used to obtain these distances).
#'
#' If dataType == "Point," users have the option of setting lonlat == TRUE 
#'    (by default lonlat == FALSE). lonlat is a logical argument that tells the
#'    function to calculate the distance between points on the WGS ellipsoid 
#'    (if lonlat == TRUE), or on a plane (lonlat == FALSE) 
#'    (see raster::pointDistance). If lonlat == TRUE, coordinates should be in 
#'    degrees. Otherwise, coordinates should represent planar ('Euclidean') 
#'    space (e.g. units of meters).This function is not currently able to 
#'    calculate distances between polygons on the WGS ellipsoid (i.e., if 
#'    dataType == "Polygon," lonlat must = FALSE). We aim to address this issue
#'    in future versions.
#'
#' Note that if inputting a separate matrix/dataframe with polygon xy 
#'    coordinates (poly.xy), coordinates must be arranged in the format of 
#'    those in referencePointToPolygon outputs (i.e., col1 = point1.x, 
#'    col2 = point1.y, col3 =point2.x, col4 = point2.y, etc., with points 
#'    listed in a clockwise (or counter-clockwise) order).
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
#' @param poly.xy Columns within x denoting polygon xy-coordinates. Polygon 
#'    coordinates must be arranged in the format of those in 
#'    referencePointToPolygon output. Defaults to NULL.
#' @param elev Vector of length nrow(data.frame(x)) or singular character data,
#'    detailing the relevant colname in x, that denotes vertical positioning of
#'    each individual in 3D space (e.g., elevation).  If argument != NULL, 
#'    relative vertical positioning will be incorporated into distance 
#'    calculations. Defaults to NULL.
#' @param parallel Logical. If TRUE, sub-functions within the dist2All wrapper 
#'    will be parallelized. Note that this can significantly speed up 
#'    processing of relatively small data sets, but may cause R to crash due 
#'    to lack of available memory when attempting to process large datasets. 
#'    Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to the maximum number of cores available
#'    (i.e., parallel::detectCores()).
#' @param dataType Character string refering to the type of real-time-location 
#'    data presented in x, taking values of "Point" or "Polygon." If 
#'    argument == "Point," individuals' locations are drawn from point.x and 
#'    point.y. If argument == "Polygon," individuals' locations are drawn from 
#'    poly.xy. Defaults to "Point."
#' @param lonlat Logical. If TRUE, point.x and point.y contain geographic 
#'    coordinates (i.e., longitude and lattitude). If FALSE, point.x and 
#'    point.y contain planar coordinates. Defaults to FALSE.
#' @param numVertices Numerical. If dataType == "Polygon," users must specify 
#'    the number of vertices contained in each polygon. Defaults to 4. Note: 
#'    all polygons must contain the same number of vertices.
#' @keywords data-processing polygon point location planar GRC
#' @return Returns a data frame (or list of data frames if \code{x} is a 
#'    list of data frames) with the following columns:
#'    
#'    \item{dateTime}{The unique date-time information corresponding to when
#'    tracked individuals were observed in \code{x}.}
#'    \item{totalIndividuals}{The total number of individuals observed at least
#'    one time within \code{x}.}
#'    \item{individualsAtTimestep}{The number of individuals in \code{x}  
#'    observed at the timepoint described in the \code{dateTime} column.}    
#'    \item{id}{The unique ID of a tracked individual for which we will 
#'    evaluate distances to all other individuals observed in \code{x}.}
#'    \item{dist.to.indiv...}{The observed distance between the individual 
#'    described in the \code{id} column to every other individual observed 
#'    at specific timepoints.}
#' @export
#' @examples
#' data(calves)
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'    time = calves$time) 
#'    
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth locations to 5-min fix intervals.  
#'
#' calves.dist<-dist2All_df(x = calves.agg, parallel = FALSE, dataType = "Point", 
#'    lonlat = FALSE) #calculate distance between all individuals at each timepoint.
#' 

dist2All_df<-function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, poly.xy = NULL, elev = NULL, parallel = FALSE, nCores = parallel::detectCores(), dataType = "Point", lonlat = FALSE, numVertices = 4){

  create.distFrame<- function(x,distMat, indivSeq, idSeq, timestepIndivSeq,time){

    dist = data.frame(matrix(ncol = (length(indivSeq) + 4), nrow = 1))
    colnames(dist) = c("dateTime","totalIndividuals","individualsAtTimestep","id", paste("dist.to.indiv_", idSeq, sep = ""))
    dist$dateTime = time
    dist$totalIndividuals = length(indivSeq)
    dist$individualsAtTimestep = length(timestepIndivSeq)
    dist$id = unlist(unname(x[1]))
    vecID = 1
    column<-as.numeric(unlist(unname(x[2]))) #The column identifier must be numeric. Without the as.numeric call here, apply functions will return the error "Error in distMat[vecID, column] : no 'dimnames' attribute for array" later on.
    col.fill = NULL

    for(j in 1:length(indivSeq)){

      if(isTRUE(indivSeq[j]%in%timestepIndivSeq) == TRUE){

        if(unlist(unname(x[3])) != indivSeq[j]){ #This if statement makes sure that individuals are not reported as being a certain distance from themselves. If the id = indivSeq[j], dist = NA
          col.fill = c(col.fill,distMat[vecID,column]) #column must be must be numeric

        }else{
          col.fill = c(col.fill,NA)
        }
        vecID = (vecID + 1)

      }else{
        col.fill = c(col.fill,NA)
      }
    }

    dist[1,5:ncol(dist)] = col.fill
    return(dist)
  }
  list.breaker1<-function(x,y,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){ #calculates dist.all or distToArea for randomized input
    input<- data.frame(y[unname(unlist(x[1]))])
    dist.all<- dist.generator1(input,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores)
    return(dist.all)
  }
  dist.generator1<-function(x,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){
    idVec1=NULL #added in the case that idVec1 isn't created when x isn't specified
    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
      if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
        originTab = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime)
      }
      if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){

        originTab = data.frame(matrix(ncol = 0, nrow = length(id)))
        originTab$id = id
        colnames(poly.xy)[seq(1,(ncol(poly.xy) - 1),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".x", sep = "")
        colnames(poly.xy)[seq(2,ncol(poly.xy),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".y", sep = "")
        dateFrame = data.frame(dateTime = dateTime)
        bindlist = list(originTab,poly.xy,dateFrame)
        originTab = data.frame(do.call("cbind", bindlist))
      }
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

      if(dataType == "point" || dataType == "Point" || dataType == "POINT"){

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
        xyFrame1<- data.frame(x = x$x, y = x$y)
      }

      if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){

        if(length(poly.xy) > 0){
          if(length(poly.xy) == (numVertices*2) & length(which(is.na(match(poly.xy, names(x)))) == TRUE) == 0){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than poly.xy being a matrix/dataframe of length(nrow(x)), it may be necessary to designate the colnames for intended coordinate values (i.e., if the xy-coordinate values in different list entries are different)
            xyFrame1<-x[,match(poly.xy,names(x))]
          }else{
            #x[,2:(1 + length(poly.xy))] = poly.xy
            xyFrame1<- data.frame(poly.xy)
          }
        }else{ #if length(poly.xy == 0)
          xyFrame1 <-  x[,(match("point1.x", names(x))):((2*numVertices) + (match("point1.x", names(x)) -1))] #if there is no poly.xy input, the code assumes the input is output from referencePointToPolygon function, and therefore, the first point of interest would be "point1.x"
        }
        colnames(xyFrame1)[seq(1,(numVertices*2),2)] = paste("point",seq(1,numVertices,1),".x", sep = "")
        colnames(xyFrame1)[seq(2,((numVertices*2) + 1),2)] = paste("point",seq(1,numVertices,1),".y", sep = "")
      }

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

    #02/01/2019 - Here I'll remove the need for integer ids by creating a separate integ.ID column in originTab
    if(length(idVec1)== 0) { #added to use idVec1 if not created above 
      idVec1 <- originTab$id
      } 
    idVec2=unique(originTab$id) #Added in the case x is not supplied
    originTab$integ.ID<-NA

    for(a in 1:length(idVec2)){
      originTab$integ.ID[which(idVec1 == idVec2[a])] <-as.integer(a)
    }

    originTab$dateTime = as.character(originTab$dateTime)
    
    elev.calc <-FALSE #logical; tells us later if we need to include elevation calculations.
    if(length(elev) > 0){
      elev.calc <-TRUE 
      if(length(elev) == 1 && is.na(match(elev[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than elev being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "elev" values (i.e., if the elevs in different list entries are different)
        originTab$elev <- x[,match(elev, names(x))]
      }else{ #if length(elev) > 1
        originTab$elev = elev
      }
    }

    if(dataType == "point" || dataType == "Point" || dataType == "POINT"){

      originTab<-originTab[order(originTab$id, originTab$dateTime),]
      indivSeq = unique(originTab$integ.ID)
      idSeq = unique(originTab$id)
      dateTimeFrame = data.frame(dateTime = unique(originTab$dateTime))
      dist.measurement = lonlat

      dist.process.point <- function(x, originTab, indivSeq, idSeq, dist.measurement, elev.calc){

        time = unlist(unname(x[1]))
        timestep = originTab[which(originTab$dateTime == time),]
        timestepIndivSeq.integ = unique(timestep$integ.ID)
        distMat = raster::pointDistance(timestep[,c(match("x",names(timestep)),match("y",names(timestep)))], timestep[,c(match("x",names(timestep)),match("y",names(timestep)))], lonlat = dist.measurement, allpairs = TRUE)

        if(is.matrix(distMat) == TRUE){ #If there's only one record in timestep, distMat will be a vector, not a matrix, and will only denote the distance between an individual and itself (which isn't useful)
          if(elev.calc == TRUE){
            elevDifference = function(x, elevVec){ #This function creates a matrix showing the difference in elevation between all observed points in the data set (timestep) 
              e1 <- unname(unlist(x[1]))
              de <- as.integer(abs(e1 - elevVec)) 
              return(de)
            }
            elevVec <- unname(unlist(timestep$elev))
            elevFrame <- data.frame(elevVec)
            elev.dif<-apply(elevFrame, 1, elevDifference, elevVec)
            distMat<-distMat + elev.dif #adds the difference in elevations to the linear distances previously calculated. Note that the elev distance units must be in meters, as that is the unit of distMat values. 
          }
          timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1), timestepIndivSeq.integ) #Note that colnum will not necessarily be the same as timestepIndivSeq.integ because originTab has been previously ordered
          distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, idSeq, timestepIndivSeq.integ,time)))
          return(distTab)
        }
      }


      if (parallel == TRUE){
        cl<-parallel::makeCluster(nCores)
        distTab<-parallel::parApply(cl, dateTimeFrame, 1, dist.process.point, originTab, indivSeq,idSeq, dist.measurement, elev.calc)
        parallel::stopCluster(cl)
      }else{
        distTab = apply(dateTimeFrame, 1, dist.process.point, originTab, indivSeq,idSeq, dist.measurement, elev.calc)
      }
      dist.all = data.frame(data.table::rbindlist(distTab))
    }

    if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){

      originTab<-originTab[order(originTab$id, originTab$dateTime),]
      naVec<-which(is.na(originTab[,match("point2.x", names(originTab))]) == TRUE) #the referencePointtoPolygon function will create some observations that are not complete polygons (i.e., only the point1 coordinates are recorded). This identifies those observations, so that they may be removed. If they are not removed, they will cause errors later.
      if(length(naVec) > 0){
        originTab <- originTab[-naVec,]
      }
      indivSeq = unique(originTab$integ.ID)
      idSeq <- unique(originTab$id)
      dateTimeFrame = data.frame(dateTime = unique(originTab$dateTime))

      dist.process.poly <- function(x, originTab, indivSeq, idSeq, numVertices, elev.calc){
        create.poly1<-function(x,y, numVertices){
          polygon = unlist(c(y[unlist(unname(x[1])),c(2:(1 + (numVertices*2)),2:3)])) #note that when making a Spatial Polygon you need to close the polygon by repeating the first point at the end (Note that this also means you need to add a "+ 1" to nrow = numVertices below)
          spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices + 1), ncol = 2, byrow = TRUE))),y$integ.ID[unlist(unname(x[1]))])
          return(spatPoly)
        }

        time = unlist(unname(x[1]))
        timestep = originTab[which(originTab$dateTime == time),]

        if(nrow(timestep) > 1){ #If there's only one record in timestep, any distMat produced will only denote the distance between an individual and itself (which isn't useful), so it's more time efficient to skip forward nrow(timestep <= 1)

          timestepIndivSeq.integ = unique(timestep$integ.ID)
          makePolyFrame<-data.frame(seq(1,length(timestepIndivSeq.integ),1))
          spatialPolygons <- apply(makePolyFrame,1,create.poly1,timestep,numVertices)
          sPolys = sp::SpatialPolygons(spatialPolygons,timestepIndivSeq.integ) #note that the second part of this argument must be an integer. Otherwise, it will return the following error: Error: is.integer(pO) is not TRUE
          distMat = rgeos::gDistance(sPolys, byid = TRUE)
          
          if(elev.calc == TRUE){
            elevDifference = function(x, elevVec){ #This function creates a matrix showing the difference in elevation between all observed points in the data set (timestep) 
              e1 <- unname(unlist(x[1]))
              de <- as.integer(abs(e1 - elevVec)) 
              return(de)
            }
            elevVec <- unname(unlist(timestep$elev))
            elevFrame <- data.frame(elevVec)
            elev.dif<-apply(elevFrame, 1, elevDifference, elevVec)
            distMat<-distMat + elev.dif #adds the difference in elevations to the linear distances previously calculated. Note that the elev distance units must be in meters, as that is the unit of distMat values. 
          }
          
          timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1), timestepIndivSeq.integ) #Note that colnum will not necessarily be the same as timestepIndivSeq.integ because originTab has been previously ordered
          distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, idSeq, timestepIndivSeq.integ,time)))
          return(distTab)
        }
      }

      if (parallel == TRUE){
        cl<-parallel::makeCluster(nCores)
        distTab<-parallel::parApply(cl, dateTimeFrame, 1, dist.process.poly, originTab, indivSeq, idSeq, numVertices, elev.calc)
        parallel::stopCluster(cl)
      }else{
        distTab = apply(dateTimeFrame, 1, dist.process.poly, originTab, indivSeq, idSeq, numVertices, elev.calc)
      }
      dist.all = data.frame(data.table::rbindlist(distTab))
    }
    return(dist.all)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.dist <- apply(breakFrame, 1, list.breaker1,y = x,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.
    return(list.dist)

  }else{ #if(is.list(x) == FALSE)

    frame.dist<- dist.generator1(x,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores)

    return(frame.dist)
  }
}
