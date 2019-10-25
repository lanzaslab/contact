#' Calculate Distances Between Individuals and Fixed Points/Polygons
#'
#' Calculate distances (either planar or great circle - see dist.all) between 
#'    each individual, reported in x, and a fixed point(s)/polygon(s), reported
#'    in y, at each timestep.
#'
#' Polygon coordinates (in both x and y inputs) must be arranged in the format 
#'    of those in referencePointToPolygon outputs (i.e., col1 = point1.x, 
#'    col2 = point1.y, col3 =point2.x, col4 = point2.y, etc., with points 
#'    listed in a clockwise (or counter-clockwise) order).
#'    
#' This variant of dist2Area requires x and y inputs to be non-shapefile data.
#' @param x Data frame or list of data frames containing real-time-location 
#'    data for individuals. 
#' @param y Data frame or list of data frames describing fixed-area 
#'    polygons/points for which we will calculate distances relative to tracked
#'    individuals at all time steps. Polygons contained within the same data 
#'    frame must have the same number of vertices.
#' @param x.id Vector of length nrow(data.frame(x)) or singular character data,
#'    detailing the relevant colname in x, that denotes what unique ids for 
#'    tracked individuals will be used. If argument == NULL, the function 
#'    assumes a column with the colname "id" exists in x. Defaults to NULL.
#' @param y.id Vector of length sum(nrow(data.frame(y[1:length(y)]))) or 
#'    singular character data, detailing the relevant colname in y, that 
#'    denotes what unique ids for fixed-area polygons/points will be used. If 
#'    argument == NULL, the function assumes a column with the colname "id" may
#'    exist in y. If such a column does exist, fixed-area polygons will be 
#'    assigned unique ids based on values in this column. If no such column 
#'    exists, fixed-area polygons/points will be assigned sequential numbers as
#'    unique identifiers. Defaults to NULL.
#' @param point.x Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what planar-x or 
#'    longitude coordinate information will be used. If argument == NULL, the 
#'    function assumes a column with the colname "x" exists in x. Defaults to 
#'    NULL.
#' @param point.y Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what planar-y 
#'    or lattitude coordinate information will be used. If argument == NULL, 
#'    the function assumes a column with the colname "y" exists in x. 
#'    Defaults to NULL.
#' @param dateTime Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what dateTime 
#'    information will be used. If argument == NULL, the function assumes a 
#'    column with the colname "dateTime" exists in x. Defaults to NULL.
#' @param poly.xy Columns within x denoting polygon xy-coordinates. Polygon 
#'    coordinates must be arranged in the format of those in 
#'    referencePointToPolygon output. Defaults to NULL.
#' @param parallel Logical. If TRUE, sub-functions within the dist2Area_df 
#'    wrapper will be parallelized. Note that this can significantly speed up 
#'    processing of relatively small data sets, but may cause R to crash due to
#'    lack of available memory when attempting to process large datasets. 
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
#'    the number of vertices contained in each polygon described in x. 
#'    Defaults to 4. Note: all polygons must contain the same number of 
#'    vertices.
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
#'    \item{dist.to...}{The observed distance between the individual 
#'    described in the \code{id} column to every each polygon/fixed location}
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
#'                   y = c(62.44815, 62.73341, 61.93864, 61.67411)) #delineate water polygon
#'
#' water_poly<-data.frame(matrix(ncol = 8, nrow = 1)) #make coordinates to dist2Area specifications
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
#'   lonlat = FALSE, numVertices = NULL) 

dist2Area_df<-function(x = NULL, y = NULL, x.id = NULL, y.id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, poly.xy = NULL, parallel = FALSE, nCores = parallel::detectCores(), dataType = "Point", lonlat = FALSE, numVertices = 4){
  
  create.distFrame<- function(x,distMat, indivSeq, timestepIndivSeq,time, origin.y){
    dist = data.frame(matrix(ncol = (nrow(origin.y) + 4), nrow = 1))
    colnames(dist) = c("dateTime","totalIndividuals","individualsAtTimestep","id", paste("dist.to.", origin.y[,1], sep = ""))
    dist$dateTime = time
    dist$totalIndividuals = length(indivSeq)
    dist$individualsAtTimestep = length(timestepIndivSeq)
    dist$id = unname(unlist(x[1]))
    #vecID = 1
    col.fill = NULL
    col.fill = unname(unlist(distMat[,unname(unlist(x[2]))])) #There's no need for the forloop used in the dist.all function b/c there will never be a distance given that represents individuals' distance to themselves.
    dist[1,5:ncol(dist)] = col.fill
    return(dist)
  }
  
  list.breaker4<-function(x, z, y, x.id, y.id, dateTime, point.x, point.y, poly.xy, parallel, dataType, lonlat, numVertices, nCores){
    input<- data.frame(z[unname(unlist(x[1]))])
    dist.all<- dist.generator2(input, y, x.id, y.id, dateTime, point.x, point.y, poly.xy, parallel, dataType, lonlat, numVertices, nCores)
    return(dist.all)
  }
  dist.generator2<-function(x, y, x.id, y.id, dateTime, point.x, point.y, poly.xy, parallel, dataType, lonlat, numVertices, nCores){
    id<-NULL #bind this variable to a local object so that R CMD check doesn't flag it.
    
    idVec1=NULL #added in the case that idVec1 isn't created when x isn't specified
    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
      if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
        originTab = data.frame(id = x.id, x = point.x, y = point.y, dateTime = dateTime)
      }
      if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
        originTab = data.frame(matrix(ncol = 0, nrow = length(id)))
        originTab$id = x.id
        colnames(poly.xy)[seq(1,(ncol(poly.xy) - 1),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".x", sep = "")
        colnames(poly.xy)[seq(2,ncol(poly.xy),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".y", sep = "")
        dateFrame = data.frame(dateTime = dateTime)
        bindlist = list(originTab,poly.xy,dateFrame)
        originTab = data.frame(do.call("cbind", bindlist))
      }
    }
    
    if(length(x) > 0){ #for some reason using an "else" statement would always result in an originTab table with 0 records...
      if(length(x.id) > 0){ 
        if(length(x.id) == 1 & is.na(match(x.id[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than id being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
          x$id <- x[,match(x.id, names(x))]
        }else{ #if length(id) > 1
          x$id = x.id
        }
      }
      idVec1 <- x$id
      
      if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
        if(length(point.x) > 0){
          if(length(point.x) == 1 & is.na(match(point.x[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
            x$x <- x[,match(point.x, names(x))]
          }else{ #if length(point.x) > 1
            x$x = point.x
          }
        }
        if(length(point.y) > 0){
          if(length(point.y) == 1 & is.na(match(point.y[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
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
            xyFrame1<- data.frame(poly.xy)
          }
        }else{ #if length(poly.xy == 0)
          xyFrame1 <-  x[,(match("point1.x", names(x))):((2*numVertices) + (match("point1.x", names(x)) -1))] #if there is no poly.xy input, the code assumes the input is output from referencePointToPolygon function, and therefore, the first point of interest would be "point1.x"
        }
        colnames(xyFrame1)[seq(1,(numVertices*2),2)] = paste("point",seq(1,numVertices,1),".x", sep = "")
        colnames(xyFrame1)[seq(2,((numVertices*2) + 1),2)] = paste("point",seq(1,numVertices,1),".y", sep = "")
      }
      
      if(length(dateTime) > 0){
        if(length(dateTime) == 1 & is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
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
    if(length(idVec1)== 0) {  idVec1 <- originTab$id} #added to use idVec1 if not created above 
    idVec2=unique(originTab$id)                      #Added in the case x is not supplied
    #idVec2<-unique(idVec1)
    
    originTab$integ.ID<-NA
    
    for(a in 1:length(idVec2)){
      originTab$integ.ID[which(idVec1 == idVec2[a])] <-as.integer(a)
    }  
    originTab$dateTime = as.character(originTab$dateTime)
    
    if(is.data.frame(y) == FALSE & is.list(y) == TRUE){ #added 02/06/2019, here we determine if y is a list.
      frameBreaker<- function(x,y){
        brokenFrame<- y[unlist(unname(x[1])),]
        return(brokenFrame)
      }
      growPoly<-function(x, vertex.colnum){
        if(is.na(match("id", names(x))) == FALSE){ #so, if there is an "id" column in x
          area.id<- x[1,match("id", names(x))]
          x<-x[,-match("id", names(x))] #we remove the id column, so it doesn't affect vertex calculations.
        }else{ #if there's no "id" column
          area.id = NULL
        }
        cols.needed<- vertex.colnum - ncol(x) #determines how many columns are needed to increase the number of vertices to reach the desired number of vertices.
        expanded.x <-x
        if(cols.needed > 0){ #so, if the polygon needs more vertices to match the maximum number detected.
          expanded.x[,(ncol(expanded.x) + 1):(ncol(expanded.x) + cols.needed)] <- expanded.x[,1:cols.needed] #this loops the polygon around the first recorded vertices. 
        }
        if(length(area.id) >0){ #if there was an id column
          cbindlist<-list(area.id,expanded.x)
          output<-do.call("cbind", cbindlist)
        }else{ #if there was no id column
          output <- expanded.x
        }
        return(output)
      }
      colLengths<-NULL
      areas<-list()
      for(a in 1:length(y)){
        areaFrame<-data.frame(y[a])
        colnum<- ifelse(is.na(match("id", names(areaFrame))) == FALSE, ncol(areaFrame) -1, ncol(areaFrame)) #if there IS an id column, it will not be counted towards the column number, which is representative of the number of xycoordinates.
        colLengths<-c(colLengths, colnum) #we compile a sequence ncols, so we can determine the maximum number of vertices the input areas have. 
        if(nrow(areaFrame) > 1){ #if multiple fixed areas are represented in a list entry
          rowSeq<-seq(1,nrow(areaFrame),1)
          rowSeqFrame<-data.frame(rowSeq)
          polys <- apply(rowSeqFrame,1,frameBreaker,areaFrame)
          for(b in 1:length(polys)){
            areas<-c(areas,polys[b])
          }
        }else{#if areaFrame only has one row (representing a single fixed area)
          areas<-c(areas, list(areaFrame))
        }
      }
      #Now we have the areas list, which we can input into the growPoly function.
      vertex.colnum<-max(colLengths) #the maximum number of columns (representative of the number of vertices) in the area set.
      new.y<-lapply(areas, growPoly, vertex.colnum)
      y <- data.frame(data.table::rbindlist(new.y)) #here we remake y as a dataframe containing polygons all with the same number of vertices. After this step, the function can move forward as if y was not a list. 
    }
    
    ycols = ncol(y)
    if(is.na(match("id", names(y))) == TRUE){ #if there is no "id" column, then all columns in y relate to positions of points
      numVertices.y = ycols/2
      if(length(y.id) > 0){
        id.y = y.id
      }else{
        id.y = seq(1,nrow(y),1)
      }
      
    }else{#if there is an "id" column
      numVertices.y = (ycols - 1)/2
      id.y = y$id
      if(length(y.id) > 0){
        id.y = y.id
      }
      y<- y[,-match("id", names(y))]
    }
    
    ids <- data.frame(orig.id = id.y, new.id = as.integer(seq(1,nrow(y),1)))
    origin.y <- do.call("cbind",list(ids,y)) #so, now we have a table that is arranged like: id1,id2,xCoord1, yCoord1..., xCoordnumVertices.y,yCoordnumVertices.y
    
    if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
      
      originTab<-originTab[order(originTab$id, originTab$dateTime),]
      indivSeq = unique(originTab$id)
      dateTimeFrame = data.frame(dateTime = unique(originTab$dateTime))
      dist.measurement = lonlat
      
      dist.process.point <- function(x, originTab, indivSeq, dist.measurement, origin.y, numVertices.y){
        create.poly2<-function(x,y, numVertices.y){
          polygon = unlist(c(y[unlist(unname(x[1])),c(3:(2 + (numVertices.y*2)),3:4)])) #note that when making a Spatial Polygon you need to close the polygon by repeating the first point at the end (Note that this also means you need to add a "+ 1" to nrow = numVertices below)
          spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices.y + 1), ncol = 2, byrow = TRUE))),y[unlist(unname(x[1])),2])
          return(spatPoly)
        }
        
        time = unlist(unname(x[1]))
        timestep = originTab[which(originTab$dateTime == time),]
        timestepIndivSeq = unique(timestep$id)
        
        if(ncol(origin.y) - 2 == 2){ #i.e., if the fixed coordinates in y only represent single points (not polygons)
          distMat = raster::pointDistance(timestep[,c(match("x",names(timestep)),match("y",names(timestep)))], origin.y[,3:ncol(origin.y)], lonlat = dist.measurement, allpairs = TRUE)
          if(is.matrix(distMat) == FALSE){ #if there's only 1 fixed point, the pointDistance function will create a named vector rather than a matrix. This turns distMat into a matrix with one column for each tracked individual
            distMat = matrix(distMat, ncol = length(timestepIndivSeq))
          }
        }else{
          xy = timestep[,c(match("x",names(timestep)), match("y",names(timestep)))]
          rownames(xy) <- timestepIndivSeq
          spatPoints =sp::SpatialPoints(xy)
          makePolyFrame<-data.frame(seq(1,nrow(origin.y),1))
          spatialPolygons <- apply(makePolyFrame,1,create.poly2,origin.y,numVertices.y)
          sPolys = sp::SpatialPolygons(spatialPolygons,as.integer(origin.y[,2])) #note that the second part of this argument must be an integer. Otherwise, it will return the following error: Error: is.integer(pO) is not TRUE
          distMat = rgeos::gDistance(spatPoints,sPolys, byid = TRUE) #columns are indivIDs, rows are fixed areas (this is how distances are presented in raster::pointDistance as well)
        }
        distMat = data.frame(distMat)
        timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(timestepIndivSeq),1))
        distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, timestepIndivSeq,time, origin.y))) 
        return(distTab)
      }
      
      if (parallel == TRUE){
        cl<-parallel::makeCluster(nCores)
        distTab<-parallel::parApply(cl, dateTimeFrame, 1, dist.process.point, originTab, indivSeq, dist.measurement, origin.y, numVertices.y)
        parallel::stopCluster(cl)
      }else{
        distTab = apply(dateTimeFrame, 1, dist.process.point, originTab, indivSeq, dist.measurement, origin.y, numVertices.y)	
      }
      dist.all = data.frame(data.table::rbindlist(distTab))
    }
    
    if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
      originTab<-originTab[order(originTab$id, originTab$dateTime),]
      naVec<-which(is.na(originTab[,match("point2.x", names(originTab))]) == TRUE) #the referencePointtoPolygon function will create some observations that are not complete polygons (i.e., only the point1 coordinates are recorded). This identifies those observations, so that they may be removed. If they are not removed, they will cause errors later.
      if(length(naVec) > 0){
        originTab <- originTab[-naVec,]
      }
      indivSeq = unique(originTab$id)
      dateTimeFrame = data.frame(dateTime = unique(originTab$dateTime))
      
      dist.process.poly <- function(x, originTab, indivSeq, origin.y, numVertices, numVertices.y){
        create.poly1<-function(x,y, numVertices){
          polygon = unlist(c(y[unlist(unname(x[1])),c(2:(1 + (numVertices*2)),2:3)])) #note that when making a Spatial Polygon you need to close the polygon by repeating the first point at the end (Note that this also means you need to add a "+ 1" to nrow = numVertices below)
          #		spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices + 1), ncol = 2, byrow = TRUE))),y$id[unlist(unname(x[1]))])
          spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices + 1), ncol = 2, byrow = TRUE))),y$integ.ID[unlist(unname(x[1]))])
          return(spatPoly)
        }
        create.poly2<-function(x,y, numVertices.y){
          polygon = unlist(c(y[unlist(unname(x[1])),c(3:(2 + (numVertices.y*2)),3:4)])) #note that when making a Spatial Polygon you need to close the polygon by repeating the first point at the end (Note that this also means you need to add a "+ 1" to nrow = numVertices below)
          spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices.y + 1), ncol = 2, byrow = TRUE))),y[unlist(unname(x[1])),2])
          return(spatPoly)
        }
        time = unlist(unname(x[1]))
        timestep = originTab[which(originTab$dateTime == time),]
        timestepIndivSeq.integ = unique(timestep$integ.ID)
        timestepIndivSeq = unique(timestep$id)
        makePolyFrame1<-data.frame(seq(1,length(timestepIndivSeq.integ),1))
        spatialPolygons1 <- apply(makePolyFrame1,1,create.poly1,timestep,numVertices)
        sPolys1 <- sp::SpatialPolygons(spatialPolygons1,as.integer(timestepIndivSeq.integ)) #note that the second part of this argument must be an integer. Otherwise, it will return the following error: Error: is.integer(pO) is not TRUE
        
        if(ncol(origin.y) - 2 == 2){ #i.e., if the fixed coordinates in y only represent single points (not polygons)
          xy = origin.y[,3:4]
          rownames(xy) <- origin.y[,2]
          spatPoints =sp::SpatialPoints(xy)
          distMat = rgeos::gDistance(sPolys1,spatPoints, byid = TRUE)
        }else{ #if y has more than only one point
          makePolyFrame2<-data.frame(seq(1,nrow(origin.y),1))
          spatialPolygons2 <- apply(makePolyFrame2,1,create.poly2,origin.y,numVertices.y)
          sPolys2 = sp::SpatialPolygons(spatialPolygons2,as.integer(origin.y[,2])) #note that the second part of this argument must be an integer. Otherwise, it will return the following error: Error: is.integer(pO) is not TRUE
          distMat = rgeos::gDistance(sPolys1,sPolys2, byid = TRUE)
        }
        distMat = data.frame(distMat)
        timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1))

        distTab <- data.frame(NULL)
        
        for(i in 1:nrow(timestepIndivSeqFrame)){ #This for-loop is pretty much just the create.distFrame function (see note immediately above).
          x<- timestepIndivSeqFrame[i,]
          dist = data.frame(matrix(ncol = (nrow(origin.y) + 4), nrow = 1))
          colnames(dist) = c("dateTime","totalIndividuals","individualsAtTimestep","id", paste("dist.to.", origin.y[,1], sep = ""))
          dist$dateTime = time
          dist$totalIndividuals = length(indivSeq)
          dist$individualsAtTimestep = length(timestepIndivSeq)
          dist$id = unname(unlist(x[1]))
          #vecID = 1
          col.fill = NULL
          col.fill = unname(unlist(distMat[,unname(unlist(x[2]))])) #There's no need for the forloop used in the dist.all function b/c there will never be a distance given that represents individuals' distance to themselves.
          dist[1,5:ncol(dist)] = col.fill
          
          distTab<-data.frame(data.table::rbindlist(list(distTab,dist)))
        }
        return(distTab)
      }
      
      if (parallel == TRUE){
        cl<-parallel::makeCluster(nCores)
        distTab<-parallel::parApply(cl, dateTimeFrame, 1, dist.process.poly, originTab, indivSeq, origin.y, numVertices, numVertices.y)
        parallel::stopCluster(cl)
      }else{
        distTab = apply(dateTimeFrame, 1, dist.process.poly, originTab, indivSeq, origin.y, numVertices, numVertices.y)	
      }
      dist.all = data.frame(data.table::rbindlist(distTab))
    }
    return(dist.all)
  }
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.dist <- apply(breakFrame, 1, list.breaker4,z = x,y, x.id, y.id, dateTime, point.x, point.y, poly.xy, parallel, dataType, lonlat, numVertices, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.
    return(list.dist)
  }else{ #if(is.list(x) == FALSE)
    frame.dist<- dist.generator2(x, y, x.id, y.id, dateTime, point.x, point.y, poly.xy, parallel, dataType, lonlat, numVertices, nCores)
    return(frame.dist)
  }
}