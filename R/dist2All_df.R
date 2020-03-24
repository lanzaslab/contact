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
#'    will be parallelized. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
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
#'    \item{dist.to.indiv_...}{The observed distance between the individual 
#'    described in the \code{id} column to every other individual observed 
#'    at specific timepoints.}
#' @import foreach
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
#' calves.dist2<-dist2All_df(x = calves.agg, parallel = FALSE, dataType = "Point",
#'    lonlat = FALSE) #calculate distance between all individuals at each timepoint.
#' 

dist2All_df<-function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, poly.xy = NULL, elev = NULL, parallel = FALSE, nCores = (parallel::detectCores()/2), dataType = "Point", lonlat = FALSE, numVertices = 4){
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  i <- NULL
  j <- NULL
  k <- NULL
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R treats dataframes as lists.
    
    listBreak_dist.generator <-function(x,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){ #this function is exactly the same as what happens when the x input to the master function is a single data frame.
      
      #write all the sub-functions first
      
      dist.generator<-function(x,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){ 
        
        create.distFrame<- function(x,distMat, indivSeq, idSeq, timestepIndivSeq,time){
          
          dist = data.frame(matrix(ncol = (length(indivSeq) + 4), nrow = 1), stringsAsFactors = TRUE)
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

        if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
          
          #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
          daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
          lub.dates = lubridate::date(x$dateTime)
          x<-x[order(x$id, lub.dates, daySecondList),] #order x 
          
          rm(list = c("daySecondList", "lub.dates")) #remove these objects because they are no longer needed.
          
          indivSeq = unique(x$integ.ID)
          idSeq = unique(x$id)
          dist.measurement = lonlat
          
          dist.process.point <- function(x, y, indivSeq, idSeq, dist.measurement, elev.calc){
            
            time = unlist(unname(x[1]))
            timestep = y[which(y$dateTime == time),]
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
                elevFrame <- data.frame(elevVec, stringsAsFactors = TRUE)
                elev.dif<-apply(elevFrame, 1, elevDifference, elevVec)
                distMat<-distMat + elev.dif #adds the difference in elevations to the linear distances previously calculated. Note that the elev distance units must be in meters, as that is the unit of distMat values. 
              }
              timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1), timestepIndivSeq.integ, stringsAsFactors = TRUE) #Note that colnum will not necessarily be the same as timestepIndivSeq.integ because y has been previously ordered
              distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, idSeq, timestepIndivSeq.integ,time)), stringsAsFactors = TRUE)
              return(distTab)
            }
          }
          
          distTab <- foreach::foreach(i = unique(x$dateTime), .packages = 'foreach') %do% dist.process.point(i, y = x, indivSeq,idSeq, dist.measurement, elev.calc)
          
          dist.all = data.frame(data.table::rbindlist(distTab), stringsAsFactors = TRUE)
        }
        
        if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
          
          #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
          daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
          lub.dates = lubridate::date(x$dateTime)
          x<-x[order(x$id, lub.dates, daySecondList),] #order x 
          
          rm(list = c("daySecondList", "lub.dates")) #remove these objects because they are no longer needed.
          
          naVec<-which(is.na(x[,match("point2.x", names(x))]) == TRUE) #the referencePointtoPolygon function will create some observations that are not complete polygons (i.e., only the point1 coordinates are recorded). This identifies those observations, so that they may be removed. If they are not removed, they will cause errors later.
          if(length(naVec) > 0){
            x <- x[-naVec,]
          }
          indivSeq = unique(x$integ.ID)
          idSeq <- unique(x$id)
          
          dist.process.poly <- function(x, y, indivSeq, idSeq, numVertices, elev.calc){
            create.poly1<-function(x,y, numVertices){
              polygon = unlist(c(y[unlist(unname(x[1])),c(2:(1 + (numVertices*2)),2:3)])) #note that when making a Spatial Polygon you need to close the polygon by repeating the first point at the end (Note that this also means you need to add a "+ 1" to nrow = numVertices below)
              spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices + 1), ncol = 2, byrow = TRUE))),y$integ.ID[unlist(unname(x[1]))])
              return(spatPoly)
            }
            
            time = unlist(unname(x[1]))
            timestep = y[which(y$dateTime == time),]
            
            if(nrow(timestep) > 1){ #If there's only one record in timestep, any distMat produced will only denote the distance between an individual and itself (which isn't useful), so it's more time efficient to skip forward nrow(timestep <= 1)
              
              timestepIndivSeq.integ = unique(timestep$integ.ID)
              makePolyFrame<-data.frame(seq(1,length(timestepIndivSeq.integ),1), stringsAsFactors = TRUE)
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
                elevFrame <- data.frame(elevVec, stringsAsFactors = TRUE)
                elev.dif<-apply(elevFrame, 1, elevDifference, elevVec)
                distMat<-distMat + elev.dif #adds the difference in elevations to the linear distances previously calculated. Note that the elev distance units must be in meters, as that is the unit of distMat values. 
              }
              
              timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1), timestepIndivSeq.integ, stringsAsFactors = TRUE) #Note that colnum will not necessarily be the same as timestepIndivSeq.integ because y has been previously ordered
              distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, idSeq, timestepIndivSeq.integ,time)), stringsAsFactors = TRUE)
              return(distTab)
            }
          }
          
          distTab <- foreach::foreach(i = unique(x$dateTime), .packages = 'foreach') %do% dist.process.poly(i, y = x, indivSeq, idSeq, numVertices, elev.calc)
          
          dist.all = data.frame(data.table::rbindlist(distTab), stringsAsFactors = TRUE)
        }
        return(dist.all)
      }
      
      date_hourSub.func<-function(x, data){ #added 02/20/2020 #This function will be used to break down data sets into hourly time blocks prior to further processing to increase speed. Admittedly, this is not a pretty fix for increasing efficiency of processing large data sets, but it's a working fix nonetheless. 
        date_hour <- droplevels(data[which(data$date_hour == unname(unlist(x[1]))),]) #subset data
        return(date_hour)
      }
      
      day_listDistance <- function(x, data.list, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){ #Because this function slows down when trying to process large data frames AND large list sets, we must concattenate both here. We did so to the former by breaking the data frame into hourly lists, and the latter by breaking these lists into daily subsets with this function.
        
        day_lists <- data.list[grep(unname(unlist(x[1])), names(data.list))] #pulls the hour lists within a given day
        names(day_lists)<-NULL #ensure that list names do not mess up column names
        list.dist <- lapply(day_lists, dist.generator, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply.
        dist.bind <- data.frame(data.table::rbindlist(list.dist, fill = TRUE), stringsAsFactors = TRUE) #bind these hours back together
        
        return(dist.bind)
      }
      
      idVec1=NULL #added in the case that idVec1 isn't created when x isn't specified
      if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
        if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
          originTab = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime, stringsAsFactors = TRUE)
        }
        if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
          
          originTab = data.frame(matrix(ncol = 0, nrow = length(id)), stringsAsFactors = TRUE)
          originTab$id = id
          colnames(poly.xy)[seq(1,(ncol(poly.xy) - 1),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".x", sep = "")
          colnames(poly.xy)[seq(2,ncol(poly.xy),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".y", sep = "")
          dateFrame = data.frame(dateTime = dateTime, stringsAsFactors = TRUE)
          bindlist = list(originTab,poly.xy,dateFrame)
          originTab = data.frame(do.call("cbind", bindlist), stringsAsFactors = TRUE)
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
          xyFrame1<- data.frame(x = x$x, y = x$y, stringsAsFactors = TRUE)
        }
        
        if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
          
          if(length(poly.xy) > 0){
            if(length(as.matrix(poly.xy)) == (numVertices*2) && length(which(is.na(match(as.matrix(poly.xy), names(x)))) == TRUE) == 0){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than poly.xy being a matrix/dataframe of length(nrow(x)), it may be necessary to designate the colnames for intended coordinate values (i.e., if the xy-coordinate values in different list entries are different)
              xyFrame1<-x[,match(poly.xy,names(x))]
            }else{
              #x[,2:(1 + length(poly.xy))] = poly.xy
              xyFrame1<- data.frame(poly.xy, stringsAsFactors = TRUE)
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
      
      rm(x) #now that originTab is created, we no longer need x
      
      #The next thing we need to do is remove any NAs in the data set 
      if(length(unique(c(which(is.na(originTab$x) == TRUE), which(is.na(originTab$y) == TRUE)))) > 0){
        
        originTab <- droplevels(originTab[- unique(c(which(is.na(originTab$x) == TRUE), which(is.na(originTab$y) == TRUE))),])
        
      }
      
      data.dates<-lubridate::date(originTab$dateTime) #now we can start concattenating the data by subsetting it into smaller lists
      
      originTab$date_hour <- paste(data.dates, lubridate::hour(originTab$dateTime), sep = "_") #create a tag for each unique date_hour combination in the data set
      date_hour.vec <- unique(originTab$date_hour)
      date.vec <- unique(data.dates)
      
      if(length(date_hour.vec) == 1){ #the processing step requires a list of data frames. If there's only a single hour represented in originTab, we can just create the list using the "list" function
        data.list <- list(originTab)
      }else{
        data.list <- foreach::foreach(i = date_hour.vec) %do% date_hourSub.func(i, originTab)
      }
      
      names(data.list)<-date_hour.vec #add names to list to pull for date lists below
      
      rm(list =  c("originTab", "data.dates", "date_hour.vec")) #remove the unneeded objects to free up memory
      
      if(parallel == TRUE){
        
        cl <- parallel::makeCluster(nCores)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))
        distances<-foreach::foreach(j = date.vec, .packages = 'foreach') %dopar% day_listDistance(j, data.list, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
        
      }else{ #if parallel == FALSE
        
        distances<-foreach::foreach(j = date.vec, .packages = 'foreach') %do% day_listDistance(j, data.list, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
        
      }
      
      frame.dist<- data.frame(data.table::rbindlist(distances, fill = TRUE), stringsAsFactors = TRUE)
      
      return(frame.dist)
      
    }
    
    list.dist<- foreach::foreach(k = 1:length(x), .packages = 'foreach') %do% listBreak_dist.generator(x[[k]],id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops
    
    return(list.dist)
    
  }else{ #if(is.data.frame(x) == TRUE)
    
    #write all the sub-functions first
    
    dist.generator<-function(x,id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){ 
      
      create.distFrame<- function(x,distMat, indivSeq, idSeq, timestepIndivSeq,time){
        
        dist = data.frame(matrix(ncol = (length(indivSeq) + 4), nrow = 1), stringsAsFactors = TRUE)
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
      
      
      if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
        
        #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
        daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
        lub.dates = lubridate::date(x$dateTime)
        x<-x[order(x$id, lub.dates, daySecondList),] #order x 
        
        rm(list = c("daySecondList", "lub.dates")) #remove these objects because they are no longer needed.
        
        indivSeq = unique(x$integ.ID)
        idSeq = unique(x$id)
        dist.measurement = lonlat
        
        dist.process.point <- function(x, y, indivSeq, idSeq, dist.measurement, elev.calc){
          
          time = unlist(unname(x[1]))
          timestep = y[which(y$dateTime == time),]
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
              elevFrame <- data.frame(elevVec, stringsAsFactors = TRUE)
              elev.dif<-apply(elevFrame, 1, elevDifference, elevVec)
              distMat<-distMat + elev.dif #adds the difference in elevations to the linear distances previously calculated. Note that the elev distance units must be in meters, as that is the unit of distMat values. 
            }
            timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1), timestepIndivSeq.integ, stringsAsFactors = TRUE) #Note that colnum will not necessarily be the same as timestepIndivSeq.integ because y has been previously ordered
            distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, idSeq, timestepIndivSeq.integ,time)), stringsAsFactors = TRUE)
            return(distTab)
          }
        }
        
        distTab <- foreach::foreach(i = unique(x$dateTime), .packages = 'foreach') %do% dist.process.point(i, y = x, indivSeq,idSeq, dist.measurement, elev.calc)
        
        dist.all = data.frame(data.table::rbindlist(distTab), stringsAsFactors = TRUE)
      }
      
      if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
        
        #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
        daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
        lub.dates = lubridate::date(x$dateTime)
        x<-x[order(x$id, lub.dates, daySecondList),] #order x 
        
        rm(list = c("daySecondList", "lub.dates")) #remove these objects because they are no longer needed.
        
        naVec<-which(is.na(x[,match("point2.x", names(x))]) == TRUE) #the referencePointtoPolygon function will create some observations that are not complete polygons (i.e., only the point1 coordinates are recorded). This identifies those observations, so that they may be removed. If they are not removed, they will cause errors later.
        if(length(naVec) > 0){
          x <- x[-naVec,]
        }
        indivSeq = unique(x$integ.ID)
        idSeq <- unique(x$id)
        
        dist.process.poly <- function(x, y, indivSeq, idSeq, numVertices, elev.calc){
          create.poly1<-function(x,y, numVertices){
            polygon = unlist(c(y[unlist(unname(x[1])),c(2:(1 + (numVertices*2)),2:3)])) #note that when making a Spatial Polygon you need to close the polygon by repeating the first point at the end (Note that this also means you need to add a "+ 1" to nrow = numVertices below)
            spatPoly = sp::Polygons(list(sp::Polygon(matrix(polygon, nrow = (numVertices + 1), ncol = 2, byrow = TRUE))),y$integ.ID[unlist(unname(x[1]))])
            return(spatPoly)
          }
          
          time = unlist(unname(x[1]))
          timestep = y[which(y$dateTime == time),]
          
          if(nrow(timestep) > 1){ #If there's only one record in timestep, any distMat produced will only denote the distance between an individual and itself (which isn't useful), so it's more time efficient to skip forward nrow(timestep <= 1)
            
            timestepIndivSeq.integ = unique(timestep$integ.ID)
            makePolyFrame<-data.frame(seq(1,length(timestepIndivSeq.integ),1), stringsAsFactors = TRUE)
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
              elevFrame <- data.frame(elevVec, stringsAsFactors = TRUE)
              elev.dif<-apply(elevFrame, 1, elevDifference, elevVec)
              distMat<-distMat + elev.dif #adds the difference in elevations to the linear distances previously calculated. Note that the elev distance units must be in meters, as that is the unit of distMat values. 
            }
            
            timestepIndivSeqFrame = data.frame(id = unique(timestep$id), colnum = seq(1,length(unique(timestep$id)),1), timestepIndivSeq.integ, stringsAsFactors = TRUE) #Note that colnum will not necessarily be the same as timestepIndivSeq.integ because y has been previously ordered
            distTab = data.frame(data.table::rbindlist(apply(timestepIndivSeqFrame,1,create.distFrame,distMat,indivSeq, idSeq, timestepIndivSeq.integ,time)), stringsAsFactors = TRUE)
            return(distTab)
          }
        }
        
        distTab <- foreach::foreach(i = unique(x$dateTime), .packages = 'foreach') %do% dist.process.poly(i, y = x, indivSeq, idSeq, numVertices, elev.calc)
        
        dist.all = data.frame(data.table::rbindlist(distTab), stringsAsFactors = TRUE)
      }
      return(dist.all)
    }
    
    date_hourSub.func<-function(x, data){ #added 02/20/2020 #This function will be used to break down data sets into hourly time blocks prior to further processing to increase speed. Admittedly, this is not a pretty fix for increasing efficiency of processing large data sets, but it's a working fix nonetheless. 
      date_hour <- droplevels(data[which(data$date_hour == unname(unlist(x[1]))),]) #subset data
      return(date_hour)
    }
    
    day_listDistance <- function(x, data.list, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores){ #Because this function slows down when trying to process large data frames AND large list sets, we must concattenate both here. We did so to the former by breaking the data frame into hourly lists, and the latter by breaking these lists into daily subsets with this function.
      
      day_lists <- data.list[grep(unname(unlist(x[1])), names(data.list))] #pulls the hour lists within a given day
      names(day_lists)<-NULL #ensure that list names do not mess up column names
      list.dist <- lapply(day_lists, dist.generator, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply.
      dist.bind <- data.frame(data.table::rbindlist(list.dist, fill = TRUE), stringsAsFactors = TRUE) #bind these hours back together
      
      return(dist.bind)
    }
    
    idVec1=NULL #added in the case that idVec1 isn't created when x isn't specified
    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
      if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
        originTab = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime, stringsAsFactors = TRUE)
      }
      if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
        
        originTab = data.frame(matrix(ncol = 0, nrow = length(id)), stringsAsFactors = TRUE)
        originTab$id = id
        colnames(poly.xy)[seq(1,(ncol(poly.xy) - 1),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".x", sep = "")
        colnames(poly.xy)[seq(2,ncol(poly.xy),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".y", sep = "")
        dateFrame = data.frame(dateTime = dateTime, stringsAsFactors = TRUE)
        bindlist = list(originTab,poly.xy,dateFrame)
        originTab = data.frame(do.call("cbind", bindlist), stringsAsFactors = TRUE)
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
        xyFrame1<- data.frame(x = x$x, y = x$y, stringsAsFactors = TRUE)
      }
      
      if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
        
        if(length(poly.xy) > 0){
          if(length(as.matrix(poly.xy)) == (numVertices*2) && length(which(is.na(match(as.matrix(poly.xy), names(x)))) == TRUE) == 0){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than poly.xy being a matrix/dataframe of length(nrow(x)), it may be necessary to designate the colnames for intended coordinate values (i.e., if the xy-coordinate values in different list entries are different)
            xyFrame1<-x[,match(poly.xy,names(x))]
          }else{
            #x[,2:(1 + length(poly.xy))] = poly.xy
            xyFrame1<- data.frame(poly.xy, stringsAsFactors = TRUE)
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
    
    rm(x) #now that originTab is created, we no longer need x
    
    #The next thing we need to do is remove any NAs in the data set 
    if(length(unique(c(which(is.na(originTab$x) == TRUE), which(is.na(originTab$y) == TRUE)))) > 0){
      
      originTab <- droplevels(originTab[- unique(c(which(is.na(originTab$x) == TRUE), which(is.na(originTab$y) == TRUE))),])
      
    }
    
    data.dates<-lubridate::date(originTab$dateTime) #now we can start concattenating the data by subsetting it into smaller lists
    
    originTab$date_hour <- paste(data.dates, lubridate::hour(originTab$dateTime), sep = "_") #create a tag for each unique date_hour combination in the data set
    date_hour.vec <- unique(originTab$date_hour)
    date.vec <- unique(data.dates)
    
    if(length(date_hour.vec) == 1){ #the processing step requires a list of data frames. If there's only a single hour represented in originTab, we can just create the list using the "list" function
      data.list <- list(originTab)
    }else{
      data.list <- foreach::foreach(i = date_hour.vec) %do% date_hourSub.func(i, originTab)
    }
    
    names(data.list)<-date_hour.vec #add names to list to pull for date lists below
    
    rm(list =  c("originTab", "data.dates", "date_hour.vec")) #remove the unneeded objects to free up memory
    
    if(parallel == TRUE){
      
      cl <- parallel::makeCluster(nCores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      distances<-foreach::foreach(j = date.vec, .packages = 'foreach') %dopar% day_listDistance(j, data.list, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
      
    }else{ #if parallel == FALSE
      
      distances<-foreach::foreach(j = date.vec, .packages = 'foreach') %do% day_listDistance(j, data.list, id, dateTime, point.x, point.y, poly.xy, elev, parallel, dataType, lonlat, numVertices, nCores) #we set the .packages argument to 'foreach' to allow us to use foreach loops within foreach loops. Note that the parallel and nCores arguments here are artifacts of previous function iterations. They do not affect anything going forward.
      
    }
    
    frame.dist<- data.frame(data.table::rbindlist(distances, fill = TRUE), stringsAsFactors = TRUE)
    
    return(frame.dist)
  }
}
