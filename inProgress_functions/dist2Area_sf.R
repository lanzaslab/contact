#' Calculate Distances Between Individuals and Fixed Points/Polygons
#'
#' Calculate distances between each individual, reported in x, and a fixed 
#'    point(s)/polygon(s), reported in y, at each timestep. 
#'
#' This variant of dist2Area requires that x and y inputs must be spatial 
#'    objects of class sf (i.e., shapefile), SpatialPointsDataFrame, or 
#'    SpatialPolygonsDataFrame.
#'    
#' If inputs are SpatialMultiPointsDataFrame, which have one-to-many 
#'    cardinality, you will receive this error: "Error in .pointsToMatrix(p) : 
#'    points should be vectors of length 2, matrices with 2 columns, or 
#'    inheriting from a SpatialPoints* object." In this case, we suggest 
#'    removing duplicate entries with the dup function and converting the input
#'    to a SpatialPointsDataFrame prior to running the function.
#' @param x Spatial object containing real-time-location data for individuals. 
#' @param y Spatial object containing fixed-area polygons/points for which we 
#'    will calculate distances relative to tracked individuals at all time 
#'    steps.
#' @param x.id Vector of length nrow(x) or singular character data, detailing 
#'    the relevant colname in x, that denotes what unique ids for tracked 
#'    individuals will be used. If argument == NULL, the function assumes a 
#'    column with the colname "id" exists in x. Defaults to NULL.
#' @param y.id Vector of length nrow(y) or singular character data, detailing 
#'    the relevant colname in y, that denotes what unique ids for fixed-area 
#'    polygons/points will be used. If argument == NULL, the function assumes a
#'    column with the colname "id" exists in y. Defaults to NULL.
#' @param dateTime Vector of length nrow(data.frame(x)) or singular character 
#'    data, detailing the relevant colname in x, that denotes what dateTime 
#'    information will be used. If argument == NULL, the function assumes a 
#'    column with the colname "dateTime" exists in x. Defaults to NULL.
#' @param parallel Logical. If TRUE, sub-functions within the dist2Area_sf 
#'    wrapper will be parallelized. Note that this can significantly speed up 
#'    processing of relatively small data sets, but may cause R to crash due to 
#'    lack of available memory when attempting to process large datasets. 
#'    Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to the maximum number of cores available
#'    (i.e., parallel::detectCores()).
#' @keywords data-processing polygon point location spatial
#' @export
#' @examples
#' #Examples imminent

dist2Area_sf<-function(x, y, x.id = NULL, y.id = NULL, dateTime = NULL, parallel = FALSE, nCores = parallel::detectCores()){
  
  if(length(x.id) > 0){ 
    if(length(x.id) == 1 & is.na(match(x.id[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than id being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
      x$id <- x[,match(x.id, names(x))]
    }else{ #if length(id) > 1
      x$id = x.id
    }
  }
  
  if(length(y.id) > 0){ 
    if(length(y.id) == 1 & is.na(match(y.id[1], names(y))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If y is a list, rather than id being a vector of length(nrow(y)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
      y$id <- y[,match(y.id, names(y))]
    }else{ #if length(id) > 1
      y$id = y.id
    }
  }
  
  if(length(dateTime) > 0){
    if(length(dateTime) == 1 & is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
      x$dateTime <- x[,match(dateTime, names(x))]
    }else{ #if length(dateTime) > 1
      x$dateTime = dateTime
    }
  }
  
  if(class(x)[1] == "sf"){ #if class(x) is a shapefile, rather than Spatial data (the 2 potential types of inputs for this function), we must convert it to Spatial data.
    x.spatial <- sf::as_Spatial(x)
  }else{ #if x is Spatial data
    x.spatial <- x
  }
  if(class(y)[1] != "sf"){ #if class(y) is Spatial data, rather than sf (the 2 potential types of inputs for this function), we must convert it sf so that we can break up the feature set into single features.
    y.sf <- methods::as(y, "sf")
  }else{ #if y is a shapefile
    y.sf <- y
  }
  
  breakFrame<-data.frame(key.id = seq(nrow(y.sf)), y.id = y.sf$id)# This frame will be used to call each specific feature in the distance function below.

  spatDistCalc<-function(x,y,z){ #here, x is breakFrame, y is y.sf, and z is x.spatial
    #spatial_object<-sf::as(y[unname(unlist(x[1])),],"Spatial") #grab the specific fixed polygon from y (the shapefile input) denoted by x[1] of breakFrame
    spatial_object<-sf::as_Spatial(y[unname(unlist(x[1])),]) #grab the specific fixed polygon from y (the shapefile input) denoted by x[1] of breakFrame
    distances<-geosphere::dist2Line(z,spatial_object) #calculate the angular distance, in meters, between features in z and the fixed polygons
    output<-data.frame(distances) #output a dataframe so that 
    output$x.id <- z$id #adds the x.spatial id column to the data frame
    output$dateTime <- z$dateTime #adds the x.spatial dateTime column to the data frame
    output$key.id <- unname(unlist(x[1])) #The geosphere::dist2Line function calculates the distance to the nearest polygon within a polygon set. The "ID" column in the output represents which feature in the set is closest to each point, p. Because we use the apply function to break the polygon set into subsets of 1 feature, the "ID" column in the output will always be "1." We fix that here, by forcing the ID column to represent the feature sequence of the full set.
    return(output)
  }
  
  if (parallel == TRUE){
    cl<-parallel::makeCluster(nCores)
    if(nrow(y.sf) > 9){ #for some reason the apply function would always return " Error in from[[1]] : subscript out of bounds " when processing rows 9-10. When run individually or run with sets <= 9, or >= 10 it worked just fine.
      spat.apply1 <-parallel::parApply(cl, breakFrame[1:9,], 1, spatDistCalc, y = y.sf, z = x.spatial)
      spat.apply2 <-parallel::parApply(cl, breakFrame[10:nrow(y.sf),], 1, spatDistCalc, y = y.sf, z = x.spatial)
      parallel::stopCluster(cl)
      bind1<-data.frame(data.table::rbindlist(spat.apply1))
      bind2<-data.frame(data.table::rbindlist(spat.apply2))
      spat.apply<-list(bind1,bind2)
    }else{ # if nrow(y.sf) <= 9
      spat.apply<-parallel::parApply(cl, breakFrame, 1, spatDistCalc, y = y.sf, z = x.spatial)
      parallel::stopCluster(cl)
    }
  }else{ #if parallel == FALSE
    if(nrow(y.sf) > 9){ #for some reason the apply function would always return " Error in from[[1]] : subscript out of bounds " when processing rows 9-10. When run individually or run with sets <= 9, or >= 10 it worked just fine.
      spat.apply1 <- apply(breakFrame[1:9,], 1, spatDistCalc, y = y.sf, z = x.spatial)
      spat.apply2 <- apply(breakFrame[10:nrow(y.sf),], 1, spatDistCalc, y = y.sf, z = x.spatial)
      bind1<-data.frame(data.table::rbindlist(spat.apply1))
      bind2<-data.frame(data.table::rbindlist(spat.apply2))
      spat.apply<-list(bind1,bind2)
    }else{ # if nrow(y.sf) <= 9
      spat.apply<-apply(breakFrame, 1, spatDistCalc, y = y.sf, z = x.spatial)
    }
  }
  sp_bind<-data.frame(data.table::rbindlist(spat.apply))
  sp_merge<-merge(sp_bind, breakFrame, by = "key.id") #this gives us a data frame showing the distance to each fixed polygon (y.id) from each individual (x.id) at every timepoint.
  
  create.distFrame<- function(x,sp_merge, indivSeq, areaSeq, indivAreaFrame){ #This function ensures that the output data frame will be formatted in such a way that data can be processed by the dur.area function.
    
    identifyDistances<-function(x, timeSub.distance, timeSub.x.id, timeSub.y.id){ #This pulls out distances between each area-individual pair in the order that they will be presented in dist
      output<-timeSub.distance[which(timeSub.x.id == unname(unlist(x[1])) & timeSub.y.id == unname(unlist(x[2])))]
      return(output)
    }
    
    timeSub<-subset(sp_merge, dateTime == unname(unlist(x[1]))) #create a subset with data only at a single timestep
    ##Here we vectorize all the timeSub components of interest
    timeSub.distance <- timeSub$distance
    timeSub.x.id <- timeSub$x.id
    timeSub.y.id <- timeSub$y.id
    
    distances.reported<-unlist(apply(indivAreaFrame,1,identifyDistances, timeSub.distance, timeSub.x.id, timeSub.y.id)) #pull out a vector of the distances for each area-individual pair

    dist = data.frame(matrix(ncol = (length(areaSeq) + 4), nrow = length(unique(timeSub.x.id))))
    colnames(dist) = c("dateTime","totalIndividuals","individualsAtTimestep","id", paste("dist.to.", areaSeq, sep = ""))
    dist$dateTime = unname(unlist(x[1])) #date taken from timestepFrame
    dist$totalIndividuals = length(indivSeq)
    dist$individualsAtTimestep = length(unique(timeSub.x.id))
    dist$id = unique(timeSub.x.id)
    dist[,5:ncol(dist)] <- distances.reported
    return(dist)
  }
  
  indivSeq <- unique(sp_merge$x.id)
  areaSeq <- unique(sp_merge$y.id)
  timeSeq <- unique(sp_merge$dateTime)
  
  timestepFrame<-data.frame(unique(sp_merge$dateTime))
  indivAreaFrame <- expand.grid(indivSeq, areaSeq)
  
  if (parallel == TRUE){
    cl<-parallel::makeCluster(nCores)
    distTab<-parallel::parApply(cl, timestepFrame, 1, create.distFrame, sp_merge, indivSeq, areaSeq, indivAreaFrame)
    parallel::stopCluster(cl)
  }else{ #if parallel == FALSE
    distTab<-apply(timestepFrame, 1, create.distFrame, sp_merge, indivSeq, areaSeq, indivAreaFrame)
  }
  distTab.agg = data.frame(data.table::rbindlist(distTab)) 
  return(distTab.agg)
}
