#' Project Geographic Coordinates onto a Plane
#'
#' This function converts lon/lat data (decimal degrees) from a geographic 
#'    coordinate system to planar coordinates using a custom azimuthal 
#'    equidistant projection, and appends these new coordinates to an input 
#'    data frame (x). By default, the function assumes longitude and latitude 
#'    coordinates were produced using the WGS84 datum, but users may change 
#'    this datum if they wish.
#'
#' Users may specify longitude and latitude coordinates to become the origin of
#'    the projection (i.e., the (0,0) coordinate). If they do not specify these
#'    values, however, the function calculates the centroid of the data and 
#'    will use this as the origin point. 
#' 
#' Note: this function does not allow any NA coordinate values in 
#'    longitude/latitude vectors. If NAs exist you will get the following 
#'    error: "Error in .local(obj, ...) : NA values in coordinates." If NAs 
#'    exist in your data, we suggest 1.) removing them, or 2.) smoothing data 
#'    using contact::tempAggregate prior to running this function.
#' @param x Data frame or matrix containing geographic data. Defaults to NULL.
#' @param x.lon Vector of length(nrow(x)) or singular character data, detailng 
#'    the relevant colname in x, that denotes what longitude information will 
#'    be used. If argument == NULL, makePlanar assumes a column with the 
#'    colname "lon" exists in x. Defaults to NULL.
#' @param x.lat Vector of length(nrow(x)) or singular character data, detailing
#'    the relevant colname in x, that denotes what latitude information will 
#'    be used. If argument == NULL, makePlanar assumes a column with the 
#'    colname "lat" exists in x. Defaults to NULL.
#' @param origin.lon Numerical. Describes the longitude will be used as the 
#'    origin-point longitude for the azimuthal-equidistant projection. If NULL,
#'    defaults to the longitude of the data set's centroid. Defaults to NULL.
#' @param origin.lat Numerical. Describes the latitude will be used as the 
#'    origin-point latitude for the azimuthal-equidistant projection. If NULL,
#'    defaults to the latitude of the data set's centroid. Defaults to NULL.
#' @param datum Character string describing the datum used to generate x.lon 
#'    and x.lat. Defaults to "WGS84."
#' @keywords data-processing point location planar GRC
#' @return Output is \code{x} appended with the following columns:
#' 
#'    \item{planar.x}{Planar x-coordinate values derived from longitude 
#'    observations.}
#'    \item{planar.y}{Planar y-coordinate values derived from latitude 
#'    observations.}
#'    \item{origin.lon}{Longitude of the origin point, either user specified 
#'    or the longitude of the data's centroid.}
#'    \item{origin.lat}{Latitude of the origin point, either user specified 
#'    or the latitude of the data's centroid.}
#'    \item{origin.distance}{Linear distance (m) between every point and the 
#'    origin point.}
#' @export
#' @examples
#' 
#' data(baboons)
#' 
#' head(baboons) #see that locations are in geographic coordinates
#' 
#' lon.na <- which(is.na(baboons$location.long) == TRUE) #pull row ids of lon NAs
#' lat.na <- which(is.na(baboons$location.lat) == TRUE) #pull row ids of lat NAs
#' 
#' baboons.naRM <- droplevels(baboons[-unique(c(lon.na, lat.na)),]) #remove NAs
#' 
#' baboons.naRM_planar <- makePlanar(x = baboons.naRM, 
#'    x.lon = baboons.naRM$location.long, x.lat = baboons.naRM$location.lat, 
#'    origin.lon = NULL, origin.lat = NULL, datum = "WGS84") #note no specified origin coords
#'    
#' head(baboons.naRM_planar) #see that planar coordinates are reported

makePlanar<-function(x = NULL, x.lon = NULL, x.lat = NULL, origin.lon = NULL, origin.lat = NULL, datum = "WGS84"){

  plane.generator<-function(x, x.lon, x.lat, origin.lon, origin.lat, datum){
    
    #columns with names proceeded by "makePlanarFunction__" are appended to the data set for analysis within this function. These columns will be removed from the output.
    
    if(length(x) == 0){ #This if statement allows users to input either a series of vectors, a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
      originTab = data.frame(lon = x.lon, lat = x.lat, makePlanarFunction__lon = x.lon, makePlanarFunction__lat = x.lat)
    }
    
    if(length(x) > 0){
      
      if(length(x.lat) > 0){ 
        if(length(x.lat) == 1 && is.na(match(x.lat[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than x.lat being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "x.lat" values (i.e., if the x.lats in different list entries are different)
          x$makePlanarFunction__lat <- x[,match(x.lat, names(x))]
        }else{ #if length(x.lat) > 1 OR length(x.lat) == 1, but is.na(match(x.lat[1], names(x))) == TRUE
          if(length(x.lat) == 1 && is.na(match(x.lat[1], names(x))) == TRUE){ #if x.lat is a vector of length ==1, but does not represent a column name in x, this one value will be the only x.lat used in analysis.
            x$makePlanarFunction__lat <- rep(x.lat, nrow(x))
          }
          if(length(x.lat) > 1){
            x$makePlanarFunction__lat <- x.lat
          }
        }
      }else{ #if length(x.lat) == 0
        x$makePlanarFunction__lat <- x$lat #if length(x.lat) == 0, then the function assumes there is a "lat" column in x
      }
      
      if(length(x.lon) > 0){ 
        if(length(x.lon) == 1 && is.na(match(x.lon[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than x.lon being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "x.lon" values (i.e., if the x.lons in different list entries are different)
          x$makePlanarFunction__lon <- x[,match(x.lon, names(x))]
        }else{ #if length(x.lon) > 1 OR length(x.lon) == 1, but is.na(match(x.lon[1], names(x))) == TRUE
          if(length(x.lon) == 1 && is.na(match(x.lon[1], names(x))) == TRUE){ #if x.lon is a vector of length ==1, but does not represent a column name in x, this one value will be the only x.lon used in analysis.
            x$makePlanarFunction__lon <- rep(x.lon, nrow(x))
          }
          if(length(x.lon) > 1){
            x$makePlanarFunction__lon <- x.lon
          }
        }
      }else{ #if length(x.lon) == 0
        x$makePlanarFunction__lon <- x$lon #if length(x.lon) == 0, then the function assumes there is a "lon" column in x
      }
      
      originTab <- x
    }
    
    if(length(origin.lon) == 0){ #if origin.lon == NULL (i.e., users did not specify a origin.lon coordinate), the origin.lon value will default to the mean x$lon value
      cntr.lon <- mean(originTab$makePlanarFunction__lon, na.rm = TRUE)
    }else{ #if origin.lon != NULL
      cntr.lon <- origin.lon 
    }
    
    if(length(origin.lat) == 0){ #if origin.lat == NULL (i.e., users did not specify a origin.lat coordinate), the origin.lat value will default to the mean x$lat value
      cntr.lat <- mean(originTab$makePlanarFunction__lat, na.rm = TRUE)
    }else{ #if origin.lat != NULL
      cntr.lat <- origin.lat
    }
    
    #identify the centroid of the data and measure the distance to it at any given timestep
    centroid.dist <- raster::pointDistance(originTab[,c(match("makePlanarFunction__lon",names(originTab)),match("makePlanarFunction__lat",names(originTab)))], c(cntr.lon, cntr.lat), lonlat = TRUE) #note that because the data will be projected onto a planar surface using an equidistant projection, all distances (and relative positions) to the origin point (in this case, the centroid) will be preserved 
    
    #Project the data using an equidistant azimuthal projection centered at the centroid.
    Coords<- sp::SpatialPoints(coords = originTab[,c(match("makePlanarFunction__lon",names(originTab)),match("makePlanarFunction__lat",names(originTab)))], proj4string = sp::CRS(paste("+proj=longlat +datum=",datum, sep = ""))) #projects the long/lat data on the elipsoid described by datum. datum default is WGS84.
    planarCoords<-sp::spTransform(Coords, sp::CRS(paste("+proj=aeqd +lon_0=",cntr.lon," +lat_0=",cntr.lat, " +units=m", sep = ""))) #creates a spatial-class feature with points projected on a planar surface using an equidistant projection centered at the centroid of the data. The "units=m" argument indicates that the coordinate unit is meters
    planarCoords.frame <- as.data.frame(planarCoords)
    names(planarCoords.frame)<-c("planar.x", "planar.y")
    
    #Now we bind this new coordinate information to originTab
    bindlist<-list(originTab, planarCoords.frame)
    outputTab<- do.call("cbind", bindlist)
    outputTab$origin.lon <- cntr.lon
    outputTab$origin.lat <- cntr.lat
    outputTab$origin.distance <- centroid.dist
    
    #Here we remove the "makePlanarFunction__" columns.
    outputTab<-outputTab[,-c(match("makePlanarFunction__lon", names(outputTab)), match("makePlanarFunction__lat", names(outputTab)))]
    
    return(outputTab)
  }
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #1/15 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    new.coords <- lapply(x, plane.generator, x.lon, x.lat, origin.lon, origin.lat, datum) 
  }else{ #if(is.list(x) == FALSE)
    new.coords<- plane.generator(x, x.lon, x.lat, origin.lon, origin.lat, datum)
  }
  return(new.coords)
}

