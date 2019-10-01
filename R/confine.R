#' Identify and Remove Data Points Outside of a Specified Area
#'
#' dup (a.k.a. Multiple instance filter) identifies and removes timepoints when
#'    tracked individuals were observed in >1 place concurrently. If avg == TRUE,
#'    duplicates are replaced by a single row describing an individuals' average
#'    location (e.g., planar xy coordinates) during the duplicated time point.
#'    If avg == FALSE, all duplicated timepoints will be removed, as there is
#'    no way for the function to determine which instance among the 
#'    duplicates should stay. If users are not actually interested in filtering 
#'    datasets, but rather, determining what observations should be filtered, 
#'    they may set filterOutput == FALSE. By doing so, this function will 
#'    append a "duplicated" column to the dataset, which reports values that 
#'    describe if any timepoints in a given individual's path are duplicated. 
#'    Values are: 0: timepoint is not duplicated, 1: timepoint is duplicated.
#'
#' Identifies and removes timepoints when tracked individuals were observed 
#'    outside of a defined polygon (note: the polygon should be described by 
#'    the vectors confinementCoord.x (x coordinates) and confinementCoord.y 
#'    (y coordinates) (note: these vectors must be the same length and the 
#'    coordinates should be listed in the clockwise or counter-clockwise 
#'    order that they are observed on the confining polygon), which identify
#'    where polygon vertices exist).
#'
#' If users are not actually interested in filtering datasets, but rather, 
#'    determining what observations should be filtered, they may set 
#'    filterOutput == FALSE. By doing so, this function will append a 
#'    "confinement_status" column to the output dataframe, which reports 
#'    the results of sp::point.in.polygon function that is used to determine
#'    if individuals are confined within a given polygon. In this column, 
#'    values are: 0: point is strictly exterior to pol; 1: point is strictly 
#'    interior to pol; 2: point lies on the relative interior of an edge of 
#'    pol; 3: point is a vertex of pol (see ?sp::point.in.polygon).
#'    
#' @param x Non-data-frame list or data frame that will be filtered.
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
#' @param confinementCoord.x Vector describing x-coordinates of 
#'    confining-polygon vertices. Each vertex should be described in clockwise 
#'    or counter-clockwise order, and ordering should be consistent with 
#'    confinementCoord.y.
#' @param confinementCoord.y Vector describing y-coordinates of 
#'    confining-polygon vertices. Each vertex should be described in clockwise 
#'    or counter-clockwise order, and ordering should be consistent with 
#'    confinementCoord.x.
#' @param filterOutput Logical. If TRUE, output will be a data frame or list of
#'    data frames (depending on whether or not x is a data frame or not) 
#'    containing only points within confinement polygons. If FALSE, no 
#'    observations are removed and a "confinement_status" column is appended to
#'    x, detailing the relationship of each point to the confinement polygon. 
#'    Defaults to TRUE.
#' @keywords filter confinement polygon
#' @export
#' @examples
#' #read in the calves data set
#' data("calves")
#' 
#' #report the x and y coordinates for a polygon. In this case the vertices of
#'    #the water trough within the feedlot pen where calves were housed.
#' water_trough.x<- c(61.43315, 61.89377, 62.37518, 61.82622)
#' water_trough.y<- c(62.44815, 62.73341, 61.93864, 61.67411)
#' 
#' #determine when calves' heads (note that real-time-location points describe
#'    #the location of radio-tracking eartags on left ears) were within the 
#'    #confines of the water trough.
#' headWater1<- confine(calves, point.x = calves$x, point.y = calves$y, 
#'    confinementCoord.x = water_trough.x, confinementCoord.y = water_trough.y,
#'    filterOutput = TRUE) #creates a data set comprised ONLY of points within 
#'    #the water polygon.
#' headWater2<- confine(calves, point.x = calves$x, point.y = calves$y, 
#'    confinementCoord.x = water_trough.x, confinementCoord.y = water_trough.y,
#'    filterOutput = FALSE) #appends the "confinement_status" column to x.

confine <- function(x, point.x = NULL, point.y = NULL, confinementCoord.x, 
                    confinementCoord.y, filterOutput = TRUE){
  
  filter2.func<-function(x, point.x, point.y, confinementCoord.x, confinementCoord.y, filterOutput){

    xVec <- NULL
    yVec <- NULL

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

    x$confinement_status = sp::point.in.polygon(point.x = xVec1, point.y = yVec1,pol.x = confinementCoord.x,pol.y = confinementCoord.y)
    if(filterOutput == TRUE){
      x<- subset(x, confinement_status >= 1)
      x <- x[,-match("confinement_status",names(x))] #confinement_filter is not necessary for future functions presented here. Researchers may want this field removed to reduce file size
    }
    if(nrow(x) > 0){
      rownames(x) <-seq(1,nrow(x),1)
    }
    return(x)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    list.conf<-lapply(x, filter2.func, point.x, point.y, confinementCoord.x, confinementCoord.y, filterOutput)
    return(list.conf)
  }else{ #if x is a dataFrame
    frame.conf<- filter2.func(x, point.x, point.y, confinementCoord.x, confinementCoord.y, filterOutput)
    return(frame.conf)
  }

}
