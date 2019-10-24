#' Identify and Remove Data Points Outside of a Specified Area
#'
#' Identifies and removes timepoints when tracked individuals were observed 
#'    outside of a defined polygon (note: the polygon should be described by 
#'    the vectors confinementCoord.x (x coordinates) and confinementCoord.y 
#'    (y coordinates). These vectors must be the same length and the 
#'    coordinates should be listed in the clockwise or counter-clockwise 
#'    order that they are observed on the confining polygon.
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
#' @param x Data frame or non-data-frame list that will be filtered.
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
#' @return If filterOutput == TRUE, returns \code{x} less observations where
#'    points were located outside of the polygon defined by points in 
#'    \code{confinementCoord.x} and \code{confinementCoord.y}.
#'    
#'    If filterOutput == FALSE, returns \code{x} appended with a 
#'    "confinement_status" column which reports the results of 
#'    sp::point.in.polygon function, which is used to determine if observed 
#'    points are confined within the polygon defined by points in 
#'    \code{confinementCoord.x} and \code{confinementCoord.y}.
#' @export
#' @examples
#' data("calves")
#' 
#' water_trough.x<- c(61.43315, 61.89377, 62.37518, 61.82622) #water polygon x-coordinates
#' water_trough.y<- c(62.44815, 62.73341, 61.93864, 61.67411) #water polygon y-coordinates
#' 
#' headWater1<- confine(calves, point.x = calves$x, point.y = calves$y, 
#'    confinementCoord.x = water_trough.x, confinementCoord.y = water_trough.y,
#'    filterOutput = TRUE) #creates a data set comprised ONLY of points within the water polygon.
#'    
#' headWater2<- confine(calves, point.x = calves$x, point.y = calves$y, 
#'    confinementCoord.x = water_trough.x, confinementCoord.y = water_trough.y,
#'    filterOutput = FALSE) #appends the "confinement_status" column to x.

confine <- function(x, point.x = NULL, point.y = NULL, confinementCoord.x, confinementCoord.y, filterOutput = TRUE){
  
  confinement_status<-NULL #bind this variable to a local object so that R CMD check doesn't flag it.
  
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
