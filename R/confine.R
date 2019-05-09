#' Identify and Remove Data Points Outside of a Specified Area
#'
#' dup (a.k.a. Multiple instance filter) identifies and removes timepoints when tracked individuals were observed in >1 place concurrently. If avg == TRUE, duplicates are replaced by a single row describing an individuals' average location (e.g., planar xy coordinates) during the duplicated time point. If avg == FALSE, all duplicated timepoints will be removed, as there is no way for the function to determine which instance among the duplicates should stay. If users are not actually interested in filtering datasets, but rather, determining what observations should be filtered, they may set filterOutput == FALSE. By doing so, this function will append a "duplicated" column to the dataset, which reports values that describe if any timepoints in a given individual's path are duplicated. Values are: 0: timepoint is not duplicated, 1: timepoint is duplicated.
#'
#' Identifies and removes timepoints when tracked individuals were observed outside of a defined polygon (note: the polygon should be described by the vectors confinementCoord.x (x coordinates) and confinementCoord.y (y coordinates) (note: these vectors must be the same length and the coordinates should be listed in the clockwise or counter-clockwise order that they are observed on the confining polygon), which identify where polygon vertices exist).
#'
#' If users are not actually interested in filtering datasets, but rather, determining what observations should be filtered, they may set filterOutput == FALSE. By doing so, this function will append a "confinement_status" column to the output dataframe, which reports the results of sp::point.in.polygon function that is used to determine if individuals are confined within a given polygon. In this column, values are: 0: point is strictly exterior to pol; 1: point is strictly interior to pol; 2: point lies on the relative interior of an edge of pol; 3: point is a vertex of pol (https://www.rdocumentation.org/packages/sp/versions/1.3-1/topics/point.in.polygon).
#' @param x List or data frame that will be filtered.
#' @param id Vector of length(nrow(data.frame(x))) or singular character data, detailng the relevant colname in x, that denotes what date information will be used. If argument == NULL, datetime.append assumes a column withe colname "id" exists in x. Defaults to NULL.
#' @param point.x Description imminent
#' @param point.y Description imminent
#' @param confinementCoord.x Description imminent
#' @param confinementCoord.y Description imminent
#' @param filterOutput Description imminent
#' @keywords filter confinement polygon
#' @export
#' @examples
#' Examples imminent

confine <- function(x, point.x = NULL, point.y = NULL, confinementCoord.x, confinementCoord.y, filterOutput = TRUE){
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
