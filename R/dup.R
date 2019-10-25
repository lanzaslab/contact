#' Identify and Remove Duplicated Data Points
#'
#' dup (a.k.a. Multiple instance filter) identifies and removes timepoints when
#'    tracked individuals were observed in >1 place concurrently. If avg == 
#'    TRUE, duplicates are replaced by a single row describing an individuals' 
#'    average location (e.g., planar xy coordinates) during the duplicated time
#'    point. If avg == FALSE, all duplicated timepoints will be removed, as 
#'    there is no way for the function to determine which instance among the 
#'    duplicates should stay. If users are not actually interested in filtering
#'    datasets, but rather, determining what observations should be filtered, 
#'    they may set filterOutput == FALSE. By doing so, this function will 
#'    append a "duplicated" column to the dataset, which reports values that 
#'    describe if any timepoints in a given individual's path are duplicated. 
#'    Values are: 0: timepoint is not duplicated, 1: timepoint is duplicated.
#'
#' If users want to remove specific duplicated observations, we suggest setting
#'    filterOutput == FALSE, reviewing what duplicated timepoints exist in 
#'    individuals' paths, and manually removing observations of interest.
#' @param x Data frame containing real-time-location data that will be filtered.
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
#' @param avg Logical. If TRUE, point.x and point.y values for duplicated 
#'    time steps will be averaged, producing a singular point for all time 
#'    steps in individuals' movement paths. If FALSE, all duplicated time 
#'    steps are removed from the data set. 
#' @param parallel Logical. If TRUE, sub-functions within the dup wrapper will 
#'    be parallelized. Note that this can significantly speed up processing of 
#'    relatively small data sets, but may cause R to crash due to lack of 
#'    available memory when attempting to process large datasets. Defaults to 
#'    FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to the maximum number of cores available
#'    (i.e., parallel::detectCores()).
#' @param filterOutput Logical. If TRUE, output will be a data frame 
#'    containing only movement paths with non-duplicated timesteps. If FALSE, 
#'    no observations are removed and a "duplicated" column is appended to x, 
#'    detailing if time steps are duplicated (column value == 1), or not 
#'    (column value == 0). Defaults to TRUE.
#' @keywords filter duplicates
#' @return If filterOutput == TRUE, returns \code{x} less observations at 
#'    duplicated timepoints.
#'    
#'    If filterOutput == FALSE, returns \code{x} appended with a 
#'    "duplicated" column which reports timepoints are duplicated (column 
#'    value == 1), or not (column value == 0).
#' @export
#' @examples
#' 
#' data(calves2018) #load the data set
#' 
#' calves_dup<- dup(calves2018, id = calves2018$calftag, 
#'    point.x = calves2018$x, point.y = calves2018$y, 
#'    dateTime = calves2018$dateTime, avg = FALSE, parallel = FALSE, 
#'    filterOutput = TRUE) #there were no duplicates to remove in the first place.

dup <- function(x, id = NULL, point.x = NULL, point.y = NULL, dateTime = NULL, avg = TRUE, parallel = FALSE, nCores = parallel::detectCores(), filterOutput = TRUE){
  filter1.func<-function(x, id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores){

    idVec <- NULL
    xVec <- NULL
    yVec <- NULL
    dateTimeVec <- NULL
    colname.x <-NULL
    colname.y <-NULL

    if(length(id) > 0){
      if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #added 1/21 Rather than id being a vector of length(nrow(x)), it may be more convenient to designate the colname for intended "id" values
        idVec <- x[,match(id, names(x))]
      }else{ #if length(id) > 1
        idVec <- id
      }
    }

    if(length(point.x) > 0){
      if(length(point.x) == 1 && is.na(match(point.x[1], names(x))) == FALSE){
        xVec <- x[,match(point.x, names(x))]
        colname.x <- point.x
      }else{ #if length(point.x) > 1
        xVec <- point.x
      }
    }
    if(length(point.y) > 0){
      if(length(point.y) == 1 && is.na(match(point.y[1], names(x))) == FALSE){
        yVec <- x[,match(point.y, names(x))]
        colname.y <- point.y
      }else{ #if length(point.y) > 1
        yVec <- point.y
      }
    }

    if(length(dateTime) > 0){

      if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){
        dateTimeVec <- x[,match(dateTime, names(x))]
      }else{ #if length(dateTime) > 1
        dateTimeVec <- dateTime
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

    if(length(idVec) > 0){
      idVec1 <- idVec
    }else{
      idVec1 <- x$id
    }

    if(length(dateTimeVec) > 0){
      dateTimeVec1 <- dateTimeVec
    }else{
      dateTimeVec1 <- x$dateTime
    }

    x$idVec1 <- idVec1
    x$dateTimeVec1 <- dateTimeVec1
    x$xVec1 <- xVec1
    x$yVec1 <- yVec1

    x<-x[order(idVec1, dateTimeVec1),] #this sorts the data so that future processes will work

    rownames(x) <-seq(1,nrow(x),1)
    x$indiv_dateTimes = paste(x$idVec1, x$dateTimeVec1, sep = " ")
    originTab = x
    a=duplicated(x$indiv_dateTimes)
    duplicates = which(a == TRUE)

    if(length(duplicates) > 0){ #This if statement prevents an error from occurring due to the lack of duplicated timepoints
      dupFrame = data.frame(unique(x$indiv_dateTimes[duplicates]))
      dupFixer1<-function(x,y){
        removeVec = which(y$indiv_dateTimes == x[1])
        return(removeVec)
      }

      dupFixer2<-function(x,y){
        oldX = y$xVec1[which(y$indiv_dateTimes == x[1])]
        oldY = y$yVec1[which(y$indiv_dateTimes == x[1])]
        newX = (sum(oldX)/length(oldX)) #calculates the average x location and adds it to the replacement row
        newY = (sum(oldY)/length(oldY)) #calculates the average y location and adds it to the replacement row
        newCoord = c(newX,newY)
        return(newCoord)
      }

      if(filterOutput == TRUE){

        if(parallel == TRUE){
          cl<-parallel::makeCluster(nCores)
          dupRemove = unlist(parallel::parApply(cl,dupFrame,1,dupFixer1, originTab))
          if(avg == TRUE){
            dupReplace = unlist(parallel::parApply(cl,dupFrame,1,dupFixer2, originTab))
          }
          parallel::stopCluster(cl)
        }else{ #if parallel == FALSE
          dupRemove = unlist(apply(dupFrame,1,dupFixer1, originTab)) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point travelled to reach it from the subsequent point.
          if(avg == TRUE){
            dupReplace = unlist(apply(dupFrame,1,dupFixer2, originTab))
          }
        }
        if(avg == TRUE){
          replaceTab = x[duplicates,] #This creates a dataframe with all the relevant data included (e.g., id, dateTime,etc.). We just need to adjust the x and y coordinates in the table.
          newXYMat = matrix(dupReplace, nrow = (nrow(replaceTab)), ncol = 2, byrow = TRUE)
          replaceTab$xVec1 = newXYMat[,1]
          replaceTab$yVec1 = newXYMat[,2]
        }
        x = x[-dupRemove,]
        if(avg == TRUE){
          x = data.frame(data.table::rbindlist(list(x,replaceTab)))

          if(length(colname.x) > 0){ #colname.x and colname.y are the modified columns, yet these columns are artifacts of this function that will ultimately be removed. This step ensures that the relevant permanent columns are modified as well.
            x[,match(colname.x,names(x))] <- x$xVec1
          }else{ #if length(colnameNum.x) == 0
            x$x <- x$xVec1
          }
          if(length(colname.y) > 0){ #colname.x and colname.y are the modified columns, yet these columns are artifacts of this function that will ultimately be removed. This step ensures that the relevant permanent columns are modified as well.
            x[,match(colname.y,names(x))] <- x$yVec1
          }else{ #if length(colnameNum.y) == 0
            x$y <- x$yVec1
          }
        }

      }else{ #i.e., if filterOutput == FALSE

        if(parallel == TRUE){
          cl<-parallel::makeCluster(nCores)
          dupRemove = unlist(parallel::parApply(cl,dupFrame,1,dupFixer1, originTab)) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point travelled to reach it from the subsequent point.
          parallel::stopCluster(cl)
        }else{
          dupRemove = unlist(apply(dupFrame,1,dupFixer1, originTab)) #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point travelled to reach it from the subsequent point.
        }
        x$duplicated = 0
        x$duplicated[dupRemove] = 1
      }

      if(nrow(x) > 0){
        rownames(x) <-seq(1,nrow(x),1)
      }
    }

    x <- x[,-match("indiv_dateTimes",names(x))]
    x <- x[,-match("idVec1",names(x))]
    x <- x[,-match("dateTimeVec1",names(x))]
    x <- x[,-match("xVec1",names(x))]
    x <- x[,-match("yVec1",names(x))]

    return(x)
  }
  list.breaker<-function(x,y,id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores){
    input<- data.frame(y[unname(unlist(x[1]))])
    dup.filter<-filter1.func(input, id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores)
    return(dup.filter)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.dup <- apply(breakFrame, 1, list.breaker,y = x, id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.

    return(list.dup)

  }else{ #if x is a dataFrame

    frame.dup<- filter1.func(x, id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores)

    return(frame.dup)
  }

}
