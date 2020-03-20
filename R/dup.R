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
#'    steps wherein individuals were observed in different locations 
#'    concurrently are removed from the data set. 
#' @param parallel Logical. If TRUE, sub-functions within the dup wrapper will 
#'    be parallelized. This is only relevant if avg == TRUE. Defaults to 
#'    FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to the maximum number of cores available
#'    (i.e., (parallel::detectCores()/2)).
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
#' @import foreach  
#' @export
#' @examples
#' 
#' data(calves2018) #load the data set
#' 
#' calves_dup<- dup(calves2018, id = calves2018$calftag, 
#'    point.x = calves2018$x, point.y = calves2018$y, 
#'    dateTime = calves2018$dateTime, avg = FALSE, parallel = FALSE, 
#'    filterOutput = TRUE) #there were no duplicates to remove in the first place.

dup <- function(x, id = NULL, point.x = NULL, point.y = NULL, dateTime = NULL, avg = TRUE, parallel = FALSE, nCores = (parallel::detectCores()/2), filterOutput = TRUE){

  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  l <- NULL
    
  #write the sub-functions
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
    
    #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
    daySecondList = lubridate::hour(dateTimeVec1) * 3600 + lubridate::minute(dateTimeVec1) * 60 + lubridate::second(dateTimeVec1) #This calculates a day-second
    lub.dates = lubridate::date(dateTimeVec1)
    x<-x[order(idVec1, lub.dates, daySecondList),] #order x 
    
    #x<-x[order(idVec1, dateTimeVec1),] #this sorts the data so that future processes will work

    rownames(x) <-seq(1,nrow(x),1)
    x$indiv_dateTimes = paste(x$idVec1, x$dateTimeVec1, sep = " ")
    a=duplicated(x$indiv_dateTimes)
    duplicates = which(a == TRUE)

    if(length(duplicates) > 0){ #This if statement prevents an error from occurring due to the lack of duplicated timepoints
      
      fullVector <- paste(x$idVec1, x$dateTimeVec1, xVec1, yVec1, sep = " ") #create a vector of all inputs
      exactDuplicates <- which(duplicated(fullVector) == TRUE) #identify which rows represent complete duplicates (i.e., all inputs are duplicated). 

      if(length(exactDuplicates) > 0){ #if there are exact duplicates, they will simply be removed later on. There's no need to calculate the average (if avg == TRUE).
        exactDup.values <- droplevels(x[exactDuplicates,]) #pull a data frame of exactDuplicates (will only be relevant if filter.output == FALSE)
        x <- droplevels(x[-exactDuplicates,]) #remove exact duplicates from x
        duplicates.adjusted<- which(duplicated(x$indiv_dateTimes) == TRUE) #re-evaluate duplicates
      }else{
        duplicates.adjusted<- duplicates #if there were no exact duplicates, then this vector need not change
      }
      
      if(length(duplicates.adjusted) > 0){
      
        dupRemove <- unique(c((duplicates.adjusted - 1), duplicates.adjusted)) #compiles a vector of rows that are duplicated

        if(filterOutput == TRUE){

          if(avg == TRUE){ 
          
            dupFixer2<-function(x,y){
              oldX = y$xVec1[which(y$indiv_dateTimes == x[1])]
              oldY = y$yVec1[which(y$indiv_dateTimes == x[1])]
              newX = (sum(oldX)/length(oldX)) #calculates the average x location and adds it to the replacement row
              newY = (sum(oldY)/length(oldY)) #calculates the average y location and adds it to the replacement row
              newCoord = c(newX,newY)
              return(newCoord)
            }
          
            dupFrame = data.frame(unique(x$indiv_dateTimes[duplicates.adjusted]), stringsAsFactors = TRUE)
            originTab = x
          
            if(parallel == TRUE){
              cl<-parallel::makeCluster(nCores)
              on.exit(parallel::stopCluster(cl))
              dupReplace = unlist(parallel::parApply(cl,dupFrame,1,dupFixer2, originTab))
            }else{ #if parallel == FALSE
              dupReplace = unlist(apply(dupFrame,1,dupFixer2, originTab))
            }

            replaceTab = x[duplicates.adjusted,] #This creates a dataframe with all the relevant data included (e.g., id, dateTime,etc.). We just need to adjust the x and y coordinates in the table.
            newXYMat = matrix(dupReplace, nrow = (nrow(replaceTab)), ncol = 2, byrow = TRUE)
            replaceTab$xVec1 = newXYMat[,1]
            replaceTab$yVec1 = newXYMat[,2]
          }
        
          x <- droplevels(x[-dupRemove,]) 
        
        if(avg == TRUE){
          x = data.frame(data.table::rbindlist(list(x,replaceTab)), stringsAsFactors = TRUE)

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

        x$duplicated = 0
        x$duplicated[dupRemove] = 1
        if(length(exactDuplicates) > 0){ #if there were exact duplicates, then they are indicated as well.
          exactDup.values$duplicated <- 1
          x <- data.frame(data.table::rbindlist(list(x, exactDup.values)), stringsAsFactors = TRUE) #bind x and exactDup.values
          x<-x[order(idVec1, lub.dates, daySecondList),] #this sorts the data, putting duplicates back into place 
          #x<-x[order(idVec1, dateTimeVec1),] #this sorts the data, putting duplicates back into place
          x$duplicated[exactDuplicates -1] <- 1
        }
      }
      }else{ #if duplicates.adjusted == 0
        if(filterOutput == FALSE){ #in this case ALL duplicates were exact duplicates
          x$duplicated = 0
          exactDup.values$duplicated <- 1
          x <- data.frame(data.table::rbindlist(list(x, exactDup.values)), stringsAsFactors = TRUE) #bind x and exactDup.values
          x<-x[order(idVec1, lub.dates, daySecondList),] #this sorts the data, putting duplicates back into place 
          #x<-x[order(idVec1, dateTimeVec1),] #this sorts the data, putting duplicates back into place
          x$duplicated[exactDuplicates -1] <- 1
        }
      }
    }

    if(nrow(x) > 0){
      rownames(x) <-seq(1,nrow(x),1)
    }
    
    x <- x[,-match("indiv_dateTimes",names(x))]
    x <- x[,-match("idVec1",names(x))]
    x <- x[,-match("dateTimeVec1",names(x))]
    x <- x[,-match("xVec1",names(x))]
    x <- x[,-match("yVec1",names(x))]

    return(x)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    
    list.dup <- foreach::foreach(l = seq(from = 1, to = length(x), by = 1)) %do% filter1.func(x[[l]], id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here.
    
    return(list.dup)

  }else{ #if x is a dataFrame

    frame.dup<- filter1.func(x, id, point.x, point.y, dateTime, avg, parallel, filterOutput, nCores)

    return(frame.dup)
  }

}
