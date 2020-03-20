#' Calculate Time Difference Between Relocations
#'
#' This function calculates the time difference between relocation events, 
#'    accounting for individuals' ids. This function has the capability to 
#'    calculate the differences between sequential timepoints related to two 
#'    different features (e.g., contactStartTime and contactEndTime) if both 
#'    dateTime1 and dateTime2 are defined, or just sequential timepoints from a
#'    single vector (e.g., contactStartTime) if only dateTime1 is defined.
#'    
#' This is a sub-function contained within contactDur variants and contactTest 
#'    functions.
#' @param x data frame containing time data. If NULL at least dateTime must be 
#'    defined. Defaults to NULL.
#' @param id Vector of length nrow(data.frame(x)) that denotes what unique ids for 
#'    tracked individuals will be used. If argument == NULL, the function 
#'    assumes a column with the colname "id" exists in x. Defaults to NULL.
#' @param dateTime1 Vector of length nrow(data.frame(x)) or singular 
#'    character data, detailing the relevant colname in x, that denotes what 
#'    dateTime information will be used. If argument == NULL, the function 
#'    assumes a column with the colname "dateTime" exists in x. Defaults to 
#'    NULL.
#' @param dateTime2 Vector of length nrow(data.frame(x)) or singular 
#'    character data, detailing the relevant colname in x, that denotes what 
#'    dateTime information will be used. If argument == NULL, the function 
#'    will calculate differences between sequential timepoints in dateTime1. If
#'    != NULL, the function will calculate differences between dateTime1 and 
#'    dateTime2 values. Defaults to NULL.
#' @param timeUnits Chracter string describing the time unit of calculated 
#'    differences. It takes the values "secs," "mins," "hours," "days," or 
#'    "weeks." Defaults to "secs."
#' @param parallel Logical. If TRUE, sub-functions within the dt.calc wrapper 
#'    will be parallelized. Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
#' @param timeStepRelation Numerical. Takes the value "1" or "2." If 
#'    argument == "1," dt values in output represent the difference between 
#'    time t and time t-1. If argument == "2," dt values in output represent 
#'    the difference between time t and time t+1. Defaults to 1.
#' @keywords data-processing contact sub-function
#' @return Output is a data frame with the following columns
#' 
#'    \item{id}{The unique ID of a tracked individual.}
#'    \item{dt}{Time difference between relocation events.}
#'    \item{units}{Temporal unit defined by \code{timeUnits} argument.}
#'    
#' @export
#' @examples
#' data(calves) #load calves data set
#' calves.datetime<-datetime.append(calves)
#' dt<-dt.calc(x = calves.datetime, id = calves.datetime$calftag, 
#'    dateTime1 = calves.datetime$dateTime, dateTime2 = NULL, 
#'    timeUnits = "secs", parallel = FALSE, timeStepRelation = 1)
#'    
#' head(dt)
   
dt.calc<-function(x = NULL, id = NULL, dateTime1 = NULL, dateTime2 = NULL, timeUnits = "secs", parallel = FALSE, nCores = (parallel::detectCores()/2), timeStepRelation = 1){
  timeDifference = function(x, timeUnits){
    t1 = unname(unlist(x[1]))
    t2 = unname(unlist(x[2]))
    dt = as.integer(difftime(time1 = t2, time2 = t1, units = timeUnits))
    return(dt)
  }
  idBreak = function(x, originTab, timeUnits){
    breakVec1 <- originTab$dateTime1[which(originTab$id == unname(unlist(x[1])))]
    breakVec2 <- originTab$dateTime2[which(originTab$id == unname(unlist(x[1])))]
    timesFrame = data.frame(breakVec1[1:(length(breakVec1) - 1)], breakVec2[2:length(breakVec2)], stringsAsFactors = TRUE)
    timedif <- apply(timesFrame,1,timeDifference, timeUnits)
    if(timeStepRelation == 1){ #dt values represent the diference between time i and time i-1
      timedif = c(NA, timedif)
    }
    if(timeStepRelation == 2){ #dt values represent the diference between time i and time i+1
      timedif = c(timedif, NA)
    }
    timeTab = data.frame(id = unlist(rep(unname(unlist(x[1])),length(timedif))), dt = timedif, stringsAsFactors = TRUE)
    return(timeTab)
  }
  
  if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime), a dataframe with columns named the same, or a combination of dataframe and vectors.
    x = data.frame(id = id, dateTime1 = dateTime1, stringsAsFactors = TRUE) #define x if the input is NULL
    if(length(dateTime2) >0){
      x$dateTime2 <- dateTime2
    }else{
      x$dateTime2 <- dateTime1
    }
  }
  
  if(length(x) > 0){ #for some reason using an "else" statement would always result in a table with 0 records...
    if(length(id) > 0){
      x$id = id
    }
    if(length(dateTime1 ) > 0){
      x$dateTime1 = dateTime1
    }
    if(length(dateTime2) >0){
      x$dateTime2 <- dateTime2
    }else{
      x$dateTime2 <- dateTime1
    }
  }
  
  #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
  daySecondList = lubridate::hour(x$dateTime1) * 3600 + lubridate::minute(x$dateTime1) * 60 + lubridate::second(x$dateTime1) #This calculates a day-second
  lub.dates = lubridate::date(x$dateTime1)
  x<-x[order(x$id, lub.dates, daySecondList),] #order x 
  
  idVecFrame<-data.frame(unique(x$id), stringsAsFactors = TRUE)
  
  if (parallel == TRUE){
    cl<-parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(cl))
    dtTime<-parallel::parApply(cl, idVecFrame, 1, idBreak,originTab = x, timeUnits)
  }else{
    dtTime = apply(idVecFrame, 1, idBreak,originTab = x, timeUnits)	
  }
  
  dt.final <- data.frame(data.table::rbindlist(dtTime), stringsAsFactors = TRUE)
  dt.final$units = timeUnits
  return(dt.final)
}
