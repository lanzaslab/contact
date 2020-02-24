#' Append TimeBlock Information to a Data Frame
#'
#' Appends "block," "block.start," "block.end," and "numBlocks" columns to an 
#'    input data frame (x) with a dateTime (see dateTime.append) column. This 
#'    allows users to "block" data into blockLength-blockUnit-long 
#'    (e.g., 10-min-long) temporal blocks. If x == NULL, the function output 
#'    will be a data frame with "dateTime" and block-related columns.
#' 
#' This is a sub-function that can be found in the contactDur functions.
#' @param x Data frame containing dateTime information, and to which block 
#'    information will be appended. if NULL, dateTime input relies solely on 
#'    the dateTime argument.
#' @param dateTime Vector of length nrow(x) or singular character data, 
#'    detailing the relevant colname in x, that denotes what dateTime 
#'    information will be used. If argument == NULL, the function assumes a 
#'    column with the colname "dateTime" exists in x. Defaults to NULL.
#' @param blockUnit Character string taking the values, "secs," "mins," 
#'    "hours," "days," or "weeks." Defaults to "hours."
#' @param blockLength Numerical. Describes the number blockUnits within each 
#'    temporal block. Defaults to 10.
#' @keywords data-processing sub-function
#' @return Appends the following columns to \code{x}.
#'    
#'    \item{block}{Integer ID describing unique blocks of time of pre-specified
#'    length.}
#'    \item{block.start}{The timepoint in \code{x} at which the \code{block}
#'    begins.}
#'    \item{block.end}{The timepoint in \code{x} at which the \code{block}
#'    ends.}
#'    \item{numBlocks}{Integer describing the total number of time blocks 
#'    observed within \code{x} at which the \code{block}}
#' @export
#' @examples
#' data("calves")
#' calves.dateTime<-contact::datetime.append(calves, date = calves$date, 
#'    time = calves$time) #add dateTime identifiers for location fixes.
#' calves.block<-contact::timeBlock.append(x = calves.dateTime, 
#'     dateTime = calves.dateTime$dateTime, blockLength = 10, 
#'     blockUnit = "mins")
#' head(calves.block) #see that block information has been appended.

timeBlock.append<-function(x = NULL, dateTime = NULL, blockLength = 10, blockUnit = "mins"){
  
  if(blockUnit == "Secs" || blockUnit == "SECS" || blockUnit == "secs"){
    blockLength1 <- blockLength
  }
  if(blockUnit == "Mins" || blockUnit == "MINS" || blockUnit == "mins"){
    blockLength1 <- blockLength*60 #num seconds in a minute
  }
  if(blockUnit == "Hours" || blockUnit == "HOURS" || blockUnit == "hours"){
    blockLength1 <- blockLength*60*60 #num seconds in an hour
  }
  if(blockUnit == "Days" || blockUnit == "DAYS" || blockUnit == "days"){
    blockLength1 <- blockLength*60*60*24 #num seconds in a day
  }
  if(blockUnit == "Weeks" || blockUnit == "WEEKS" || blockUnit == "weeks"){
    blockLength1 <- blockLength*60*60*24*7 #num seconds in a week
  }
  
  if(length(x) == 0){ #if there is no x input (i.e., x == NULL), #assumes that if x == NULL, dateTime does not.
    x<- data.frame(dateTime = dateTime)
  }else{ # length(x) > 0
    if(length(dateTime) > 0){ #dateTime == NULL, the function assumes that there is a "dateTime" column in x.
      if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
        x$dateTime <- x[,match(dateTime, names(x))]
      }else{ #if length(dateTime) > 1
        x$dateTime = dateTime
      }
    }
  }

  daySecondList = lubridate::hour(x$dateTime) * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
  lub.dates = lubridate::date(x$dateTime)
  x<-x[order(lub.dates, daySecondList),] #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data
  rm(list = c("daySecondList", "lub.dates")) #remove these objects because they are no longer needed.
  x$totalSecond<- difftime(x$dateTime ,x$dateTime[1] , units = c("secs")) #adds the total second column to the dataframe

  studySecond <- (x$totalSecond -min(x$totalSecond)) + 1
  x<-x[,-match("totalSecond", names(x))]
  numblocks <- ceiling((max(studySecond) - 1)/blockLength1)
  block <-rep(0,length(studySecond))
  for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
    block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
  }
  if(length(which(block == 0)) > 0){ #identifies the last block
    block[which(block == 0)] = numblocks
  }
  blockVec<-unique(block)
  dateTimeVec2<-x$dateTime #dateTime after x is sorted
  minBlockTimeSeq <- rep(0, length(block)) #Added 2/4/2019. This vector will identify the minimum timepoint in each block.
  maxBlockTimeSeq <- rep(0, length(block)) #Added 2/4/2019. This vector will identify the maximum timepoint in each block.
  for(f in 1:length(blockVec)){
    minBlockTime<-as.character(dateTimeVec2[min(which(block == blockVec[f]))])
    minBlockTimeSeq[which(block == blockVec[f])] <- minBlockTime
    maxBlockTime<-as.character(dateTimeVec2[max(which(block == blockVec[f]))])
    maxBlockTimeSeq[which(block == blockVec[f])] <- maxBlockTime
  }
  x$block <- block
  x$block.start <- minBlockTimeSeq
  x$block.end <- maxBlockTimeSeq
  x$numBlocks <- max(blockVec) #the contactTest function will require thus information (i.e. the number of blocks in the dataset)
  return(x)
}
