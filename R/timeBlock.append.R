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
#' @param blockLength Integer. Describes the number blockUnits within each 
#'    temporal block. Defaults to 1.
#' @param blockingStartTime Character string or date object describing the date
#'     OR dateTime starting point of the first time block. For example, if 
#'     blockingStartTime = "2016-05-01" OR "2016-05-01 00:00:00", the first 
#'     timeblock would begin at "2016-05-01 00:00:00." If NULL, the 
#'     blockingStartTime defaults to the minimum dateTime point in x. Note: 
#'     any blockingStartTime MUST precede or be equivalent to the minimum 
#'     timepoint in x. Additional note: If blockingStartTime is a character 
#'     string, it must be in the format ymd OR ymd hms.
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
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'    time = calves$time) #add dateTime identifiers for location fixes.
#' calves.block<-timeBlock.append(x = calves.dateTime, 
#'     dateTime = calves.dateTime$dateTime, blockLength = 10, 
#'     blockUnit = "mins")
#' head(calves.block) #see that block information has been appended.

timeBlock.append<-function(x = NULL, dateTime = NULL, blockLength = 1, blockUnit = "hours", blockingStartTime = NULL){
  
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
    x<- data.frame(dateTime = dateTime, stringsAsFactors = TRUE)
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
  x<-x[order(lub.dates, daySecondList),] #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data. It takes a bit longer (especially with larger data sets), but that's the price of accuracy
  rm(list = c("daySecondList", "lub.dates")) #remove these objects because they are no longer needed.
  
  ##totalSecond<- difftime(x$dateTime ,x$dateTime[1] , units = c("secs")) #calculate total seconds
  ##studySecond <- (totalSecond -min(totalSecond)) + 1
  
  if(length(blockingStartTime) == 1){ #if the blockingStartTime argument is defined, we calculate how far it is away (in seconds) from the minimum timepoint in x
    
    blockTimeAdjustment <- difftime(x$dateTime[1], blockingStartTime, units = c("secs"))
    
  }else{ #if the blockingStartTime argument is NOT defined, the adjustment is 0
    
    blockTimeAdjustment <- 0
    
  }
  
  #for some odd reason, difftime will output mostly zeroes (incorrectly) if there are > 1 correct 0 at the beginning. We use a crude fix here to address this. Basically, we create the zeroes first and combine it with other values afterwards
  totSecond <- rep(0, length(which(x$dateTime == x$dateTime[1])))
  if(nrow(x) > length(totSecond)){
    totSecond2<-as.integer(difftime(x$dateTime[(length(totSecond) +1): nrow(x)] ,x$dateTime[1], units = c("secs")))
  }else{
    totSecond2 <- NULL
  }
  studySecond <- as.integer((c(totSecond, totSecond2) -min(c(totSecond, totSecond2))) + 1) + blockTimeAdjustment
  
  numblocks <- as.integer(ceiling(max(studySecond)/blockLength1))
  block<- ceiling(studySecond/blockLength1)
  
  #numblocks <- as.integer(ceiling((max(studySecond) - 1)/blockLength1))
  #block <-rep(0,length(studySecond))
  #for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
  #  block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
  #}
  #if(length(which(block == 0)) > 0){ #identifies the last block
  #  block[which(block == 0)] = numblocks
  #}
  
  block.start<-as.character((as.POSIXct(x$dateTime[1]) - blockTimeAdjustment) + ((block - 1)*blockLength1)) #identify the timepoint where each block starts (down to the second resolution)
  block.end<-as.character((as.POSIXct(x$dateTime[1]) - blockTimeAdjustment) + ((block - 1)*blockLength1) + (blockLength1 -1)) #identify the timepoint where each block ends (down to the second resolution)
  
  x$block <- block
  x$block.start <- block.start
  x$block.end <- block.end
  x$numBlocks <- max(block) #the contactTest function will require this information (i.e. the number of blocks in the dataset)
  
  return(x)
}
