#' Create Fake Date Information
#'
#' This function assigns fake date information, beginning 01/01/startYear, to 
#'    each empirical timestamp. Users can control what format the output vector
#'    is in by changing the dateFormat argument 
#'    (format: "mdy" =  month-day-year, "ymd" =  year-month-day, 
#'    "dmy" =  day-month-year, or "ydm" =  year-day-month).
#'    
#' This is a sub-function that can be found within datetime.append.
#' 
#' Note that the timestamp argument should be a vector of all relevant 
#'    timepoints. Additionally, timepoints should be in hms ("hour, minute, 
#'    second") format.
#' @param timestamp Vector of time information with format 
#'    "hour:minute:second."
#' @param dateFormat Character string. Defines how date information will be 
#'    presented in output. Takes values "mdy" (i.e., month/day/year), "ymd" 
#'    (i.e., year/month/day), "dmy" (i.e., day/month/year), or "ydm" 
#'    (i.e., year/day/month). Defaults to "mdy."  
#' @param startYear Numerical. Denotes what year fake date information will 
#' begin if dateFake == TRUE. Defaults to 2000.
#' @keywords data-processing sub-function
#' @return Output is a vector of date values (e.g., "01-1-2000") with length 
#'    length(\code{timestamp}).
#' @export
#' @examples
#' 
#' data("calves")
#' dateFake(calves$time, dateFormat = "mdy", startYear = 2000) 

dateFake<-function(timestamp, dateFormat ="mdy", startYear = 2000){
  monthSeq<-c("01","02","03","04","05","06","07","08","09","10","11","12")
  numDays1<-NULL
  daySeqFrame<-NULL
  for(a in 1:12){
    monthDays<- lubridate::days_in_month(a)
    numDays1<-c(numDays1, unname(unlist(monthDays)))
    Seq1<-seq(1, unname(unlist(monthDays)),1)
    seqFrame<- expand.grid(monthSeq[a], Seq1, stringsAsFactors = TRUE)
    daySeqFrame<-data.frame(data.table::rbindlist(list(daySeqFrame,seqFrame)), stringsAsFactors = TRUE)
  }
  dataTime<- lubridate::hms(timestamp)
  dataHours<- lubridate::hour(dataTime) + 1 #assuming timestamp system is hours 0-23, not 1-24
  numYears<- ceiling(max(dataHours)/(24*365)) #number hours in the dataset / (number of hours in a day * number days in a year)
  endYear<-startYear + numYears - 1 #minus 1, because the first year with data will be the startYear.
  yearFrame<-NULL
  for(b in startYear:endYear){ #added 02/06/2019. No longer will this function always start on year 2001
    year <- b
    daySeqFrame$year <- year
    yearFrame<-data.frame(data.table::rbindlist(list(yearFrame,daySeqFrame)), stringsAsFactors = TRUE)
  }    
  rowCall<- ceiling(dataHours/24)
  if(dateFormat == "mdy"){
    outSeq<- paste(yearFrame[rowCall,1],yearFrame[rowCall,2],yearFrame[rowCall,3],sep = "-")
  }
  if(dateFormat == "ymd"){
    outSeq<- paste(yearFrame[rowCall,3],yearFrame[rowCall,1],yearFrame[rowCall,2],sep = "-")
  }
  if(dateFormat == "dmy"){
    outSeq<- paste(yearFrame[rowCall,2],yearFrame[rowCall,1],yearFrame[rowCall,3],sep = "-")
  }
  if(dateFormat == "ydm"){
    outSeq<- paste(yearFrame[rowCall,3],yearFrame[rowCall,2],yearFrame[rowCall,1],sep = "-")
  }
  
  return(outSeq)
}