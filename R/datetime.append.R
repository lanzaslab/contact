#' Append Date-Time Information to a Dataset
#'
#' This function appends date-time information to a dataset in POSIXct date_time format. It also uses functions from the lubridate package and minor calculations to parse out month, day, hour, minute, second, daySecond (the sequentially ordered second of a day), and totalSecond (sequentially ordered second over the course of the study period) of observations in a given dataset with date (format: "mdy" =  month/day/year, "ymd" =  year/month/day, "dmy" =  day/month/year, or "ydm" =  year/day/month (note: no preceding zeroes should be included before numbers <10)) and time (format: hour:minute:second (note:preceding zeroes must be included before numbers < 10, ex. 00:00:01)) information, appends this metadata to the dataset,and can assign each day a unique ID.
#' @param x List or data frame to which new information will be appended.
#' @param date Vector of length(nrow(data.frame(x))) or singular character data, detailng the relevant colname in x, that denotes what date information will be used. If argument == NULL, datetime.append assumes a column withe colname "date" exists in x, or that the dateTime argument != NULL. Defaults to NULL.
#' @param time Vector of length(nrow(data.frame(x))) or singular character data, detailng the relevant colname in x, that denotes what time information will be used. If argument == NULL, datetime.append assumes a column withe colname "time" exists in x, or that the dateTime argument != NULL. Defaults to NULL.
#' @param dateTime Description imminent
#' @param dateFormat Description imminent
#' @param dateFake Logical. If dateFake == TRUE, the function will assign fake date information, beginning 01/01/startYear, to each of the timestamps. Defaults to FALSE.
#' @param startYear Description imminent
#' @param changeTimezone Description imminent
#' @param tz.in Character. Identifies the timezone associated with the time/dateTime argument input. Defaults to "UTC." Timezone names often take the form "Country/City." See the listing of timezones at: http://en.wikipedia.org/wiki/List_of_tz_database_time_zones.
#' @param tz.out Character. Identifies the timezone that the output dateTime information will be converted to. If NULL, tz.out will be identical to tz.in. Defaults to NULL. Timezone names often take the form "Country/City." See the listing of timezones at: http://en.wikipedia.org/wiki/List_of_tz_database_time_zones.
#' @param month Description imminent
#' @param day Description imminent
#' @param year Description imminent
#' @param hour Description imminent
#' @param minute Description imminent
#' @param second Description imminent
#' @param daySecond Description imminent
#' @param dayID Description imminent
#' @param startID Description imminent
#' @param totalSecond Description imminent
#' @keywords date time date-time
#' @export
#' @examples
#' Examples imminent

datetime.append <- function(x, date = NULL, time = NULL, dateTime = NULL, dateFormat = "mdy", dateFake = FALSE, startYear = 2000, tz.in = "UTC", tz.out = NULL, month = FALSE, day = FALSE, year = FALSE, hour = FALSE, minute = FALSE, second = FALSE, daySecond = FALSE, dayID = FALSE, startID = 1, totalSecond = FALSE){
  
  generate.append<-function(x, date, time, dateTime, dateFormat, dateFake, startYear, changeTimezone, timezone, month, day, year, hour, minute, second, daySecond, dayID, startID, totalSecond){
    dateFake.func<-function(timestamp, dateFormat, startYear){
      monthSeq<-c("01","02","03","04","05","06","07","08","09","10","11","12")
      numDays1<-NULL
      daySeqFrame<-NULL
      for(a in 1:12){
        monthDays<- lubridate::days_in_month(a)
        numDays1<-c(numDays1, unname(unlist(monthDays)))
        Seq1<-seq(1, unname(unlist(monthDays)),1)
        seqFrame<- expand.grid(monthSeq[a], Seq1)
        daySeqFrame<-data.frame(data.table::rbindlist(list(daySeqFrame,seqFrame)))
      }
      dataTime<- lubridate::hms(timestamp)
      dataHours<- lubridate::hour(dataTime) + 1 #assuming timestamp system is hours 0-23, not 1-24
      numYears<- ceiling(max(dataHours)/(24*365)) #number hours in the dataset / (number of hours in a day * number days in a year)
      endYear<-startYear + numYears - 1 #minus 1, because the first year with data will be the startYear.
      yearFrame<-NULL
      for(b in startYear:endYear){ #added 02/06/2019. No longer will this function always start on year 2001
        year <- b
        daySeqFrame$year <- year
        yearFrame<-data.frame(data.table::rbindlist(list(yearFrame,daySeqFrame)))
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

    
    if(length(dateTime) > 0){
      if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){
        timevec <- x[,match(dateTime, names(x))]
      }else{ #if length(dateTime) > 1
        timevec <- dateTime
      }
      rm.dateTime <- TRUE #if there is already dateTime information present, new dateTime information will not be appended.
      
      #convert timevec to POSIXct format
      if(dateFormat == "mdy"){timevec = lubridate::mdy_hms(timevec)}
      if(dateFormat == "ymd"){timevec = lubridate::ymd_hms(timevec)}
      if(dateFormat == "dmy"){timevec = lubridate::dmy_hms(timevec)}
      if(dateFormat == "ydm"){timevec = lubridate::ydm_hms(timevec)}
      
    }else{ #if length(dateTime) == 0
      
      rm.dateTime <- FALSE
      
      if(length(time) > 0){
        if(length(time) == 1 && is.na(match(time[1], names(x))) == FALSE){ #added 2/02 Rather than time being a vector of length(nrow(x)), it may be more convenient to designate the colname for intended "time" values
          time.append <- x[,match(time, names(x))]
        }else{ #if length(time) > 1
          time.append <- time
        }
      }else{ #if length(time) == 0
        time.append = x$time #if time == NULL, the function assumes there is a column in x named "time"
      }
      
      time.append<-as.character(time.append)
      
      if(dateFake == TRUE){
        date.append<- dateFake.func(time.append, dateFormat, startYear)
      }else{ #if dateFake == FALSE
        
        if(length(date) > 0){
          if(length(date) == 1 && is.na(match(date[1], names(x))) == FALSE){
            date.append <- x[,match(date, names(x))]
          }else{ #if length(date) > 1
            date.append <- date
          }
        }else{ #if length(date) == 0
          date.append = x$date #if date == NULL, the function assumes there is a column in x named "date"
        }
      }
      
      date.time = paste(date.append, time.append, sep = " ")
      
      if(dateFormat == "mdy"){timevec = lubridate::mdy_hms(date.time)}
      if(dateFormat == "ymd"){timevec = lubridate::ymd_hms(date.time)}
      if(dateFormat == "dmy"){timevec = lubridate::dmy_hms(date.time)}
      if(dateFormat == "ydm"){timevec = lubridate::ydm_hms(date.time)}
    }


    lubridate::tz(timevec) <- tz.in #sets the appropriate input timezone. Defaults to UTC.
    
    if(is.null(tz.out) == FALSE){ #If tz.out is anything other than NULL, the output timezone will be changed to tz.out
      timevec =	lubridate::with_tz(timevec,tz=tz.out)
    }

    cbindTab<- data.frame(matrix(ncol = 0, nrow = length(timevec)))

    if(rm.dateTime == FALSE){
      cbindTab$dateTime <- timevec
    }

    if(month == TRUE){
      cbindTab$month = lubridate::month(timevec)
    }
    if(day == TRUE){
      cbindTab$day = lubridate::day(timevec)
    }
    if(year == TRUE){
      cbindTab$year = lubridate::year(timevec)
    }
    if(hour == TRUE){
      cbindTab$hour = lubridate::hour(timevec)
    }
    if(minute == TRUE){
      cbindTab$minute = lubridate::minute(timevec)
    }
    if(second == TRUE){
      cbindTab$second = lubridate::second(timevec)
    }
    if(daySecond == TRUE | totalSecond == TRUE){
      daySecondVec = lubridate::hour(timevec) * 3600 + lubridate::minute(timevec) * 60 + lubridate::second(timevec) #This calculates a day-second
      if(daySecond == TRUE){
        cbindTab$daySecond <- daySecondVec
      }
    }
    if(dayID == TRUE | totalSecond == TRUE){ #in some cases, calculating specific dayIDs may not be worth the time it takes to process it. In the event that specific dayIDs are not needded for future analyses, you can opt to skip evaluating it.
      x<-x[order(timevec),] #Just in case the data wasn't already ordered in this way.
      cbindTab<-cbindTab[order(timevec),]
      daySecondVec <-daySecondVec[order(timevec)] #added 02/05/2019 after Dan pointed out that that daySecondVec was not sorted properly and was therefore not creating accurate totalSeconds. I chose to sort this vector rather than direct the function to cbindTab$daySecond, because doing the latter would require users to have set daySecond == TRUE. 
      timevec <-timevec[order(timevec)]
      lub.dates = lubridate::date(timevec)
      dateseq = unique(lub.dates)
      dayIDVec = NULL
      dayIDseq = seq(startID,(length(dateseq) + startID - 1),1)
      for(b in dayIDseq){
        ID = rep(b,length(which(lub.dates == dateseq[which(dayIDseq == b)])))
        dayIDVec = c(dayIDVec, ID)
      } #This part of the function takes awhile (especially for large datasets), but may be useful in the future for subsetting and viewing data.
      if(dayID == TRUE){
        cbindTab$dayID = dayIDVec
      }
    }

    if(totalSecond == TRUE){ #assumes dayIDs are sequential
      cbindTab$totalSecond = ((dayIDVec - min(dayIDVec))*86400) + daySecondVec #This calculates the total second (the cumulative second across the span of the study's timeframe)
    }
    bindlist<-list(x,cbindTab)
    appendedTab<-do.call("cbind", bindlist)
    return(appendedTab)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    list.dateTime <- lapply(x, generate.append, date, time, dateTime, dateFormat, dateFake, startYear, changeTimezone, timezone, month, day, year, hour, minute, second, daySecond, dayID, startID, totalSecond)
    return(list.dateTime)
  }else{ #if x is a dataFrame
    frame.dateTime<- generate.append(x, date, time, dateTime, dateFormat, dateFake, startYear, changeTimezone, timezone, month, day, year, hour, minute, second, daySecond, dayID, startID, totalSecond)
    return(frame.dateTime)
  }
}
