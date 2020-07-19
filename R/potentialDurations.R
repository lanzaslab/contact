#' Identify Potential Contact Durations
#'
#' This function uses the output from dist2... functions to determine the 
#'    potential maximum number of direct-contact durations between 
#'    individuals in a data set. The max number of durations potentially 
#'    observed is the number of TSWs both individuals (or an individual and 
#'    fixed area) were simulataneously observed at the same time over the 
#'    study period/temporal block.
#'    
#'  Please note that this function assumes the desired minimum contact 
#'     duration (MCD), as defined by Dawson et al. (2019), is 1 (i.e., a 
#'     "contact" occurs when individuals are within a specified distance 
#'     threshold for a single timestep). In a future version of this 
#'     function we will aim to increase flexability by allowing for variable 
#'     MCD values. For further clarification on the MCD definition and various
#'     contact-determination assumptions, please see:
#'     
#'    Dawson, D.E., Farthing, T.S., Sanderson, M.W., and Lanzas, C. 
#'    2019. Transmission on empirical dynamic contact networks is influenced by
#'    data processing decisions. Epidemics 26:32-42. 
#'    https://doi.org/10.1016/j.epidem.2018.08.003/
#'     
#' @param x Output from the dist2All or dist2Area function. Can be either a 
#'     data frame or non-data-frame list.
#' @param blocking Logical. If TRUE, contacts will be evaluated for temporal 
#'     blocks spanning blockLength blockUnit (e.g., 6 hours) within the data 
#'     set. Defaults to FALSE.
#' @param blockLength Integer. Describes the number blockUnits within each 
#'     temporal block. Defaults to 1.
#' @param blockUnit Character string taking the values: "secs," "mins," 
#'     "hours," "days," or "weeks." Describes the temporal unit associated with
#'     each block. Defaults to "hours."
#' @param blockingStartTime Character string or date object describing the date
#'     OR dateTime starting point of the first time block. For example, if 
#'     blockingStartTime = "2016-05-01" OR "2016-05-01 00:00:00", the first 
#'     timeblock would begin at "2016-05-01 00:00:00." If NULL, the 
#'     blockingStartTime defaults to the minimum dateTime point in x. Note: 
#'     any blockingStartTime MUST precede or be equivalent to the minimum 
#'     timepoint in x. Additional note: If blockingStartTime is a character 
#'     string, it must be in the format ymd OR ymd hms.
#' @param distFunction Character string taking the values: "dist2All_df",
#'     or "dist2Area_df." Describes the contact-package function used to
#'     generate x.
#' 
#' @keywords data-processing contact
#' @return Returns a data frame (or list of data frames if \code{x} is a 
#'    list of data frames) with the following columns:
#'    
#'    \item{id}{The unique ID of an individual observed in the data set.}
#'    \item{potenDegree}{The maximum degree possible for individual \code{id} 
#'    based on the number of other individuals observed during the time 
#'    period.}
#'    \item{potenTotalContactDurations}{The maximum number of contact durations 
#'    individual \code{id} may experience during the time period.}    
#'    \item{potenContactDurations_...}{The maximum number of contact durations 
#'    individual \code{id} may experience with each specific individual/fixed 
#'    area during the time period.}
#'    
#'    If blocking == TRUE, the following columns are appended to the output
#'    data frame described above:
#'    
#'    \item{block}{Integer ID describing unique blocks of time during which 
#'    contacts may occur.}
#'    \item{block.start}{The timepoint in \code{x} at which the \code{block}
#'    begins.}
#'    \item{block.end}{The timepoint in \code{x} at which the \code{block}
#'    ends.}
#' @import foreach
#' @export
#' @examples
#' 
#' data(calves)
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, time =
#'     calves$time) #create a dataframe with dateTime identifiers for location foxes
#'     
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag,
#'     dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x,
#'     point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE,
#'     extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE,
#'     na.rm = TRUE, smooth.type = 1) #smooth locations to 5-min fix intervals.
#' 
#' calves.dist<-dist2All_df(x = calves.agg, parallel = FALSE, dataType = "Point",
#'     lonlat = FALSE) #calculate distance between all individuals at each timepoint
#'     
#' calves.potentialContacts<-potentialDurations(x = calves.dist, blocking = FALSE)

potentialDurations<-function(x, blocking = FALSE, blockLength = 1, blockUnit = "hours", blockingStartTime = NULL, distFunction = "dist2All_df"){ 
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  listBreak <- NULL
  breakBlock <- NULL
  j <- NULL
  l <- NULL
  m <- NULL
  
  thisEnvironment<-environment() #tag the environment of the parent function so that sub-functions may work within it.
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #if the x input is a list of data frames
  
    y<-x #rename x to prevent function confusion later on
    rm(x) #remove x to restore memory
    
    listBreakFrame<-data.frame(seq(1:length(y)), stringsAsFactors = TRUE)
    
    listBreak.function <- function(m, environmentTag = thisEnvironment){
      #subEnvir1 <- environment() #tag the environment of the apply function
      assign("listBreak", m, envir = environmentTag) #assign a "listBreak" object in the parent environment that takes the same value as x in apply function
      
      eval(expr = {
        
        x<- y[[listBreak]] #break the original input into a single data frame
        y[[listBreak]] <- NULL #replace data frame of interest in the original input with NULL to free up memory
        
        if(blocking == TRUE){
    
    #the code immediately below comes from the contact::timeBlock.append function. We just refrain from calling it here to prevent the inputs being cloned within the function and because we don't need ALL the info that particular function would append to x. Thus, we save time and memory.
    
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
    
    lub.hours = lubridate::hour(x$dateTime)
    lub.dates = lubridate::date(x$dateTime)
    
    daySecondList = lub.hours * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
    x<-x[order(lub.dates, daySecondList),] #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data

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
    
    block.start<-as.character((as.POSIXct(x$dateTime[1]) - blockTimeAdjustment) + ((block - 1)*blockLength1)) #identify the timepoint where each block starts (down to the second resolution)
    block.end<-as.character((as.POSIXct(x$dateTime[1]) - blockTimeAdjustment) + ((block - 1)*blockLength1) + (blockLength1 -1)) #identify the timepoint where each block ends (down to the second resolution)
    
    x$block <- block
    x$block.start <- block.start
    x$block.end <- block.end

    rm(list = c("daySecondList", "block.start","block.end", "lub.dates", "lub.hours", "numblocks", "totSecond", "totSecond2", "studySecond", "block")) #remove these objects because they are no longer needed.
    
    #now that we've assigined the block information, we can begin processing the data in earnest
    
    #The first thing we need to do is remove any NAs in the data set (it's unlikely)
    
    
    idSeq<- unique(x$id)
    blockSeq <- unique(x$block)
    dist.colnames <-substring(colnames(x)[grep("dist.to.", colnames(x))],9) #Note specific id in "dist.to." colnames start at character 9. This will be different depending on if dist2all or dist2area was used
    dist.colnames_modified <- gsub("indiv_","",dist.colnames) #remove "indiv_" from dist2all output
    
    blockBreak.function <-function(l, environmentTag = thisEnvironment, columnNames = dist.colnames_modified){
      
      #browser()
      #Because the input data files can be quite large, we do not want to clone it excessively by running subfunctions. So, we relate the following expression to the previously-generated environment (where the data already exist)
      assign("breakBlock", unname(unlist(l)), envir = environmentTag)
      
      eval(expr = { 
        blockSub <- droplevels(x[which(x$block == breakBlock),]) #subset the data set by block
          
          fe1 <- foreach::foreach(j = idSeq, .noexport = c("idSub", "outputMat")) %do% {
            
            idSub <- droplevels(blockSub[which(blockSub$id == j),]) #subset dist.input to only contain the id-value of interest
            outputMat<- matrix(nrow = 1, ncol = (3 + length(grep("dist.to.", colnames(blockSub))))) #set up the empty matrix to hold the output data.
            colnames(outputMat)<- c("id", "potenDegree", "potenTotalContactDurations", paste("potenContactDurations_", columnNames, sep =""))  
            
            #if the individual was not present (i.e., nrow(idSub) == 0), then 0s well be put in the matrix. However, how many potential nodes present in the block takes a bit more effort to calculate because the distance input object may have come from dist2All_df or dist2Area_df.
            
            potentialDegreeIdentifier <-NULL #create an empty vector to describe when nodes are present in the data set.
            
            for(i in 1:length(grep("dist.to.", colnames(blockSub)))){ #loops through the potentially-contactable nodes.
              
              presenceTest<- is.na(blockSub[,grep("dist.to.", colnames(blockSub))[i]]) #if nodes were not present all entries will be TRUE.
              
              if(length(which(presenceTest == FALSE)) > 0){ #identify when the reported value is FALSE (i.e., when nodes WERE present)
                
                potentialDegreeIdentifier <-c(potentialDegreeIdentifier, 1) #if nodes are present add 1 value to potentialDegreeIdentifier
                
              }else{next}
              
            }
            if(distFunction == "dist2All_df"){ #if the dist2All_df function was used to generate x
              potentialDegree <- ifelse(nrow(idSub) > 0, (length(potentialDegreeIdentifier) - 1), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed minus 1 (because individual j could not be in contact with itself).
            }
            if(distFunction == "dist2Area_df"){ #if the dist2Area_df function was used to generate x, then 1 does not need to be subtracted here 
              potentialDegree <- ifelse(nrow(idSub) > 0, length(potentialDegreeIdentifier), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed.
            }            
            if(nrow(idSub) > 0){ #for some reason, if and else statements kept returning an error, so I was forced to use 2 if statements instead
              potentialIndivDurations <- unname(apply(idSub[,grep("dist.to.", colnames(idSub))], 2, function(x){length(which(is.na(x) == FALSE))})) #This means if the individual WAS present the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. 
            }
            if(nrow(idSub) == 0){
              potentialIndivDurations <- rep(0, length(4:ncol(outputMat)))
            }
            
            potentialDurations <- sum(as.integer(potentialIndivDurations)) #if the individual WAS present the maximum potential contact durations is the sum of all individuals observed at each time step during the time period/block
            
            outputMat[1,(2:ncol(outputMat))] <- c(potentialDegree, potentialDurations, potentialIndivDurations)
            
            outputFrame <- data.frame(outputMat, stringsAsFactors = TRUE) #convert to a data frame to allow storage of multiple data types (the ids will be character strings)
            outputFrame$id <- j #define the id
            return(outputFrame)
            
          }
          
          potentialDurationsBlock<-data.frame(data.table::rbindlist(fe1), stringsAsFactors = TRUE) #bind the data together
          
          potentialDurationsBlock$block <- unique(blockSub$block) #define block info
          potentialDurationsBlock$block.start <- unique(blockSub$block.start) #define block info
          potentialDurationsBlock$block.end <- unique(blockSub$block.end) #define block info
      }, envir = environmentTag)
      
      return(potentialDurationsBlock) #note that all the other processes took place in the master-function frame, so we can just return NULL here.
      
    }
    
    potentialDurationsFrame<-data.frame(data.table::rbindlist(foreach::foreach(l = unique(x$block)) %do% blockBreak.function(l)), stringsAsFactors = TRUE)
    potentialDurationsFrame<- potentialDurationsFrame[order(potentialDurationsFrame$block, potentialDurationsFrame$id),] #order the final output by block and id
    
  }else{ #if blocking == FALSE
  
    idSeq<- unique(x$id)
    dist.colnames <-substring(colnames(x)[grep("dist.to.", colnames(x))],9) #Note specific id in "dist.to." colnames start at character 9. This will be different depending on if dist2all or dist2area was used
    dist.colnames_modified <- gsub("indiv_","",dist.colnames) #remove "indiv_" from dist2all output
      
      fe1 <- foreach::foreach(j = idSeq, .noexport = c("idSub", "outputMat")) %do% {
        
        idSub <- droplevels(x[which(x$id == j),]) #subset dist.input to only contain the id-value of interest
        outputMat<- matrix(nrow = 1, ncol = (3 + length(grep("dist.to.", colnames(x))))) #set up the empty matrix to hold the output data.
        colnames(outputMat)<- c("id", "potenDegree", "potenTotalContactDurations", paste("potenContactDurations_", dist.colnames_modified, sep =""))  
        
        #if the individual was not present (i.e., nrow(idSub) == 0), then 0s well be put in the matrix.  However, how many potential nodes present in the block takes a bit more effort to calculate because the distance input object may have come from dist2All_df or dist2Area_df.
        
        potentialDegreeIdentifier <-NULL #create an empty vector to describe when nodes are present in the data set.
        
        for(i in 1:length(grep("dist.to.", colnames(x)))){ #loops through the potentially-contactable nodes.
          
          presenceTest<- is.na(x[,grep("dist.to.", colnames(x))[i]]) #if nodes were not present all entries will be TRUE.
          
          if(length(which(presenceTest == FALSE)) > 0){ #identify when the reported value is FALSE (i.e., when nodes WERE present)
            
            potentialDegreeIdentifier <-c(potentialDegreeIdentifier, 1) #if nodes are present add 1 value to potentialDegreeIdentifier
            
          }else{next}
          
        }
        if(distFunction == "dist2All_df"){ #if the dist2All_df function was used to generate x
          potentialDegree <- ifelse(nrow(idSub) > 0, (length(potentialDegreeIdentifier) - 1), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed minus 1 (because individual j could not be in contact with itself).
        }
        if(distFunction == "dist2Area_df"){ #if the dist2Area_df function was used to generate x, then 1 does not need to be subtracted here 
          potentialDegree <- ifelse(nrow(idSub) > 0, length(potentialDegreeIdentifier), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed.
        }            

        if(nrow(idSub) > 0){ #for some reason, if and else statements kept returning an error, so I was forced to use 2 if statements instead
          potentialIndivDurations <- unname(apply(idSub[,grep("dist.to.", colnames(idSub))], 2, function(x){length(which(is.na(x) == FALSE))})) #This means if the individual WAS present the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. 
        }
        if(nrow(idSub) == 0){
          potentialIndivDurations <- rep(0, length(4:ncol(outputMat))) #just zeroes
        }
        
        potentialDurations <- sum(as.integer(potentialIndivDurations)) #if the individual WAS present the maximum potential contact durations is the sum of all individuals observed at each time step during the time period/block
        
        outputMat[1,(2:ncol(outputMat))] <- c(potentialDegree, potentialDurations, potentialIndivDurations)
        
        outputFrame <- data.frame(outputMat, stringsAsFactors = TRUE) #convert to a data frame to allow storage of multiple data types (the ids will be character strings)
        outputFrame$id <- j #define the id
        return(outputFrame)
        
      }
      
      potentialDurationsFrame<-data.frame(data.table::rbindlist(fe1), stringsAsFactors = TRUE) #bind the data together
      potentialDurationsFrame<- potentialDurationsFrame[order(potentialDurationsFrame$id),] #order the final output by id
      rm(fe1) #remove to free up memory
    
        }
        
        #assign("potentialDurationsFrame", potentialDurations, envir = subEnvir1) #bring the function output back into the apply environment
        
      }, envir = environmentTag)
      return(potentialDurationsFrame)
    }
    
    potentialDurationsList<- foreach::foreach(m = seq(1:length(y))) %do% listBreak.function(m)
    #potentialDurationsList<- apply(listBreakFrame, 1, listBreak.function, environmentTag = thisEnvironment)
    return(potentialDurationsList)
    
  }else{ #if x input is a single data frame
    
    if(blocking == TRUE){
      
      #the code immediately below comes from the contact::timeBlock.append function. We just refrain from calling it here to prevent the inputs being cloned within the function and because we don't need ALL the info that particular function would append to x. Thus, we save time and memory.
      
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
      
      lub.hours = lubridate::hour(x$dateTime)
      lub.dates = lubridate::date(x$dateTime)
      
      daySecondList = lub.hours * 3600 + lubridate::minute(x$dateTime) * 60 + lubridate::second(x$dateTime) #This calculates a day-second
      x<-x[order(lub.dates, daySecondList),] #in case this wasn't already done, we order by date and second. Note that we must order it in this round-about way (using the date and daySecond vectors) to prevent ordering errors that sometimes occurs with dateTime data
      
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
      
      numblocks <- as.integer(ceiling((max(studySecond) - 1)/blockLength1))
      block <-rep(0,length(studySecond))
      for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
        block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
      }
      if(length(which(block == 0)) > 0){ #identifies the last block
        block[which(block == 0)] = numblocks
      }
      
      block.start<-as.character((as.POSIXct(x$dateTime[1]) - blockTimeAdjustment) + ((block - 1)*blockLength1)) #identify the timepoint where each block starts (down to the second resolution)
      block.end<-as.character((as.POSIXct(x$dateTime[1]) - blockTimeAdjustment) + ((block - 1)*blockLength1) + (blockLength1 -1)) #identify the timepoint where each block ends (down to the second resolution)
      
      x$block <- block
      x$block.start <- block.start
      x$block.end <- block.end

      rm(list = c("daySecondList", "block.start","block.end", "lub.dates", "lub.hours", "numblocks", "totSecond", "totSecond2", "studySecond", "block")) #remove these objects because they are no longer needed.
      
      #now that we've assigined the block information, we can begin processing the data in earnest
      
      idSeq<- unique(x$id)
      blockSeq <- unique(x$block)
      dist.colnames <-substring(colnames(x)[grep("dist.to.", colnames(x))],9) #Note specific id in "dist.to." colnames start at character 9. This will be different depending on if dist2all or dist2area was used
      dist.colnames_modified <- gsub("indiv_","",dist.colnames) #remove "indiv_" from dist2all output
      
      blockBreak.function <-function(l, environmentTag = thisEnvironment, columnNames = dist.colnames_modified){
        
        #browser()
        #Because the input data files can be quite large, we do not want to clone it excessively by running subfunctions. So, we relate the following expression to the previously-generated environment (where the data already exist)
        assign("breakBlock", unname(unlist(l)), envir = environmentTag)
        
        eval(expr = { 
          blockSub <- droplevels(x[which(x$block == breakBlock),]) #subset the data set by block
            
            fe1 <- foreach::foreach(j = idSeq, .noexport = c("idSub", "outputMat")) %do% {
              
              idSub <- droplevels(blockSub[which(blockSub$id == j),]) #subset dist.input to only contain the id-value of interest
              outputMat<- matrix(nrow = 1, ncol = (3 + length(grep("dist.to.", colnames(blockSub))))) #set up the empty matrix to hold the output data.
              colnames(outputMat)<- c("id", "potenDegree", "potenTotalContactDurations", paste("potenContactDurations_", columnNames, sep =""))  
              
              #if the individual was not present (i.e., nrow(idSub) == 0), then 0s well be put in the matrix
              potentialDegreeIdentifier <-NULL #create an empty vector to describe when nodes are present in the data set.
              
              for(i in 1:length(grep("dist.to.", colnames(blockSub)))){ #loops through the potentially-contactable nodes.
                
                presenceTest<- is.na(blockSub[,grep("dist.to.", colnames(blockSub))[i]]) #if nodes were not present all entries will be TRUE.
                
                if(length(which(presenceTest == FALSE)) > 0){ #identify when the reported value is FALSE (i.e., when nodes WERE present)
                  
                  potentialDegreeIdentifier <-c(potentialDegreeIdentifier, 1) #if nodes are present add 1 value to potentialDegreeIdentifier
                  
                }else{next}
                
              }
              if(distFunction == "dist2All_df"){ #if the dist2All_df function was used to generate x
                potentialDegree <- ifelse(nrow(idSub) > 0, (length(potentialDegreeIdentifier) - 1), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed minus 1 (because individual j could not be in contact with itself).
              }
              if(distFunction == "dist2Area_df"){ #if the dist2Area_df function was used to generate x, then 1 does not need to be subtracted here 
                potentialDegree <- ifelse(nrow(idSub) > 0, length(potentialDegreeIdentifier), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed.
              }     
              if(nrow(idSub) > 0){ #for some reason, if and else statements kept returning an error, so I was forced to use 2 if statements instead
                potentialIndivDurations <- unname(apply(idSub[,grep("dist.to.", colnames(idSub))], 2, function(x){length(which(is.na(x) == FALSE))})) #This means if the individual WAS present the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. 
              }
              if(nrow(idSub) == 0){
                potentialIndivDurations <- rep(0, length(4:ncol(outputMat)))
              }
              
              potentialDurations <- sum(as.integer(potentialIndivDurations)) #if the individual WAS present the maximum potential contact durations is the sum of all individuals observed at each time step during the time period/block
              
              outputMat[1,(2:ncol(outputMat))] <- c(potentialDegree, potentialDurations, potentialIndivDurations)
              
              outputFrame <- data.frame(outputMat, stringsAsFactors = TRUE) #convert to a data frame to allow storage of multiple data types (the ids will be character strings)
              outputFrame$id <- j #define the id
              return(outputFrame)
              
            }
            
            potentialDurationsBlock<-data.frame(data.table::rbindlist(fe1), stringsAsFactors = TRUE) #bind the data together
            
            potentialDurationsBlock$block <- unique(blockSub$block) #define block info
            potentialDurationsBlock$block.start <- unique(blockSub$block.start) #define block info
            potentialDurationsBlock$block.end <- unique(blockSub$block.end) #define block info
            
        }, envir = environmentTag)
        
        return(potentialDurationsBlock) #note that all the other processes took place in the master-function frame, so we can just return NULL here.
        
      }
      
      potentialDurationsFrame<-data.frame(data.table::rbindlist(foreach::foreach(l = unique(x$block)) %do% blockBreak.function(l)), stringsAsFactors = TRUE)
      potentialDurationsFrame<- potentialDurationsFrame[order(potentialDurationsFrame$block, potentialDurationsFrame$id),] #order the final output by block and id
      
    }
    else{ #if blocking == FALSE
      
      idSeq<- unique(x$id)
      dist.colnames <-substring(colnames(x)[grep("dist.to.", colnames(x))],9) #Note specific id in "dist.to." colnames start at character 9. This will be different depending on if dist2all or dist2area was used
      dist.colnames_modified <- gsub("indiv_","",dist.colnames) #remove "indiv_" from dist2all output
        
        fe1 <- foreach::foreach(j = idSeq) %do% {
          
          idSub <- droplevels(x[which(x$id == j),]) #subset dist.input to only contain the id-value of interest
          outputMat<- matrix(nrow = 1, ncol = (3 + length(grep("dist.to.", colnames(x))))) #set up the empty matrix to hold the output data.
          colnames(outputMat)<- c("id", "potenDegree", "potenTotalContactDurations", paste("potenContactDurations_", dist.colnames_modified, sep =""))  

          #if the individual was not present (i.e., nrow(idSub) == 0), then 0s well be put in the matrix
          potentialDegreeIdentifier <-NULL #create an empty vector to describe when nodes are present in the data set.
          
          for(i in 1:length(grep("dist.to.", colnames(x)))){ #loops through the potentially-contactable nodes.
            
            presenceTest<- is.na(x[,grep("dist.to.", colnames(x))[i]]) #if nodes were not present all entries will be TRUE.
            
            if(length(which(presenceTest == FALSE)) > 0){ #identify when the reported value is FALSE (i.e., when nodes WERE present)
              
              potentialDegreeIdentifier <-c(potentialDegreeIdentifier, 1) #if nodes are present add 1 value to potentialDegreeIdentifier
              
            }else{next}
            
          }
          if(distFunction == "dist2All_df"){ #if the dist2All_df function was used to generate x
            potentialDegree <- ifelse(nrow(idSub) > 0, (length(potentialDegreeIdentifier) - 1), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed minus 1 (because individual j could not be in contact with itself).
          }
          if(distFunction == "dist2Area_df"){ #if the dist2Area_df function was used to generate x, then 1 does not need to be subtracted here 
            potentialDegree <- ifelse(nrow(idSub) > 0, length(potentialDegreeIdentifier), 0) #if the individual WAS present the maximum potential degree is the number of individuals observed over the course of the time period/block when individual j was also observed.
          }   
          
          if(nrow(idSub) > 0){ #for some reason, if and else statements kept returning an error, so I was forced to use 2 if statements instead
            potentialIndivDurations <- unname(apply(idSub[,grep("dist.to.", colnames(idSub))], 2, function(x){length(which(is.na(x) == FALSE))})) #This means if the individual WAS present the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. 
          }
          if(nrow(idSub) == 0){
            potentialIndivDurations <- rep(0, length(4:ncol(outputMat)))
          }
          
          potentialDurations <- sum(as.integer(potentialIndivDurations)) #if the individual WAS present the maximum potential contact durations is the sum of all individuals observed at each time step during the time period/block
          
          outputMat[1,(2:ncol(outputMat))] <- c(potentialDegree, potentialDurations, potentialIndivDurations)
          outputFrame <- data.frame(outputMat, stringsAsFactors = TRUE) #convert to a data frame to allow storage of multiple data types (the ids will be character strings)
          outputFrame$id <- j #define the id
          return(outputFrame)
        }
        
        potentialDurationsFrame<-data.frame(data.table::rbindlist(fe1), stringsAsFactors = TRUE) #bind the data together
        potentialDurationsFrame<- potentialDurationsFrame[order(potentialDurationsFrame$id),] #order the final output by id
        rm(fe1) #remove to free up memory
      
    }
    return(potentialDurationsFrame)
  }
}
