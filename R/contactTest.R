#' Determine if Observed Contacts are More or Less Frequent than in a Random Distribution
#'
#' This function is used to determine if tracked individuals in an empirical dataset had more or fewer contacts with other tracked individuals/specified locations than would be expected at random. The function works by comparing an empirically-based summarizeContacts function output (x) to the contactReport.avg generated from randomized data (y). Users can specify if they want to look at individual or overall (population-level) comparisons.
#' 
#' The default tested column (i.e., categorical data column from which data is drawn to be compared to randomized sets herein) is "id." This means that contacts involving each individual (defined by a unique "id") will be compared to randomized sets. Users may not use any data column for analysis other than "id." If users want to use another categorical data column in analyses rather than "id," we recommend re-processing data (starting from the dist.all/distToArea functions), while specifying this new data as an "id." For example, users may annotate an illness status column to the empirical input, wherein they describe if the tracked individual displayed gastrointestinal ("gastr"), respiratory ("respr"), both ("both"), illness symptoms, or were consistently healthy ("hel") over the course of the tracking period. Users could set this information as the "id," and carry it forward as such through the data-processing pipeline. Ultimately, they could determine if each of these disease states affected contact rates, relative to what would be expected at random.
#'
#' Note: Currently, this function requires defined TSWs (temporal sampling windows; see the tempAggregate function) (i.e., contact durations must represent time periods that are equidistant).
#'
#' Note: The current functionality is limited to comparisons using the X-squared test of independence. This works by comparing the number of TSWs in a given dataset/block between individuals/fixed locations and the number of TSWs during which individuals were observed, but were not in contact with specific individuals/places, for empirical (x) and randomized (y) datasets. The dist (z ; from dist.all or distToArea output) input is used here to determine how frequently each individual was observed in the empirical dataset/block of interest, allowing us to calculate the number of TSWs each individual was present but not involved in contacts. Note here that if X-squared expected values will be very small, approximations of p may not be right (and in fact, all estimates will be poor). It may be best to weight these tests differently. To address this, I've added the "warning" column to the output which notifies users when the chi-sq function reported that results may be inaccurate.
#'
#' Note: if blocking == TRUE, blockLength, and blockUnit should be consistent with previous function inputs. 
#' @param emp.input List or data frame containing summarizeContact-function output refering to the empirical data.
#' @param rand.input List or data frame containing averaged summarizeContact-function output (i.e., argument "avg" was set to "TRUE" when the summarizeContacts function was run) refering to the randomized-path data.
#' @param dist.input List or data frame containing dist.all-/distToArea-function output refering to the empirical data.
#' @param test Describes the statistical test used to evaluate differences. Currently only takes the value "chisq," more tests will be added in later versions.
#' @param indivComparison Logical. If TRUE output will return a data frame showing individual-level comparsions. Defaults to TRUE.
#' @param overallComparison Logical. If TRUE output will return a data frame showing population-level comparisons. Defaults to TRUE.
#' @param blocking Description imminent
#' @param blockUnit Description imminent
#' @param blockLength Description imminent
#' @param shuffle.type Description imminent
#' @param parallel Description imminent
#' @keywords data-processing smoothing location point
#' @export
#' @examples
#' Examples imminent

contactTest<-function(emp.input, rand.input, dist.input, test = "chisq", indivComparison = TRUE, overallComparison = TRUE, blocking = TRUE, blockUnit = "mins", blockLength = 10, shuffle.type = 0, parallel = TRUE){
  
  timeBlock.append<-function(x, dateTime = NULL, blockLength, blockUnit){
    datetime.append1 = function(x){
      
      timevec = x$dateTime
      daySecondList = lubridate::hour(timevec) * 3600 + lubridate::minute(timevec) * 60 + lubridate::second(timevec) #This calculates a day-second
      lub.dates = lubridate::date(x$dateTime)
      dateseq = unique(lub.dates)
      dayIDList = NULL
      dayIDseq = seq(1,(length(dateseq)),1)
      for(b in dayIDseq){
        ID = rep(b,length(which(lub.dates == dateseq[which(dayIDseq == b)])))
        dayIDList = c(dayIDList, ID)
      } 
      x$totalSecond = ((dayIDList - min(dayIDList))*86400) + daySecondList #This calculates the total second (the cumulative second across the span of the study's timeframe)
      
      return(x)
    }
    
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
    
    x<-x[order(x$dateTime),] #Just in case the data wasn't already ordered in this way.
    x<-datetime.append1(x) #adds the total second column to the dataframe
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
    x$block = block
    return(x)
  }
 # chiFunc.noBlock <- function(x, empirical, randomized, dist){
#    
#    empDurations<- empirical[which(empirical$id == unlist(unname(x[1]))), unlist(unname(x[2]))]
#    
#    if(is.na(empDurations) == FALSE){ #if is.na == TRUE, nothing will happen
#      
#      if(unlist(unname(x[2])) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in dist
#        maxDurations = length(which(dist$id == unlist(unname(x[1])) & is.na(dist[,(unlist(unname(x[2])) + 1)]) == FALSE)) #This means the the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. For example, if one of the individuals disappears because we have no recorded time points for them at a given time, we can't say that they have the potential to be in contact with others during this period.
#      }else{ #i.e., x[2] == 2 or 3
#        maxDurations = length(which(dist$id == unlist(unname(x[1]))))
#      }
#      
#      randDurations<-randomized[which(randomized$id == unlist(unname(x[1]))), unlist(unname(x[2]))]
#      
#      compareFrame<- data.frame(matrix(ncol = 2, nrow = 2*maxDurations))
#      colnames(compareFrame)<-c("contact", "data")
#      compareFrame[,1] <- c(rep(1,empDurations), rep(0,maxDurations-empDurations),rep(1,randDurations), rep(0,maxDurations-randDurations))
#      compareFrame[,2] <- c(rep("empirical",maxDurations), rep("randomized",maxDurations))
#      
#      compareTable <- table(compareFrame$data,compareFrame$contact)
#      test<- chisq.test(compareTable)
#      
#      summaryFrame <- data.frame(id1 = unlist(unname(x[1])), id2 = NA, method = test[4], X.squared = test[1], df = test[2], p.val = test[3], empiricalContactDurations = empDurations, randContactDurations.mean = randDurations, empiricalNoContactDurations = (maxDurations - empDurations), randNoContactDurations.mean = (maxDurations - randDurations), difference = abs((maxDurations - empDurations) - (maxDurations - randDurations)))
#      
#      summaryFrame[1,2] <- ifelse(grep("contactDuration_",names(empirical)[unlist(unname(x[2]))]) != 1, names(empirical)[unlist(unname(x[2]))], substring(names(empirical)[unlist(unname(x[2]))],17)) #This identifies whom the individual noted in summaryFrame[1,1] is in contact with. Remember that the input for the contactTest function is the output of the summarizeContact function. Columns 2-3 of this input represent totalDegree and totalContactDurations. This line of code puts those names in the summaryFrame$id2 column. Otherwise the entry will be taken from the column names such as "contactDuration_Indiv.1", beginning with the 17th character, for example "Indiv.1".
#      
#      return(summaryFrame)
#    }
#  }
#  chiFunc.Block <- function(x, empirical, randomized, dist){
#    
#    empDurations<- empirical[which(empirical$id == unlist(unname(x[1])) & empirical$block == unlist(unname(x[3]))), unlist(unname(x[2]))] #pull the value of a given column for a specific id in a specific block
#    summaryFrame<-NULL
#    
#    if(length(empDurations) > 0 && is.na(empDurations) == FALSE){ #if there is no entry OR is.na == TRUE, nothing will happen
#      block <- unlist(unname(x[3]))
#      warn = ""
#      if(unlist(unname(x[2])) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in the dist.input.
#        maxDurations = length(which(dist$id == unlist(unname(x[1])) & is.na(dist[,(unlist(unname(x[2])) + 1)]) == FALSE & dist$block == unlist(unname(x[3])))) #This means the the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. For example, if one of the individuals disappears because we have no recorded time points for them at a given time, we can't say that they have the potential to be in contact with others during this period.
#      }else{ #i.e., x[2] == 2 or 3
#        maxDurations = length(which(dist$id == unlist(unname(x[1])) & dist$block == unlist(unname(x[3]))))
#      }
#      
#      if(unlist(unname(x[4])) == 2){ #recall that shuffle.type 2 (from the randomizeLocations function) produces only 1 shuffle.unit's worth of data, rather than a dataset with the same length of x. As such, there may be a different number of blocks in y compared to x. Here we assume that the mean randomized durations per block, are representative of mean randomized durations per block across each shuffle unit (e.g., day)
#        if(unlist(unname(x[3])) > max(randomized$blocks)){
#          block = unlist(unname(x[3])) - (ceiling((unlist(unname(x[3])) - max(randomized$blocks)) /max(randomized$blocks))*max(randomized$blocks))
#          randDurations<-randomized[which(randomized$id == unlist(unname(x[1])) & randomized$block == block), unlist(unname(x[2]))]
#        }
#        
#      }else{
#        randDurations<-randomized[which(randomized$id == unlist(unname(x[1])) & randomized$block == unlist(unname(x[3]))), unlist(unname(x[2]))]
#      }
#      
#      compareFrame<- data.frame(matrix(ncol = 2, nrow = 2*maxDurations))
#      colnames(compareFrame)<-c("contact", "data")
#      compareFrame[,1] <- c(rep(1,empDurations), rep(0,maxDurations-empDurations),rep(1,randDurations), rep(0,maxDurations-randDurations))
#      compareFrame[,2] <- c(rep("empirical",maxDurations), rep("randomized",maxDurations)) 
#      
#      compareTable <- table(compareFrame$data,compareFrame$contact) #generates a table showing the total number of contacts (1) and non-contacts (0) for the empirical and randomized sets
#      warning() #clears the warnings
#      test<- chisq.test(compareTable)
#      warn<- names(warnings())[1] #used to identify when the chi-sq.test produces the "Chi-squared approximation may be incorrect" warning, indicating that that expected X-squared values are very small, and generated estimates will be poor.
#      warning() #clears the warnings
#      summaryFrame <- data.frame(id1 = unlist(unname(x[1])), id2 = NA, method = test[4], X.squared = test[1], df = test[2], p.val = test[3], empiricalContactDurations = empDurations, randContactDurations.mean = randDurations, empiricalNoContactDurations = (maxDurations - empDurations), randNoContactDurations.mean = (maxDurations - randDurations), difference = abs((maxDurations - empDurations) - (maxDurations - randDurations)), empBlock = unlist(unname(x[3])), randBlock = block, warning = warn)
#      rownames(summaryFrame) <-1
#      summaryFrame[1,2] <- ifelse(length(grep("contactDuration_",names(empirical)[unlist(unname(x[2]))])) < 1, names(empirical)[unlist(unname(x[2]))], substring(names(empirical)[unlist(unname(x[2]))],17)) #This identifies whom the individual noted in summaryFrame[1,1] is in contact with. Remember that the input for the contactTest function is the output of the summarizeContact function. Columns 2-3 of this input represent totalDegree and totalContactDurations. This line of code puts those names in the summaryFrame$id2 column. Otherwise the entry will be taken from the column names such as "contactDuration_Indiv.1", beginning with the 17th character, for example "Indiv.1".
#    }
#    
#    return(summaryFrame)
#  }
  
  chisq.forLoop<-function(x, empirical = x, randomized = y, dist = dist.input, blocking){ #I hate that I have to do this in a for-loop, but I couldn't get the apply functions (above) to work.
    output<-NULL
    
    #vectorize the input
    ids<-x[,1]
    cols<-x[,2]
    shuffle.type<-unique(x[,ncol(x)])
    
    if(blocking == TRUE){
      blocks<-x[,3]
      for(i in 1:nrow(x)){
        
        empDurations<- empirical[which(empirical$id == unlist(unname(ids[i])) & empirical$block == unlist(unname(blocks[i]))), unlist(unname(cols[i]))] #pull the value of a given column for a specific id in a specific block
        empDurations<-empDurations[is.na(empDurations) == FALSE] #remove the NAs
        summaryFrame<-NULL
        
        if(length(empDurations) > 0){ #if there is no entry OR is.na == TRUE, nothing will happen
          
          block <- unlist(unname(blocks[i]))
          warn = ""
          if(unlist(unname(cols[i])) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in the dist.input.
            maxDurations = length(which(dist$id == unlist(unname(ids[i])) & is.na(dist[,(unlist(unname(cols[i])) + 1)]) == FALSE & dist$block == unlist(unname(blocks[i])))) #This means the the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. For example, if one of the individuals disappears because we have no recorded time points for them at a given time, we can't say that they have the potential to be in contact with others during this period. Also, note that the "+1" here is because column 4 in empirical refers to contacts with indiv1, but in the dist data, column 5 (i.e., 4 + 1) refers to distance to indiv1, while column 4 is individuals' "id."
          }else{ #i.e., x[2] == 2 or 3
            maxDurations = length(which(dist$id == unlist(unname(ids[i])) & dist$block == unlist(unname(blocks[i]))))
          }
          
          if(maxDurations == 0){
            next
          }
          
          if(shuffle.type == 2){ #recall that shuffle.type 2 (from the randomizeLocations function) produces only 1 shuffle.unit's worth of data, rather than a dataset with the same length of x. As such, there may be a different number of blocks in y compared to x. Here we assume that the mean randomized durations per block, are representative of mean randomized durations per block across each shuffle unit (e.g., day)
            if(unlist(unname(blocks[i])) > max(randomized$blocks)){
              block = unlist(unname(blocks[i])) - (ceiling((unlist(unname(blocks[i])) - max(randomized$blocks)) /max(randomized$blocks))*max(randomized$blocks))
              randDurations<-randomized[which(randomized$id == unlist(unname(ids[i])) & randomized$block == block), unlist(unname(cols[i]))]
            }
            
          }else{
            randDurations<-randomized[which(randomized$id == unlist(unname(ids[i])) & randomized$block == unlist(unname(blocks[i]))), unlist(unname(cols[i]))]
          }
          randDurations<-randDurations[is.na(randDurations) == FALSE] #remove the NAs
          compareFrame<- data.frame(matrix(ncol = 2, nrow = 2*maxDurations))
          colnames(compareFrame)<-c("contact", "data")
          compareFrame[,1] <- c(rep(1,sum(empDurations)), rep(0,maxDurations-sum(empDurations)),rep(1,sum(randDurations)), rep(0,maxDurations-sum(randDurations)))
          compareFrame[,2] <- c(rep("empirical",maxDurations), rep("randomized",maxDurations)) 
          
          compareTable <- table(compareFrame$data,compareFrame$contact) #generates a table showing the total number of contacts (1) and non-contacts (0) for the empirical and randomized sets
          warning() #clears the warnings
          test<- chisq.test(compareTable)
          warn<- names(warnings())[1] #used to identify when the chi-sq.test produces the "Chi-squared approximation may be incorrect" warning, indicating that that expected X-squared values are very small, and generated estimates will be poor.
          warning() #clears the warnings
          summaryFrame <- data.frame(id1 = unlist(unname(ids[i])), id2 = NA, method = test[4], X.squared = test[1], df = test[2], p.val = test[3], empiricalContactDurations = sum(empDurations), randContactDurations.mean = sum(randDurations), empiricalNoContactDurations = (maxDurations - sum(empDurations)), randNoContactDurations.mean = (maxDurations - sum(randDurations)), difference = abs((maxDurations - sum(empDurations)) - (maxDurations - sum(randDurations))), empBlock = unlist(unname(blocks[i])), randBlock = block, warning = warn)
          rownames(summaryFrame) <-1
          summaryFrame[1,2] <- ifelse(length(grep("contactDuration_",names(empirical)[unlist(unname(cols[i]))])) < 1, names(empirical)[unlist(unname(cols[i]))], substring(names(empirical)[unlist(unname(cols[i]))],17)) #This identifies whom the individual noted in summaryFrame[1,1] is in contact with. Remember that the input for the contactTest function is the output of the summarizeContact function. Columns 2-3 of this input represent totalDegree and totalContactDurations. This line of code puts those names in the summaryFrame$id2 column. Otherwise the entry will be taken from the column names such as "contactDuration_Indiv.1", beginning with the 17th character, for example "Indiv.1".
          
          bindlist<-list(output, summaryFrame)
          output<- data.table::rbindlist(bindlist)
        }else{ #if there is no entry OR is.na == TRUE, nothing will happen
          next
        }
      }
    }else{ #if blocking == False
      for(i in 1:nrow(x)){
        
        empDurations<- empirical[which(empirical$id == unlist(unname(ids[i]))), unlist(unname(cols[i]))] #pull the value of a given column for a specific id
        empDurations<-empDurations[is.na(empDurations) == FALSE] #remove the NAs
        summaryFrame<-NULL
        
        if(length(empDurations) > 0 && is.na(empDurations) == FALSE){ #if there is no entry OR is.na == TRUE, nothing will happen
          warn = ""
          if(unlist(unname(cols[i])) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in the dist.input.
            maxDurations = length(which(dist$id == unlist(unname(ids[i])) & is.na(dist[,(unlist(unname(cols[i])) + 1)]) == FALSE)) #This means the the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. For example, if one of the individuals disappears because we have no recorded time points for them at a given time, we can't say that they have the potential to be in contact with others during this period. Also, note that the "+1" here is because column 4 in empirical refers to contacts with indiv1, but in the dist data, column 5 (i.e., 4 + 1) refers to distance to indiv1, while column 4 is individuals' "id."
          }else{ #i.e., x[2] == 2 or 3
            maxDurations = length(which(dist$id == unlist(unname(ids[i]))))
          }
          if(maxDurations == 0){
            next
          }
          
          randDurations<-randomized[which(randomized$id == unlist(unname(ids[i]))), unlist(unname(cols[i]))]
          randDurations<-randDurations[is.na(randDurations) == FALSE] #remove the NAs
          compareFrame<- data.frame(matrix(ncol = 2, nrow = 2*maxDurations))
          colnames(compareFrame)<-c("contact", "data")
          compareFrame[,1] <- c(rep(1,sum(empDurations)), rep(0,maxDurations-sum(empDurations)),rep(1,sum(randDurations)), rep(0,maxDurations-sum(randDurations)))
          compareFrame[,2] <- c(rep("empirical",maxDurations), rep("randomized",maxDurations)) 
          
          compareTable <- table(compareFrame$data,compareFrame$contact) #generates a table showing the total number of contacts (1) and non-contacts (0) for the empirical and randomized sets
          warning() #clears the warnings
          test<- chisq.test(compareTable)
          warn<- names(warnings())[1] #used to identify when the chi-sq.test produces the "Chi-squared approximation may be incorrect" warning, indicating that that expected X-squared values are very small, and generated estimates will be poor.
          warning() #clears the warnings
          summaryFrame <- data.frame(id1 = unlist(unname(ids[i])), id2 = NA, method = test[4], X.squared = test[1], df = test[2], p.val = test[3], empiricalContactDurations = sum(empDurations), randContactDurations.mean = sum(randDurations), empiricalNoContactDurations = (maxDurations - sum(empDurations)), randNoContactDurations.mean = (maxDurations - sum(randDurations)), difference = abs((maxDurations - sum(empDurations)) - (maxDurations - sum(randDurations))), warning = warn)
          rownames(summaryFrame) <-1
          summaryFrame[1,2] <- ifelse(length(grep("contactDuration_",names(empirical)[unlist(unname(cols[i]))])) < 1, names(empirical)[unlist(unname(cols[i]))], substring(names(empirical)[unlist(unname(cols[i]))],17)) #This identifies whom the individual noted in summaryFrame[1,1] is in contact with. Remember that the input for the contactTest function is the output of the summarizeContact function. Columns 2-3 of this input represent totalDegree and totalContactDurations. This line of code puts those names in the summaryFrame$id2 column. Otherwise the entry will be taken from the column names such as "contactDuration_Indiv.1", beginning with the 17th character, for example "Indiv.1".
          
          bindlist<-list(output, summaryFrame)
          output<- data.table::rbindlist(bindlist)
        }else{ #if there is no entry OR is.na == TRUE, nothing will happen
          next
        }
      }         
    }
    final.out<-data.frame(data.frame(output))
    return(final.out)
  }
  
  x <- emp.input
  y <- rand.input
  
  x$id<-as.character(x$id)
  y$id<-as.character(y$id)
  dist.input$id<-as.character(dist.input$id)
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #if the x input is a list (not a data frame), the function assumes only the first list entry is relevant to our purposes.
    x = data.frame(x[1])
  }
  
  if(is.data.frame(y) == FALSE & is.list(y) == TRUE){ #if the y input is a list (not a data frame), the function assumes only the first list entry is relevant to our purposes.
    y = data.frame(y[1]) #in the case of y (a.k.a. the randomized data set), the first list entry is likely to be the average summary report (i.e., users set avg = TRUE in the summarizeContacts function)
  }
  
  if(is.data.frame(dist.input) == FALSE & is.list(dist.input) == TRUE){ #if the dist.input input is a list (not a data frame), the function assumes only the first list entry is relevant to our purposes.
    dist.input = data.frame(dist.input[1]) 
  }
  
  #if y represents the average summary report, all colnames will be preceded by "avg..". This code block removes that tag from all relevant columns.
  y.colnameSubstring<- substring(names(y),1,5) 
  alterNames.y<-which(y.colnameSubstring == "avg..")
  if(length(alterNames.y) > 0){
    names(y)<- substring(names(y)[alterNames.y],6,10000)
  }
  
  if(test == "chisq" || test == "CHISQ" || test == "Chisq" || test == "chisquare" || test == "CHISQUARE" || test == "Chisquare" || test == "chi-square" || test == "CHI-SQUARE" || test == "Chi-square"){
    
    idSeq <- unique(x$id)
    ColSeq <- 2:max(grep("contactDuration_", names(x)))
    
    if(blocking == TRUE){
      
      #Because we need to know the maximum number of potential blocks, we pull this information from the empirical numBlocks column
      blockSeq <- seq(1,as.numeric(as.character(unique(x$numBlocks))),1)
      
      #We need to use the dist.input file to determine the number of timesteps individuals were actually observed in each block (if individuals were not present, then by definition no contact could have occurred).
      #So, we need to ensure that timeblocks are appended to the dist.input file.
      dist.input <- timeBlock.append(dist.input, dateTime = NULL, blockLength, blockUnit) #dateTime = NULL assumes a dateTime column exists in dist.input
      
    }
      
      if(indivComparison == TRUE){
        if(blocking == TRUE){
          testFrame1 <- expand.grid(idSeq, ColSeq, blockSeq)
        }else{
          testFrame1 <- expand.grid(idSeq, ColSeq)
        }
        testFrame1$shuffle.type <-shuffle.type
        
        ###for some reason the apply functions below would only ever return NULL values (despite haveing removed all factor information from input data), so I had to switch to a for-loop. If we can get teh applies to work, this would be much more efficient
        
       # if (parallel == TRUE){
      #    cl<-parallel::makeCluster(parallel::detectCores())
      #    chiSqTests.indiv<-parallel::parApply(cl,testFrame, 1, chiFunc.Block, empirical = x, randomized = y, dist = dist.input)
      #    parallel::stopCluster(cl)
      #  }else{ #if parallel == FALSE
      #    chiSqTests.indiv<- apply(testFrame, 1, chiFunc.Block, empirical = x, randomized = y, dist = dist.input)
      #  }
        
      #  if(is.data.frame(chiSqTests.indiv) == FALSE){
      #    testResultsFrame.indiv <- data.frame(data.table::rbindlist(chiSqTests.indiv))
      #  }else{ #if there's only one comparison, there may not be multiple entries in chiSqTests. In this case, chiSqTests will be a single dataFrame.
      #    testResultsFrame.indiv <- chiSqTests.indiv	
      #  }
      #
      
        testResultsFrame.indiv<- chisq.forLoop(testFrame1, empirical = x, randomized = y, dist = dist.input, blocking)
      }
        
      if(overallComparison == TRUE){
        x$id = "."
        y$id = "."
        dist.input$id = "."
        idSeq = "."
        if(blocking == TRUE){
          testFrame2 <- expand.grid(idSeq, ColSeq, blockSeq)
        }else{
          testFrame2 <- expand.grid(idSeq, ColSeq)
        }
        testFrame2$shuffle.type <-shuffle.type
        
        #if (parallel == TRUE){
        #  cl<-parallel::makeCluster(parallel::detectCores())
        #  #chiSqTests<-parallel::parApply(cl,testFrame, 1, chiFunc.Block, x.Agg, y.Agg, durations.max)
        #  chiSqTests.overall<-parallel::parApply(cl,testFrame, 1, chiFunc.Block, x, y, dist.input)
        #  parallel::stopCluster(cl)
        #}else{ #if parallel == FALSE
        #  #chiSqTests<- apply(testFrame, 1, chiFunc.Block, x.Agg, y.Agg, durations.max)
        #  chiSqTests.overall<- apply(testFrame, 1, chiFunc.Block, x, y, dist.input)
        #}
        #
        #if(is.list(chiSqTests.overall) == TRUE){
        #  testResultsFrame.overall <- data.frame(data.table::rbindlist(chiSqTests.overall))
        #}else{ #if there's only one comparison, there may not be multiple entries in chiSqTests. In this case, chiSqTests will be a single dataFrame.
        #  testResultsFrame.overall <- chiSqTests.overall	
        #}
        
        testResultsFrame.overall<- chisq.forLoop(testFrame2, empirical = x, randomized = y, dist = dist.input, blocking)
      }
      
      
  #  }else{ #if blocking != TRUE
  #    if(indivComparison == TRUE){
  #      testFrame <- expand.grid(idSeq, ColSeq)
  #      if (parallel == TRUE){
  #        cl<-parallel::makeCluster(parallel::detectCores())
  #        #chiSqTests<-parallel::parApply(cl,testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, durations.max)
  #        chiSqTests.indiv<-parallel::parApply(cl,testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, dist.all)
  #        parallel::stopCluster(cl)
  #      }else{ #if parallel == FALSE
  #        #chiSqTests<- apply(testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, durations.max)
  #        chiSqTests.indiv<- apply(testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, dist.all)
  #      }
  #      if(is.list(chiSqTests.indiv) == TRUE){
  #        testResultsFrame.indiv <- data.frame(data.table::rbindlist(chiSqTests.indiv))
  #      }else{ #if there's only one comparison, there may not be multiple entries in chiSqTests. In this case, chiSqTests will be a single dataFrame.
   #       testResultsFrame.indiv <- chiSqTests.indiv	
  #      }
  #    }
  #    if(overallComparison == TRUE){
  #      x.Agg$id = "."
  #      idSeq = "."
  #      testFrame <- expand.grid(idSeq, ColSeq)
   #     if (parallel == TRUE){
  #        cl<-parallel::makeCluster(parallel::detectCores())
  #        #chiSqTests<-parallel::parApply(cl,testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, durations.max)
  #        chiSqTests.overall<-parallel::parApply(cl,testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, dist.all)
  #        parallel::stopCluster(cl)
  #      }else{ #if parallel == FALSE
   #       #chiSqTests<- apply(testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, durations.max)
  #        chiSqTests.overall<- apply(testFrame, 1, chiFunc.noBlock, x.Agg, y.Agg, dist.all)
  #      }
  #      if(is.list(chiSqTests.overall) == TRUE){
  #        testResultsFrame.overall <- data.frame(data.table::rbindlist(chiSqTests.overall))
  #      }else{ #if there's only one comparison, there may not be multiple entries in chiSqTests. In this case, chiSqTests will be a single dataFrame.
  #        testResultsFrame.overall <- chiSqTests.overall	
   #     }
  #    }
 #   }
    
  }
  
  if(indivComparison == TRUE & overallComparison == TRUE){
    finalReturn = list(testResultsFrame.overall, testResultsFrame.indiv)
    names(finalReturn)<- c("overallResults", "individualResults")
  }
  if(indivComparison == TRUE & overallComparison == FALSE){
    finalReturn = list(testResultsFrame.indiv)
    names(finalReturn) <- "individualResults"
  }
  if(indivComparison == FALSE & overallComparison == TRUE){
    finalReturn = list(testResultsFrame.overall)
    names(finalReturn) <- "overallResults"
  }
  
  return(finalReturn)
  
}
