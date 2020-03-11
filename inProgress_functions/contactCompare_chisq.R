



#' @param importBlocks Logical. If true, each block in x.summary will be 
#'    analyzed separately. Defaults to FALSE. Note that the "block" column must
#'    exist in x.summary AND x.potential, otherwise an error will be returned.
#' @param shuffle.type Integer. Describes which shuffle.type (from the 
#'    randomizePaths function) was used to randomize the y.summary data 
#'    set(s). Takes the values "0," "1," or "2." For tests other than 
#'    "chisq" this value is irrelevant.



contactCompare_chisq<-function(x.summary, y.summary, x.potential, y.potential = NULL, importBlocks = FALSE, shuffleType = 1){
  
  
  tryCatch.W.E <- function(expr) #this function comes from https://stat.ethz.ch/pipermail/r-help/2010-December/262626.html
  {
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
  }
  
  summaryAgg.block1<-function(x,y){ #calculates the mean potential contacts by id and block. Using this apply function is faster than simply aggregating the data set by id and block
    sumTable<-y[which(y$block == unname(unlist(x[1]))),]
    
    if(nrow(sumTable) == 0){output <- NULL #if there's nothing in the subset, the function will not proceed any further.
    
    }else{
      
      blockStart<- unique(lubridate::as_datetime(sumTable$block.start)) #added 02/05/2019 - had to keep track of this new information ; updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
      blockEnd<- unique(lubridate::as_datetime(sumTable$block.end)) #added 02/05/2019 - had to keep track of this new information ;  updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
      sumTable.redac<-sumTable[,-c(match("id", names(sumTable)), match("block", names(sumTable)), match("block.start", names(sumTable)), match("block.end", names(sumTable)))]  #Remove the columns that cannot/shoud not be averaged.
      output<-aggregate(sumTable.redac, list(id = sumTable$id), mean) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set.
      output$block = unname(unlist(x[1])) #add this information back into the table
      output$block.start = blockStart #add this information back into the table
      output$block.end = blockEnd #add this information back into the table
      
    }
    return(output)
  }
  summaryAgg.block2<-function(x,y){ #calculates the sum potential contacts by id and replicateID. Using this apply function is faster than simply aggregating the data set by id and replicateID
    sumTable<-y[which(y$replicateID == unname(unlist(x[1]))),]
    
    if(nrow(sumTable) == 0){output <- NULL #if there's nothing in the subset, the function will not proceed any further.
    
    }else{
      
      sumTable.redac<-sumTable[,-c(match("id", names(sumTable)), match("block", names(sumTable)), match("block.start", names(sumTable)), match("block.end", names(sumTable)), match("replicateID", names(sumTable)))]  #Remove the columns that cannot/shoud not be averaged.
      output<-aggregate(sumTable.redac, list(id = sumTable$id), sum) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set.
      
    }
    return(output)
  }
  
  #first we reduce the inputs
  
  ##summary inputs (i.e., outputs from the summarizeContacts function)
  
  if(is.data.frame(x.summary) == FALSE & is.list(x.summary) == TRUE){ #Because R treats dataframes as lists, we assess here if the input was only a single data frame or a list of data frames.
    x.summary <- x.summary[[1]] #if x.summary is a list of data frames, only the first one is used in processing.
  }
  
  colnames(x.summary)<- gsub("avg..","",colnames(x.summary)) #if x.summary represents the average summary report (from the summarizeContacts function), all colnames may be preceded by "avg..". This line removes that tag from all relevant columns.

  if(importBlocks == TRUE){ #if importBlocks == TRUE, but there is no block information in x.summary, the function stops and returns an error.
    
    if(length(x.summary$block) == 0){
      
      stop("importBlocks set to TRUE, but no block column exists in x.summary")
      
    }
    
  }
  
  if(length(x.summary$block) > 0){ #Alternatively, if there are blocked time sets in the empirical set, we need to dictate how to handle them.
    
    if(importBlocks == FALSE){ #if there ARE blocks in x.summary, but importBlocks == FALSE, we have to recreate x.summary without the bloking information.
      
      x.summaryRedac <- droplevels(x.summary[, - c(which(colnames(x.summary) == "block"):ncol(x.summary))]) #remove the block information (four columns: block, block.start, block.end, and numBlocks)
      x.summary<- aggregate(x.summaryRedac[,-match("id", colnames(x.summary))], list(id = x.summaryRedac$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set.
      ###so now we have an x.summary object with no block information
      rm(x.summaryRedac) #remove x.summaryRedac to free up local memory
    }
    
  }
  
  #we do the same thing for y.summary as we did for x.summary above.
  if(is.data.frame(y.summary) == FALSE & is.list(y.summary) == TRUE){ #Because R treats dataframes as lists, we assess here if the input was only a single data frame or a list of data frames.
    y.summary <- y.summary[[1]] #if y.summary is a list of data frames, only the first one is used in processing.
  }
  
  colnames(y.summary)<- gsub("avg..","",colnames(y.summary)) #if y.summary represents the average summary report (from the summarizeContacts function), all colnames may be preceded by "avg..". This line removes that tag from all relevant columns.
  
  if(importBlocks == TRUE){ #if importBlocks == TRUE, but there is no block information in y.summary, the function returns a warning, but precedes as if the y input is relevant to all columns.
    
    if(length(y.summary$block) == 0){
      
      warning("importBlocks set to TRUE, but no block column exists in y.summary. Proceding as if y.summary values are stable and relevant to EVERY block.")
      
      y.summaryBlock<-NULL #create an empty object to contain new block information
      
      for(i in unique(x$block)){
        
        y.summary$block <- i #add a block column containing only the i value to y.summary
        y.summaryBlock <- data.table::rbindlist(list(y.summaryBlock, y.summary), fill = TRUE) #rbind the data frames together
        
      }
      
      y.summary <- data.frame(y.summaryBlock, stringsAsFactors = TRUE) #redefine y.summary as the object containing block information
      rm(y.summaryBlock) #remove y.summaryBlock to free up local memory
    }
    
  }
  
  if(length(y.summary$block) > 0){ #Alternatively, if there are blocked time sets in the empirical set, we need to dictate how to handle them.
    
    if(importBlocks == FALSE){ #if there ARE blocks in y.summary, but importBlocks == FALSE, we have to recreate y.summary without the bloking information.
      
      y.summaryRedac <- droplevels(y.summary[, - c(which(colnames(y.summary) == "block"):ncol(y.summary))]) #remove the block information (four columns: block, block.start, block.end, and numBlocks)
      y.summary<- aggregate(y.summaryRedac[,-match("id", colnames(y.summary))], list(id = y.summaryRedac$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set.
      ###so now we have an y.summary object with no block information
      rm(y.summaryRedac) #remove y.summaryRedac to free up local memory
      
    }
    
  }
  
  ##potential inputs (i.e., outputs from potentialContacts function)
  
  if(is.null(y.potential) == TRUE){ #if there was no y.potential input, the function assumes that x.potential is relevant to the random set as well. This is a completely fair assumption whenever there is no blocking or when the shuffleUnit == 0. In cases when teh shuffle unit is 1 or 2, however, this can lead to erroneous results and/or errors in the function.
    
    y.potential<-x.potential
    
  }
  
  if(is.data.frame(x.potential) == FALSE & is.list(x.potential) == TRUE){ #Because R treats dataframes as lists, we assess here if the input was only a single data frame or a list of data frames.
    
    x.potentialAgg<- data.frame(data.table::rbindlist(x.potential, fill = TRUE), stringsAsFactors = TRUE) #bind the lists together. 
    
    if(importBlocks == TRUE){
      
      if(is.na(match("block", colnames(x.potentialAgg))) == TRUE){ #if there is NOT a "block" column present, the function will return an error. 
        
        stop("importBlocks set to TRUE, but no block column exists in x.potential")
        
      }else{ #If there IS a block column, we average observations by id AND block
        
        blockSeq1<-unique(x.potentialAgg$block) #there may be differences in blocks containing contact info between the empirical and randomized data sets
        meanTab <- apply(data.frame(blockSeq1, stringsAsFactors = TRUE), 1, summaryAgg.block1, y = x.potentialAgg) #Note that block information MUST be included in the x.potentialAgg input
        x.potential <- data.frame(data.table::rbindlist(meanTab, fill = TRUE), stringsAsFactors = TRUE) #We keep the same name for simplicity's sake below. 
        rm(list = c("x.potentialAgg", "meanTab")) #remove x.potentialAgg and meanTab to free up local memory
        
      }
      
    }else{ #if importBlocks == FALSE we only average the data by id.
      
      if(length(x.potentialAgg$block) > 0){ #Alternatively, if there are blocked time sets in the empirical input, we need to dictate how to handle them.
        
        if(importBlocks == FALSE){ #if there ARE blocks in x.potential, but importBlocks == FALSE, we have to recreate x.potential without the bloking information.
          
          for(i in 1:length(x.potential)){ #add replicateID to the input. This was not required when processing x.summary because only the first entry of the list is returned.
            x.potential[[i]]$replicateID == i}
          
          x.potentialAgg<- data.frame(data.table::rbindlist(x.potential, fill = TRUE), stringsAsFactors = TRUE) #remake x.potentialAgg
          
          blockSeq1<-unique(x.potentialAgg$block) #there may be differences in blocks containing contact info between the empirical and randomized data sets
          meanTab <- apply(data.frame(blockSeq1, stringsAsFactors = TRUE), 1, summaryAgg.block1, y = x.potentialAgg) #Note that block information MUST be included in the x.potentialAgg input
          x.potentialAggBlockMean <- data.frame(data.table::rbindlist(meanTab, fill = TRUE), stringsAsFactors = TRUE) #We keep the same name for simplicity's sake below. 
          #now we have an object that represents average contact durations in each block
          repSeq<-unique(x.potentialAggBlockMean$replicateID)
          sumTab <- apply(data.frame(repSeq, stringsAsFactors = TRUE), 1, summaryAgg.block2, y = x.potentialAggBlockMean) #sum across columns by id and replicateID
          x.potentialAgg <- data.frame(data.table::rbindlist(sumTab, fill = TRUE), stringsAsFactors = TRUE) #We keep the same name for simplicity's sake below. 
          ##and now we have an object contained summed averages by id, replicateID
          rm(list = c("meanTab", "sumTab", "x.potentialAggBlockMean")) #remove meanTab to free up local memory
          
        }
        
      }
      ##now we take the average across all replicateIDs
      x.potential<-aggregate(x.potentialAgg[,-match("id", colnames(x.potentialAgg))], list(id = x.potentialAgg$id), mean) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 
      rm(x.potentialAgg) #remove x.potentialAgg to free up local memory
      
    }
    
  }else{ #if x.potential is only a single data frame
    
    if(importBlocks == TRUE){
      
      if(is.na(match("block", colnames(x.potential))) == TRUE){ #if there is NOT a "block" column present, the function will return an error. 
        
        stop("importBlocks set to TRUE, but no block column exists in x.potential")
        
      }
      
    }else{ #if importBlocks == FALSE we must sum the data by id if blocks exist in the input data. If blocks do not exist in the data, nothing needs to be done.
      
      if(length(x.potential$block) > 0){ #Alternatively, if there are blocked time sets in the empirical input, we need to dictate how to handle them.
        
        x.potential<-aggregate(x.potential[,-match("id", colnames(x.potential))], list(id = x.potential$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 
        
      }

    }
    
  }
  
  
  #we do the same thing for y.potential as we did for x.potential above.
  
  
  
  
  
  indivSummaryTest<- ifelse(length(grep("contactDuration_Indiv", colnames(x.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
  
  
  
}



















