#' Summarize Contact Events
#'
#' This function takes the output from contactDur.all or contactDur.area and 
#'    reports the number of durations when tracked individuals are in "contact"
#'    with one another (contactDur.all) or with specified fixed points/polygons
#'    (contactDur.area).
#' 
#' If x is a list, and avg == TRUE, this function will produce an extra data 
#'    frame containing the mean column values for each id (per block if 
#'    importBlocks == TRUE).
#' 
#' This is a sub-function found within the contactTest and ntrkEdges function.
#' @param x Output from the contactDur.all or contactDur.area functions. Can 
#'    be either a data frame or list of data frames.
#' @param importBlocks Logical. If true, each block in x will be analyzed 
#'    separately. Defaults to FALSE. Note that the "block" column must exist 
#'    in x.
#' @param avg Logical. If TRUE, summary output from all data frames contained 
#'    in x will be averaged together. Output will produce an extra data frame 
#'    containing the mean column values for each id (per block if 
#'    importBlocks == TRUE). Defaults to FALSE.
#' @keywords data-processing contact sub-function
#' @return Returns a data frame (or list of data frames if \code{x} is a 
#'    list of data frames) with the following columns:
#'    
#'    \item{id}{The unique ID of a tracked individual for which we will 
#'    summarize to all other individuals/fixed locations observed in \code{x}.}
#'    \item{id}{Sum number of individuals/fixed locations observed in contact 
#'    specific individuals.}
#'    \item{id}{Sum number of contacts associated with specific individuals.}
#'    \item{contactDuration_...}{Number of contacts between specific dyads.} 
#'    
#'    If importBlocks == TRUE, the following columns are appended to the output
#'    data frame described above:
#'    
#'    \item{block}{Integer ID describing unique blocks of time during which 
#'    contacts occur.}
#'    \item{block.start}{The timepoint in \code{x} at which the \code{block}
#'    begins.}
#'    \item{block.end}{The timepoint in \code{x} at which the \code{block}
#'    ends.}
#'    \item{numBlocks}{Integer describing the total number of time blocks 
#'    observed within \code{x} at which the \code{block}}   
#' @export
#' @examples
#' data(calves)
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'    time = calves$time) #create a dataframe with dateTime identifiers for location fixes
#'
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth to 5-min fix intervals.
#'
#' calves.dist<-dist2All_df(x = calves.agg, parallel = FALSE, 
#'    dataType = "Point", lonlat = FALSE) 
#' calves.contact.block<-contactDur.all(x = calves.dist, dist.threshold=1, 
#'    sec.threshold=10, blocking = TRUE, blockUnit = "hours", blockLength = 1, 
#'    equidistant.time = FALSE, parallel = FALSE, reportParameters = TRUE) 
#'    
#' calves.contactSumm.NOblock <- summarizeContacts(calves.contact.block)
#' head(calves.contactSumm.NOblock)
#' 
#' calves.contactSumm.block <- summarizeContacts(calves.contact.block, 
#'    importBlocks = TRUE)
#' head(calves.contactSumm.block)

summarizeContacts<- function(x, importBlocks = FALSE, avg = FALSE){
  
  distributeContacts1<- function(x,y, me){
    if(unname(unlist(x[1])) == me){
      spec.durations = 0
    }else{
      contact1 <- y[c(which(as.character(y$dyadMember1) == unname(unlist(x[1])))),]
      contact2 <- y[c(which(as.character(y$dyadMember2) == unname(unlist(x[1])))),]
      if((nrow(contact1) >= 1) & (nrow(contact2) >= 1)){
        contact.full <- data.frame(data.table::rbindlist(list(contact1,contact2)))
      }
      if((nrow(contact1) >= 1) & (nrow(contact2) == 0)){
        contact.full <- contact1
      }
      if((nrow(contact2) >= 1) & (nrow(contact1) == 0)){
        contact.full <- contact2
      }
      if((nrow(contact2) == 0) & (nrow(contact1) == 0)){
        contact.full <- contact1 #if neither contact1 or contact2 have any rows, contact.full won't have any rows either.
      }
      spec.durations <- ifelse(nrow(contact.full) >= 1, sum(contact.full$contactDuration),0)
    }
    return(spec.durations)
  }
  distributeContacts2<- function(x,y){
    contact.full <- y[c(which(y$area.id == unname(unlist(x[1])))),]
    spec.durations <- ifelse(nrow(contact.full) >= 1, sum(contact.full$contactDuration),0)
    return(spec.durations)
  }
  contSum <-function(x,y, indivSeq, areaSeq){
    me = (unname(unlist(x[1])))
    if(length(y$dyadMember1) > 0){ #This essentially determines if the input was created with dist.all or distToArea. If length(dyadMember1) >0, it was created with dist.all
      indivContact1 <- y[c(which(as.character(y$dyadMember1) == me)),] #had to make this as.character b/c "me" is a factor with levels that may be different than y$dyadMember1
      indivContact2 <- y[c(which(as.character(y$dyadMember2) == me)),] #had to make this as.character b/c "me" is a factor with levels that may be different than y$dyadMember2
    }else{
      indivContact1 <- y[c(which(as.character(y$indiv.id) == me)),] 
      indivContact2 <- matrix(nrow=0,ncol=0)
    }
    
    #Here identify the number of contact durations individuals had with others. How we do this is determined by the number of times individuals appear in y's dyadMember1 and dyadMember2 columns. Note that the only option if the function input originated from contactDur.area is (nrow(indivContact1) >= 1) & (nrow(indivContact2) == 0)
    if((nrow(indivContact1) >= 1) & (nrow(indivContact2) >= 1)){
      indivContact.full <- data.frame(data.table::rbindlist(list(indivContact1,indivContact2)))
      specIndivSeq = unique(c(as.character(indivContact.full$dyadMember1),as.character(indivContact.full$dyadMember2))) #had to add as.character call because contactDur functions now produce factor data. 02/05/2019
      specIndivSeq1 = specIndivSeq[-which(specIndivSeq == me)]
    }
    if((nrow(indivContact1) >= 1) & (nrow(indivContact2) == 0)){
      indivContact.full <- indivContact1
      if(length(y$dyadMember1) > 0){ #This essentially determines if the input was created with dist.all or distToArea. If length(dyadMember1) >0, it was created with dist.all
        specIndivSeq = unique(c(as.character(indivContact.full$dyadMember1),as.character(indivContact.full$dyadMember2))) #had to add as.character call because contactDur functions now produce factor data. 02/05/2019
        specIndivSeq1 = specIndivSeq[-which(specIndivSeq == me)]
      }else{
        specIndivSeq1 = unique(as.character(indivContact.full$area.id)) #had to add as.character call because contactDur functions now produce factor data. 02/05/2019
      }
    }
    if((nrow(indivContact2) >= 1) & (nrow(indivContact1) == 0)){
      indivContact.full <- indivContact2
      specIndivSeq = unique(c(as.character(indivContact.full$dyadMember1),as.character(indivContact.full$dyadMember2))) #had to add as.character call because contactDur functions now produce factor data. 02/05/2019
      specIndivSeq1 = specIndivSeq[-which(specIndivSeq == me)]
    }
    if((nrow(indivContact2) == 0) & (nrow(indivContact1) == 0)){
      indivContact.full <- indivContact1 #if neither indivContact1 or indivContact2 have any rows, indivContact.full won't have any rows either.
      specIndivSeq1 = 0
    }
    if(length(y$dyadMember1) > 0){ #This essentially determines if the input was created with dist.all or distToArea. If length(dyadMember1) >0, it was created with dist.all
      if(nrow(indivContact.full) > 1){
        indivSeqFrame1 <-data.frame(indivSeq)
        contactSum<-apply(indivSeqFrame1, 1, distributeContacts1, indivContact.full, me)
        sumTable <- data.frame(matrix(ncol = (3+length(indivSeq)), nrow = 1))
        colnames(sumTable) <- c("id","totalDegree","totalContactDurations", paste("contactDuration_Indiv",indivSeq, sep = ""))
        sumTable$id = me
        sumTable$totalDegree <- length(specIndivSeq1)
        sumTable$totalContactDurations = sum(indivContact.full$contactDuration)
        sumTable[1,4:ncol(sumTable)] <- contactSum
        sumTable[,match(paste("contactDuration_Indiv",me, sep = ""), names(sumTable))] = NA
      }else{ #if nrow !>0
        if(nrow(indivContact.full) == 1){
          sumTable <- data.frame(matrix(ncol = (3+length(indivSeq)), nrow = 1))
          colnames(sumTable) <- c("id","totalDegree","totalContactDurations", paste("contactDuration_Indiv",indivSeq, sep = ""))
          sumTable$id = me
          sumTable$totalDegree <- 1
          sumTable$totalContactDurations = indivContact.full$contactDuration
          sumTable[1,4:ncol(sumTable)] <- 0
          sumTable[,match(paste("contactDuration_Indiv",specIndivSeq1, sep = ""), names(sumTable))] = indivContact.full$contactDuration
          sumTable[,match(paste("contactDuration_Indiv",me, sep = ""), names(sumTable))] = NA
        }
        if(nrow(indivContact.full) == 0){
          sumTable <- data.frame(matrix(ncol = (3+length(indivSeq)), nrow = 1))
          colnames(sumTable) <- c("id","totalDegree","totalContactDurations", paste("contactDuration_Indiv",indivSeq, sep = ""))
          sumTable$id = me
          sumTable[1,2:ncol(sumTable)] <- 0
          sumTable[,match(paste("contactDuration_Indiv",me, sep = ""), names(sumTable))] = NA
        }
      }
    }else{ # length(y$dyadMember1) == 0
      if(nrow(indivContact.full) > 1){
        areaSeqFrame <- data.frame(areaSeq)
        contactSum<-apply(areaSeqFrame, 1, distributeContacts2, indivContact.full)
        sumTable <- data.frame(matrix(ncol = (3+length(areaSeq)), nrow = 1))
        colnames(sumTable) <- c("id","totalDegree","totalContactDurations", paste("contactDuration_Area_",areaSeq, sep = ""))
        sumTable$id = me
        sumTable$totalDegree <- length(specIndivSeq1)
        sumTable$totalContactDurations = sum(indivContact.full$contactDuration)
        sumTable[1,4:ncol(sumTable)] <- contactSum
      }else{ #if nrow !>1
        if(nrow(indivContact.full) == 1){
          areaVec <- unique(y$area.id)
          sumTable <- data.frame(matrix(ncol = (3+length(areaSeq)), nrow = 1))
          colnames(sumTable) <- c("id","totalDegree","totalContactDurations", paste("contactDuration_Area_",areaSeq, sep = ""))
          sumTable$id = me
          sumTable$totalDegree <- 1
          sumTable$totalContactDurations = indivContact.full$contactDuration
          sumTable[1,4:ncol(sumTable)] <- 0
          sumTable[,match(paste("contactDuration_Area_",areaVec, sep = ""), names(sumTable))] = indivContact.full$contactDuration
        }
        if(nrow(indivContact.full) == 0){
          sumTable <- data.frame(matrix(ncol = (3+length(areaSeq)), nrow = 1))
          colnames(sumTable) <- c("id","totalDegree","totalContactDurations", paste("contactDuration_Area_",areaSeq, sep = ""))
          sumTable$id = me
          sumTable[1,2:ncol(sumTable)] <- 0
        }
      }			
    }
    return(sumTable)
  }
  blockSum <-function(x,y, indivSeq, areaSeq){
    blockDurFrame<-y[which(y$block == unname(unlist(x[1]))),]
    indivSeqFrame <- data.frame(indivSeq)
    summary.contacts<-apply(indivSeqFrame, 1, contSum, blockDurFrame, indivSeq, areaSeq)
    indivSum.full<- data.frame(data.table::rbindlist(summary.contacts))
    indivSum.full$block <- unname(unlist(x[1]))
    
    #added 02/05/2019 - to maintain this new information created in the newest version of the contactDur functions.
    indivSum.full$block.start <- unique(lubridate::as_datetime(blockDurFrame$block.start)) # updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
    indivSum.full$block.end <- unique(lubridate::as_datetime(blockDurFrame$block.end)) # updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
    indivSum.full$numBlocks <- unique(blockDurFrame$numBlocks)
    return(indivSum.full)
  }
  summaryAgg.block<-function(x,y){ #calculates the mean contacts from multiple summarizeContacts outputs (i.e., only applicable if avg == TRUE)
    sumTable<-y[which(y$id == unname(unlist(x[1])) & y$block == unname(unlist(x[2]))),]
    blockStart<- unique(lubridate::as_datetime(sumTable$block.start)) #added 02/05/2019 - had to keep track of this new information ; updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
    blockEnd<- unique(lubridate::as_datetime(sumTable$block.end)) #added 02/05/2019 - had to keep track of this new information ;  updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
    blockNum<- unique(sumTable$numBlocks) #added 02/05/2019 - had to keep track of this new information
    sumTable.redac<-sumTable[,-c(match("id", names(sumTable)),match("block", names(sumTable)), match("block.start", names(sumTable)), match("block.end", names(sumTable)), match("numBlocks", names(sumTable)))]  #Remove the columns that cannot/shoud not be averaged.
    contact.mean <- apply(sumTable.redac,2,mean, na.rm = TRUE)
    output = sumTable[1,]
    output[1,match("id", names(sumTable))] = unname(unlist(x[1])) #add this information back into the table
    output[1,match("block", names(sumTable))] = unname(unlist(x[2])) #add this information back into the table
    output[1,match("block.start", names(sumTable))] = blockStart #add this information back into the table
    output[1,match("block.end", names(sumTable))] = blockEnd #add this information back into the table
    output[1,match("numBlocks", names(sumTable))] = blockNum #add this information back into the table
    output[1,match(names(sumTable.redac), names(output))] = contact.mean
    return(output)
  }
  summaryAgg.NoBlock<-function(x,y){
    sumTable<-y[which(y$id == unname(unlist(x[1]))),]
    sumTable.redac<-sumTable[,-match("id", names(sumTable))] #Remove the columns that cannot/shoud not be averaged.
    contact.mean <- apply(sumTable.redac,2,mean, na.rm = TRUE)
    output = sumTable[1,]
    output[1,match("id", names(sumTable))] = unname(unlist(x[1])) #add this information back into the table
    output[1,match(names(sumTable.redac), names(output))] = contact.mean
    return(output)
  }
  summary.generator<-function(x, importBlocks){
    
    if(importBlocks == TRUE){
      
      if(length(x$dyadMember1) > 0){ #This essentially determines if the input was created with dist.all or distToArea. If length(dyadMember1) >0, it was created with dist.all
        x<-x[order(x$block,x$dyadMember1,x$dyadMember2),]
        indivVec <- c(as.character(x[,match("dyadMember1", names(x))]), as.character(x[,match("dyadMember2", names(x))]))
        areaSeq = NULL
      }else{
        x<-x[order(x$block,x$indiv.id),]
        indivVec <- x[,match("indiv.id", names(x))]
        areaVec <- x[,match("area.id", names(x))]
        areaVec <- areaVec[order(areaVec)] #forces the data type to become character so that there will be no issues with apply functions later.
        areaSeq<-as.character(unique(areaVec))
      }
      indivSeq <- unique(indivVec)
      indivSeq<-indivSeq[order(indivSeq)]
      indivSeq<-as.character(indivSeq) #forces the data type to become character so that there will be no issues with apply functions later.
      blockVecFrame <- data.frame(unique(as.character(x$block)))
      summary.block <- apply(blockVecFrame, 1, blockSum, x, indivSeq, areaSeq) #according to Dan, this apply function is faster than parApply, so I've removed the parApply option 1/17
      summaryTable<- data.frame(data.table::rbindlist(summary.block))
      summaryTable<-summaryTable[order(as.numeric(as.character(summaryTable$block)),summaryTable$id),]
      
    }else{ #importBlocks == FALSE
      if(length(x$dyadMember1) > 0){ #This essentially determines if the input was created with dist.all or distToArea. If length(dyadMember1) >0, it was created with dist.all
        x<-x[order(x$dyadMember1,x$dyadMember2),]
        indivVec <- c(as.character(x[,match("dyadMember1", names(x))]), as.character(x[,match("dyadMember2", names(x))])) #as.character forces the data type to become character so that there will be no issues with apply functions later.
        areaSeq = NULL
      }else{
        x<-x[order(x$indiv.id),]
        indivVec <- x[,match("indiv.id", names(x))]
        areaVec <- x[,match("area.id", names(x))]
        areaVec <- areaVec[order(areaVec)] #forces the data type to become character so that there will be no issues with apply functions later.
        areaSeq<-as.character(unique(areaVec))
      }
      indivSeq <- unique(indivVec)
      indivSeq<-indivSeq[order(indivSeq)]
      indivSeq<-as.character(indivSeq) #forces the data type to become character so that there will be no issues with apply functions later.
      indivSeqFrame <- data.frame(indivSeq)
      summary.contacts <- apply(indivSeqFrame, 1, contSum, x, indivSeq, areaSeq) #according to Dan, this apply function is faster than parApply, so I've removed the parApply option 1/17
      summaryTable<- data.frame(data.table::rbindlist(summary.contacts))
      summaryTable<-summaryTable[order(summaryTable$id),]
    }
    return(summaryTable)
  }
  
  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){
    summaryList<-lapply(x, summary.generator, importBlocks) #changed to lapply 02/02/2019
    
    if(avg == TRUE){
      full.summary<- data.frame(data.table::rbindlist(summaryList, fill = TRUE)) #Now we need to average the number of contacts by id and block
      idSeq<-unique(full.summary$id)
      if(importBlocks == TRUE){
        blockSeq<-unique(full.summary$block)
        aggTab<- expand.grid(as.character(idSeq),as.character(blockSeq))
        sumTab <- apply(aggTab, 1, summaryAgg.block, y = full.summary)
      }else{ #if importBlocks == FALSE
        aggTab<-data.frame(idSeq)
        sumTab <- apply(aggTab, 1, summaryAgg.NoBlock, y = full.summary)
      }
      sumTab.agg <- data.frame(data.table::rbindlist(sumTab))
      summary.output<-list(sumTab.agg, summaryList)
      names(summary.output)<-c("avg.","contactSummaries.")
    }else{ #if avg == FALSE
      summary.output<- summaryList
    }
  }else{ #if x is NOT a list
    summary.output <- summary.generator(x, importBlocks)
  }
  return(summary.output)
}
