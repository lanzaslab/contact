#' Determine if Observed Contacts are More or Less Frequent than in a Random Distribution
#'
#' This function is used to determine if tracked individuals in an empirical dataset had more or fewer contacts with other tracked individuals/specified locations than would be expected at random. The function works by comparing an empirically-based contactDur.all or contactDur.area function output (emp.input) to the contactDur.all or contactDur.area output generated from randomized data (rand.input).
#' 
#' The output of this function is a list contataining two data frames. The first data frame contains pairwise analyses of node degree and total edge weight (i.e., the sum of all observed contacts involving each individual). The second data frame contains results of pairwise analyses specific dyadic relationships (e.g., contacts between individuals 1 and 2). Each data frame contains the following columns: 
#' id1 - the id of the first individual involved in the contact.
#' id2 - designation of what is being compared (e.g., totalDegree, totalContactDurations, individual 2, etc.). Content will change depending on which data frame is being observed.
#' method - statistical test used to determine significance.
#' statistic - test statistic associated with the specific method.
#' p.value - p.values associated with each comparison.
#' df - degrees of freedom associated with the statistical test.
#' block (if applicable) - denotes the relevant time block for each analysis.
#' warning - denotes if any specific warning occurred during analysis.
#' Additional columns may exist depending on the statistical method employed. For example, if method == "chisq," additional columns include: 
#' empiricalContactDurations - Describes the number of observed events in emp.input.
#' randContactDurations.mean - Describes the average number of observed events in rand.input.
#' empiricalNoContactDurations - Describes the number of events that were not observed given the total number of potential events in emp.input.
#' randNoContactDurations.mean - Describes the average number of events that were not observed given the total number of potential events in rand.input.
#' difference - The value given by subtracting randContactDurations.mean from empiricalContactDurations.
#' 
#' This function was inspired by the methods described by Spiegel et al. 2016. Who determined individuals to be expressing social behavior when nodes had greater degree values than would be expected at random, with randomized contact networks derived from movement paths randomized according to their novel methodology (i.e., shuffle.type == 2). Here, however, by specifying a p-value threshold, users can also identify when more or fewer (demonstrated by the sign of values in the "difference" column) contacts with specific individuals than would be expected at random. Such relationships suggest social affinities or aversions, respectively, may exist between specific individuals.
#' 
#' Note:The default tested column (i.e., categorical data column from which data is drawn to be compared to randomized sets herein) is "id." This means that contacts involving each individual (defined by a unique "id") will be compared to randomized sets. Users may not use any data column for analysis other than "id." If users want to use another categorical data column in analyses rather than "id," we recommend re-processing data (starting from the dist.all/distToArea functions), while specifying this new data as an "id." For example, users may annotate an illness status column to the empirical input, wherein they describe if the tracked individual displayed gastrointestinal ("gastr"), respiratory ("respr"), both ("both"), illness symptoms, or were consistently healthy ("hel") over the course of the tracking period. Users could set this information as the "id," and carry it forward as such through the data-processing pipeline. Ultimately, they could determine if each of these disease states affected contact rates, relative to what would be expected at random.
#' 
#' Note: Currently, this function requires defined TSWs (temporal sampling windows; see the tempAggregate function) (i.e., contact durations must represent time periods that are equidistant).
#' 
#' Note: The current functionality is limited to comparisons using the X-squared test of independence. This works by comparing the number of TSWs in a given dataset/block between individuals/fixed locations and the number of TSWs during which individuals were observed, but were not in contact with specific individuals/places, for empirical (x) and randomized (y) datasets. The dist (z ; from dist.all or distToArea output) input is used here to determine how frequently each individual was observed in the empirical dataset/block of interest, allowing us to calculate the number of TSWs each individual was present but not involved in contacts. Note here that if X-squared expected values will be very small, approximations of p may not be right (and in fact, all estimates will be poor). It may be best to weight these tests differently. To address this, I've added the "warning" column to the output which notifies users when the chi-sq function reported that results may be inaccurate.
#' 
#' Note: if importBlocks == TRUE, a "block" column MUST exist in emp.input. However, if a "block" column need not exist in rand.contact. If no "block" column exists in rand.input, empirical values in all emp.input blocks will be compared to the overall average values in rand.input. 
#' 
#' @param emp.input List or data frame containing contactDur.all or contactDur.area output refering to the empirical data.
#' @param rand.input List or data frame containing contactDur.all or contactDur.area output refering to the randomized-path data.
#' @param dist.input List or data frame containing dist.all-/distToArea-function output refering to the empirical data.
#' @param test Character string. Describes the statistical test used to evaluate differences. Currently only takes the value "chisq," more tests will be added in later versions.
#' @param importBlocks Logical. If true, each block in emp.input will be analyzed separately. Defaults to FALSE. Note that the "block" column must exist in emp.input.
#' @param shuffle.type Numerical. describing which shuffle.type was used to randomize the rand.input data set(s). Takes the values "0," "1," or "2."
#' @keywords network-analysis social-network
#' @references Spiegel, O., Leu, S.T., Sih, A., and C.M. Bull. 2016. Socially interacting or indifferent neighbors? Randomization of movement paths to tease apart social preference and spatial constraints. Methods in Ecology and Evolution 7:971-979. https://doi.org/10.1111/2041-210X.12553.
#' @export
#' @examples
#' #load the calves data set
#' data(calves)
#' 
#' #pre-process the data
#' calves.dateTime<-datetime.append(calves, date = calves$date, time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, point.y = calves.dateTime$y, secondAgg = 10, extrapolate.left = FALSE, extrapolate.right = FALSE, resolutionLevel = "Full", parallel = TRUE, na.rm = FALSE, smooth.type = 1) #smooth locations to 10-second fix intervals. Note that na.rm was set to "FALSE" because randomizing this data set according to Spiegel et al.'s method (see below) requires equidistant time points.
#'
#' #generate empirical time-ordered network edges.
#' calves.dist<-dist2All_df(x = calves.agg, parallel = TRUE, dataType = "Point", lonlat = FALSE) #calculate distance between all individuals at each timepoint.
#' calves.contact.block<-contactDur.all(x = calves.dist, dist.threshold=1, sec.threshold=10, blocking = TRUE, blockUnit = "hours", blockLength = 1, equidistant.time = FALSE, parallel = TRUE, reportParameters = TRUE) #compile inter-calf contacts with 1-hr blocking. Contacts are defined here as occurring when calves were within 1 m of one another.
#' 
#' #generate randomized time-ordered network edges. (2 examples)
#' calves.agg.rand1<-randomizePaths(x = calves.agg, id = "id", dateTime = "dateTime", point.x = "x", point.y = "y", poly.xy = NULL, parallel = TRUE, dataType = "Point", numVertices = 1, blocking = TRUE, blockUnit = "mins", blockLength = 10, shuffle.type = 2, shuffleUnit = "Hours", indivPaths = TRUE, numRandomizations = 1) #create 1 replicate of an hour (shuffleUnit) with 10-minute blocks shuffled around. This randomization process was first described by Spiegel et al. 2016.
#' calves.agg.rand2<-randomizePaths(x = calves.agg, id = "id", dateTime = "dateTime", point.x = "x", point.y = "y", poly.xy = NULL, parallel = TRUE, dataType = "Point", numVertices = 1, blocking = TRUE, blockUnit = "mins", blockLength = 10, shuffle.type = 0, shuffleUnit = NA, indivPaths = TRUE, numRandomizations = 1) #create 1 replicate of the calves.agg data set with calves' xy coordinates within 10-minute blocks pseudo-randomized.
#' 
#' calves.dist.rand1<-dist2All(x = calves.agg.rand1, point.x = "x.rand", point.y = "y.rand", parallel = TRUE, dataType = "Point", lonlat = FALSE) #calculate distance between all individuals at each timepoint in calves.rand1. Note that point.x and point.y must be specified as "x.rand" and "y.rand," respectively.
#' calves.dist.rand2<-dist2All(x = calves.agg.rand2, point.x = "x.rand", point.y = "y.rand", parallel = TRUE, dataType = "Point", lonlat = FALSE) #calculate distance between all individuals at each timepoint in calves.rand2. Note that point.x and point.y must be specified as "x.rand" and "y.rand," respectively.
#' calves.contact.rand1<-contactDur.all(x = calves.dist.rand1, dist.threshold=1, sec.threshold=10, blocking = FALSE, equidistant.time = FALSE, parallel = TRUE, reportParameters = TRUE) #Because the data sets contained in calves.dist.rand1 cover timepoints within a single hour, this is equivalent to randomized contacts within a one-hour block. Contacts are defined here as occurring when calves were within 1 m of one another.
#' calves.contact.rand2<-contactDur.all(x = calves.dist.rand2, dist.threshold=1, sec.threshold=10, blocking = TRUE, blockUnit = "hours", blockLength = 1, equidistant.time = FALSE, parallel = TRUE, reportParameters = TRUE) #compile randomized inter-calf contacts with 1-hr blocking. Contacts are defined here as occurring when calves were within 1 m of one another.
#' 
#' #No block
#' contactTest1<- contactTest(emp.input = calves.contact.block, rand.input = calves.contact.rand2, dist.input = calves.dist, importBlocks = FALSE, shuffle.type = 0)
#' #Block
#' ##Both of these examples show the comparisons at each hour in the dataset, however, they have different associated assumptions because rand.input was calculated differently
#' contactTest2<- contactTest(emp.input = calves.contact.block, rand.input = calves.contact.rand1, dist.input = calves.dist, importBlocks = TRUE, shuffle.type = 2)
#' contactTest3<- contactTest(emp.input = calves.contact.block, rand.input = calves.contact.rand2, dist.input = calves.dist, importBlocks = TRUE, shuffle.type = 0)

contactTest<-function(emp.input, rand.input, dist.input, test = "chisq", importBlocks = FALSE, shuffle.type = 0){
  
  summarizeContacts<- function(x, importBlocks, avg){
    
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
    summaryAgg.block<-function(x,y){ #calculates the mean contacts from multiple summarizeContacts outputs
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
  
  chisq.forLoop<-function(x, empirical = x, randomized = y, dist = dist.input, x.blocking, y.blocking){ #I hate that I have to do this in a for-loop, but I couldn't get the apply functions (above) to work.
    output<-NULL
    
    #vectorize the input
    ids<-x[,1]
    cols<-x[,2]
    shuffle.type<-unique(x[,ncol(x)])
    
    if(x.blocking == TRUE){
      blocks<-x[,3]
      for(i in 1:nrow(x)){
        
        empDurations<- empirical[which(empirical$id == unlist(unname(ids[i])) & empirical$block == unlist(unname(blocks[i]))), grep(unlist(unname(cols[i])),colnames(empirical))] #pull the value of a given column for a specific id in a specific block
        empDurations<-empDurations[is.na(empDurations) == FALSE] #remove the NAs
        summaryFrame<-NULL
        
        if(length(empDurations) > 0){ #if there is no entry OR is.na == TRUE, nothing will happen
          
          block <- unlist(unname(blocks[i]))
          warn = ""
          if(grep(unlist(unname(cols[i])),colnames(empirical)) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in the dist.input.
            maxDurations = length(which(dist$id == unlist(unname(ids[i])) & is.na(dist[,(grep(unlist(unname(cols[i])),colnames(empirical)) + 1)]) == FALSE & dist$block == unlist(unname(blocks[i])))) #This means the the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. For example, if one of the individuals disappears because we have no recorded time points for them at a given time, we can't say that they have the potential to be in contact with others during this period. Also, note that the "+1" here is because column 4 in empirical refers to contacts with indiv1, but in the dist data, column 5 (i.e., 4 + 1) refers to distance to indiv1, while column 4 is individuals' "id."
          }else{ #i.e., x[2] == 2 or 3
            if(grep(unlist(unname(cols[i])),colnames(empirical)) == 2){ #if cols = totalDegree, the maximum degree possible would be the number of individuals observed during the time period
              maxDurations<- max(dist[which(dist$id == unlist(unname(ids[i])) & dist$block == unlist(unname(blocks[i]))), match("individualsAtTimestep", names(dist))])
            }
            if(grep(unlist(unname(cols[i])),colnames(empirical)) == 3){ #if cols = totalcontactDurations, the maximum possible number of contacts would be would be the sum of all individuals observed at each time step during the time period, excluding individual i (hence the "-1" below).
              maxDurations<- sum((dist[which(dist$id == unlist(unname(ids[i])) & dist$block == unlist(unname(blocks[i]))), match("individualsAtTimestep", names(dist))] - 1))
            }
          }
          
          if(maxDurations == 0){ #if there was no potential for contacts to occur, the loop moves on
            next
          }
          
          if(length(grep(unlist(unname(cols[i])),colnames(randomized))) == 0){ #if i is trying to reference a specific object not present in the random set, rand duration will equal "0."
            randDurations <- 0
          }else{ # if the relavent column DOES exist in randomized
          
            if(y.blocking == TRUE){ #If there are blocks in the randomized input
          
              if(shuffle.type == 2){ #recall that shuffle.type 2 (from the randomizeLocations function) produces only 1 shuffle.unit's worth of data, rather than a dataset with the same length of x. As such, there may be a different number of blocks in y compared to x. Here we assume that the mean randomized durations per block, are representative of mean randomized durations per block across each shuffle unit (e.g., day)
                if(unlist(unname(blocks[i])) > max(randomized$blocks)){
                  block = unlist(unname(blocks[i])) - (ceiling((unlist(unname(blocks[i])) - max(randomized$blocks)) /max(randomized$blocks))*max(randomized$blocks))
                  randDurations<-randomized[which(randomized$id == unlist(unname(ids[i])) & randomized$block == block), grep(unlist(unname(cols[i])),colnames(randomized))]
                }
              }else{ #if there is any other shuffle type than 2
                randDurations<-randomized[which(randomized$id == unlist(unname(ids[i])) & randomized$block == unlist(unname(blocks[i]))), grep(unlist(unname(cols[i])),colnames(randomized))]
              }
            }else{ #if y.blocking == FALSE (i.e., if there are no blocks in the randomized input)
              randDurations<-randomized[which(randomized$id == unlist(unname(ids[i]))), grep(unlist(unname(cols[i])),colnames(randomized))]
            }
          }
          
          randDurations<-randDurations[is.na(randDurations) == FALSE] #remove the NAs
          if(length(randDurations) == 0){ #if removing the NAs removed the only observations, then randDurations reverts to 0
            randDurations <- 0
          }
          randDurations <- round(randDurations) #because randDurations likely represents the averaged durations from randomized data sets, there is potential for it to be a non-integer. We need it to be an integer to prevent errors from arising due to data frame columns of differing lengths below.
          compareFrame<- data.frame(matrix(ncol = 2, nrow = 2*maxDurations))
          colnames(compareFrame)<-c("contact", "data")
          compareFrame[,1] <- c(rep(1,sum(empDurations)), rep(0,maxDurations-sum(empDurations)),rep(1,sum(randDurations)), rep(0,maxDurations-sum(randDurations)))
          compareFrame[,2] <- c(rep("empirical",maxDurations), rep("randomized",maxDurations)) 
          
          compareTable <- table(compareFrame$data,compareFrame$contact) #generates a table showing the total number of contacts (1) and non-contacts (0) for the empirical and randomized sets
          warning() #clears the warnings
          test<- chisq.test(compareTable)
          warn<- names(warnings())[1] #used to identify when the chi-sq.test produces the "Chi-squared approximation may be incorrect" warning, indicating that that expected X-squared values are very small, and generated estimates will be poor.
          warning() #clears the warnings
          summaryFrame <- data.frame(id1 = unlist(unname(ids[i])), id2 = unlist(unname(cols[i])), method = test[4], X.squared = test[1], df = test[2], p.val = test[3], empiricalContactDurations = sum(empDurations), randContactDurations.mean = sum(randDurations), empiricalNoContactDurations = (maxDurations - sum(empDurations)), randNoContactDurations.mean = (maxDurations - sum(randDurations)), difference = abs((maxDurations - sum(empDurations)) - (maxDurations - sum(randDurations))), empBlock = unlist(unname(blocks[i])), randBlock = block, warning = warn)
          rownames(summaryFrame) <-1

          bindlist<-list(output, summaryFrame)
          output<- data.table::rbindlist(bindlist)
        }else{ #if there is no entry OR is.na == TRUE, nothing will happen
          next
        }
      }
    }else{ #if blocking == False
      for(i in 1:nrow(x)){
        
        empDurations<- empirical[which(empirical$id == unlist(unname(ids[i]))),grep(unlist(unname(cols[i])),colnames(empirical))] #pull the value of a given column for a specific id
        empDurations<-empDurations[is.na(empDurations) == FALSE] #remove the NAs
        summaryFrame<-NULL
        
        if(length(empDurations) > 0 && is.na(empDurations) == FALSE){ #if there is no entry OR is.na == TRUE, nothing will happen
          warn = ""
          if(grep(unlist(unname(cols[i])),colnames(empirical)) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in the dist.input.
            maxDurations = length(which(dist$id == unlist(unname(ids[i])) & is.na(dist[,(grep(unlist(unname(cols[i])),colnames(empirical)) + 1)]) == FALSE)) #This means the the max number of durations potentially observed is the number of TSWs both individuals (or an individual and fixed area) were observed at the same time. For example, if one of the individuals disappears because we have no recorded time points for them at a given time, we can't say that they have the potential to be in contact with others during this period. Also, note that the "+1" here is because column 4 in empirical refers to contacts with indiv1, but in the dist data, column 5 (i.e., 4 + 1) refers to distance to indiv1, while column 4 is individuals' "id."
          }else{ #i.e., x[2] == 2 or 3
            if(grep(unlist(unname(cols[i])),colnames(empirical)) == 2){ #if cols = totalDegree, the maximum degree possible would be the number of individuals observed during the time period
              maxDurations<- max(dist[which(dist$id == unlist(unname(ids[i]))), match("individualsAtTimestep", names(dist))])
            }
            if(grep(unlist(unname(cols[i])),colnames(empirical)) == 3){ #if cols = totalcontactDurations, the maximum possible number of contacts would be would be the sum of all individuals observed at each time step during the time period, excluding individual i (hence the "-1" below).
              maxDurations<- sum((dist[which(dist$id == unlist(unname(ids[i]))), match("individualsAtTimestep", names(dist))] - 1))
            }
          }
          if(maxDurations == 0){  #if there was no potential for contacts to occur, the loop moves on
            next
          }
          
          if(length(grep(unlist(unname(cols[i])),colnames(randomized))) == 0){ #if i is trying to reference a specific object not present in the random set, rand duration will equal "0."
            randDurations <- 0
          }else{ # if the relavent column DOES exist in randomized
            randDurations<-randomized[which(randomized$id == unlist(unname(ids[i]))), grep(unlist(unname(cols[i])),colnames(randomized))]
          }
          randDurations<-randDurations[is.na(randDurations) == FALSE] #remove the NAs
          if(length(randDurations) == 0){ #if removing the NAs removed the only observations, then randDurations reverts to 0
            randDurations <- 0
          }
          randDurations<-round(randDurations) #because randDurations likely represents the averaged durations from randomized data sets, there is potential for it to be a non-integer. We need it to be an integer to prevent errors from arising due to data frame columns of differing lengths below.
          compareFrame<- data.frame(matrix(ncol = 2, nrow = 2*maxDurations))
          colnames(compareFrame)<-c("contact", "data")
          compareFrame[,1] <- c(rep(1,sum(empDurations)), rep(0,maxDurations-sum(empDurations)),rep(1,sum(randDurations)), rep(0,maxDurations-sum(randDurations)))
          compareFrame[,2] <- c(rep("empirical",maxDurations), rep("randomized",maxDurations)) 
          
          compareTable <- table(compareFrame$data,compareFrame$contact) #generates a table showing the total number of contacts (1) and non-contacts (0) for the empirical and randomized sets
          warning() #clears the warnings
          test<- chisq.test(compareTable)
          warn<- names(warnings())[1] #used to identify when the chi-sq.test produces the "Chi-squared approximation may be incorrect" warning, indicating that that expected X-squared values are very small, and generated estimates will be poor.
          warning() #clears the warnings
          summaryFrame <- data.frame(id1 = unlist(unname(ids[i])), id2 = unlist(unname(cols[i])), method = test[4], X.squared = test[1], df = test[2], p.val = test[3], empiricalContactDurations = sum(empDurations), randContactDurations.mean = sum(randDurations), empiricalNoContactDurations = (maxDurations - sum(empDurations)), randNoContactDurations.mean = (maxDurations - sum(randDurations)), difference = abs((maxDurations - sum(empDurations)) - (maxDurations - sum(randDurations))), warning = warn)
          rownames(summaryFrame) <-1

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
  
  x.blocking<-importBlocks

  
  if(is.data.frame(emp.input) == FALSE & is.list(emp.input) == TRUE){ #if the emp.input is a list (not a data frame), the function assumes only the first list entry is relevant to our purposes.
    emp.input = data.frame(emp.input[1])
  }
  
  if(is.data.frame(dist.input) == FALSE & is.list(dist.input) == TRUE){ #if the dist.input input is a list (not a data frame), the function assumes only the first list entry is relevant to our purposes.
    dist.input = data.frame(dist.input[1]) 
  }
  

  #summarize the contacts in each input
  x<- summarizeContacts(emp.input, importBlocks, avg = FALSE)
  if(is.data.frame(rand.input) == FALSE & is.list(rand.input) == TRUE){ #if the rand.input is a list (not a data frame),
    rand.frame<-data.frame(rand.input[1])#we need to know if a "block" column exists in the data frames compiled into the rand.input list.
    
    if(x.blocking == TRUE & length(rand.frame$block) > 0){ #if blocking is true and there are blocks denoted in rand.input
      y<- summarizeContacts(rand.input, importBlocks = TRUE, avg = TRUE)
      y.blocking <-TRUE
    }
    if(x.blocking == TRUE & length(rand.frame$block) == 0){#if blocking is true, but there are no blocks in rand.input
      y<- summarizeContacts(rand.input, importBlocks = FALSE, avg = TRUE)
      y.blocking <-FALSE
    }
    if(x.blocking == FALSE){#if blocking is false
      y<- summarizeContacts(rand.input, importBlocks = FALSE, avg = TRUE)
      y.blocking <-FALSE
    }
  }else{ #if rand.input is a data frame
    if(x.blocking == TRUE & length(rand.input$block) > 0){ #if blocking is true and there are blocks denoted in rand.input
      y<- summarizeContacts(rand.input, importBlocks = TRUE, avg = TRUE)
      y.blocking <-TRUE
    }
    if(x.blocking == TRUE & length(rand.input$block) == 0){#if blocking is true, but there are no blocks in rand.input
      y<- summarizeContacts(rand.input, importBlocks = FALSE, avg = TRUE)
      y.blocking <-FALSE
    }
    if(x.blocking == FALSE){#if blocking is false
      y<- summarizeContacts(rand.input, importBlocks = FALSE, avg = FALSE)
      y.blocking <- FALSE
    }
  }
  
  if(is.data.frame(y) == FALSE & is.list(y) == TRUE){ #if the y input is a list (not a data frame), the function assumes only the first list entry is relevant to our purposes.
    y = data.frame(y[1]) #in the case of y (a.k.a. the randomized data set), the first list entry is likely to be the average summary report (because rand.input is likely a list of data sets, and the summarizeContacts subfunction processes y with "avg" ==  TRUE). 
  }
  #if y represents the average summary report, all colnames will be preceded by "avg..". This code block removes that tag from all relevant columns.
  y.colnameSubstring<- substring(names(y),1,5) 
  alterNames.y<-which(y.colnameSubstring == "avg..")
  if(length(alterNames.y) > 0){
    names(y)<- substring(names(y)[alterNames.y],6,10000)
  }
  
  x$id<-as.character(x$id)
  y$id<-as.character(y$id)
  dist.input$id<-as.character(dist.input$id)
  
  if(test == "chisq" || test == "CHISQ" || test == "Chisq" || test == "chisquare" || test == "CHISQUARE" || test == "Chisquare" || test == "chi-square" || test == "CHI-SQUARE" || test == "Chi-square"){
    
    idSeq <- unique(x$id) #pulls the unique ids for each individual
    #empColSeq <- 2:max(grep("contactDuration_", names(x))) #pulls the column number for columns of interest for our analysis. (i.e., total degree, total durations, and durations in contact with each specific object)
    empColnames<- c("totalDegree", "totalContactDurations", substring(names(x[,4:max(grep("contactDuration_", names(x)))]), 22, 1000000)) #identifies which object each column relates to
    #randColnames<- c("totalDegree", "totalContactDurations", substring(names(x[,4:max(grep("contactDuration_", names(y)))]), 22, 1000000)) #identifies which object each column relates to
    
    if(x.blocking == TRUE){
      
      #Because we need to know the maximum number of potential blocks, we pull this information from the empirical numBlocks column
      blockSeq <- seq(1,as.numeric(as.character(unique(x$numBlocks))),1)
      #We also need to know how large each block is (e.g., how many seconds are in each block)
      blockLength<-round(as.numeric(difftime(time1 = unique(x$block.start)[2], time2 = unique(x$block.start)[1], units = "secs"))) #calculates the number of seconds in the first block of the empirical data set. Because all blocks are the same length, we know 
      
      #We need to use the dist.input file to determine the number of timesteps individuals were actually observed in each block (if individuals were not present, then by definition no contact could have occurred).
      #So, we need to ensure that timeblocks are appended to the dist.input file.
      dist.input <- timeBlock.append(dist.input, dateTime = NULL, blockLength, blockUnit = "secs") #dateTime = NULL assumes a dateTime column exists in dist.input
    }

      if(x.blocking == TRUE){
          testFrame1 <- expand.grid(idSeq, empColnames, blockSeq)
        }else{
          testFrame1 <- expand.grid(idSeq, empColnames)
        }
        testFrame1$shuffle.type <-shuffle.type
        
        ###for some reason the apply functions below would only ever return NULL values (despite having removed all factor information from input data), so I had to switch to a for-loop. If we can get the applies to work, this would be much more efficient
        
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
      testResultsFrame.indiv<- chisq.forLoop(testFrame1, empirical = x, randomized = y, dist = dist.input, x.blocking, y.blocking)
  
  #  }else{ #if blocking != TRUE
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
 #   }
    
  }
  totalDegree_Durations<-subset(testResultsFrame.indiv, id2 == "totalDegree" | id2 == "totalContactDurations")
  specific_contacts <- subset(testResultsFrame.indiv, id2 != "totalDegree" & id2 != "totalContactDurations")
  final_out<-list(totalDegree_Durations, specific_contacts)
  return(final_out)
}
