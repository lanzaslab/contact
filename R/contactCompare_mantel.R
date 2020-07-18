#' Statistically Compare Two Contact Matrices
#'
#' Tests for similarity of the x.summary input to y.summary. Please note that 
#'    this is a function of convience that is essentially a wrapper for the 
#'    ape::mantel.test function, that allows users to easily compare contact
#'    networks created using our pipeline of contact:: functions.
#'    Please understand that ape::mantel.test does not allow for missing
#'    values in matrices, so all NAs will be treated as zeroes.
#'
#' @param x.summary List or single-data frame output from the summarizeContacts
#'    function refering to the empirical data. Note that if x.summary is a list
#'    of data frames, only the first data frame will be used in the function.
#' @param y.summary List or single-data frame output from the summarizeContacts
#'    function refering to the randomized data (i.e., NULL model 
#'    contact-network edge weights). Note that if y.summary is a list
#'    of data frames, only the first data frame will be used in the function.
#' @param numPermutations Integer. Number of times to permute the data.
#'    Defaults to 1000.
#' @param alternative.hyp Character string. Describes the nature of the 
#'    alternative hypothesis being tested when test == "mantel." Takes the 
#'    values "two.sided," "less," or "greater." Defaults to "two.sided."
#' @param importBlocks Logical. If true, each block in x.summary will be 
#'    analyzed separately. Defaults to FALSE. Note that the "block" column must
#'    exist in .summary objects AND values must be identical (i.e., if block 
#'    100 exists in x.summary, it must also exist in y.summary), otherwise an 
#'    error will be returned.
#'    
#' @keywords network-analysis 
#' @return Output format is a single data frame with the following columns.
#' 
#'    \item{method}{Statistical test used to determine significance.}
#'    \item{z.val}{z statistic.}
#'    \item{p.value}{p.values associated with each comparison.}
#'    \item{x.mean}{mean contacts in x.summary overall or by block (if 
#'    applicable). Note that these means are calculated BEFORE any NAs are 
#'    converted to zeroes (see note above)}
#'    \item{y.mean}{mean contacts in y.summary overall or by block (if
#'    applicable). Note that these means are calculated BEFORE any NAs are 
#'    converted to zeroes (see note above)} 
#'    \item{alternative.hyp}{The nature of the alternative hypothesis being 
#'    tested.}
#'    \item{nperm}{Number of permutations used to generate p value.}
#'    \item{warning}{Denotes if any specific warning occurred during analysis.}
#'    
#' @references Mantel, N. 1967. The detection of disease clustering and a 
#'    generalized regression approach. Cancer Research, 27:209–220.
#'
#' @export
#' @examples
#' \donttest{
#' data(calves) #load data
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date,
#'                                  time = calves$time) #add dateTime column
#' 
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag,
#'                        dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x,
#'                        point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE,
#'                        extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE,
#'                        na.rm = TRUE, smooth.type = 1) #aggregate to 5-min timepoints
#' 
#' calves.dist<-dist2All_df(x = calves.agg, parallel = FALSE,
#'                        dataType = "Point", lonlat = FALSE) #calculate  inter-calf distances
#' 
#' calves.contact.block<-contactDur.all(x = calves.dist, dist.threshold=1,
#'                        sec.threshold=10, blocking = TRUE, blockUnit = "hours", blockLength = 1,
#'                        equidistant.time = FALSE, parallel = FALSE, reportParameters = TRUE)
#' 
#' emp.summary <- summarizeContacts(calves.contact.block, 
#'                                  importBlocks = TRUE) #empirical contact summ.
#' 
#' 
#'
#' calves.agg.rand<-randomizePaths(x = calves.agg, id = "id",
#'                        dateTime = "dateTime", point.x = "x", point.y = "y", poly.xy = NULL,
#'                        parallel = FALSE, dataType = "Point", numVertices = 1, blocking = TRUE,
#'                        blockUnit = "mins", blockLength = 20, shuffle.type = 0, shuffleUnit = NA,
#'                        indivPaths = TRUE, numRandomizations = 2) #randomize calves.agg
#' 
#' calves.dist.rand<-dist2All_df(x = calves.agg.rand, point.x = "x.rand",
#'                        point.y = "y.rand", parallel = FALSE, dataType = "Point", lonlat = FALSE)
#' 
#' calves.contact.rand<-contactDur.all(x = calves.dist.rand,
#'                        dist.threshold=1, sec.threshold=10, blocking = TRUE, blockUnit = "hours",
#'                        blockLength = 1, equidistant.time = FALSE, parallel = FALSE,
#'                        reportParameters = TRUE) #NULL model contacts (list of 2)
#' 
#' rand.summary <- summarizeContacts(calves.contact.rand, avg = TRUE,
#'                                   importBlocks = TRUE) #NULL contact summary
#' 
#' 
#' contactCompare_mantel(x.summary = emp.summary, y.summary = rand.summary,
#'                     importBlocks = FALSE, numPermutations = 100,
#'                     alternative.hyp = "two.sided") #no blocking
#' 
#' contactCompare_mantel(x.summary = emp.summary, y.summary = rand.summary,
#'                     importBlocks = TRUE, numPermutations = 100,
#'                     alternative.hyp = "two.sided") #blocking
#'    }

contactCompare_mantel<-function(x.summary, y.summary, numPermutations = 1000, alternative.hyp = "two.sided", importBlocks = FALSE){
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  block <- NULL
  
  #define sub-function
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
      x.summary<- stats::aggregate(x.summaryRedac[,-match("id", colnames(x.summary))], list(id = x.summaryRedac$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set.
      ###so now we have an x.summary object with no block information. However, we must recognize that the "totalDegree" column cannot be summed like other columns can (it is the number of individuals/fixed areas that each unique individual was observed in contact with.)
      ###below we recalculate the totalDegree for each individual. Note that this will be a different value depending on whether or not x.summary was based on dist2All or dist2Area outputs.
      indivSummaryTest.x<- ifelse(length(grep("contactDuration_Indiv", colnames(x.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
      
      if(indivSummaryTest.x == TRUE){ #if x.summary was derived from dist2All
        x.summary$totalDegree <- length(unique(x.summary$id)) - 1 #the potential degree is the maximum number of edges that can extend from any one node in the network
      }
      if(indivSummaryTest.x == FALSE){ #if x.summary was derived from dist2Area
        x.summary$totalDegree <- length(grep("contactDuration_", colnames(x.summary))) #the potential degree is the maximum number of edges that can extend from any one node in the network
      }
    }
    
  }
  
  #we do the same thing for y.summary as we did for x.summary above.
  if(is.data.frame(y.summary) == FALSE & is.list(y.summary) == TRUE){ #Because R treats dataframes as lists, we assess here if the input was only a single data frame or a list of data frames.
    y.summary <- y.summary[[1]] #if y.summary is a list of data frames, only the first one is used in processing.
  }
  
  colnames(y.summary)<- gsub("avg..","",colnames(y.summary)) #if y.summary represents the average summary report (from the summarizeContacts function), all colnames may be preceded by "avg..". This line removes that tag from all relevant columns.
  
  if(importBlocks == TRUE){ #if importBlocks == TRUE, but there is no block information in y.summary, the function returns a warning, but precedes as if the y input is relevant to all columns.
    
    if(length(y.summary$block) == 0){
      
      warning("importBlocks set to TRUE, but no block column exists in y.summary. Proceeding as if y.summary values are stable across time and relevant to EVERY block.", immediate. = TRUE)
      
      y.summaryBlock<-NULL #create an empty object to contain new block information
      
      for(i in unique(x.summary$block)){
        
        y.summary$block <- i #add a block column containing only the i value to y.summary
        y.summaryBlock <- data.table::rbindlist(list(y.summaryBlock, y.summary), fill = TRUE) #rbind the data frames together
        
      }
      
      y.summary <- data.frame(y.summaryBlock, stringsAsFactors = TRUE) #redefine y.summary as the object containing block information
    }
    
  }
  
  if(length(y.summary$block) > 0){ #Alternatively, if there are blocked time sets in the empirical set, we need to dictate how to handle them.
    
    if(importBlocks == FALSE){ #if there ARE blocks in y.summary, but importBlocks == FALSE, we have to recreate y.summary without the bloking information.
      
      y.summaryRedac <- droplevels(y.summary[, - c(which(colnames(y.summary) == "block"):ncol(y.summary))]) #remove the block information (four columns: block, block.start, block.end, and numBlocks)
      y.summary<- stats::aggregate(y.summaryRedac[,-match("id", colnames(y.summary))], list(id = y.summaryRedac$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set.
      ###so now we have an y.summary object with no block information. However, we must recognize that the "totalDegree" column cannot be summed like other columns can (it is the number of individuals/fixed areas that each unique individual was observed in contact with.)     
      ###below we recalculate the totalDegree for each individual.Note that this will be a different value depending on whether or not y.summary was based on dist2All or dist2Area outputs.
      indivSummaryTest.y<- ifelse(length(grep("contactDuration_Indiv", colnames(y.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
      
      if(indivSummaryTest.x == TRUE){ #if x.summary was derived from dist2All
        y.summary$totalDegree <- length(unique(y.summary$id)) - 1 #the potential degree is the maximum number of edges that can extend from any one node in the network
      }
      if(indivSummaryTest.x == FALSE){ #if x.summary was derived from dist2Area
        y.summary$totalDegree <- length(grep("contactDuration_", colnames(y.summary))) #the potential degree is the maximum number of edges that can extend from any one node in the network
      }
    }
    
  }
  
  #before executing the mantel tests we must ensure that warnings for sub-level functions will not be silenced (so that we can record them as they occur, allowing users to pinpoint what part of their output might be erroneous).
  oldw <- getOption("warn") #pull the current warn setting 
  options(warn = 1) #make warnings appear as they occur, rather than be stored until the top-level function returns
  on.exit(options(warn = oldw)) #ensure that when the function ends, users' options are reset.
  
  #now we run the mantel test
  
  if(importBlocks == FALSE){ #if there is no blocking
    
    #make them matrices for x & y summaries with equal dimensions
    emp.matrix<-as.matrix(x.summary[,4:(3+nrow(x.summary))]) #create the matrix detailing observed contacts. Note that this matrix excludes id, totalDegree, and totalContactDuration columns
    rownames(emp.matrix)<-x.summary$id #rename the rows of the matrix the ids

    randCols<-which(colnames(y.summary)%in%colnames(emp.matrix)) #pull the colnumbers of columns in y.summary that appear in the empirical matrix as well.
    randRows<-which(y.summary$id%in%x.summary$id) #pull the ids in y.summary that also exist in the x.summary set
    y.redac<-droplevels(y.summary[randRows, c(1,randCols)]) #reduce y.summary to only include individuals present in the x.summary (note that id column is also appended here)
    
    rand.matrix<-matrix(ncol = ncol(emp.matrix), nrow = nrow(emp.matrix)) #create the empty random-set matrix of appropriate length.
    rownames(rand.matrix)<-x.summary$id #name the rows appropriately
    colnames(rand.matrix)<-colnames(emp.matrix) #name the columns appropriately
    for(i in 1:ncol(rand.matrix)){ #loop through the randRows and randCols values to add data to the rand.matrix
      if(is.na(match(colnames(rand.matrix)[i], colnames(y.redac))) == TRUE){next} #if the column was not found in y.redac, the loop moves on
      rand.matrix[match(y.redac$id, rownames(rand.matrix)),i]<-y.redac[,match(colnames(rand.matrix)[i], colnames(y.redac))] #pull the values from y.redac into rand.matrix
    }
    
    emp.mean<-mean(emp.matrix, na.rm = TRUE) #calculate the mean number of contacts in the empirical set
    rand.mean <- mean(rand.matrix, na.rm = TRUE) #calculate the mean number of contacts in the random set
    
    #the ape::mantel.test function below DOES NOT accept NAs in the input data frames, so have to replace NAs with 0
    emp.matrix[which(is.na(emp.matrix) == TRUE)]<-0 
    rand.matrix[which(is.na(rand.matrix) == TRUE)]<-0
    
    test.tryCatch<-tryCatch.W.E(ape::mantel.test(emp.matrix, rand.matrix, nperm = numPermutations, alternative = alternative.hyp)) #perform the mantel test
    output.test<-test.tryCatch$value #pull the value
    warn1<-test.tryCatch$warning$message #used to identify when the mantel.test produces any errors.
    
    if(length(warn1) == 0){
      warn1 = ""
    }
    
    summaryFrame <- data.frame(method = "Mantel's permutation test for similarity of two matrices", z.val = output.test[1], alternative.hyp = output.test[3], nperm = numPermutations, p.val = output.test[2], x.mean = emp.mean, y.mean = rand.mean, warning = warn1, stringsAsFactors = TRUE)
    final_out <- summaryFrame #redefine summaryFrame as the final function output.
  }
  
  if(importBlocks == TRUE){ #if there is blocking in the empirical set AND the random data set, all empirical blocks will be compared to the paired random-summary block.
    blockSeq.x <- unique(x.summary$block) #pull the unique block ids in x.summary
    final_out <-NULL #create the empty final output that we will bind loop items to.
    for(j in blockSeq.x){ #use a for-loop to subset x.summary and y.summary by each block id
      block.x <- droplevels(subset(x.summary, block == j)) #subset x
      block.y <- droplevels(subset(y.summary, block == j)) #subset y. Note that blocks in x and y must be identical. 
      #make them matrices for x & y summaries with equal dimensions
      emp.matrix<-as.matrix(block.x[,4:(3+nrow(block.x))]) #create the matrix detailing observed contacts. Note that this matrix excludes id, totalDegree, and totalContactDuration columns
      rownames(emp.matrix)<-block.x$id #rename the rows of the matrix the ids
      
      randCols<-which(colnames(block.y)%in%colnames(emp.matrix)) #pull the colnumbers of columns in block.y that appear in the empirical matrix as well.
      randRows<-which(block.y$id%in%block.x$id) #pull the ids in block.y that also exist in the block.x set 
      y.redac<-droplevels(block.y[randRows, c(1,randCols)]) #reduce y.summary to only include individuals present in the x.summary (note that id column is also appended here)
      
      rand.matrix<-matrix(ncol = ncol(emp.matrix), nrow = nrow(emp.matrix)) #create the empty random-set matrix of appropriate length.
      rownames(rand.matrix)<-block.x$id #name the rows appropriately
      colnames(rand.matrix)<-colnames(emp.matrix) #name the columns appropriately
      for(i in 1:ncol(rand.matrix)){ #loop through the randRows and randCols values to add data to the rand.matrix
        if(is.na(match(colnames(rand.matrix)[i], colnames(y.redac))) == TRUE){next} #if the column was not found in y.redac, the loop moves on
        rand.matrix[match(y.redac$id, rownames(rand.matrix)),i]<-y.redac[,match(colnames(rand.matrix)[i], colnames(y.redac))] #pull the values from y.redac into rand.matrix
      }
      
      emp.mean<-mean(emp.matrix, na.rm = TRUE) #calculate the mean number of contacts in the empirical set
      rand.mean <- mean(rand.matrix, na.rm = TRUE) #calculate the mean number of contacts in the random set
      
      #the ape::mantel.test function below DOES NOT accept NAs in the input data frames, so have to replace NAs with 0
      emp.matrix[which(is.na(emp.matrix) == TRUE)]<-0 
      rand.matrix[which(is.na(rand.matrix) == TRUE)]<-0
      
      assign("last.warning", NULL, envir = baseenv()) #clears any warnings
      
      test.tryCatch<-tryCatch.W.E(ape::mantel.test(emp.matrix, rand.matrix, nperm = numPermutations, alternative = alternative.hyp)) #perform the mantel test
      output.test<-test.tryCatch$value #pull the value
      warn1<-test.tryCatch$warning$message #used to identify when the mantel.test produces any errors.
      
      if(length(warn1) == 0){
        warn1 = ""
      }
      
      summaryFrame <- data.frame(method = "Mantel's permutation test for similarity of two matrices", z.val = output.test[1], alternative.hyp = output.test[3], nperm = numPermutations, 
                                 p.val = output.test[2], x.mean = emp.mean, y.mean = rand.mean, block = unique(block.x$block), block.start = unique(block.x$block.start), 
                                 block.end = unique(block.x$block.end) , warning = warn1, stringsAsFactors = TRUE)
      
      bindlist<-list(final_out, summaryFrame) #create the list describing what should be bound together
      final_out <- data.frame(data.table::rbindlist(bindlist), stringsAsFactors = TRUE) #rbind summaryFrame to final function output.
    }
  }
  
  return(final_out)
}
