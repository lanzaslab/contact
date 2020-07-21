#' Compare Observed Contacts to a Random Distribution Using Chi-Square GoF
#' 
#' This function is used to determine if tracked individuals in an 
#'    empirical dataset had more or fewer contacts with other tracked 
#'    individuals/specified locations than would be expected at random. The
#'    function works by comparing an empirical contact distribution (generated 
#'    using x.summary and x.potential) to a NULL distribution (generated using 
#'    y.summary and y.potential) using a X-square goodness-of-fit test. Note
#'    that here, the NULL hypothesis is that empirical data are consistent with
#'    the NULL distribution, and the alternative hypothesis is that the data
#'    are NOT consistent. This function SHOULD NOT be used to compare two 
#'    empirical networks using Chi-squared tests, as the function assumes 
#'    x.summary and y.summary represent observed and expected values, 
#'    respectively. Please note that this is a function of convience that is 
#'    essentially a wrapper for the chisq.test function, that allows users to 
#'    easily compare contact networks created using our pipeline of contact:: 
#'    functions.
#'    
#' This function was inspired by the methods described by Spiegel et al. 2016. 
#'    They determined individuals to be expressing social behavior when nodes 
#'    had greater degree values than would be expected at random, with 
#'    randomized contact networks derived from movement paths randomized 
#'    according to their novel methodology (that can be implemented using our 
#'    randomizePaths function). Here, users can also identify when more or 
#'    fewer contacts (demonstrated by the sign of values in the "difference" 
#'    column in the output) with specific individuals than would be expected 
#'    at random, given a pre-determined p-value threshold. Such relationships 
#'    suggest social affinities or aversions, respectively, may exist between 
#'    specific individuals.
#'    
#' Note:The default tested column (i.e., categorical data column from which 
#'    data is drawn to be compared to randomized sets herein) is "id." This 
#'    means that contacts involving each individual (defined by a unique "id") 
#'    will be compared to randomized sets. Users may not use any data column 
#'    for analysis other than "id." If users want to use another categorical 
#'    data column in analyses rather than "id," we recommend re-processing 
#'    data (starting from the dist.all/distToArea functions), while specifying 
#'    this new data as an "id." For example, users may annotate an illness 
#'    status column to the empirical input, wherein they describe if the 
#'    tracked individual displayed gastrointestinal ("gastr"), respiratory 
#'    ("respr"), both ("both"), illness symptoms, or were consistently healthy 
#'    ("hel") over the course of the tracking period. Users could set this 
#'    information as the "id," and carry it forward as such through the 
#'    data-processing pipeline. Ultimately, they could determine if each of 
#'    these disease states affected contact rates, relative to what would be 
#'    expected at random.    
#'    
#' Take care to ensure that the same shuffle.type is denoted as was originally 
#'    used to randomize individuals' locations (assuming the randomizePaths 
#'    function was used to do so). This is important for two reasons: 1.) If 
#'    there was no y.potential input, the function assumes that x.potential is 
#'    relevant to the random set as well. This is a completely fair assumption 
#'    when importBlocks == FALSE or when the shuffleUnit == 0. In cases when 
#'    the shuffle.type is 1 or 2, however, this assumption can lead to 
#'    erroneous results and/or errors in the function. 2.) In the 
#'    randomizePaths function, setting shuffle.type == 2 produces only 1 
#'    shuffle.unit's worth of data (e.g., 1 day), rather than a dataset with 
#'    the same length of x. As such, there may be a different number of blocks 
#'    in y compared to x. Here we assume that the mean randomized durations 
#'    per block in y.summary and y.potential, are representative of mean 
#'    randomized durations per block across each shuffle unit (e.g., day 1 is 
#'    represntative of day 3, etc.).
#'    
#' Finally, if X-square expected values will be very small, 
#'    approximations of p may not be correct (and in fact, all estimates will 
#'    be poor). It may be best to weight these tests differently. In the event 
#'    that this is the case, \code{\link{contactCompare_binom}} may be used to 
#'    obtain more-accurate estimates.
#'
#' @param x.summary List or single-data frame output from the summarizeContacts
#'    function refering to the empirical data. Note that if x.summary is a list
#'    of data frames, only the first data frame will be used in the function.
#' @param y.summary List or single-data frame output from the summarizeContacts
#'    function refering to the randomized data (i.e., NULL model 
#'    contact-network edge weights). Note that if y.summary is a list
#'    of data frames, only the first data frame will be used in the function.
#' @param x.potential List or single-data frame output from the 
#'    potentialDurations function refering to the empirical data. Note that if 
#'    x.potential is a list of data frames, potential contact durations used in
#'    the function will be determined by averaging those reported in each list 
#'    entry. 
#' @param y.potential List or single-data frame output from the 
#'    potentialDurations function refering to the randomized data. Note that if 
#'    y.potential is a list of data frames, potential contact durations used in
#'    the function will be determined by averaging those reported in each list 
#'    entry. If NULL, reverts to x.potential. Defaults to NULL.
#' @param importBlocks Logical. If true, each block in x.summary will be 
#'    analyzed separately. Defaults to FALSE. Note that the "block" column must
#'    exist in .summary AND .potential objects, and values must be identical 
#'    (i.e., if block 100 exists in x inputs, it must also exist in y inputs), 
#'    otherwise an error will be returned.
#' @param shuffle.type Integer. Describes which shuffle.type (from the 
#'    randomizePaths function) was used to randomize the y.summary data 
#'    set(s). Takes the values "0," "1," or "2." This is important because 
#'    there are different assumptions associated with each shuffle.type.
#' @param pairContacts Logical. If TRUE individual id columns from x.summary 
#'    and y.summary inputs will be included in analyses. Defaults to TRUE.
#' @param totalContacts Logical. If TRUE totalDegree and totalContactDurations
#'    columns from x.summary and y.summary inputs will be included in analyses.
#'    Defaults to TRUE.
#' @param popLevelOutput Logical. If TRUE a secondary output describing 
#'    population-level comparisons will be appended to the standard, 
#'    individual-level function output.
#' @param parallel Logical. If TRUE, sub-functions within the summarizeContacts
#'    wrapper will be parallelized. Note that the only sub-function 
#'    parallelized here is called ONLY when importBlocks == TRUE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to half of the maximum number of cores 
#'    available (i.e., (parallel::detectCores()/2)).
#' @param ... Other arguments to be passed to the chisq.test function.
#' @keywords network-analysis social-network
#' @return Output format is dependent on \code{popLevelOutput} value.
#' 
#'    If \code{popLevelOut} == FALSE output will be a single two data frame 
#'    containing individual-level pairwise analyses of node degree, total 
#'    edge weight (i.e., the sum of all observed contacts involving each 
#'    individual), and specific dyad weights (e.g., contacts between 
#'    individuals 1 and 2). The data frame contains the following columns: 
#'    
#'    \item{id}{the id of the specific individual.}
#'    \item{metric}{designation of what is being compared (e.g., totalDegree, 
#'    totalContactDurations, individual 2, etc.). Content will 
#'    change depending on which data frame is being observed.}
#'    \item{method}{Statistical test used to determine significance.}
#'    \item{X.squared}{Test statistic associated with the comparison.}
#'    \item{p.val}{p.values associated with each comparison.}
#'    \item{df}{Degrees of freedom associated with the statistical test.}
#'    \item{contactDurations.x}{Describes the number of observed events
#'    in x.summary.}
#'    \item{contactDurations.y}{Describes the number of observed events in 
#'    y.summary.}
#'    \item{noContactDurations.x}{Describes the number of empirical events that
#'    were not observed given the total number of potential events in 
#'    x.potential.}
#'    \item{noContactDurations.y}{Describes the number of random events that
#'    were not observed given the total number of potential events in 
#'    y.potential.}
#'    \item{difference}{The absolute value given by subtracting 
#'    contactDurations.y from contactDurations.x.}
#'    \item{warning}{Denotes if any specific warning occurred during analysis.}
#'    \item{block.x}{Denotes the specific time block from x.(Only if 
#'    \code{importBlocks} == TRUE)}
#'    \item{block.start.x}{Denotes the specific timepoint at the beginning of 
#'    each time block. (Only if \code{importBlocks} == TRUE)}
#'    \item{block.end.x}{Denotes the specific timepoint at the end of each time
#'    block. (Only if \code{importBlocks} == TRUE)}
#'    \item{block.y}{Denotes the specific time block from y.(Only if 
#'    \code{importBlocks} == TRUE)}
#'    \item{block.start.y}{Denotes the specific timepoint at the beginning of 
#'    each time block. (Only if \code{importBlocks} == TRUE)}
#'    \item{block.end.y}{Denotes the specific timepoint at the end of each time
#'    block. (Only if \code{importBlocks} == TRUE)}
#'    
#'    If \code{popLevelOutput} == TRUE, output will be a list of two data 
#'    frames: The one described above, and second describing the 
#'    population-level comparisons. Columns in each data frame are identical.
#'    
#' @references Agresti, A. 2007. An introduction to categorical data analysis, 
#'    2nd ed. New York: John Wiley & Sons. 38.
#' 
#'    Farine, D.R., 2017. A guide to null models for animal social 
#'    network analysis. Methods in Ecology and Evolution 8:1309-1320.
#'    https://doi.org/10.1111/2041-210X.12772.
#'    
#'    Spiegel, O., Leu, S.T., Sih, A., and C.M. Bull. 2016. Socially 
#'    interacting or indifferent neighbors? Randomization of movement paths to 
#'    tease apart social preference and spatial constraints. Methods in Ecology
#'    and Evolution 7:971-979. https://doi.org/10.1111/2041-210X.12553.
#'  
#' @import foreach  
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
#' emp.potential <- potentialDurations(calves.dist, blocking = TRUE, 
#'                                     blockUnit = "hours", blockLength = 1, 
#'                                     distFunction = "dist2All_df") 
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
#' rand.potential <- potentialDurations(calves.dist.rand, blocking = TRUE, 
#'                                      blockUnit = "hours", blockLength = 1, 
#'                                      distFunction = "dist2All_df") 
#' 
#' 
#' contactCompare_chisq(x.summary = emp.summary, y.summary = rand.summary, 
#'                      x.potential = emp.potential, y.potential = rand.potential,
#'                      importBlocks = FALSE, shuffle.type = 0, 
#'                      popLevelOut = TRUE, parallel = FALSE) #no blocking
#' 
#' contactCompare_chisq(x.summary = emp.summary, y.summary = rand.summary, 
#'                      x.potential = emp.potential, y.potential = rand.potential,
#'                      importBlocks = TRUE, shuffle.type = 0, 
#'                      popLevelOut = TRUE, parallel = FALSE) #blocking
#'    }

contactCompare_chisq<-function(x.summary, y.summary, x.potential, y.potential = NULL, importBlocks = FALSE, shuffle.type = 1, pairContacts = TRUE, totalContacts = TRUE, popLevelOutput = FALSE, parallel = FALSE, nCores = (parallel::detectCores()/2), ...){
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  block <- NULL
  id<-NULL 
  output<-NULL
  j <- NULL
  metric <-NULL
  block.x<-NULL
  
  
  chisq.forLoop<-function(x, empirical, randomized, emp.potential, rand.potential, ...){ #I hate that I have to do this in a for-loop, but I couldn't get the apply functions to work. Note "x" here is not x.summary or x.potential. Those are represented by the empirical and emp.potential arguments, respectively.
    
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

    output<-NULL #create empty object
    
    indivSummaryTest<- ifelse(length(grep("contactDuration_Indiv", colnames(empirical))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
    
    #browser()
    
    #vectorize the input
    ids<-x[,1]
    cols<-x[,2]
    
    #pull the id colnames (as they relate to cols values) for empirical and random sets
    emp.IdCols<-c(colnames(empirical)[1:3], substring(names(x.summary[,4:max(grep("contactDuration_", names(empirical)))]), 22)) #pull the empirical colnames for the 1st 3 columns, then the unique IDs for later columns
    rand.IdCols<-c(colnames(randomized)[1:3], substring(names(y.summary[,4:max(grep("contactDuration_", names(randomized)))]), 22)) #pull the empirical colnames for the 1st 3 columns, then the unique IDs for later columns
    
    for(i in 1:nrow(x)){
      
      if(match(cols[i], emp.IdCols) >= 4){ #empirical[,1:3] do not represent contacts derived from singular columns in the contact::dist2... function outputs.
        
        if(indivSummaryTest == TRUE){ #if the summarizeContacts output represented contacts with individuals
          
          empDurations<- empirical[which(empirical$id == unlist(unname(ids[i]))), match(paste("contactDuration_Indiv",unlist(unname(cols[i])), sep = ""),colnames(empirical))] #pull the value of a given column for a specific id 
          
        }
        
        if(indivSummaryTest == FALSE){ #if the summarizeContacts output represented contacts with individuals
          
          empDurations<- empirical[which(empirical$id == unlist(unname(ids[i]))), match(paste("contactDuration_Area_",unlist(unname(cols[i])), sep = ""),colnames(empirical))] #pull the value of a given column for a specific id 
          
        }
      }else{ #i.e., x[2] == empirical[,2] or empirical[,3]
        
        empDurations<- empirical[which(empirical$id == unlist(unname(ids[i]))), grep(unlist(unname(cols[i])),colnames(empirical))] #pull the value of a given column for a specific id 
        
      }
      
      
      empDurations<-empDurations[is.na(empDurations) == FALSE] #remove the NAs
      summaryFrame<-NULL
      
      if(length(empDurations) > 0){ #if there is no entry OR is.na == TRUE, nothing will happen
        
        if(match(cols[i], emp.IdCols) >= 4){ #empirical[,2:3] do not represent contacts derived from singular columns in the contact::dist2... function outputs.
          
          maxDurations<- emp.potential[which(emp.potential$id == unlist(unname(ids[i]))), match(paste("potenContactDurations_", unlist(unname(cols[i])), sep = ""), names(emp.potential))]
          
        }else{ #i.e., x[2] == 2 or 3
          if(match(cols[i], emp.IdCols) == 2){ #if cols = totalDegree, the maximum degree possible would be the number of individuals observed during the time period (i.e., potenDegree in emp.potential)
            maxDurations<- emp.potential[which(emp.potential$id == unlist(unname(ids[i]))), match("potenDegree", names(emp.potential))]
          }
          if(match(cols[i], emp.IdCols) == 3){ #if cols = totalcontactDurations, the maximum possible number of contacts would be would be the sum of all individuals observed at each time step during the time period, excluding individual i (i.e., potenTotalContactDurations in emp.potential).
            maxDurations<- emp.potential[which(emp.potential$id == unlist(unname(ids[i]))), match("potenTotalContactDurations", names(emp.potential))]
          }
        }
        
        if(maxDurations == 0 | is.infinite(maxDurations) == TRUE){ #if there was no potential for contacts to occur, the loop moves on
          next
        }
        
        if(is.na(match(cols[i], rand.IdCols)) == TRUE ){ #if i is trying to reference a specific object not present in the random set, rand duration will equal "0."
          randDurations <- 0
        }else{ # if the relavent column DOES exist in randomized
          
          
          if(match(cols[i], rand.IdCols) >= 4){ #randomized[,2:3] do not represent contacts derived from singular columns in the contact::dist2... function outputs.
            
            
            if(indivSummaryTest == TRUE){ #if the summarizeContacts output represented contacts with individuals
              
              randDurations<- randomized[which(randomized$id == unlist(unname(ids[i]))), 
                                         match(paste("contactDuration_Indiv",unlist(unname(cols[i])), sep = ""),colnames(randomized))] #pull the value of a given column for a specific id in a specific block
              
            }
            
            if(indivSummaryTest == FALSE){ #if the summarizeContacts output represented contacts with individuals
              
              randDurations<- randomized[which(randomized$id == unlist(unname(ids[i]))), 
                                         match(paste("contactDuration_Area_",unlist(unname(cols[i])), sep = ""),colnames(randomized))] #pull the value of a given column for a specific id in a specific block
              
            }
          }else{ #i.e., x[2] == randomized[,2] or randomized[,3]
            
            randDurations<- randomized[which(randomized$id == unlist(unname(ids[i]))), 
                                       grep(unlist(unname(cols[i])),colnames(randomized))] #pull the value of a given column for a specific id in a specific block
            
          }
          
        }
        
        randDurations<-randDurations[is.na(randDurations) == FALSE] #remove the NAs
        if(length(randDurations) == 0){ #if removing the NAs removed the only observations, then randDurations reverts to 0
          randDurations <- 0
        }
        
        
        if(is.na(match(cols[i], rand.IdCols)) == TRUE | match(cols[i], rand.IdCols) >= 4){ #randomized[,1:3] do not represent contacts derived from singular columns in the contact::dist2... function outputs.
          
          maxDurations.rand<- rand.potential[which(rand.potential$id == unlist(unname(ids[i]))), match(paste("potenContactDurations_", unlist(unname(cols[i])), sep = ""), names(rand.potential))]
          
        }else{ #i.e., x[2] == 2 or 3
          if(grep(unlist(unname(cols[i])),colnames(randomized))[1] == 2){ #if cols = totalDegree, the maximum degree possible would be the number of individuals observed during the time period (i.e., potenDegree in rand.potential)
            maxDurations.rand<- rand.potential[which(rand.potential$id == unlist(unname(ids[i]))), match("potenDegree", names(rand.potential))]
          }
          if(grep(unlist(unname(cols[i])),colnames(randomized))[1] == 3){ #if cols = totalcontactDurations, the maximum possible number of contacts would be would be the sum of all individuals observed at each time step during the time period, excluding individual i (i.e., potenTotalContactDurations in rand.potential).
            maxDurations.rand<- rand.potential[which(rand.potential$id == unlist(unname(ids[i]))), match("potenTotalContactDurations", names(rand.potential))]
          }
        }
        
        
        maxDurations.rand<-maxDurations.rand[is.na(maxDurations.rand) == FALSE] #remove the NAs
        if(length(maxDurations.rand) == 0){ #if removing the NAs removed the only observations, then maxDurations.rand reverts to 0
          maxDurations.rand <- 0
        }
        
        observedVec <- c(sum(empDurations), maxDurations-sum(empDurations)) #create a vector describing observed counts

        if(maxDurations.rand > 0){ #we must prevent dividing by zero here. 
          expectedProb <- c((sum(randDurations)/maxDurations.rand), ((maxDurations.rand-sum(randDurations))/maxDurations.rand)) #convert the observed random durations to probabilities.

        }else{ #if maxDurations.rand == 0, expectedProb will always be c(0,0) and therefore chisq.test will return an error
          
          expectedProb <- c(0,1) #here we force expected probability to be (0, 1) (i.e., 100% probability of no contact), which would be the case if the individual was not observed.
          
        }
        
        assign("last.warning", NULL, envir = baseenv()) #clears the warnings
        
        test<-tryCatch.W.E(stats::chisq.test(x = observedVec, p = expectedProb, ...))$value
        warn1<-tryCatch.W.E(stats::chisq.test(x = observedVec, p = expectedProb, ...))$warning$message
        
        if(length(warn1) == 0){
          warn1 = ""
        }
        
        summaryFrame <- data.frame(id = unlist(unname(ids[i])), metric = unlist(unname(cols[i])), method = unname(test[4]), X.squared = unname(test[1]), 
                                   df = unname(test[2]), p.val = unname(test[3]), empiricalContactDurations = sum(empDurations), 
                                   randContactDurations.mean = sum(randDurations), empiricalNoContactDurations = (maxDurations - sum(empDurations)), 
                                   randNoContactDurations.mean = (maxDurations.rand - sum(randDurations)), difference = abs((maxDurations - sum(empDurations)) - (maxDurations.rand - sum(randDurations))), 
                                   warning = warn1, stringsAsFactors = TRUE)

        colnames(summaryFrame)<-c("id", "metric", "method", "X.squared", "df", "p.val", "contactDurations.x", "contactDurations.y", "noContactDurations.x", 
                                  "noContactDurations.y", "difference", "warning") #for some reason the data.frame command above kept producing incorrect colNames.
        
        
        
        rownames(summaryFrame) <-1
        
        bindlist<-list(output, summaryFrame)
        output<- data.table::rbindlist(bindlist, fill = TRUE)
      }else{ #if there is no entry OR is.na == TRUE, nothing will happen
        next
      }
    }
    
    final.out<-data.frame(output, stringsAsFactors = TRUE)
    return(final.out)
  }
  
  summaryAgg.block1<-function(x,y){ #calculates the mean potential contacts by id and block. Using this apply function is faster than simply aggregating the data set by id and block
    sumTable<-y[which(y$block == unname(unlist(x[1]))),]
    
    if(nrow(sumTable) == 0){output <- NULL #if there's nothing in the subset, the function will not proceed any further.
    
    }else{
      
      blockStart<- unique(lubridate::as_datetime(sumTable$block.start)) #added 02/05/2019 - had to keep track of this new information ; updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
      blockEnd<- unique(lubridate::as_datetime(sumTable$block.end)) #added 02/05/2019 - had to keep track of this new information ;  updated 06/02/2019 - converted the factor data to POSIXct format in order to avoid a "length is too large for hashing" error.
      sumTable.redac<-sumTable[,-c(match("id", names(sumTable)), match("block", names(sumTable)), match("block.start", names(sumTable)), match("block.end", names(sumTable)))]  #Remove the columns that cannot/shoud not be averaged.
      output<-stats::aggregate(sumTable.redac, list(id = sumTable$id), mean) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set.
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
      output<-stats::aggregate(sumTable.redac, list(id = sumTable$id), sum) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set.
      output$potenDegree <- length(grep("potenContactDurations_", colnames(output))) #redefine the potential degree as the total number potentially-contactable entities present at in the data. Note that this is only accurate if y is derived from contact::dist2Area_df output. If y was instead derived from dist2All_df, then these values must have 1 subtracted from them. If this is the case, the subtraction will take place outside of this function. 
      
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
  
  if(importBlocks == TRUE){ #if importBlocks == TRUE, but there is no block information in y.summary, the function returns a warning, but proceeds as if the y input is relevant to all columns.
    
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
      }
      
    }else{ #if importBlocks == FALSE we only average the data by id.
      
      if(length(x.potentialAgg$block) > 0){ #Alternatively, if there are blocked time sets in the empirical input, we need to dictate how to handle them.
        
        if(importBlocks == FALSE){ #if there ARE blocks in x.potential, but importBlocks == FALSE, we have to recreate x.potential without the bloking information.
          
          for(i in 1:length(x.potential)){ #add replicateID to the input. This was not required when processing x.summary because only the first entry of the list is returned.
            x.potential[[i]]$replicateID = i}
          
          x.potentialAgg<- data.frame(data.table::rbindlist(x.potential, fill = TRUE), stringsAsFactors = TRUE) #remake x.potentialAgg
          
          #now we have an object that represents average contact durations in each block
          repSeq<-unique(x.potentialAgg$replicateID)
          sumTab <- apply(data.frame(repSeq, stringsAsFactors = TRUE), 1, summaryAgg.block2, y = x.potentialAgg) #sum across columns by id and replicateID
          x.potentialAgg <- data.frame(data.table::rbindlist(sumTab, fill = TRUE), stringsAsFactors = TRUE) #We keep the same name for simplicity's sake below. 
          ##and now we have an object contained summed averages by id, replicateID. However, we may need to retractively adjust potenDegree values (see summaryAgg.block2)
          indivSummaryTest.x<- ifelse(length(grep("contactDuration_Indiv", colnames(x.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
          
          if(indivSummaryTest.x == TRUE){ #if x.summary was derived from dist2All (Note that this assumes x.summary and x.potential are derived from the same distance function)
            x.potentialAgg$potenDegree <- x.potentialAgg$potenDegree - 1 #the potential degree is the maximum number of edges that can extend from any one node in the network
          }
          
        }
        
      }
      ##now we take the average across all replicateIDs
      x.potential<-stats::aggregate(x.potentialAgg[,-match("id", colnames(x.potentialAgg))], list(id = x.potentialAgg$id), mean) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 
    }
    
  }else{ #if x.potential is only a single data frame
    
    if(importBlocks == TRUE){
      
      if(is.na(match("block", colnames(x.potential))) == TRUE){ #if there is NOT a "block" column present, the function will return an error. 
        
        stop("importBlocks set to TRUE, but no block column exists in x.potential")
        
      }
      
    }else{ #if importBlocks == FALSE we must sum the data by id if blocks exist in the input data. If blocks do not exist in the data, nothing needs to be done.
      
      if(length(x.potential$block) > 0){ #Alternatively, if there are blocked time sets in the empirical input, we need to dictate how to handle them.
        
        x.potential<-stats::aggregate(x.potential[,- c(match("id", colnames(x.potential)),match("block", colnames(x.potential)):ncol(x.potential))], list(id = x.potential$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 
        ### We must recognize, however, that the "totalDegree" column cannot be summed like other columns can (it is the number of individuals/fixed areas that each unique individual was observed in contact with.)
        ###below we recalculate the totalDegree for each individual.Note that this will be a different value depending on whether or not x.summary was based on dist2All or dist2Area outputs.
        indivSummaryTest.x<- ifelse(length(grep("contactDuration_Indiv", colnames(x.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
        
        if(indivSummaryTest.x == TRUE){ #if x.summary was derived from dist2All (Note that this assumes x.summary and x.potential are derived from the same distance function)
          x.potential$potenDegree <- length(unique(x.potential$id)) - 1 #the potential degree is the maximum number of edges that can extend from any one node in the network
        }
        if(indivSummaryTest.x == FALSE){ #if x.summary was derived from dist2Area
          x.potential$potenDegree <- length(grep("potentialContactDurations_", colnames(x.potential))) #the potential degree is the maximum number of edges that can extend from any one node in the network
        }        
      }

    }
    
  }
  
  
  #we do the same thing for y.potential as we did for x.potential above.
  
  if(is.data.frame(y.potential) == FALSE & is.list(y.potential) == TRUE){ #Because R treats dataframes as lists, we assess here if the input was only a single data frame or a list of data frames.
    
    y.potentialAgg<- data.frame(data.table::rbindlist(y.potential, fill = TRUE), stringsAsFactors = TRUE) #bind the lists together. 
    
    if(importBlocks == TRUE){
      
      if(is.na(match("block", colnames(y.potentialAgg))) == TRUE){ #if there is NOT a "block" column present, the function will return a warning and proceed as if importBlocks == FALSE. 
        
        warning("importBlocks set to TRUE, but no block column exists in y.potential. Proceeding as if y.potential values are stable across time and relevant to EVERY block.", immediate. = TRUE)
        
        #return y.potential as if importBlocks == FALSE and no block column exists in y.potential then repeat it for each block
        y.potential<-stats::aggregate(y.potentialAgg[,-match("id", colnames(y.potentialAgg))], list(id = y.potentialAgg$id), mean) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 

        y.potentialBlock<-NULL #create an empty object to contain new block information
        
        for(i in unique(x.potential$block)){
          
          y.potential$block <- i #add a block column containing only the i value to y.potential
          y.potentialBlock <- data.table::rbindlist(list(y.potentialBlock, y.potential), fill = TRUE) #rbind the data frames together
          
        }
        
        y.potential <- data.frame(y.potentialBlock, stringsAsFactors = TRUE) #redefine y.potential as the object containing block information

      }else{ #If there IS a block column, we average observations by id AND block
        
        blockSeq1<-unique(y.potentialAgg$block) #there may be differences in blocks containing contact info between the empirical and randomized data sets
        meanTab <- apply(data.frame(blockSeq1, stringsAsFactors = TRUE), 1, summaryAgg.block1, y = y.potentialAgg) #Note that block information MUST be included in the y.potentialAgg input
        y.potential <- data.frame(data.table::rbindlist(meanTab, fill = TRUE), stringsAsFactors = TRUE) #We keep the same name for simplicity's sake below. 

      }
      
    }else{ #if importBlocks == FALSE we only average the data by id.
      
      if(length(y.potentialAgg$block) > 0){ #Alternatively, #if there ARE blocks in y.potential, but importBlocks == FALSE, we have to recreate y.potential without the blocking information.
        
          for(i in 1:length(y.potential)){ #add replicateID to the input. This was not required when processing x.summary because only the first entry of the list is returned.
            y.potential[[i]]$replicateID = i}
        
          y.potentialAgg<- data.frame(data.table::rbindlist(y.potential, fill = TRUE), stringsAsFactors = TRUE) #remake y.potentialAgg
          repSeq<-unique(y.potentialAgg$replicateID) #make a vector of the unique replicateIDs
          sumTab <- apply(data.frame(repSeq, stringsAsFactors = TRUE), 1, summaryAgg.block2, y = y.potentialAgg) #sum across columns by id and replicateID
          y.potentialAgg <- data.frame(data.table::rbindlist(sumTab, fill = TRUE), stringsAsFactors = TRUE) #We keep the same name for simplicity's sake below. 
          ##and now we have an object contained summed averages by id, replicateID. However, we may need to retractively adjust potenDegree values (see summaryAgg.block2)
          indivSummaryTest.y<- ifelse(length(grep("contactDuration_Indiv", colnames(y.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
          
          if(indivSummaryTest.y == TRUE){ #if x.summary was derived from dist2All (Note that this assumes x.summary and x.potential are derived from the same distance function)
            y.potentialAgg$potenDegree <- y.potentialAgg$potenDegree - 1 #the potential degree is the maximum number of edges that can extend from any one node in the network
          }
          
      }
      ##now we take the average across all replicateIDs
      y.potential<-stats::aggregate(y.potentialAgg[,-match("id", colnames(y.potentialAgg))], list(id = y.potentialAgg$id), mean) #this not only calculates the mean of each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 

    }
    
  }else{ #if y.potential is only a single data frame
    
    if(importBlocks == TRUE){
      
      if(is.na(match("block", colnames(y.potential))) == TRUE){ #if there is NOT a "block" column present.
        
        warning("importBlocks set to TRUE, but no block column exists in y.potential. Proceeding as if y.potential values are stable across time and relevant to EVERY block.")
        
        #return y.potential as if importBlocks == FALSE and no block column exists in y.potential then repeat it for each block
        
        y.potentialBlock<-NULL #create an empty object to contain new block information
        
        for(i in unique(x.potential$block)){
          
          y.potential$block <- i #add a block column containing only the i value to y.potential
          y.potentialBlock <- data.table::rbindlist(list(y.potentialBlock, y.potential), fill = TRUE) #rbind the data frames together
          
        }
        
        y.potential <- data.frame(y.potentialBlock, stringsAsFactors = TRUE) #redefine y.potential as the object containing block information
      }
      
    }else{ #if importBlocks == FALSE we must sum the data by id if blocks exist in the input data. If blocks do not exist in the data, nothing needs to be done.
      
      if(length(y.potential$block) > 0){ #Alternatively, if there are blocked time sets in the empirical input, we need to dictate how to handle them.
        
        y.potential<-stats::aggregate(y.potential[,-c(match("id", colnames(y.potential)),match("block", colnames(y.potential)):ncol(y.potential))], list(id = y.potential$id), sum) #this not only sums each column by id, but also adds the "id" column back into the data set. #We keep the same name for simplicity's sake below. 
        ### We must recognize, however, that the "totalDegree" column cannot be summed like other columns can (it is the number of individuals/fixed areas that each unique individual was observed in contact with.)
        ###below we recalculate the totalDegree for each individual.Note that this will be a different value depending on whether or not x.summary was based on dist2All or dist2Area outputs.
        indivSummaryTest.y<- ifelse(length(grep("contactDuration_Indiv", colnames(y.summary))) >0, TRUE, FALSE) #the summarizeContacts function can either represent contacts with individuals OR fixed areas. We need to confirm which it is here. 
        
        if(indivSummaryTest.y == TRUE){ #if y.summary was derived from dist2All (Note that this assumes y.summary and y.potential are derived from the same distance function)
          y.potential$potenDegree <- length(unique(y.potential$id)) - 1 #the potential degree is the maximum number of edges that can extend from any one node in the network
        }
        if(indivSummaryTest.y == FALSE){ #if y.summary was derived from dist2Area
          y.potential$potenDegree <- length(grep("potentialContactDurations_", colnames(y.potential))) #the potential degree is the maximum number of edges that can extend from any one node in the network
        }     
      }
    }
  }
  
  #OK. Now we have x.summary, y.summary, x.potential, and y.potential objects that are gauranteed to work with the chisq.forLoop function.
  
  #define the loop objects.
  
  idSeq <- unique(c(x.summary$id, y.summary$id)) #pulls the unique ids for each individual
  empiricalColNames <- c("totalDegree", "totalContactDurations", substring(names(x.summary[,4:max(grep("contactDuration_", names(x.summary)))]), 22)) #identifies which field each column relates to in x.summary
  loopFrame <- expand.grid(idSeq, empiricalColNames, stringsAsFactors = TRUE) #create the loopFrame to run through the chisqLoop function by combining the 2 previously-created vectors
  
  #based on the pairContacts and totalContacts arguments, users can choose to exclude certain comparisons to speed up processing
  if(pairContacts == FALSE){
    loopFrame<- droplevels(loopFrame[c(which(loopFrame$Var2 == "totalDegree"), which(loopFrame$Var2 == "totalContactDurations")),]) #remove othe observations
  }
  
  if(totalContacts == FALSE){
    loopFrame<- droplevels(loopFrame[-c(which(loopFrame$Var2 == "totalDegree"), which(loopFrame$Var2 == "totalContactDurations")),]) #keep only these observations
  }
  
  if(pairContacts == FALSE & totalContacts == FALSE){
    
    stop("All observed contacts have been removed prior to analyses because both pairContacts AND totalContacts are set to FALSE.")
    
  }
  
  #ensure that warnings in the chisq.forLoop function are able to be recorded in output
  ##before executing the tests we must ensure that warnings for sub-level functions will not be silenced (so that we can record them as they occur, allowing users to pinpoint what part of their output might be erroneous).
  oldw <- getOption("warn") #pull the current warn setting 
  options(warn = 1) #make warnings appear as they occur, rather than be stored until the top-level function returns
  on.exit(options(warn = oldw)) #ensure that when the function ends, users' options are reset.
  
  #run the chisq.forLoop function
  
  if(length(x.summary$block) > 0){ #importBlocks == TRUE and blocks exist in the inputs (note: if importBlocks == FALSE, code above ensures that no "block" column will exist at this point.)
    
    blockSeq <- unique(x.summary$block) #pull unique blocks from x.summary
    
    if(parallel == TRUE){
      
      cl <- parallel::makeCluster(nCores) #set up cluster
      doParallel::registerDoParallel(cl) #register the cluster
      on.exit(parallel::stopCluster(cl), add = TRUE) #ensure the cluster is closed out when the function ends. This is added to the previous on.exit call that resets user's warning settings.
      
      chisqOut.list <- foreach::foreach(j = blockSeq) %dopar% {
        x.summaryBlock <- droplevels(subset(x.summary, block == j)) #subset x.summary by block
        x.potentialBlock <- droplevels(subset(x.potential, block == j)) #subset x.potential by block

        if(shuffle.type == 2){ #recall that shuffle.type 2 (from the randomizeLocations function) produces only 1 shuffle.unit's worth of data, rather than a dataset with the same length of x. As such, there may be a different number of blocks in y compared to x. Here we assume that the mean randomized durations per block, are representative of mean randomized durations per block across each shuffle unit (e.g., day)
          blockSeq.integ <- as.integer(as.character(blockSeq)) #ensure that the blocks in blockSeq are integers so that they can be put through mathematical operations
          blockSeqY.integ <- as.integer(as.character(unique(y.summary$block))) #do the same thing for the y input
          block.y = as.integer(as.character(j)) - (ceiling((as.integer(as.character(j)) - max(blockSeqY.integ)) /max(blockSeqY.integ))*max(blockSeqY.integ)) #identify what block in y to compare to block j in x
        }else{ #if shuffle.type is anything other than 2, block.y just equals j
          block.y = j
        }
        
        y.summaryBlock <- droplevels(subset(y.summary, block == block.y)) #subset y.summary by block
        y.potentialBlock <- droplevels(subset(y.potential, block == block.y)) #subset x.potential by block
        chisqOut<-chisq.forLoop(loopFrame, empirical = x.summaryBlock, randomized = y.summaryBlock, emp.potential = x.potentialBlock, rand.potential = y.potentialBlock, ...) #generate the individual-level chisq output.
        ##add block information to chisqOut (Note: we add the information for both x AND y even though the y information will be redundant unless shuffle.type == 2)
        chisqOut$block.x <- j
        chisqOut$block.start.x <- unique(x.summaryBlock$block.start)
        chisqOut$block.end.x <- unique(x.summaryBlock$block.end)
        chisqOut$block.y <- block.y
        chisqOut$block.start.y <- unique(y.summaryBlock$block.start)
        chisqOut$block.end.y <- unique(y.summaryBlock$block.end)
        
        return(chisqOut)
        
      }
      
    }else{ #if parallel == FALSE
      
      chisqOut.list <- foreach::foreach(j = blockSeq) %do% {
        
        x.summaryBlock <- droplevels(subset(x.summary, block == j)) #subset x.summary by block
        x.potentialBlock <- droplevels(subset(x.potential, block == j)) #subset x.potential by block
        
        if(shuffle.type == 2){ #recall that shuffle.type 2 (from the randomizeLocations function) produces only 1 shuffle.unit's worth of data, rather than a dataset with the same length of x. As such, there may be a different number of blocks in y compared to x. Here we assume that the mean randomized durations per block, are representative of mean randomized durations per block across each shuffle unit (e.g., day)
          blockSeq.integ <- as.integer(as.character(blockSeq)) #ensure that the blocks in blockSeq are integers so that they can be put through mathematical operations
          blockSeqY.integ <- as.integer(as.character(unique(y.summary$block))) #do the same thing for the y input
          block.y = as.integer(as.character(j)) - (ceiling((as.integer(as.character(j)) - max(blockSeqY.integ)) /max(blockSeqY.integ))*max(blockSeqY.integ)) #identify what block in y to compare to block j in x
        }else{ #if shuffle.type is anything other than 2, block.y just equals j
          block.y = j
        }
        
        y.summaryBlock <- droplevels(subset(y.summary, block == block.y)) #subset y.summary by block
        y.potentialBlock <- droplevels(subset(y.potential, block == block.y)) #subset x.potential by block
        chisqOut<-chisq.forLoop(loopFrame, empirical = x.summaryBlock, randomized = y.summaryBlock, emp.potential = x.potentialBlock, rand.potential = y.potentialBlock, ...) #generate the individual-level chisq output.
        ##add block information to chisqOut (Note: we add the information for both x AND y even though the y information will be redundant unless shuffle.type == 2)
        chisqOut$block.x <- j
        chisqOut$block.start.x <- unique(x.summaryBlock$block.start)
        chisqOut$block.end.x <- unique(x.summaryBlock$block.end)
        chisqOut$block.y <- block.y
        chisqOut$block.start.y <- unique(y.summaryBlock$block.start)
        chisqOut$block.end.y <- unique(y.summaryBlock$block.end)
        
        return(chisqOut)
        
      }
      
    }
    
    chisqOut <- data.frame(data.table::rbindlist(chisqOut.list, fill = TRUE), stringsAsFactors = TRUE) #bind the lists together
    
    if(popLevelOutput == TRUE){
      
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
      
      uniqueMetrics <- unique(chisqOut$metric) #pull out the unique metrics for the data set
      blocks <- unique(chisqOut$block.x) #pull out the unique empirical blocks in chisqOut
      
      processFrame<- expand.grid(uniqueMetrics, blocks, stringsAsFactors = TRUE) #create the frame to loop through below
      
      popLevelOut<-foreach::foreach(i = seq(from = 1, to = nrow(processFrame), by = 1)) %do% { #no need to make this parallel. Doing so would only slow things down.
        
        summaryFrame <- data.frame(matrix(ncol = 11, nrow = 1), stringsAsFactors = TRUE) #just in case there are no entries in summaryFrame, all entries in this summaryFrame will be NA
        summaryFrame[1,] <- NA
        
        popLevelMetric <- droplevels(subset(chisqOut, metric == processFrame[i,1] & block.x == processFrame[i,2])) # pull the unique value of interest
        
        if(nrow(popLevelMetric) > 0){
          
          observedVec <- c(sum(popLevelMetric[,match("contactDurations.x", colnames(chisqOut))]), sum(popLevelMetric[,match("noContactDurations.x", colnames(chisqOut))])) #create a vector describing total observed counts of empirical contacts and non-contacting timepoints
          maxDurations.x <- sum(c(popLevelMetric[,match("contactDurations.x", colnames(chisqOut))], popLevelMetric[,match("noContactDurations.x", colnames(chisqOut))])) #total number of random temporal sampling-windows observed. 
          expectedVec <- c(sum(popLevelMetric[,match("contactDurations.y", colnames(chisqOut))]), sum(popLevelMetric[,match("noContactDurations.y", colnames(chisqOut))])) #create a vector describing total observed counts of empirical contacts and non-contacting timepoints
          maxDurations.y <- sum(c(popLevelMetric[,match("contactDurations.y", colnames(chisqOut))], popLevelMetric[,match("noContactDurations.y", colnames(chisqOut))])) #total number of random temporal sampling-windows observed. This will be used to created the probability distribution for the NULL model
        
          maxDurations.y<-maxDurations.y[is.na(maxDurations.y) == FALSE] #remove the NAs
          if(length(maxDurations.y) == 0){ #if removing the NAs removed the only observations, then maxDurations.y reverts to 0
            maxDurations.y <- 0
          }
          
          if(maxDurations.y > 0){ #we must prevent dividing by zero here. 
            expectedProb <- c((expectedVec[1]/maxDurations.y), (expectedVec[2]/maxDurations.y)) #convert the observed random durations to probabilities.
            
          }else{ #if maxDurations.y == 0, expectedProb will always be c(0,0) and therefore chisq.test will return an error
            
            expectedProb <- c(0,1) #here we force expected probability to be (0, 1) (i.e., 100% probability of no contact), which would be the case if the individual was not observed.
            
          }
          
          assign("last.warning", NULL, envir = baseenv()) #clears previous warnings so that we may record any new warnings
        
          test<-tryCatch.W.E(stats::chisq.test(x = observedVec, p = expectedProb))$value #get the chisq values
          warn1<-tryCatch.W.E(stats::chisq.test(x = observedVec, p = expectedProb))$warning$message #get any warning message that the chisq may have triggered
        
          if(length(warn1) == 0){
            warn1 = ""
          }
        
          summaryFrame <- data.frame(metric = processFrame[i,1], method = unname(test[4]), X.squared = unname(test[1]), 
                                   df = unname(test[2]), p.val = unname(test[3]), empiricalContactDurations = observedVec[1], 
                                   randContactDurations.mean = expectedVec[1], empiricalNoContactDurations = observedVec[2], 
                                   randNoContactDurations.mean = expectedVec[2], difference = abs((maxDurations.x - observedVec[1]) - (maxDurations.y - expectedVec[1])), 
                                   warning = warn1, stringsAsFactors = TRUE)
        
        }
          colnames(summaryFrame)<-c("metric",  "method", "X.squared", "df", "p.val", "contactDurations.x", "contactDurations.y", "noContactDurations.x", 
                                  "noContactDurations.y", "difference", "warning") #for some reason the data.frame command above kept producing incorrect colNames.
        
          ##add block information to summaryFrame (Note: we add the information for both x AND y even though the y information will be redundant unless shuffle.type == 2)
          
          summaryFrame$block.x <- unique(droplevels(processFrame[i,2]))
          summaryFrame$block.start.x <- unique(chisqOut$block.start.x[which(chisqOut$block.x == processFrame[i,2])])
          summaryFrame$block.end.x <- unique(chisqOut$block.end.x[which(chisqOut$block.x == processFrame[i,2])])
          summaryFrame$block.y <- unique(chisqOut$block.y[which(chisqOut$block.x == processFrame[i,2])])
          summaryFrame$block.start.y <- unique(chisqOut$block.start.y[which(chisqOut$block.x == processFrame[i,2])])
          summaryFrame$block.end.y <- unique(chisqOut$block.end.y[which(chisqOut$block.x == processFrame[i,2])])
        
          rownames(summaryFrame) <-1
          
        return(summaryFrame)
        
        }
      
      popLevelOut<- data.frame(data.table::rbindlist(popLevelOut), stringsAsFactors = TRUE) #bind these frames together
      
      chisqOut <- list(chisqOut, popLevelOut) #group the outputs together
      names(chisqOut) <- c("individualLevel", "populationLevel") #rename the objects in the list
    }
    
  }else{ #if importBlocks == FALSE
    
    chisqOut<-chisq.forLoop(loopFrame, empirical = x.summary, randomized = y.summary, emp.potential = x.potential, rand.potential = y.potential, ...) #generate the individual-level chisq output.
    
    if(popLevelOutput == TRUE){
      
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
      
      uniqueMetrics <- unique(chisqOut$metric) #pull out the unique metrics for the data set
      
      popLevelOut<-foreach::foreach(i = uniqueMetrics) %do% { #no need to make this parallel. Doing so would only slow things down.
        
        summaryFrame <- data.frame(matrix(ncol = 11, nrow = 1), stringsAsFactors = TRUE) #just in case there are no entries in summaryFrame, all entries in this summaryFrame will be NA
        summaryFrame[1,] <- NA
         
        popLevelMetric <- droplevels(subset(chisqOut, metric == i)) # pull the unique value of interest
        
        if(nrow(popLevelMetric) > 0){
        
          observedVec <- c(sum(popLevelMetric[,match("contactDurations.x", colnames(chisqOut))]), sum(popLevelMetric[,match("noContactDurations.x", colnames(chisqOut))])) #create a vector describing total observed counts of empirical contacts and non-contacting timepoints
          maxDurations.x <- sum(c(popLevelMetric[,match("contactDurations.x", colnames(chisqOut))], popLevelMetric[,match("noContactDurations.x", colnames(chisqOut))])) #total number of random temporal sampling-windows observed. 
          expectedVec <- c(sum(popLevelMetric[,match("contactDurations.y", colnames(chisqOut))]), sum(popLevelMetric[,match("noContactDurations.y", colnames(chisqOut))])) #create a vector describing total observed counts of empirical contacts and non-contacting timepoints
          maxDurations.y <- sum(c(popLevelMetric[,match("contactDurations.y", colnames(chisqOut))], popLevelMetric[,match("noContactDurations.y", colnames(chisqOut))])) #total number of random temporal sampling-windows observed. This will be used to created the probability distribution for the NULL model
      
          maxDurations.y<-maxDurations.y[is.na(maxDurations.y) == FALSE] #remove the NAs
          if(length(maxDurations.y) == 0){ #if removing the NAs removed the only observations, then maxDurations.y reverts to 0
            maxDurations.y <- 0
          }
          
          if(maxDurations.y > 0){ #we must prevent dividing by zero here. 
            expectedProb <- c((expectedVec[1]/maxDurations.y), (expectedVec[2]/maxDurations.y)) #convert the observed random durations to probabilities.
            
          }else{ #if maxDurations.y == 0, expectedProb will always be c(0,0) and therefore chisq.test will return an error
            
            expectedProb <- c(0,1) #here we force expected probability to be (0, 1) (i.e., 100% probability of no contact), which would be the case if the individual was not observed.
            
          }
        
          assign("last.warning", NULL, envir = baseenv()) #clears previous warnings so that we may record any new warnings
        
          test<-tryCatch.W.E(stats::chisq.test(x = observedVec, p = expectedProb))$value #get the chisq values
          warn1<-tryCatch.W.E(stats::chisq.test(x = observedVec, p = expectedProb))$warning$message #get any warning message that the chisq may have triggered
        
          if(length(warn1) == 0){
            warn1 = ""
          }
        
          summaryFrame <- data.frame(metric = i, method = unname(test[4]), X.squared = unname(test[1]), 
                                   df = unname(test[2]), p.val = unname(test[3]), empiricalContactDurations = observedVec[1], 
                                   randContactDurations.mean = expectedVec[1], empiricalNoContactDurations = observedVec[2], 
                                   randNoContactDurations.mean = expectedVec[2], difference = abs((maxDurations.x - observedVec[1]) - (maxDurations.y - expectedVec[1])), 
                                   warning = warn1, stringsAsFactors = TRUE)
        
        }
        
        colnames(summaryFrame)<-c("metric",  "method", "X.squared", "df", "p.val", "contactDurations.x", "contactDurations.y", "noContactDurations.x", 
                                  "noContactDurations.y", "difference", "warning") #for some reason the data.frame command above kept producing incorrect colNames.
        
        rownames(summaryFrame) <-1
        
        return(summaryFrame)
        
      }
      
       popLevelOut<- data.frame(data.table::rbindlist(popLevelOut), stringsAsFactors = TRUE) #bind these frames together
       
       chisqOut <- list(chisqOut, popLevelOut) #group the outputs together
       names(chisqOut) <- c("individualLevel", "populationLevel") #rename the objects in the list
    }
    
  }

  return(chisqOut) #outputs the chisqOut object (Note that this can be either a single data frame or a list of 2 data frames)
  
}
