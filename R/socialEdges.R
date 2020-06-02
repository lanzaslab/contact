#' Identify Edges in Social Networks
#' 
#' This function identifies edges in social networks representative of 
#'    instances where there are greater or fewer contacts than would
#'    be expected at random, given a pre-determined p-value threshold for 
#'    significance (i.e., alpha level).
#' 
#' This function will automatically import defined time blocks if applicable. 
#'    Furthermore, because this function is intended describe social 
#'    relationships between individuals, any "totalDegree" and 
#'    "totalContactDurations" metrics are not included in function output, 
#'    even if they are present in x.
#' 
#' @param x A data frame created by a contactCompare function (e.g., 
#'    contactCompare_chisq).
#' @param alpha Numerical threshold for determining social significance given 
#'    p-values reported in x. Observations in x with p.values >= alpha will be 
#'    returned by this function.
#' @param weight Vector of length nrow(data.frame(x)) denoting what 
#'    information should be carried over from x to the function output (e.g., 
#'    number of observed contacts). If the weight is not specified, the 
#'    "weight" in function output is presented as the proportion of total 
#'    potential contact durations that nodes were observed in contact with one 
#'    another (in each separate timeblock if applicable).
#' @param removeDuplicates Logical. If removeDuplicates == true, duplicated 
#'    edges are removed are removed from the output. Defaults to TRUE.
#' @keywords network-analysis social-network
#' @return Returns a list with three objects
#'    
#'    \item{Greater}{Data frame of dyads with more contacts than would be 
#'       expected at random given the chosen alpha level.}
#'    \item{Fewer}{Data frame of dyads with fewer contacts than would be 
#'       expected at random given the chosen alpha level.}
#'    \item{p.val_threshold}{Reports the chosen alpha level.}
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
#' CC1 <-contactCompare_chisq(x.summary = emp.summary, y.summary = rand.summary, 
#'                      x.potential = emp.potential, y.potential = rand.potential,
#'                      importBlocks = FALSE, shuffle.type = 0, 
#'                      popLevelOut = TRUE, parallel = FALSE) #no blocking
#' 
#' socEdges <- socialEdges(x = CC1[[1]], alpha = 0.05, weight = NULL, 
#'                      removeDuplicates = TRUE)
#' 
#'    }

socialEdges<-function(x, alpha = 0.05, weight = NULL, removeDuplicates = TRUE){
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  i <- NULL
  j <- NULL
  block <- NULL
  
  #we first remove any instances of totalDegree and totalContactDurations if they exist
  
  if(length(grep("totalDegree", x$metric)) > 0){
    
    td.rm <- which(x$metric == "totalDegree") #identify which observations to remove
    
    x<-droplevels(x[- td.rm,]) #remove observations from x
    
    if(length(weight) > 0){ #if weight is defined, we must also remove observations from here to ensure everything is of the proper length
      
      weight <- weight[- td.rm] #remove observations from weight
      
    }
    
  }
  
  if(length(grep("totalContactDurations", x$metric)) > 0){
    
    tcd.rm <- which(x$metric == "totalContactDurations") #identify which observations to remove
    
    x<-droplevels(x[- tcd.rm,]) #remove observations from x
    
    if(length(weight) > 0){ #if weight is defined, we must also remove observations from here to ensure everything is of the proper length
      
      weight <- weight[- tcd.rm] #remove observations from weight
      
    }
    
  }
  
  
  originFrame<- data.frame(from = x$id, to = x$metric, p.val = x$p.val, dif = (x$contactDurations.x - x$contactDurations.y)) #put these vectors together in a single data frame
  

  
  if(length(weight) > 0){ #if there is a non-zero length weight vector, add the weight column to originFrame
    
    originFrame$weight <- weight
    
  }else{ #if the weight is not specified, the weight is presented as the proportion of total potential contact durations that nodes were observed in contact with one another (in each separate timeblock if applicable). 
    
    originFrame$weight <- x$contactDurations.x/(x$contactDurations.x + x$noContactDurations.x)
    
  }
  
  if(length(grep("block", colnames(x))) > 0){ #if there is a non-zero length block column in x, add the relevant empirical block columns to originFrame
    
    originFrame$block <- x$block.x
    originFrame$block.start <- x$block.start.x
    originFrame$block.end <- x$block.end.x
  }
  
  sigVec <- which(originFrame$p.val >= alpha) #identify which values in pval are statistically significant
  
  if(length(sigVec) == 0){ #if there were no significant observations, we stop the function and report this fact to users
    
    stop("no significant social relationships given the alpha value")
    
  }
  
  originFrame.sig<-droplevels(originFrame[sigVec,-match("p.val", colnames(originFrame))]) #leave only the significant interactions (Note we remove the pval column because it's no longer needed and smaller data frames are processed more quickly than their larger counterparts)
  
  if(removeDuplicates == TRUE){ #we're going to remove duplicate edges for each timeblock (if applicable)
    
    if(length(grep("block", colnames(x))) == 0){ #if there's no block column
      
      rm.vec<- unlist(foreach::foreach(i = unique(originFrame.sig$to)) %do% { #creates a vector of rows describing duplicated edges.
        last_edge<-max(which(originFrame.sig$to == i)) #identifies the last observation of an edge attached to node i
        rm<-(which(originFrame.sig$from[(last_edge + 1):nrow(originFrame.sig)] == i) + last_edge) #identifies which rows after max describe edges attached to node i
        return(rm)
      })
      
      if(length(rm.vec) > 0){ #if there are any identified rows in the rm.vec, we remove them
        originFrame.sig<-droplevels(originFrame.sig[-rm.vec,]) #remove duplicates
      }
      
    }else{ #if there IS a block column (Note the only difference here is that block information is present in the dyad_ids)
      
      blocks.adjusted<- foreach::foreach(j = unique(originFrame.sig$block), .packages = "foreach") %do% { #creates a vector of rows describing duplicated edges.
        
        blockSub<- droplevels(subset(originFrame.sig, block == j)) #subset the data by block
        
        rm_block<-unlist(foreach::foreach(i = unique(blockSub$to)) %do% {
          
          last_edge<-max(which(blockSub$to == i)) #identifies the last observation of an edge attached to node i
          rm<-(which(blockSub$from[(last_edge + 1):nrow(blockSub)] == i) + last_edge) #identifies which rows should be removed to ensure duplicates are taken out of each block
          return(rm)
        })
        
        if(length(rm_block) > 0){ #if there are any identified rows in the rm.vec, we remove them from blockSub
          blockSub<- droplevels(blockSub[-rm_block,]) #remove duplicates
        }
        
        return(blockSub) #return the blockSub data frame
      }
      
      originFrame.sig <- data.frame(data.table::rbindlist(blocks.adjusted)) #bind the adjusted block sets together to redefine originFrame.sig
      
    }
  }
  
  #finally, we identify which significant social realtionships are associated with a GREATER number of contacts and which are associated with FEWER contacts
  greaterContacts<- droplevels(originFrame.sig[which(originFrame.sig$dif > 0), -match("dif", colnames(originFrame.sig))]) #pull significant interactions where empirical contacts were greater than randomized ones.
  
  if(nrow(greaterContacts) == 0){ #if no significantly greater contact events were present, we make sure to inform users.
    
    warning("No instances of significantly more empirical contact events exist given this alpha value. The 'Greater' list object will be empty.")
    
  }
  
  fewerContacts<- droplevels(originFrame.sig[which(originFrame.sig$dif < 0), -match("dif", colnames(originFrame.sig))]) #pull significant interactions where empirical contacts were greater than randomized ones.
  
  if(nrow(fewerContacts) == 0){ #if no significantly fewer contact events were present, we make sure to inform users.
    
    warning("No instances of significantly fewer empirical contact events exist given this alpha value. The 'Fewer' list object will be empty.")
    
  }
  
  out.list<- list(greaterContacts, fewerContacts, alpha) #create the function output
  names(out.list) <- c("Greater", "Fewer", "p.val_threshold") #ensure list objects have the proper names.
  return(out.list)
}
