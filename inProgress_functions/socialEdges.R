#' Identify Edges in Social Networks
#' 
#' This function 
#' 
#' 



socialEdges<-function(x, pThresh = 0.05, weight = NULL, removeDuplicates = TRUE){
  
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
  
  sigVec <- which(originFrame$p.val >= pThresh) #identify which values in pval are statistically significant
  
  if(length(sigVec) == 0){ #if there were no significant observations, we stop the function and report this fact to users
    
    stop("no significant social relationships given the pThresh value")
    
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
    
    warning("No instances of significantly more empirical contact events exist given this pThresh value. The 'Greater' list object will be empty.")
    
  }
  
  fewerContacts<- droplevels(originFrame.sig[which(originFrame.sig$dif < 0), -match("dif", colnames(originFrame.sig))]) #pull significant interactions where empirical contacts were greater than randomized ones.
  
  if(nrow(fewerContacts) == 0){ #if no significantly fewer contact events were present, we make sure to inform users.
    
    warning("No instances of significantly fewer empirical contact events exist given this pThresh value. The 'Fewer' list object will be empty.")
    
  }
  
  out.list<- list(greaterContacts, fewerContacts, pThresh) #create the function output
  names(out.list) <- c("Greater", "Fewer", "p.val_threshold") #ensure list objects have the proper names.
  return(out.list)
}