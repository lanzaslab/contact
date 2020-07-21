RTLS_SensSpec<-function(inContact.n = 1000, outContact.n = NULL, acc.Dist1 = 0.5, acc.Dist2 = NULL, pWithin1 = 90, pWithin2 = NULL, spTh = 0.666, outContact.range = c(0.667, 1.167)){
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  i <- NULL
  
  euc=function(x) { #distance calculation function
    point1 = x.cor=unlist(c(unname(unlist(x[1])),unname(unlist(x[2]))))
    point2 = x.cor=unlist(c(unname(unlist(x[3])),unname(unlist(x[4]))))
    euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
    return(euc.dis)
  }
  
  if(is.null(pWithin2) == TRUE){pWithin2 = pWithin1} #If == NULL, defaults to pWithin1 value.
  if(is.null(acc.Dist2) == TRUE){acc.Dist2 = acc.Dist1} #If == NULL, defaults to acc.Dist1 value.
  if(is.null(outContact.n) == TRUE){outContact.n = inContact.n} #If == NULL, defaults to inContact.n value.
  
  if(is.na(outContact.range[1]) == TRUE){ #if the minimum out-of-contact distance is not specified, the minimum distance updates to spTh + 0.01 and users are informed.
    
    outContact.range[1] = spTh + 0.01
    warning("Minimum distance for out-of-contact points is unspecified. Defaulting to spTh + 0.01.", immediate. = TRUE)
    
  }
  
  if(is.na(outContact.range[2]) == TRUE){ #if the maximum out-of-contact distance is not specified, the maximum distance updates to acc.Dist1 + spTh + 0.01 and users are informed.
    
    outContact.range[2] = spTh + acc.Dist1 + 0.01 #we chose this value because it is the maximum value that two out-of-contact points could be away from one another while plausibly being observed "in-contact" with one another due to positioning errors (NOTE that this is only the case if acc.Dist1 = acc.Dist2).
    warning("Maximum distance for out-of-contact points is unspecified. Defaulting to spTh + acc.Dist1 + 0.01.", immediate. = TRUE)
    
  }
  
  if(outContact.range[1] <= spTh | outContact.range[2] <= spTh){ #if either the minimum or maximum out-of-contact distance values are less than the spTh, the process stops
    
    stop("Out-of-contact distance value(s) is less than or equal to the designated spatial threshold for contact.")
    
  }
  
  if(outContact.range[1] > outContact.range[2]){ #if the minimum out-of-contact distance value is greater than the maximum, they are swapped.
    
    newMin <- outContact.range[2] #pull the original (smaller) maximum value
    outContact.range[2] <- outContact.range[1] #overwrite the maximum value with the original minimum one.
    outContact.range[1] <- newMin #overwrite the minimum value with the original maximum one.
    
    warning("Minimum and maximum out-of-contact distances are reversed.", immediate. = TRUE)
    
  }
  
  #calculate point-location standard deviations
  
  conf1.1 <- (1 - unname(unlist(pWithin1))/100)/2 #calculate alpha/2 for the accuracy distribution
  conf1.2<-unlist(ifelse(conf1.1 >0, conf1.1, ((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
  if(conf1.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
    
    warning("pWithin1 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist1 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
    
  }
  zscore1<-abs(stats::qnorm(conf1.2)) #calculate z-score
  
  conf2.1<- (1 - unname(unlist(pWithin2))/100)/2 
  conf2.2<-unlist(ifelse(conf2.1 >0, conf2.1,((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
  if(conf2.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
    
    warning("pWithin2 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist2 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
    
  }
  
  zscore2<-abs(stats::qnorm(conf2.2)) #calculate z-score
  
  #generate distance distribution for in-contact points
  x1<-data.frame(p1.x = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p1.y = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p2.x = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist2))/zscore2),p2.y = stats::rnorm(unname(unlist(inContact.n)), mean = 1+unname(unlist(spTh)), sd = unname(unlist(acc.Dist2))/zscore2), stringsAsFactors = TRUE) #note that one of the points is spTh-m away from the other
  dist.distr1<-apply(x1,1,euc) #calculate inter-point distances
  dist.mean1<-mean(dist.distr1) #calculate mean of distance distribution
  dist.sd1<-stats::sd(dist.distr1) #calculate standard deviation of the distance distribution
  dist.min1 <- min(dist.distr1) #calculate min of distance distribution
  dist.max1 <- max(dist.distr1) #calculate max of distance distribution
  TruePositive <- length(which(dist.distr1 <= as.numeric(spTh))) #number of contacts correctly identified.
  FalseNegative <- length(which(dist.distr1 > as.numeric(spTh))) #number of contacts incorrectly identified.
  TPR<- TruePositive/as.integer(inContact.n) #True positive rate (i.e., sensitivity). This is the same as calculating TP/ (TP + FN).
  FNR <- 1 - TPR #False negative rate
  
  #generate distance distribution for out-of-contact points
  x2.list <- foreach::foreach(i = 1:outContact.n) %do% data.frame(p1.x = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p1.y = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p2.x = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist2))/zscore2),p2.y = stats::rnorm(n = 1, mean = 1+runif(n = 1, min = unname(unlist(outContact.range[1])),max = unname(unlist(outContact.range[2]))), sd = unname(unlist(acc.Dist2))/zscore2), stringsAsFactors = TRUE) #note that one of the points is a random number within the outContact.range range-m away from the other
  x2 <- data.frame(data.table::rbindlist(x2.list)) #rbind the list together
  dist.distr2<-apply(x2,1,euc) #calculate inter-point distances
  dist.mean2<-mean(dist.distr2) #calculate mean of distance distribution
  dist.sd2<-stats::sd(dist.distr2) #calculate standard deviation of the distance distribution
  dist.min2 <- min(dist.distr2) #calculate min of distance distribution
  dist.max2 <- max(dist.distr2) #calculate max of distance distribution
  FalsePositive <- length(which(dist.distr2 <= as.numeric(spTh))) #number of non-contacts incorrectly identified.
  TrueNegative <- length(which(dist.distr2 > as.numeric(spTh))) #number of non-contacts correctly identified.
  TNR<- TrueNegative/as.integer(outContact.n) #True negative rate (i.e., specificity). This is the same as calculating TN/ (FP + TN).
  FPR <- 1 - TNR #False positive rate
  
  #compile output list
  
  distr1.summary <- c(dist.mean1, dist.sd1, dist.min1, dist.max1) #string together summary statistics for output
  names(distr1.summary) <- c("mean", "sd", "min", "max")
  
  distr2.summary <- c(dist.mean2, dist.sd2, dist.min2, dist.max2) #string together summary statistics for output
  names(distr2.summary) <- c("mean", "sd", "min", "max")
  
  dist.counts <- c(TruePositive, TrueNegative, FalsePositive, FalseNegative) #string together observed counts
  names(dist.counts) <- c("TruePositive", "TrueNegative", "FalsePositive", "FalseNegative")
  
  dist.rates <- c(TPR, TNR, FPR, FNR) #string together observed rates
  names(dist.rates) <- c("sensitivity", "specificity", "falsePosRate", "falseNegRate")
  
  out.list <- list(distr1.summary, distr2.summary , dist.counts, dist.rates)
  names(out.list) <- c("inContactDistr.summ", "outContactDistr.summ", "counts", "rates")
  
  return(out.list) 
}

  
estim.TrueContacts <- function(x, TPR, FPR, FNR){
  
  #x is the sum of observed contacts
  #y is a representation of the sum TRUE (i.e., real-world) contacts that is likely unobserved due to RTLS/GPS accuracy
  
  #x = TPR*(y) + FNR*(y) - FPR*(y)  # observed contacts equals the sum of True contacts and false negatives taken from the theoretical distribution, minus the false positives. 
  ##Because everything on the right-hand side is multiplied by y, this formula can be written as:
  #x = ((TPR + FNR) - FPR)*(y) #note we add the extra parantheses just to ensure proper order of operations in the code.
  #x / ((TPR + FNR) - FPR) = y #So, we can estimate TRUE contacts given that we know sensitivity and specificity of the system.
  
  y <- (x / ((TPR + FNR) - FPR)) #solve for y
  
  return(y) #return the number of True contacts
  
}  

#maximize sensitivity & specificity (technically, the function minimizes FPR and FNR, but these are inversely proportional to sensitivity and specificity) to find the ideal contact threshold given 
##Note that we can't just maximize one or the other because maximizing sensitivity would result in the largest posible threshold, while maximizing specificity would result in the smallest possible threshold.
#because optim can only minimize an input variable (spTh) based on a single output, we're technically minimizing (FPR*FNR), and thus maximizing sensitivity*specificity 


spTh_finder<-function(spTH.init = 0, inContact.n = 1000, outContact.n = NULL, acc.Dist1 = 0.5, acc.Dist2 = NULL, pWithin1 = 90, pWithin2 = NULL, outContact.range = c(NA, NA)){
  
  RTLS_SensSpec<-function(inContact.n = 1000, outContact.n = NULL, acc.Dist1 = 0.5, acc.Dist2 = NULL, pWithin1 = 90, pWithin2 = NULL, spTh = 0.666, outContact.range = c(0.667, 1.167)){
    
    #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
    i <- NULL
    
    euc=function(x) { #distance calculation function
      point1 = x.cor=unlist(c(unname(unlist(x[1])),unname(unlist(x[2]))))
      point2 = x.cor=unlist(c(unname(unlist(x[3])),unname(unlist(x[4]))))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    
    if(is.null(pWithin2) == TRUE){pWithin2 = pWithin1} #If == NULL, defaults to pWithin1 value.
    if(is.null(acc.Dist2) == TRUE){acc.Dist2 = acc.Dist1} #If == NULL, defaults to acc.Dist1 value.
    if(is.null(outContact.n) == TRUE){outContact.n = inContact.n} #If == NULL, defaults to inContact.n value.
    
    if(is.na(outContact.range[1]) == TRUE){ #if the minimum out-of-contact distance is not specified, the minimum distance updates to spTh + 0.01 and users are informed.
      
      outContact.range[1] = spTh + 0.01
      warning("Minimum distance for out-of-contact points is unspecified. Defaulting to spTh + 0.01.", immediate. = TRUE)
      
    }
    
    if(is.na(outContact.range[2]) == TRUE){ #if the maximum out-of-contact distance is not specified, the maximum distance updates to acc.Dist1 + spTh + 0.01 and users are informed.
      
      outContact.range[2] = spTh + acc.Dist1 + 0.01 #we chose this value because it is the maximum value that two out-of-contact points could be away from one another while plausibly being observed "in-contact" with one another due to positioning errors (NOTE that this is only the case if acc.Dist1 = acc.Dist2).
      warning("Maximum distance for out-of-contact points is unspecified. Defaulting to spTh + acc.Dist1 + 0.01.", immediate. = TRUE)
      
    }
    
    if(outContact.range[1] <= spTh | outContact.range[2] <= spTh){ #if either the minimum or maximum out-of-contact distance values are less than the spTh, the process stops
      
      stop("Out-of-contact distance value(s) is less than or equal to the designated spatial threshold for contact.")
      
    }
    
    if(outContact.range[1] > outContact.range[2]){ #if the minimum out-of-contact distance value is greater than the maximum, they are swapped.
      
      newMin <- outContact.range[2] #pull the original (smaller) maximum value
      outContact.range[2] <- outContact.range[1] #overwrite the maximum value with the original minimum one.
      outContact.range[1] <- newMin #overwrite the minimum value with the original maximum one.
      
      warning("Minimum and maximum out-of-contact distances are reversed.", immediate. = TRUE)
      
    }
    
    #calculate point-location standard deviations
    
    conf1.1 <- (1 - unname(unlist(pWithin1))/100)/2 #calculate alpha/2 for the accuracy distribution
    conf1.2<-unlist(ifelse(conf1.1 >0, conf1.1, ((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
    if(conf1.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
      
      warning("pWithin1 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist1 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
      
    }
    zscore1<-abs(stats::qnorm(conf1.2)) #calculate z-score
    
    conf2.1<- (1 - unname(unlist(pWithin2))/100)/2 
    conf2.2<-unlist(ifelse(conf2.1 >0, conf2.1,((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
    if(conf2.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
      
      warning("pWithin2 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist2 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
      
    }
    
    zscore2<-abs(stats::qnorm(conf2.2)) #calculate z-score
    
    #generate distance distribution for in-contact points
    x1<-data.frame(p1.x = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p1.y = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p2.x = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist2))/zscore2),p2.y = stats::rnorm(unname(unlist(inContact.n)), mean = 1+unname(unlist(spTh)), sd = unname(unlist(acc.Dist2))/zscore2), stringsAsFactors = TRUE) #note that one of the points is spTh-m away from the other
    dist.distr1<-apply(x1,1,euc) #calculate inter-point distances
    dist.mean1<-mean(dist.distr1) #calculate mean of distance distribution
    dist.sd1<-stats::sd(dist.distr1) #calculate standard deviation of the distance distribution
    dist.min1 <- min(dist.distr1) #calculate min of distance distribution
    dist.max1 <- max(dist.distr1) #calculate max of distance distribution
    TruePositive <- length(which(dist.distr1 <= as.numeric(spTh))) #number of contacts correctly identified.
    FalseNegative <- length(which(dist.distr1 > as.numeric(spTh))) #number of contacts incorrectly identified.
    TPR<- TruePositive/as.integer(inContact.n) #True positive rate (i.e., sensitivity). This is the same as calculating TP/ (TP + FN).
    FNR <- 1 - TPR #False negative rate
    
    #generate distance distribution for out-of-contact points
    x2.list <- foreach::foreach(i = 1:outContact.n) %do% data.frame(p1.x = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p1.y = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p2.x = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist2))/zscore2),p2.y = stats::rnorm(n = 1, mean = 1+runif(n = 1, min = unname(unlist(outContact.range[1])),max = unname(unlist(outContact.range[2]))), sd = unname(unlist(acc.Dist2))/zscore2), stringsAsFactors = TRUE) #note that one of the points is a random number within the outContact.range range-m away from the other
    x2 <- data.frame(data.table::rbindlist(x2.list)) #rbind the list together
    dist.distr2<-apply(x2,1,euc) #calculate inter-point distances
    dist.mean2<-mean(dist.distr2) #calculate mean of distance distribution
    dist.sd2<-stats::sd(dist.distr2) #calculate standard deviation of the distance distribution
    dist.min2 <- min(dist.distr2) #calculate min of distance distribution
    dist.max2 <- max(dist.distr2) #calculate max of distance distribution
    FalsePositive <- length(which(dist.distr2 <= as.numeric(spTh))) #number of non-contacts incorrectly identified.
    TrueNegative <- length(which(dist.distr2 > as.numeric(spTh))) #number of non-contacts correctly identified.
    TNR<- TrueNegative/as.integer(outContact.n) #True negative rate (i.e., specificity). This is the same as calculating TN/ (FP + TN).
    FPR <- 1 - TNR #False positive rate
    
    #compile output list
    
    distr1.summary <- c(dist.mean1, dist.sd1, dist.min1, dist.max1) #string together summary statistics for output
    names(distr1.summary) <- c("mean", "sd", "min", "max")
    
    distr2.summary <- c(dist.mean2, dist.sd2, dist.min2, dist.max2) #string together summary statistics for output
    names(distr2.summary) <- c("mean", "sd", "min", "max")
    
    dist.counts <- c(TruePositive, TrueNegative, FalsePositive, FalseNegative) #string together observed counts
    names(dist.counts) <- c("TruePositive", "TrueNegative", "FalsePositive", "FalseNegative")
    
    dist.rates <- c(TPR, TNR, FPR, FNR) #string together observed rates
    names(dist.rates) <- c("sensitivity", "specificity", "falsePosRate", "falseNegRate")
    
    out.list <- list(distr1.summary, distr2.summary , dist.counts, dist.rates)
    names(out.list) <- c("inContactDistr.summ", "outContactDistr.summ", "counts", "rates")
    
    return(out.list) 
  }
  
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  i <- NULL
  
  if(is.null(pWithin2) == TRUE){pWithin2 = pWithin1} #If == NULL, defaults to pWithin1 value.
  if(is.null(acc.Dist2) == TRUE){acc.Dist2 = acc.Dist1} #If == NULL, defaults to acc.Dist1 value.
  if(is.null(outContact.n) == TRUE){outContact.n = inContact.n} #If == NULL, defaults to inContact.n value.
  
  if(is.na(outContact.range[1]) == TRUE){ #if the minimum out-of-contact distance is not specified, the minimum distance updates to spTh + 0.01 and users are informed.
    
    #outContact.range[1] = spTh + 0.01
    warning("Minimum distance for out-of-contact points is unspecified. Defaulting to spTh + 0.01.", immediate. = TRUE)
    
  }
  
  if(is.na(outContact.range[2]) == TRUE){ #if the maximum out-of-contact distance is not specified, the maximum distance updates to acc.Dist1 + spTh + 0.01 and users are informed.
    
    #outContact.range[2] = spTh + acc.Dist1 + 0.01 #we chose this value because it is the maximum value that two out-of-contact points could be away from one another while plausibly being observed "in-contact" with one another due to positioning errors (NOTE that this is only the case if acc.Dist1 = acc.Dist2).
    warning("Maximum distance for out-of-contact points is unspecified. Defaulting to spTh + acc.Dist1 + 0.01.", immediate. = TRUE)
    
  }
  
  #calculate point-location standard deviations
  
  conf1.1 <- (1 - unname(unlist(pWithin1))/100)/2 #calculate alpha/2 for the accuracy distribution
  conf1.2<-unlist(ifelse(conf1.1 >0, conf1.1, ((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
  if(conf1.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
    
    warning("pWithin1 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist1 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
    
  }
  zscore1<-abs(stats::qnorm(conf1.2)) #calculate z-score
  
  conf2.1<- (1 - unname(unlist(pWithin2))/100)/2 
  conf2.2<-unlist(ifelse(conf2.1 >0, conf2.1,((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
  if(conf2.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
    
    warning("pWithin2 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist2 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
    
  }
  
  zscore2<-abs(stats::qnorm(conf2.2)) #calculate z-score
  
  optim.vec <- c(spTH.init, inContact.n, outContact.n, acc.Dist1, acc.Dist2, outContact.range, zscore1, zscore2)
  
  optimizeFunc <- function(x, inContact.n = optim.vec[2], outContact.n = optim.vec[3], acc.Dist1 = optim.vec[4], acc.Dist2 = optim.vec[5], outContact.range = c(optim.vec[6],optim.vec[7]) , zscore1 = optim.vec[8], zscore2 = optim.vec[9]){
    
    euc=function(x) { #distance calculation function
      point1 = x.cor=unlist(c(unname(unlist(x[1])),unname(unlist(x[2]))))
      point2 = x.cor=unlist(c(unname(unlist(x[3])),unname(unlist(x[4]))))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    
    spTh <- x #define the spatial threshold to be used
    
    if(is.na(outContact.range[1]) == TRUE){ #if the minimum out-of-contact distance is not specified, the minimum distance updates to spTh + 0.01 and users are informed.
      
      outContact.range[1] = spTh + 0.01
      
    }
    
    if(is.na(outContact.range[2]) == TRUE){ #if the maximum out-of-contact distance is not specified, the maximum distance updates to acc.Dist1 + spTh + 0.01 and users are informed.
      
      outContact.range[2] = spTh + acc.Dist1 + 0.01 #we chose this value because it is the maximum value that two out-of-contact points could be away from one another while plausibly being observed "in-contact" with one another due to positioning errors (NOTE that this is only the case if acc.Dist1 = acc.Dist2).
      
    }
    
    #browser()
    
    #generate distance distribution for in-contact points
    x1<-data.frame(p1.x = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p1.y = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p2.x = stats::rnorm(unname(unlist(inContact.n)), mean = 1, sd = unname(unlist(acc.Dist2))/zscore2),p2.y = stats::rnorm(unname(unlist(inContact.n)), mean = 1+unname(unlist(spTh)), sd = unname(unlist(acc.Dist2))/zscore2), stringsAsFactors = TRUE) #note that one of the points is spTh-m away from the other
    dist.distr1<-apply(x1,1,euc) #calculate inter-point distances
    TruePositive <- length(which(dist.distr1 <= as.numeric(spTh))) #number of contacts correctly identified.
    FalseNegative <- length(which(dist.distr1 > as.numeric(spTh))) #number of contacts incorrectly identified.
    TPR<- TruePositive/as.integer(inContact.n) #True positive rate (i.e., sensitivity). This is the same as calculating TP/ (TP + FN).
    FNR <- 1 - TPR #False negative rate
    
    #generate distance distribution for out-of-contact points
    x2.list <- foreach::foreach(i = 1:outContact.n) %do% data.frame(p1.x = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p1.y = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist1))/zscore1),p2.x = stats::rnorm(n = 1, mean = 1, sd = unname(unlist(acc.Dist2))/zscore2),p2.y = stats::rnorm(n = 1, mean = 1+runif(n = 1, min = unname(unlist(outContact.range[1])),max = unname(unlist(outContact.range[2]))), sd = unname(unlist(acc.Dist2))/zscore2), stringsAsFactors = TRUE) #note that one of the points is a random number within the outContact.range range-m away from the other
    x2 <- data.frame(data.table::rbindlist(x2.list)) #rbind the list together
    dist.distr2<-apply(x2,1,euc) #calculate inter-point distances
    FalsePositive <- length(which(dist.distr2 <= as.numeric(spTh))) #number of non-contacts incorrectly identified.
    TrueNegative <- length(which(dist.distr2 > as.numeric(spTh))) #number of non-contacts correctly identified.
    TNR<- TrueNegative/as.integer(outContact.n) #True negative rate (i.e., specificity). This is the same as calculating TN/ (FP + TN).
    FPR <- 1 - TNR #False positive rate
    
    output <- FPR*FNR #generate output
    #if(output = 0){output = 1000000} #if output is 0, because either FPR or FNR is 0, we adjust it, so that optimization will skip this value.This value would indicate that either TPR or TNR was 100% or 0%, and we already know how to get there.
    
    return(output)
    
  }
  
  spTh.optim <- stats::optim(par = optim.vec[1], fn = optimizeFunc, method = "Brent", lower = optim.vec[1], upper = (3*acc.Dist1))
  
  spTh.optim$par #once we have the optimized value, let's run it through the sensitivity tester (essentially recreating the optim function) to pull the sennsitivity and specificity
  
  sensiTest <-RTLS_SensSpec(inContact.n, outContact.n, acc.Dist1, acc.Dist2, pWithin1, pWithin2, spTh = spTh.optim$par, outContact.range)
  
  out.list<-list(spTh.optim$par, sensiTest)
  
  return(out.list) 
}
