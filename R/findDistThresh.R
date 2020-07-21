#' Identify Point-Based Distance Threshold for Contact
#'
#' Sample from a multivariate normal distribution to create "in-contact" n 
#'    point pairs based on real-time-location systems accuracy, and 
#'    generate a distribution describing observed distances between point
#'    pairs. 
#' 
#' This function is for adjusting contact-distance thresholds (spTh) to account
#'    for positional accuracy of real-time-location systems, assuming random 
#'    (non-biased) error in location-fix positions relative to true locations. 
#'    Essentially this function can be used to determine an adjusted spTh value
#'    that likely includes >= 99-percent of true contacts defined using the 
#'    initial spTh.
#'    
#'    NOTE that the description given in our paper, Farthing et al. (2020) 
#'    (full reference below), for how this function works is slightly 
#'    incorrect. In the paper, we state that contact definitions should be 
#'    updated according to an upper 99% confidence interval calculated from a 
#'    distance distribution generated from a multivariate in-contact point 
#'    distribution. THIS IS NOT CORRECT, however, as this confidence interval
#'    relates to determination of the mean of the distance distribution (which
#'    is not what we want). Instead, the desired spatial threshold should be 
#'    updated based on standard deviations away from the mean. Because of the
#'    empirical rule, we can capture the spatial threshold values capable of 
#'    capturing roughly 84, 98 and 100% of observed contacts. We are in the 
#'    process of writing another paper addressing this mistake and improving
#'    the method.
#'    
#' @param n Integer. Number of "in-contact" point-pairs used in the 
#'    expected-distance distribution(s). Defaults to 1000.
#' @param acc.Dist1 Numerical. Accuracy distance for point 1.
#' @param acc.Dist2 Numerical. Accuracy distance for point 2. If == NULL, 
#'    defaults to acc.Dist1 value.
#' @param pWithin1 Numerical. Percentage of data points within acc.Dist of true
#'    locations for point 1. 
#' @param pWithin2 Numerical. Percentage of data points within acc.Dist of true
#'    locations for point 2. If == NULL, defaults to pWithin1 value.
#' @param spTh Numerical. Pre-determined distance representing biological 
#'    threshold for contact.
#' @keywords contact location point
#' @references Farthing, T.S., Dawson, D.E., Sanderson, M.W., and Lanzas, 
#'    C. 2020. Accounting for space and uncertainty in real-time-location-
#'    system-derived contact networks. Ecology and Evolution 10(11):4702-4715.
#' @import foreach
#' @export
#' @return Output is a list containing 3 named vectors. The first vector 
#'    describes summary statistics of the simulated distance distribution. The
#'    second vector described adjusted spTh values that will capture 
#'    approximately 84, 98, and 100% of true contacts given the pre-determined
#'    spTh value (all calculated using the Empirical rule). Finally, the third 
#'    vector describes the observed frequency of contact observation given the
#'    spTh adjustments listed in the second vector. 
#' @examples
#' findDistThresh(n = 10,  acc.Dist1 = 0.5, acc.Dist2 = NULL, 
#'    pWithin1 = 90, pWithin2 = NULL, spTh = 0.5) 
#' 

findDistThresh<-function(n = 1000, acc.Dist1 = 0.5, acc.Dist2 = NULL, pWithin1 = 90, pWithin2 = NULL, spTh = 0.666){
  
  #bind the following variables to the global environment so that the CRAN check doesn't flag them as potential problems
  i <- NULL
  
  dist.distributionFunc<-function(x, simFPR, n.outOfContact, range.outOfContact){ #generate an expected-distance distribution from two points pulled from normal distributions with mean values at point1 and point2 x- and y-values, and standard deviations derived from the acc.Dist and pWithin values. 
  
    euc=function(x) {
      point1 = x.cor=unlist(c(unname(unlist(x[1])),unname(unlist(x[2]))))
      point2 = x.cor=unlist(c(unname(unlist(x[3])),unname(unlist(x[4]))))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    
    conf1.1 <- (1 - unname(unlist(x[4]))/100)/2 #calculate alpha/2 for the accuracy distribution
    conf1.2<-unlist(ifelse(conf1.1 >0, conf1.1, ((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
    if(conf1.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
      
      warning("pWithin1 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist1 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
      
    }
    zscore1<-abs(stats::qnorm(conf1.2)) #calculate z-score
  
    conf2.1<- (1 - unname(unlist(x[5]))/100)/2 
    conf2.2<-unlist(ifelse(conf2.1 >0, conf2.1,((1 - 99.99/100)/2))) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
    if(conf2.1 == 0){ #Ensure that users know that users know if conf1.1 was changed
      
      warning("pWithin2 == 100%. To prevent an error here, probability of true locations falling outside of an acc.Dist2 radius around reported point locations changed from 0 to 0.01 (i.e., pWithin1 = 99.99). As a result, spTh estimates may be inflated.", immediate. = TRUE)
      
    }
    
    zscore2<-abs(stats::qnorm(conf2.2)) #calculate z-score
  
    x1<-data.frame(p1.x = stats::rnorm(unname(unlist(x[1])), mean = 1, sd = unname(unlist(x[2]))/zscore1),p1.y = stats::rnorm(unname(unlist(x[1])), mean = 1, sd = unname(unlist(x[2]))/zscore1),p2.x = stats::rnorm(unname(unlist(x[1])), mean = 1, sd = unname(unlist(x[3]))/zscore2),p2.y = stats::rnorm(unname(unlist(x[1])), mean = 1+unname(unlist(x[6])), sd = unname(unlist(x[3]))/zscore2), stringsAsFactors = TRUE) #note that one of the points is spTh-m away from the other
    dist.distr<-apply(x1,1,euc)
    dist.mean<-mean(dist.distr)
    dist.sd<-stats::sd(dist.distr)

    dataFall84_98_100 <- c((dist.mean + dist.sd), (dist.mean + (dist.sd*2)), (dist.mean + (dist.sd*3))) #due to the empirical rule, approximately we can calculate where 68, 95, and 99% of data fall (i.e., the distances required to capture contacts in these cases).
    names(dataFall84_98_100 ) <- c("spTh_84%Capture", "spTh_98%Capture", "spTh_100%Capture") #names
    
    TruePositive <- length(which(dist.distr <= as.numeric(x[6]))) #number of contacts correctly identified.
    TPR<- TruePositive/as.integer(x[1]) #True positive rate (i.e., sensitivity)
    
    distr.summary <- c(dist.mean, dist.sd, min(dist.distr), max(dist.distr), TPR)
    names(distr.summary) <- c("mean", "sd", "min", "max", "TPR")
    
    #dataFall.noise <- (unlist(foreach::foreach(i = 1:length(dataFall84_98_100)) %do% length(which(dist.distr <= as.numeric(dataFall84_98_100[i])))) / as.integer(x[1]))/TPR #this gives you the relative number of observed contact durations for each sd interval, relative to the SpTh value. It's a measure of noise, as we expect the TPR to be constant, but observed positive-contact rates reported here inflate this value. Divide any weighted contacts generated using the chosen CI value by the constant given here, to revert back to the TRUE positive rate associated with the original spTh value.
    dataFall.noise <- (unlist(foreach::foreach(i = 1:length(dataFall84_98_100)) %do% length(which(dist.distr <= as.numeric(dataFall84_98_100[i])))) / as.integer(x[1])) #this gives you the proportion of observed contact durations for each sd interval, relative to the SpTh value. It's a measure of noise, as we expect the TPR to be constant, but observed positive-contact rates reported here inflate this value.
    names(dataFall.noise) <- paste("freq_",c("84%Capture", "98%Capture", "100%Capture"), sep ="") 
    
    out.list <- list(distr.summary, dataFall84_98_100 , dataFall.noise)
    names(out.list) <- c("distribution.summary", "spTh.adjustments", "contact.frequency")
    
    return(out.list)
  }
  
  if(is.null(pWithin2) == TRUE){pWithin2 = pWithin1} #If == NULL, defaults to pWithin1 value.
  if(is.null(acc.Dist2) == TRUE){acc.Dist2 = acc.Dist1} #If == NULL, defaults to acc.Dist1 value.

  inputFrame<-data.frame(matrix(ncol =6, nrow = 1), stringsAsFactors = TRUE)
  inputFrame[,1]<-n
  inputFrame[,2]<-acc.Dist1
  inputFrame[,3]<-acc.Dist2
  inputFrame[,4]<-pWithin1
  inputFrame[,5]<-pWithin2
  inputFrame[,6]<-spTh
  distributions<-apply(inputFrame,1,dist.distributionFunc) #creates a matrix of distance-distribution CI values
  return(distributions[[1]]) #for some reason distributions is a nested list, so we have to designate that we want the further nested objects.
}
