#' Identify Distance Threshold for Contact
#'
#' Sample from a multivariate normal distribution to create "in-contact" point 
#'    pairs (n = n1) based on real-time-location systems accuracy, and 
#'    generate distribution of length n2 describing average distances between 
#'    point pairs. This function outputs upper confidence intervals associated 
#'    with average distances between "in-contact" points.
#' 
#' This function is for adjusting contact-distance thresholds (spTh) to account
#'    for positional accuracty of real-time-location systems, assuming random 
#'    (non-biased) error in location-fix positions relative to true locations. 
#'    Essentially this function can be used to determine an adjusted spTh value
#'    that likely includes >= 99-percent of true contacts defined using the 
#'    initial spTh.
#' @param n1 Numerical. Number of points used in the expected-distance 
#'    distribution(s). Defaults to 1000.
#' @param n2 Numerical. Number of expected-distance distribution iterations to 
#'    be averaaged. Defaults to 1000.
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
#'    C. in Review. Accounting for space and uncertainty in real-time-location-
#'    system-derived contact networks. Ecology and Evolution.
#' @export
#' @return Output is a named vector with 22 observations describing the mean, 
#'    max, and upper 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,
#'    80, 85, 90, 95, and 99-percent CI values calculated from the 
#'    contact-distance distribution.
#' @examples
#' findDistThresh(n1 = 10, n2 = 10, acc.Dist1 = 0.5, acc.Dist2 = NULL, 
#'    pWithin1 = 90, pWithin2 = NULL, spTh = 0.5) 
#' 
#' findDistThresh(n1 = 10, n2 = 10, acc.Dist1 = 0.5, acc.Dist2 = NULL, 
#'    pWithin1 = 90, pWithin2 = NULL, spTh = 0)

findDistThresh<-function(n1 = 1000, n2 = 1000, acc.Dist1 = 0.5, acc.Dist2 = NULL, pWithin1 = 90, pWithin2 = NULL, spTh = 0.666){
  
  dist.distributionFunc<-function(x){ #generate an expected-distance distribution from two points pulled from normal distributions with mean values at point1 and point2 x- and y-values, and standard deviations derived from the acc.Dist and pWithin values. 
  
    euc=function(x) {
      point1 = x.cor=unlist(c(unname(unlist(x[1])),unname(unlist(x[2]))))
      point2 = x.cor=unlist(c(unname(unlist(x[3])),unname(unlist(x[4]))))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    
    conf1.1 <- (1 - unname(unlist(x[4]))/100)/2 #calculate alpha/2 for the accuracy distribution
    conf1.2<-unlist(ifelse(conf1.1 >0, conf1.1,0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
    zscore1<-abs(stats::qnorm(conf1.2)) #calculate z-score
  
    conf2.1<- (1 - unname(unlist(x[5]))/100)/2 
    conf2.2<-unlist(ifelse(conf2.1 >0, conf2.1,0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)) #if conf = 0, it is replaced with this value which is extremely close to 0, so that qnorm doesn't return an inf value
    zscore2<-abs(stats::qnorm(conf2.2)) #calculate z-score
  
    x1<-data.frame(p1.x = stats::rnorm(unname(unlist(x[1])), mean = 1, sd = unname(unlist(x[2]))/zscore1),p1.y = stats::rnorm(unname(unlist(x[1])), mean = 1, sd = unname(unlist(x[2]))/zscore1),p2.x = stats::rnorm(unname(unlist(x[1])), mean = 1, sd = unname(unlist(x[3]))/zscore2),p2.y = stats::rnorm(unname(unlist(x[1])), mean = 1+unname(unlist(x[6])), sd = unname(unlist(x[3]))/zscore2)) #note that one of the points is spTh-m away from the other
    dist.distr<-apply(x1,1,euc)
    dist.mean<-mean(dist.distr)
    dist.sd<-stats::sd(dist.distr)
    CIseq<-seq(5,95,5)
    CIseq<-c(CIseq, 99)
    CIvec<-NULL #create a vector to hold the upper 5-95% & 99% CI limits 
    for(i in CIseq){
      alphaOver2<-(1 - i/100)/2
      marginOfError<-abs(stats::qnorm(alphaOver2))*(dist.sd/sqrt(unname(unlist(x[1])))) #calculate the margin of error. note that n = number of points used in the expected-distance distribution (i.e., n1).
      CIvec<-c(CIvec, dist.mean + marginOfError) # Note that we're only interested in taking the upper limit, as this gaurentees that the true mean is covered.
    }
    CIvec<-c(dist.mean, CIvec,max(dist.distr)) #add the 0% (i.e., just the mean) and upper 100% (i.e., the max) CI limit
    names(CIvec)<- c("mean", paste(CIseq, "%-CI", sep =""), "max")
    return(CIvec)
  }
  
  if(is.null(pWithin2) == TRUE){pWithin2 = pWithin1} #If == NULL, defaults to pWithin1 value.
  if(is.null(acc.Dist2) == TRUE){acc.Dist2 = acc.Dist1} #If == NULL, defaults to acc.Dist1 value.
  
  inputFrame<-data.frame(matrix(ncol =6, nrow = n2))
  inputFrame[,1]<-n1
  inputFrame[,2]<-acc.Dist1
  inputFrame[,3]<-acc.Dist2
  inputFrame[,4]<-pWithin1
  inputFrame[,5]<-pWithin2
  inputFrame[,6]<-spTh
  distributions<-apply(inputFrame,1,dist.distributionFunc) #creates a matrix of distance-distribution CI values
  dist.means <- apply(distributions, 1,mean)
  return(dist.means)
}
