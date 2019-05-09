#' Randomize or Pseudorandomize Categorical Variables Associated With Individuals
#'
#' This function randomizes the values in a given column (or set of columns (i.e., c(colname(x)[1], colname(x)[2]))), identified by "feature," in a dataset(x). If maintainDistr == TRUE, the number of each unique value in the column will be maintained in the function output. Otherwise, the function will draw on the initial distribution to assign randomized values, but the specific number of each unique value may not be maintained. If shuffle == TRUE, unique values will be replaced with another, random unique value from the same distribution with 100% certainty. For example if the values in a "dose" column c(0mg,100mg,300mg) were shuffled, one possible outcome would be: x$dose.shuff[which(x$dose == "0mg")] <- 300mg, x$dose.shuff[which(x$dose == "100mg")] <- 0mg, and x$dose.shuff[which(x$dose == "300mg")] <- 100mg.
#'
#' Note: the only subfunction utilized in this function (for which parallel is relevant) is used ONLY if shuffle == TRUE.
#' @param x Description imminent
#' @param feature Description imminent
#' @param shuffle Description imminent
#' @param maintainDistr Description imminent
#' @param numRandomizations Description imminent
#' @keywords randomize analysis
#' @export
#' @examples
#' Examples imminent

randomizeFeature<-function(x, feature = NULL, shuffle = FALSE, maintainDistr = TRUE, numRandomizations = 10, parallel = TRUE){
  shuffle.func<-function(x,y){
    shufVec <- unlist(rep(x[2],length(which(y == x[1]))))
    return(shufVec)
  }
  randomization.func <-function(x, feature, shuffle, maintainDistr, numRandomizations, parallel){
    replicate<- unique(x$rep)
    x<-x[,-match(rep,names(x))] #we want to position this column at the end of the output dataframe.
    
    for(i in 1:length(feature)){
      
      if(shuffle == TRUE){
        x <- x[order(x[,match(feature[i],colnames(x))]),]
        y <- as.character(x[,match(feature[i],colnames(x))])
        emp.feat <-unique(y)
        rand.feat <- sample(emp.feat,length(emp.feat), replace = FALSE)
        randTab <- data.frame(emp.feat,rand.feat)
        
        if (parallel == TRUE){
          cl<-parallel::makeCluster(parallel::detectCores())
          randVec<-parallel::parApply(cl, randTab, 1, shuffle.func, y)
          parallel::stopCluster(cl)
        }else{ #if parallel == FALSE
          randVec <- apply(randTab, 1, shuffle.func, y)	
        }
        
        x[,(length(x) + 1)] <- randVec
        colnames(x)[length(x)] <- paste(feature[i],".shuff",sep = "")
        
      }else{ #if shuffle == FALSE
        if(maintainDistr == TRUE){
          emp.feat <- unname(unlist(x[,match(feature[i],colnames(x))]))
          rand.feat <- sample(emp.feat,length(emp.feat), replace = FALSE)
          x[,(length(x) + 1)] <- rand.feat
          colnames(x)[length(x)] <- paste(feature[i],".rand",sep = "")
        }else{ #if maintainDitr == FALSE
          emp.feat <- unname(unlist(x[,match(feature[i],colnames(x))]))
          rand.feat <- sample(emp.feat,length(emp.feat), replace = TRUE)
          x[,(length(x) + 1)] <- rand.feat
          colnames(x)[length(x)] <- paste(feature[i],".rand",sep = "")
        }
      }
    }
    x$randomRep <- replicate
    return(x)
  }
  input<-x
  inputList<-list()
  for(a in 1:numRandomizations){
    input$rep<-a
    inputList<-c(inputList, input)
  }
  output<-lapply(inputList, randomization.func, feature, shuffle, maintainDistr, numRandomizations, parallel) #I don't know how lapply will interact witht the parApply called in the subfunction... My gut tells me it won't do well because the lapply also parallelizes the function. We may need to go in and remove the parallel argument. 
  return(output)
}