#' Randomize or Pseudorandomize Categorical Variables
#'
#' This function randomizes the values in a given column (or set of columns 
#'    (i.e., c(colname(x)[1], colname(x)[2]))), identified by the "feature" 
#'    argument in a dataset (x). 
#'
#' Note: the shuffle argument supercedes the maintainDistr argument. Therefore,
#'    if shuffle == TRUE, the maintainDistr argument is irrelevant. 
#' @param x Data frame containing real-time-location data.
#' @param feature Vector of 1 or more column names describing variables
#'    in x to be randomized.
#' @param shuffle Logical. If TRUE, unique values will be replaced with 
#'    another, random unique value from the same distribution with 100% 
#'    certainty. For example if the values in a "dose" column 
#'    c(0mg,100mg,300mg) were shuffled, one possible outcome would be: 
#'    x$dose.shuff[which(x$dose == "0mg")] <- 300mg, 
#'    x$dose.shuff[which(x$dose == "100mg")] <- 0mg, 
#'    and x$dose.shuff[which(x$dose == "300mg")] <- 100mg.
#'    Defaults to FALSE. 
#' @param maintainDistr Logical. If TRUE, the number of each unique value in 
#'    the column will be maintained in the function output. Otherwise, the 
#'    function will draw on the initial distribution to assign randomized 
#'    values, but the specific number of each unique value may not be 
#'    maintained. Defaults to TRUE.
#' @param numRandomizations Description imminent
#' @keywords randomize analysis
#' @return Output is \code{x} appended with columns described below. 
#'    
#'    \item{...shuff}{Randomized value of specified variables.}
#'    \item{randomRep}{Randomization replicate.}
#'    
#' @references Farine, D.R., 2017. A guide to null models for animal social 
#'    network analysis. Methods in Ecology and Evolution 8:1309-1320.
#'    https://doi.org/10.1111/2041-210X.12772.
#' @export
#' @examples
#' 
#' data(calves)
#' 
#' system.time(randomizedValues<-contact::randomizeFeature(x = calves, 
#'    feature = c("calftag", "date"), shuffle = TRUE, maintainDistr = TRUE, 
#'    numRandomizations = 3)) 
#'    
#' randomizedFrame<-data.frame(randomizedValues[[1]])
#' 
#' head(randomizedFrame) #see that randomized-value columns have been appended.

randomizeFeature<-function(x, feature = NULL, shuffle = FALSE, maintainDistr = TRUE, numRandomizations = 10){
  shuffle.func<-function(x,y){
    shufVec <- unlist(rep(x[2],length(which(y == x[1]))))
    return(shufVec)
  }
  randomization.func <-function(x, feature, shuffle, maintainDistr, numRandomizations){
    replicate<- unique(x$rep)
    x<-x[,-match("rep",names(x))] #we want to position this column at the end of the output dataframe.
    
    for(i in 1:length(feature)){
      
      if(shuffle == TRUE){
        x <- x[order(x[,match(feature[i],colnames(x))]),]
        y <- as.character(x[,match(feature[i],colnames(x))])
        emp.feat <-unique(y)
        rand.feat <- sample(emp.feat,length(emp.feat), replace = FALSE)
        randTab <- data.frame(emp.feat,rand.feat)
        randList <- apply(randTab, 1, shuffle.func, y)	
        randVec<-unname(unlist(randList))
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
    inputList[[a]]<-input
  }
  output<-lapply(inputList, randomization.func, feature, shuffle, maintainDistr, numRandomizations) 
  if(length(output) == 1){ #if there's only one list object in the output, the output will be a data frame, not a list.
    output <- data.frame(output[[1]])
  }
    
  return(output)
}
