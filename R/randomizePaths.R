#' Randomize or Pseudorandomize Individuals' Relocation Events
#'
#' Randomizes or pseudorandomizes individuals' spatial locations. Randomized 
#'    datasets can later be compared to empirical ones to determine if 
#'    individuals' space use differ from what would be expected at random 
#'    (using the contactTest function).
#' 
#' Paths can be randomized, or pseudorandomized differently according to what 
#'    logical arguments are set to TRUE. 
#' 
#' Detailed shuffle.type description:
#'    If shuffle.type == 0, within-block timepoints will be randomized by 
#'    sampling from observed timepoints only within the relevant block 
#'    (e.g., points in block 1 may be redistributed). Block order in the 
#'    data set does not change, and no timepoints will be represented more than
#'    once in the randomized set. If shuffle.type == 1, blocks are shuffled 
#'    around, but points within blocks are not redistributed 
#'    (e.g., block 1 <- block 5, block 3 <- block 2, block 5 <- block 4). If 
#'    shuffle.type == 2, blocks are shuffled, but their relative temporal 
#'    location in the shuffleUnit is maintained. For example, there are 144 
#'    10-min blocks in a 24-hour day. Say our data set contains 3 days worth of
#'    data (i.e., 432 blocks). During the randomization, because blocks 
#'    maintain their relative locations in shuffleUnits (e.g., Days), block 1 
#'    in the random set will be determined by sampling from a distribution of 
#'    blocks 1,145,and 289, which each representing the first block of a given 
#'    shuffleUnit (e.g., Day 1, Day 2, Day 3). All blocks in the randomized set
#'    are decided in this fashion (e.g., block 2 of the randomized set is 
#'    identified by sampling from a distribution of blocks 2, 146, and 290). 
#'    Therefore, blocklength-units may never exceed 1 shuffleUnit 
#'    (e.g., 25-hour blocks cannot be shuffled using shuffleUnit == "Days," 
#'    but 1:24-hour blocks work just fine). Points within blocks are not 
#'    redistributed. This particular shuffle.type (i.e., 2) is based off of 
#'    the methodology described by Spiegel et al. 2016.
#' 
#' Note that, if shuffle.type == 2, all dateTime values in individuals movement
#'    paths described in x must be equidistant (e.g., relocations for 
#'    individual i occur every 10 seconds), or erroneous dateTimes will be 
#'    reported. If raw dateTime values are not equidistant, we recommend using 
#'    our tempAggregate function, with na.rm == FALSE, to make it so. 
#' 
#' @param x Data frame containing real-time-location data.  
#' @param id Vector of length nrow(x) or singular character data, detailing the
#'    relevant colname in x, that denotes what unique ids for tracked 
#'    individuals will be used. If argument == NULL, the function assumes a 
#'    column with the colname "id" exists in x. Defaults to NULL.
#' @param point.x Vector of length nrow(x) or singular character data, 
#'    detailing the relevant colname in x, that denotes what planar-x or 
#'    longitude coordinate information will be used. If argument == NULL, the 
#'    function assumes a column with the colname "x" exists in x. Defaults to 
#'    NULL.
#' @param point.y Vector of length nrow(x) or singular character data, 
#'    detailing the relevant colname in x, that denotes what planar-y or 
#'    lattitude coordinate information will be used. If argument == NULL, the 
#'    function assumes a column with the colname "y" exists in x. Defaults to 
#'    NULL.
#' @param dateTime Vector of length nrow(x) or singular character data, 
#'    detailing the relevant colname in x, that denotes what dateTime 
#'    information will be used. If argument == NULL, the function assumes a 
#'    column with the colname "dateTime" exists in x. Defaults to NULL.
#' @param poly.xy Columns within x denoting polygon xy-coordinates. Polygon 
#'    coordinates must be arranged in the format of those in 
#'    referencePointToPolygon output. Defaults to NULL.
#' @param parallel Logical. If TRUE, sub-functions within the randomizePaths 
#'    wrapper will be parallelized. Note that this can significantly speed up 
#'    processing of relatively small data sets, but may cause R to crash due 
#'    to lack of available memory when attempting to process large datasets. 
#'    Defaults to FALSE.
#' @param nCores Integer. Describes the number of cores to be dedicated to 
#'    parallel processes. Defaults to the maximum number of cores available
#'    (i.e., parallel::detectCores()).
#' @param dataType Character string refering to the type of real-time-location 
#'    data presented in x, taking values of "Point" or "Polygon." If 
#'    dataType == "Point," individuals' locations are drawn from point.x and 
#'    point.y. If argument == "Polygon," individuals' locations are drawn from 
#'    poly.xy. Defaults to "Point."
#' @param numVertices Numerical. If dataType == "Polygon," users must specify 
#'    the number of vertices contained in each polygon. Defaults to 4. Note: 
#'    all polygons must contain the same number of vertices.
#' @param blocking Logical. If TRUE, prior to randomization, timepoints will be
#'    categorized into a series of temporal blocks of blockLength-blockUnit 
#'    length (e.g., 10 mins). After generating blocks, the spatial-location 
#'    randomization methodology will follow shuffle.type. If FALSE, paths will 
#'    be randomized by sampling from observed timepoints. No timepoints will be
#'    represented more than once in the randomized set. Defaults to TRUE.
#' @param blockUnit Numerical. Describes the number blockUnits within each 
#'    temporal block. Defaults to 1.
#' @param blockLength Character string taking the values, "secs," "mins," 
#'    "hours," "days," or "weeks." Describes the temporal unit associated with 
#'    each block. Defaults to "hours."
#' @param shuffle.type Numerical. Describes which shuffle.type is used to 
#'    randomize the rand.input data set(s), given that blocking == TRUE 
#'    (Note: this value is irrelevant if blocking == FALSE). Takes the values 
#'    "0," "1," or "2," and defaults to 0. Descriptions of each shuffle.type 
#'    value can be found under Details.
#' @param shuffleUnit Character string taking the values, "secs," "mins," 
#'    "hours," "days," or "weeks." Defaults to "days." Describes what temporal 
#'    unit blocks will be shuffled across given shuffle.type == 2. 
#'    Blocklength-units may never exceed 1 shuffleUnit 
#'    (e.g., 25-hour blocks cannot be shuffled using shuffleUnit == "Days," but
#'    1:24-hour blocks work just fine).
#' @param indivPaths Logical. If TRUE, paths will be randomized with no 
#'    location switching between ids (e.g., randomized xy locations for 
#'    individual 1 will be generated by sampling only from individual 1's 
#'    location distribution). If FALSE, paths will be randomized with potential
#'    location switching between ids (e.g., randomized xy locations for 
#'    individual 1 will be generated by sampling from the entire dataset's 
#'    location distribution). Defaults to TRUE.
#' @param numRandomizations Numerical. The number of replicate data frames 
#'    produced in output. Defaults to 1.
#' @keywords data-processing contact
#' @return Output is \code{x} appended with columns described below.
#'    
#'    \item{x.rand}{Randomized values taken from the \code{point.x} 
#'    argument.}
#'    \item{y.rand}{Randomized values taken from the \code{point.y} 
#'    argument.}
#'    \item{shuffle.type}{The value specified by the \code{shuffle.type} 
#'    argument.}
#'    \item{rand.rep}{Randomization replicate.}
#'    
#'    If blocking == TRUE, the following codes are appended to the output
#'    data frame described above:
#'    
#'    \item{block}{Integer ID describing unique blocks of time during which 
#'    contacts occur.}
#'    \item{block.start}{The timepoint in \code{x} at which the \code{block}
#'    begins.}
#'    \item{block.end}{The timepoint in \code{x} at which the \code{block}
#'    ends.}
#'    \item{numBlocks}{Integer describing the total number of time blocks 
#'    observed within \code{x} at which the \code{block}}
#'    \item{blockLength}{Character string describing the length of blocks 
#'    described by \code{blockLength} and \code{blockUnit} arguments.}
#'    
#' @references Spiegel, O., Leu, S.T., Sih, A., and C.M. Bull. 2016. Socially 
#'    interacting or indifferent neighbors? Randomization of movement paths to 
#'    tease apart social preference and spatial constraints. Methods in 
#'    Ecology and Evolution 7:971-979. https://doi.org/10.1111/2041-210X.12553.
#'    
#'    Farine, D.R., 2017. A guide to null models for animal social network 
#'    analysis. Methods in Ecology and Evolution 8:1309-1320.
#'    https://doi.org/10.1111/2041-210X.12772.
#' @export
#' @examples
#' data(calves)
#' 
#' calves.dateTime<-datetime.append(calves, date = calves$date, 
#'    time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#'    
#' calves.agg<-tempAggregate(calves.dateTime, id = calves.dateTime$calftag, 
#'    dateTime = calves.dateTime$dateTime, point.x = calves.dateTime$x, 
#'    point.y = calves.dateTime$y, secondAgg = 300, extrapolate.left = FALSE, 
#'    extrapolate.right = FALSE, resolutionLevel = "reduced", parallel = FALSE, 
#'    na.rm = TRUE, smooth.type = 1) #smooth to 5-min fix intervals. 
#'
#' calves.agg.rand<-randomizePaths(x = calves.agg, id = "id", 
#'    dateTime = "dateTime", point.x = "x", point.y = "y", poly.xy = NULL, 
#'    parallel = FALSE, dataType = "Point", numVertices = 1, blocking = TRUE, 
#'    blockUnit = "mins", blockLength = 10, shuffle.type = 0, shuffleUnit = NA,
#'    indivPaths = TRUE, numRandomizations = 1) 
#'    

randomizePaths<-function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, poly.xy = NULL, parallel = FALSE, nCores = parallel::detectCores(), dataType = "Point", numVertices = 4, blocking = TRUE, blockUnit = "hours", blockLength = 1, shuffle.type = 0, shuffleUnit = "days", indivPaths = TRUE, numRandomizations = 1){
  
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
  randomization<-function(x,y){
    
    dataType = unlist(unname(x[2]))
    numVertices = unlist(unname(x[3]))
    shuffle.type = unlist(unname(x[4]))
    indivPaths = unlist(unname(x[5]))
    blocking = unlist(unname(x[6]))
    shuffleUnit = unlist(unname(x[7]))
    blockUnit = unlist(unname(x[8]))
    blockLength = as.integer(unlist(unname(x[9]))) #1/31 had to ass the as.integer call to avoid a "Error in (24)/blockLength : non-numeric argument to binary operator" error in the apply function.
    y <- y[order(y$id, y$dateTime),]
    
    if(indivPaths == TRUE){
      idVecSeq <- unique(y$id)
    }else{
      idVecSeq <- 1
    }
    numids <- length(idVecSeq)
    idTab<- data.frame(idVecSeq)
    randomness<- apply(idTab,1,datasub.func, y, numids, blocking, shuffle.type, dataType, shuffleUnit, blockUnit, blockLength) #prior to 01/31, this was missing the "blockUnit" argument, resulting in an error. Obviously, this is now fixed.
    if(is.data.frame(randomness) == TRUE){
      randomxyCoords<-randomness
    }else{
      randomxyCoords<-data.frame(data.table::rbindlist(randomness))
    }
    
    if(blocking == TRUE & shuffle.type == 2){ #To account for potential length differences in blocks if shuffle.type ==2, the outputTab is essentially compiled earlier on in the function (in the datasub.func subfunction) 
      outputTab <- randomxyCoords
    }else{
      outputTab<- cbind(y, randomxyCoords)
    }
    outputTab$shuffle.type = shuffle.type #added 1/30/2019, to ensure that people remember what shuffle.type they set.
    outputTab$rand.rep = unlist(unname(x[1]))
    return(outputTab)
  }
  datasub.func<-function(x,y,z, blocking, shuffle.type, dataType, shuffleUnit, blockUnit, blockLength){
    if(z > 1){
      idFrame <- y[which(y$id == unlist(unname(x[1]))),]
    }else{
      idFrame <- y
    }
    if(blocking == TRUE){
      if(shuffle.type == 0){ #This is the same as shuffle == FALSE. Blocks maintain their order in the dataset, but the xy locations within each block are randomized.
        blockSeq <-unique(idFrame$block)
        blockTab <-data.frame(blockSeq)
        rand <-apply(blockTab,1,xy.assignNoShuff, idFrame, dataType, numVertices)
        randomxy <- data.frame(data.table::rbindlist(rand))
      }
      if(shuffle.type == 1){ #This rearranges the order of blocks in the dataset. xycoordinates within each block remain unchanged.
        blockSeq <-unique(idFrame$block)
        rand.block <- sample(blockSeq,length(blockSeq), replace = FALSE)
        rand.blockTab <- data.frame(rand.block)
        rand <-apply(rand.blockTab,1,xy.assignShuff, idFrame, dataType, numVertices)
        randomxy <- data.frame(data.table::rbindlist(rand))
      }
      if(shuffle.type == 2){ #This potentially replaces each block in the dataset with one of the blocks from a distribution described by shuffleUnit. Each block has the potential to appear in the dataset up to numShuffUnit times.
        
        if(shuffleUnit == "Mins" || shuffleUnit == "MINS" || shuffleUnit == "mins"){
          if(blockUnit == "Secs" || blockUnit == "SECS" || blockUnit == "secs"){
            blocksPerUnit <- ceiling(60/blockLength)
          }
        }
        if(shuffleUnit == "Hours" || shuffleUnit == "HOURS" || shuffleUnit == "hours"){
          if(blockUnit == "Secs" || blockUnit == "SECS" || blockUnit == "secs"){
            blocksPerUnit <- ceiling((60*60)/blockLength)
          }
          if(blockUnit == "Mins" || blockUnit == "MINS" || blockUnit == "mins"){
            blocksPerUnit <- ceiling(60/blockLength)
          }
        }
        if(shuffleUnit == "Days" || shuffleUnit == "DAYS" || shuffleUnit == "days"){
          if(blockUnit == "Secs" || blockUnit == "SECS" || blockUnit == "secs"){
            blocksPerUnit <- ceiling((60*60*24)/blockLength)
          }
          if(blockUnit == "Mins" || blockUnit == "MINS" || blockUnit == "mins"){
            blocksPerUnit <- ceiling((24*60)/blockLength)
          }
          if(blockUnit == "Hours" || blockUnit == "HOURS" || blockUnit == "hours"){
            blocksPerUnit <- ceiling((24)/blockLength)
          }
        }
        if(shuffleUnit == "Weeks" || shuffleUnit == "WEEKS" || shuffleUnit == "weeks"){
          if(blockUnit == "Secs" || blockUnit == "SECS" || blockUnit == "secs"){
            blocksPerUnit <- ceiling((60*60*24*7)/blockLength)
          }
          if(blockUnit == "Mins" || blockUnit == "MINS" || blockUnit == "mins"){
            blocksPerUnit <- ceiling((60*24*7)/blockLength)
          }
          if(blockUnit == "Hours" || blockUnit == "HOURS" || blockUnit == "hours"){
            blocksPerUnit <- ceiling((24*7)/blockLength)
          }
          if(blockUnit == "Days" || blockUnit == "DAYS" || blockUnit == "days"){
            blocksPerUnit <- ceiling((7)/blockLength)
          }
        }
        blocksPerUnit<- unname(unlist(blocksPerUnit))
        blockSeq <-unique(idFrame$block)
        randVec1 <- NULL
        
        for(d in 1:blocksPerUnit){
          if(d > max(blockSeq)){next} #to prevent returning an error in the extremely unlikely event that someone tried randomize using a shuffle unit that incorporate all blocksPerUnit (which would not actually randomize anything).
          shufUnitSeq <- seq(d,max(blockSeq), blocksPerUnit)
          if(length(shufUnitSeq) == 1){ #because of the nature of the sample function, if the length of shufUnitSeq ==1, the function will draw from a distribution of 1:shufUnitSeq. This is not what we want....
            shufUnitSeq <- c(shufUnitSeq, shufUnitSeq)
          }
          rand.block <- sample(shufUnitSeq,1, replace = TRUE)
          randVec1 <- c(randVec1, rand.block)
        }
        
        if(length(randVec1) > length(blockSeq)){ #if blockLength does not divide evenly into blocksPerUnit, length(randVec) may be greater than the actual number of blocks
          randVec1<-randVec1[1:length(blockSeq)]
        }
        blockAbsence<-which(randVec1%in%unique(idFrame$block) == FALSE) #determine if individuals were actually present during randomly-pulled blocks. So, now we know which blocks are missing and can build the output data frame accordingly.
        if(length(blockAbsence) == length(randVec1)){ #if all blocks were empty, then there's no need to go any further
          randomxy <-NULL
        }else{ #if there was at least one populated block
        
        rand.blockTab <- data.frame(randVec1)
        rand <-apply(rand.blockTab,1,xy.assignShuff, idFrame, dataType, numVertices)
        randCoord<- data.frame(data.table::rbindlist(rand))
        
        blockFrame.full<-subset(y, block <= blocksPerUnit)
        
        if(length(blockAbsence > 0)){ #So, if there is a missing block
          timestep.per.block <- length(unique(blockFrame.full$dateTime))/length(unique(blockFrame.full$block)) #pulls the number of timesteps per block. In the event that randVec pulls a starting block(s) that does not exist for the individual, NAs will be returned. Assuming that all timesteps are equidistant (as they should be when shuffleType == 2), this can only be the case if a block is pulled from a time before the individual began reporting tracks.
          for(a in blockAbsence){ #Put NAs into randCoord where needed
            randCoord.brkPt <-timestep.per.block*(a-1) #Identifies the row immediately prior to where NAs will be placed
            naBlock<-data.frame(x.rand<-rep(NA, timestep.per.block), y.rand<-rep(NA, timestep.per.block)) #add NAs to columns when applicable
            if(randCoord.brkPt == 0){ #if the first block is missing
              randCoord<-data.frame(data.table::rbindlist(list(naBlock, randCoord)))
            }else{ #if a block other than the 1st block is missing
              if(a < length(randVec1)){ #if the missing block is not the last recorded one
                preNACoord<- randCoord[1:randCoord.brkPt,]
                postNACoord<- randCoord[(randCoord.brkPt+1):ncol(randCoord),]
                randCoord<-data.frame(data.table::rbindlist(list(preNACoord, naBlock, postNACoord)))
              }
              if(a == length(randVec1)){ #if the missing block IS the last recorded one
                randCoord<-data.frame(data.table::rbindlist(list(randCoord, naBlock)))
              }              
            }
          }
        }
        #outTable<- idFrame[1:nrow(randCoord),] #If the final block of the dataset has fewer observations than other blocks, and the randomized locations used this final block, this reduces the number of rows in idFrame to fit the number of rows in randCoord
        outTable<-data.frame(id = unname(unlist(x[1])), dateTime = unique(blockFrame.full$dateTime), blockLength = unique(blockFrame.full$blockLength)) #create a dataframe with the number of rows corresponding to the number of times individuals were observed. This data frame will contain essentially the same info. as idFrame, but will ensure that blocks are consistent for all individuals.
        outTable<- outTable[1:nrow(randCoord),] #If the final block of the dataset has fewer observations than other blocks, and the randomized locations used this final block, this reduces the number of rows in idFrame to fit the number of rows in randCoord
        
        #now we have to redefine block information for the single shuffleUnit (so that users can recall the block information they set)
        bindlist <- list(outTable, randCoord)
        randomxy <- do.call("cbind", bindlist)
        randomxy<-randomxy[order(randomxy$id, randomxy$dateTime),]
        }
      }
    }else{ #if blocking == FALSE
      if(dataType == "Point" || dataType == "POINT" || dataType == "point"){
        rand.noblock <- sample(seq(1,nrow(idFrame),1), nrow(idFrame), replace = FALSE)
        rand.x <- idFrame$x[rand.noblock]
        rand.y <- idFrame$y[rand.noblock]
        randomxy <- data.frame(x = rand.x, y = rand.y)
        colnames(randomxy) <- c("x.rand","y.rand")
      }
      if(dataType == "Polygon" || dataType == "POLYGON" || dataType == "polygon"){
        rand.noblock <- sample(seq(1,nrow(idFrame),1), nrow(idFrame), replace = FALSE)
        rand.poly <- data.frame(matrix(ncol = 0, nrow = length(rand.noblock)))
        vertex <-1 #This is the base column number that, if added to "a" in the following loop will call the appropriate xy coordinates.
        for(a in 1:numVertices){
          rand.x <- idFrame[c(rand.noblock),(vertex+a)]
          rand.y <- idFrame[c(rand.noblock),(vertex+a+1)]
          vertex<- vertex+1
          rand.point <-data.frame(rand.x, rand.y)
          colnames(rand.point)<- c(paste("point",a,".x", sep = ""), paste("point",a,".y", sep = ""))
          bindlist <-list(rand.poly, rand.point)
          rand.poly <- do.call("cbind", bindlist)
        }
        randomxy <- rand.poly
      }
    }
    return(randomxy)
  }
  xy.assignShuff<-function(x,y, dataType,numVertices){
    if(dataType == "Point" || dataType == "POINT" || dataType == "point"){
      xVec <-y$x[which(y$block == unlist(unname(x[1])))]
      yVec <-y$y[which(y$block == unlist(unname(x[1])))]
      xyrand <- data.frame(x.rand = xVec, y.rand = yVec)
    }
    if(dataType == "Polygon" || dataType == "POLYGON" || dataType == "polygon"){
      xTab <-y[which(y$block == x[1]),c(seq(2,(numVertices*2),2))]
      yTab <-y[which(y$block == x[1]),c(seq(3,(numVertices*2)+1,2))]
      polygon <-data.frame(NULL)
      for(i in 1:numVertices){
        vert.x <-xTab[,i]
        vert.y <-yTab[,i]
        bindlist <-list(polygon,vert.x,vert.y)
        polygon<-do.call("cbind", bindlist)
      }
      xyrand <- data.frame(polygon)
      colnames(xyrand)[c(seq(1,(numVertices*2)-1,2))]<-paste("point",seq(1,numVertices,1),".x",sep ="")
      colnames(xyrand)[c(seq(2,(numVertices*2),2))]<-paste("point",seq(1,numVertices,1),".y",sep ="")
    }
    return(xyrand)
  }
  xy.assignNoShuff<-function(x,y, dataType,numVertices){
    blockFrame <- y[which(y$block == unlist(unname(x[1]))),]
    randSeq <- sample(seq(1,nrow(blockFrame),1),nrow(blockFrame), replace = FALSE)
    
    if(dataType == "Point" || dataType == "POINT" || dataType == "point"){
      xVec <-blockFrame$x[randSeq]
      yVec <-blockFrame$y[randSeq]
      xyrand <- data.frame(x = xVec, y = yVec)
      colnames(xyrand) <- c("x.rand","y.rand")
    }
    if(dataType == "Polygon" || dataType == "POLYGON" || dataType == "polygon"){
      xTab <-blockFrame[randSeq,c(seq(2,(numVertices*2),2))]
      yTab <-blockFrame[randSeq,c(seq(3,(numVertices*2)+1,2))]
      polygon <-data.frame(NULL)
      for(i in 1:numVertices){
        vert.x <-xTab[,i]
        vert.y <-yTab[,i]
        bindlist <-list(polygon,vert.x,vert.y)
        polygon<-do.call("cbind", bindlist)
      }
      xyrand <- data.frame(polygon)
      colnames(xyrand)[c(seq(1,(numVertices*2)-1,2))]<-paste("point",seq(1,numVertices,1),".x",sep ="")
      colnames(xyrand)[c(seq(2,(numVertices*2),2))]<-paste("point",seq(1,numVertices,1),".y",sep ="")
    }
    return(xyrand)
  }  		
  
  if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors. No matter the input format, a table called "originTab" will be created.
    if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
      originTab = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime)
    }
    if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
      originTab = data.frame(matrix(ncol = 0, nrow = length(id)))
      originTab$id = id
      colnames(poly.xy)[seq(1,(ncol(poly.xy) - 1),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".x", sep = "")
      colnames(poly.xy)[seq(2,ncol(poly.xy),2)] = paste("point",seq(1,(ncol(poly.xy)/2),1),".y", sep = "")
      dateFrame = data.frame(dateTime = dateTime)
      bindlist = list(originTab,poly.xy,dateFrame)
      originTab = data.frame(do.call("cbind", bindlist))
    }
  }
  
  if(length(x) > 0){ #for some reason using an "else" statement would always result in an originTab table with 0 records...
    if(length(id) > 0){ 
      if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ 
        x$id <- x[,match(id, names(x))]
      }else{ #if length(id) > 1
        x$id = id
      }
    }
    idVec1 <- x$id
    
    if(dataType == "point" || dataType == "Point" || dataType == "POINT"){
      if(length(point.x) > 0){
        if(length(point.x) == 1 && is.na(match(point.x[1], names(x))) == FALSE){ 
          x$x <- x[,match(point.x, names(x))]
        }else{ #if length(point.x) > 1
          x$x = point.x
        }
      }
      if(length(point.y) > 0){
        if(length(point.y) == 1 && is.na(match(point.y[1], names(x))) == FALSE){ 
          x$y <- x[,match(point.y, names(x))]
        }else{ #if length(point.y) > 1
          x$y = point.y
        }
      }
      xyFrame1<- data.frame(x = x$x, y = x$y)
    }
    
    if(dataType == "polygon" || dataType == "Polygon" || dataType == "POLYGON"){
      
      if(length(poly.xy) > 0){
        if(length(poly.xy) == (numVertices*2) & length(which(is.na(match(poly.xy, names(x)))) == TRUE) == 0){ 
          xyFrame1<-x[,match(poly.xy,names(x))] 
        }else{
          xyFrame1<- data.frame(poly.xy)
        }
      }else{ #if length(poly.xy == 0)
        xyFrame1 <-  x[,(match("point1.x", names(x))):((2*numVertices) + (match("point1.x", names(x)) -1))] #if there is no poly.xy input, the code assumes the input is output from referencePointToPolygon function, and therefore, the first point of interest would be "point1.x"
      }
      colnames(xyFrame1)[seq(1,(numVertices*2),2)] = paste("point",seq(1,numVertices,1),".x", sep = "")
      colnames(xyFrame1)[seq(2,((numVertices*2) + 1),2)] = paste("point",seq(1,numVertices,1),".y", sep = "")
    }
    
    if(length(dateTime) > 0){
      if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
        x$dateTime <- x[,match(dateTime, names(x))]
      }else{ #if length(dateTime) > 1
        x$dateTime = dateTime
      }
    }
    dateTimeVec<-x$dateTime
    bindlist1<-list(idVec1, xyFrame1, dateTimeVec)
    originTab <- do.call("cbind", bindlist1)
    names(originTab)[c(1,ncol(originTab))]<-c("id", "dateTime")
  }
  originTab$dateTime = as.character(originTab$dateTime)
  
  if(blocking == TRUE){
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
    originTab<-originTab[order(originTab$dateTime),] #Just in case the data wasn't already ordered in this way.
    originTab<-datetime.append1(originTab) #adds the total second column to the dataframe
    studySecond <- (originTab$totalSecond -min(originTab$totalSecond)) + 1
    originTab<-originTab[,-match("totalSecond", names(originTab))]
    numblocks <- ceiling((max(studySecond) - 1)/blockLength1)
    block <-rep(0,length(studySecond))
    for(g in 1:(numblocks -1)){ #numblocks - 1 because the last block in the dataset may be smaller than previous blocks (if blockLength1 does not divide evenly into timedif)
      block[which(studySecond >= ((g-1)*blockLength1 + 1) & studySecond <= (g*blockLength1))] = g
    }
    if(length(which(block == 0)) > 0){ #identifies the last block
      block[which(block == 0)] = numblocks
    }
    blockVec<-unique(block)
    dateTimeVec2<-originTab$dateTime
    minBlockTimeSeq <- rep(0, length(block)) #Added 2/4/2019. This vector will identify the minimum timepoint in each block.
    maxBlockTimeSeq <- rep(0, length(block)) #Added 2/4/2019. This vector will identify the maximum timepoint in each block.
    for(f in 1:length(blockVec)){
      minBlockTime<-dateTimeVec2[min(which(block == blockVec[f]))]
      minBlockTimeSeq[which(block == blockVec[f])] <- minBlockTime
      maxBlockTime<-dateTimeVec2[max(which(block == blockVec[f]))]
      maxBlockTimeSeq[which(block == blockVec[f])] <- maxBlockTime
    }
    originTab$block <- block
    originTab$block.start <- minBlockTimeSeq
    originTab$block.end <- maxBlockTimeSeq
    originTab$numBlocks <- max(blockVec) 
    originTab$blockLength <- paste(blockLength, blockUnit, sep = " ")
  }
  originTab <- originTab[order(originTab$id, originTab$dateTime),] #order by id, then dateTime
  
  if(numRandomizations > 1){
    repeatTab<- data.frame(seq(1,numRandomizations, 1))
    repeatTab$dataType = dataType
    repeatTab$numVertices = numVertices
    repeatTab$shuffle.type = shuffle.type
    repeatTab$indivPaths = indivPaths
    repeatTab$blocking = blocking
    repeatTab$subUnit = shuffleUnit
    repeatTab$blockUnit = blockUnit
    repeatTab$blockLength = blockLength
    
    if(parallel == TRUE){
      cl<-parallel::makeCluster(nCores)
      randomness<-parallel::parApply(cl, repeatTab, 1, randomization, originTab)
      parallel::stopCluster(cl)
    }else{
      randomness <- apply(repeatTab, 1, randomization, originTab)
    }
    return(randomness)
    
  }else{ #if numRandomizations == 1
    if(indivPaths == TRUE){
      idVecSeq <- unique(originTab$id)
    }else{
      idVecSeq <- 1
    }
    numids <- length(idVecSeq)
    idTab<- data.frame(idVecSeq)
    
    if(parallel == TRUE){
      cl<-parallel::makeCluster(nCores)
      randomness<-parallel::parApply(cl, idTab,1,datasub.func, originTab, numids,blocking, shuffle.type, dataType, shuffleUnit, blockUnit, blockLength)
      parallel::stopCluster(cl)
    }else{
      randomness <- apply(idTab,1,datasub.func, originTab, numids,blocking, shuffle.type, dataType, shuffleUnit, blockUnit, blockLength)
    }
    
    if(is.data.frame(randomness) == TRUE){
      randomxyCoords<-randomness
    }else{
      randomxyCoords<-data.frame(data.table::rbindlist(randomness))
    }
    
    if(blocking == TRUE & shuffle.type == 2){ #To account for potential length differences in blocks if shuffle.type ==2, the outputTab is compiled earlier on in the function (in the datasub.func subfunction) 
      outputTab <- randomxyCoords
    }else{
      outputTab<- cbind(originTab, randomxyCoords)
    }
    if(blocking == TRUE){
      outputTab$shuffle.type = shuffle.type #added 1/30/2019, to ensure that people remember what shuffle.type they set.
      if(shuffle.type == 2){
        outputTab$shuffleUnit <- shuffleUnit
      }
    }
    outputTab$rand.rep = 1    
    return(outputTab)
  }
}
