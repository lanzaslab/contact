#' Move Data Point a Specified Distance
#'
#' Translates locations of a single rfid tag/gps transmitter to a different location on the equipped animal (i.e., move the coordinate up, down, left, right, or a combination thereof (i.e., up & right, up & left, down & right, down & left)) . For example, calves in our study (the test dataset) were equiped with RFID tags on their left ear. With this code, we can move this reference point somewhere else on the body of each individual. This might be done for a number of reasons, but is very useful for use in the referencePointToPolygon function later on (for delineating polygons representing entire individuals). Currently, this function only supports input data with coordinates representing planar ('Euclidean') space (e.g. units of meters).
#'
#' In this function, if the distance individuals moved was less than/equal to the noted immobThreshold, individuals are said to be "immob," and their position will not change relative to their previous one. (i.e., you assume that any observed movement less than immobThreshold was due to errors or miniscule bodily movements (e.g., head shaking) that are not indicative of actual movement.)
#'
#' Because this function needs information (dist, dx, dy) from 2 points on an individual's path to work, at least the first point in each individual's path will be removed (the function will report NAs for adjusted locations). Also note that if the distance between an individual's first point in their path and the second one is 0, the function will also report NAs for the second point's adjusted coordinates. The first non-NA values will only be reported for the instance where dist > 0.
#'
#' In the output, if input was previously processed using tempAggregate with resolutionLevel == "reduced," dt > secondAgg indicates that tracked individuals were missing in the original dataset for a period of time. In this case, the assumption that individuals are facing a given direction because they moved from the previous timepoint may not be accurate. Consider removing these rows (rows following one with dt > secondAgg; remember that dt indicates the time between recording xy coordinates in row i to row i + 1) from your dataset.
#' @param x Data frame or list of data frames containing real-time-location point data.  
#' @param id Vector of length nrow(data.frame(x)) or singular character data, detailing the relevant colname in x, that denotes what unique ids for tracked individuals will be used. If argument == NULL, the function assumes a column with the colname "id" exists in x. Defaults to NULL.
#' @param point.x Vector of length nrow(data.frame(x)) or singular character data, detailing the relevant colname in x, that denotes what planar-x or longitude coordinate information will be used. If argument == NULL, the function assumes a column with the colname "x" exists in x. Defaults to NULL.
#' @param point.y Vector of length nrow(data.frame(x)) or singular character data, detailing the relevant colname in x, that denotes what planar-y or lattitude coordinate information will be used. If argument == NULL, the function assumes a column with the colname "y" exists in x. Defaults to NULL.
#' @param dateTime Vector of length nrow(data.frame(x)) or singular character data, detailing the relevant colname in x, that denotes what dateTime information will be used. If argument == NULL, the function assumes a column with the colname "dateTime" exists in x. Defaults to NULL.
#' @param RepositionDir Character string taking the values "UP," "DOWN," "LEFT," "RIGHT," "UP&RIGHT," "UP&LEFT," "DOWN&RIGHT," or "DOWN&LEFT." Describes the direction(s) that points will be moved. 
#' @param UpDownRepositionLen Numerical. Describes the height, in planar units (e.g., meters) of the output polygon. Planar units are inherent in the real-time-location input. Defaults to 1.
#' @param LeftRightRepositionLen Numerical. Describes the width, in planar units (e.g., meters) of the output polygon. Planar units are inherent in the real-time-location input. Defaults to 1.
#' @param immobThreshold Numerical. Describes what we call, the immobility threshold, which is a movement distance (in planar units) within which we assume individuals’ physical locations and orientations remain unchanged. This immobility threshold allows us to discount observed movements so miniscule that the majority of animals’ physical-space usage is likely unaffected (e.g., head shaking). Defaults to 0.
#' @param parallel Logical. If TRUE, sub-functions within the repositionReferencePoint wrapper will be parallelized. Note that this can significantly speed up processing of relatively small data sets, but may cause R to crash due to lack of available memory when attempting to process large datasets. Defaults to TRUE.
#' @keywords data-processing location point planar
#' @export
#' @examples
#' #read in the calves data set
#' data("calves")
#' calves.dateTime<-datetime.append(calves, date = calves$date, time = calves$time) #create a dataframe with dateTime identifiers for location fixes.
#' 
#' #create our data set that shows calves average position every 10 seconds
#' system.time(calves.10secSmoothed <- tempAggregate(x = calves, id = calves$calftag, point.x = calves$x, point.y = calves$y, dateTime = calves$dateTime, secondAgg = 10, extrapolate.left = TRUE, resolutionLevel = "Full", extrapolate.right = FALSE, na.rm = TRUE, smooth.type = 2))
#' 
#' ##Create 0.333 m X 0.333 m calf head polygons.
#' #Note that this is done using the original reference points which denote the locations of RFID tags on individuals' left ears.
#' system.time(calf_heads <- referencePointToPolygon(x = calves.10secSmoothed, id = calves.10secSmoothed$id, dateTime = calves.10secSmoothed$dateTime, point.x = calves.10secSmoothed$x, point.y = calves.10secSmoothed$y, StartLocation = "DL", UpDownRepositionLen = 0.333, LeftRightRepositionLen = 0.333, CenterPoint = FALSE, MidPoints = FALSE, immobThreshold  = 0.1, parallel = TRUE))
#' 
#' #Because the head is not the same width of the body and is assumed to be centered at the front of the body, before creating body polygons, we must move reference points (on the left ear) to the left by 0.3335 m to reposition them at the upper-left corner of calves bodies. Note that we are assuming ears are parallel to shoulder-tips. 
#' system.time(leftShoulder.point<-repositionReferencePoint(x = calves.10secSmoothed, id = calves.10secSmoothed$id, dateTime = calves.10secSmoothed$dateTime, point.x = calves.10secSmoothed$x, point.y = calves.10secSmoothed$y, RepositionDir = "L", UpDownRepositionLen = 0, LeftRightRepositionLen = 0.3335, immobThreshold  = 0, parallel = TRUE)) #Note that we do not specify a standing threshold here. Rather, we will do so when we create the polygon.
#' 
#' #Now we can create generate the vertices for anterior- and posterior-body polygons. Rather than running the referencePointToPolygon function twice, we instead set MidPoints = TRUE, which will effectively identify vertices for the bottom of anterior bodies/top of posterior ones. 
#' system.time(calf_bods <- referencePointToPolygon(x = leftShoulder.point, id = leftShoulder.point$id, dateTime = leftShoulder.point$dateTime, point.x = leftShoulder.point$x.adjusted, point.y = leftShoulder.point$y.adjusted, StartLocation = "UL", UpDownRepositionLen = 2, LeftRightRepositionLen = 1, CenterPoint = FALSE, MidPoints = TRUE, immobThreshold  = 0.1, parallel = TRUE))

repositionReferencePoint <- function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, RepositionDir = "Up", UpDownRepositionLen = 1, LeftRightRepositionLen = 1, immobThreshold = 0, parallel = TRUE){

  reposition.generator<-function(x, id, dateTime, point.x, point.y, RepositionDir, UpDownRepositionLen, LeftRightRepositionLen, immobThreshold, parallel){
    euc=function(x) {
      point1 = x.cor=unlist(c(x[1],x[2]))
      point2 = x.cor=unlist(c(x[3],x[4]))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    translation <-function(x){

      dx = as.numeric(x[4])
      dy = as.numeric(x[5])
      dist = as.numeric(x[3])

      if (dist == 0 || is.na(x[1]) == TRUE || is.na(x[2]) == TRUE || is.na(dx) == TRUE || is.na(dy) == TRUE || is.na(dist) == TRUE){ #if individuals did not move or any of the necessary variables are undefined, it reports an NA.
        newX = NA
        newY = NA}else{

          if (x[6] == "Up" || x[6] == "UP" || x[6] == "up" || x[6] == "u" || x[6] == "U"){ #If you want to reposition the reference point at a part of the animal/object ahead of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              newX = as.numeric(x[1]) + as.numeric(x[7])
              newY = as.numeric(x[2])
            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              newX = as.numeric(x[1]) - as.numeric(x[7])
              newY = as.numeric(x[2])
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) + as.numeric(x[7])
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) - as.numeric(x[7])
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 2 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              newX = (as.numeric(x[1]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              newY = (as.numeric(x[2]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              newX = (as.numeric(x[1]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              newY = (as.numeric(x[2]) - (adj2 - adj1)) #We know that the animal is moving downward. This reduces as.numeric(x[2]) by the difference between adj2 and adj 1.
            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              newX = (as.numeric(x[1]) - (opp2 - opp1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between opp2 and opp1.
              newY = (as.numeric(x[2]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.
            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              newX = (as.numeric(x[1]) - (adj2 - adj1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between adj2 and adj1.
              newY = (as.numeric(x[2]) - (opp2 - opp1)) #We know that the individual is moving downwards (towards the negative side of the y axis). This reduces as.numeric(x[2]) by the difference between opp2 and opp1.
            }

          }

          if (x[6] == "Down" || x[6] == "DOWN" || x[6] == "down" || x[6] == "d" || x[6] == "D"){ #If you want to reposition the reference point at a part of the animal/object behind of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              newX = as.numeric(x[1]) - as.numeric(x[7])
              newY = as.numeric(x[2])
            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              newX = as.numeric(x[1]) + as.numeric(x[7])
              newY = as.numeric(x[2])
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) - as.numeric(x[7])
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) + as.numeric(x[7])
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 2 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              newX = (as.numeric(x[1]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              newY = (as.numeric(x[2]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              newX = (as.numeric(x[1]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              newY = (as.numeric(x[2]) - (adj2 - adj1)) #We know that the animal is moving downward. This reduces as.numeric(x[2]) by the difference between adj2 and adj 1.
            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              newX = (as.numeric(x[1]) - (opp2 - opp1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between opp2 and opp1.
              newY = (as.numeric(x[2]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.
            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              newX = (as.numeric(x[1]) - (adj2 - adj1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between adj2 and adj1.
              newY = (as.numeric(x[2]) - (opp2 - opp1)) #We know that the individual is moving downwards (towards the negative side of the y axis). This reduces as.numeric(x[2]) by the difference between opp2 and opp1.
            }
          }

          if (x[6] == "Left" || x[6] == "LEFT" || x[6] == "left" || x[6] == "l" || x[6] == "L"){ #If you want to reposition the reference point at a part of the animal/object to the left of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) + as.numeric(x[8])
            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) - as.numeric(x[8])
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              newX = as.numeric(x[1]) - as.numeric(x[8])
              newY = as.numeric(x[2])
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              newX = as.numeric(x[1]) + as.numeric(x[8])
              newY = as.numeric(x[2])
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 2 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned higher than the one described above, but shares a vertex at the original reference coordinates. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) - opp2)
              newY = (as.numeric(x[2]) + adj2)
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the right of the one described above, but shares a vertex at the original reference coordinates. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) + adj2)
              newY = (as.numeric(x[2]) + opp2)
            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the x axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the y axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the left of the one described above, but shares a vertex at the original reference coordinates. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) - adj2)
              newY = (as.numeric(x[2]) - opp2)
            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned below the one described above, but shares a vertex at the original reference coordinates. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) + opp2)
              newY = (as.numeric(x[2]) - adj2)
            }
          }

          if (x[6] == "Right" || x[6] == "RIGHT" || x[6] == "right" || x[6] == "r" || x[6] == "R"){ #If you want to reposition the reference point at a part of the animal/object to the right of the RFID tag/GPS locator #Note that the calculations are the same as for left translation, but the effects of adjustments on xy coordinates are reversed.

            if(dx > 0 & dy == 0){ #if individuals moved right
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) - as.numeric(x[8])
            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              newX = as.numeric(x[1])
              newY = as.numeric(x[2]) + as.numeric(x[8])
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              newX = as.numeric(x[1]) + as.numeric(x[8])
              newY = as.numeric(x[2])
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              newX = as.numeric(x[1]) - as.numeric(x[8])
              newY = as.numeric(x[2])
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 2 right triangles.
            if(dx > 0 & dy > 0){ #if animals moved up & right
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the right of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) + opp2)
              newY = (as.numeric(x[2]) - adj2)
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the below the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) - adj2)
              newY = (as.numeric(x[2]) - opp2)

            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned above the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) + adj2)
              newY = (as.numeric(x[2]) + opp2)
            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the left of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (as.numeric(x[1]) - opp2)
              newY = (as.numeric(x[2]) + adj2)
            }
          }

          if (x[6] == "Up&Right" || x[6] == "UP&RIGHT" || x[6] == "up&right" || x[6] == "Up & Right" || x[6] == "UP & RIGHT" || x[6] == "up & right" || x[6] == "u&r" || x[6] == "ur" || x[6] == "U&R" || x[6] == "UR"){ #If you want to reposition the reference point at a part of the animal/object ahead of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              intermediateX = as.numeric(x[1]) + as.numeric(x[7]) #moves point up
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY - as.numeric(x[8]) #moves point right
            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              intermediateX = as.numeric(x[1]) - as.numeric(x[7])
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY + as.numeric(x[8]) #moves point right

            }

            if(dx == 0 & dy > 0){ #if animals moved up
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) + as.numeric(x[7])
              newX = intermediateX + as.numeric(x[8]) #moves point right
              newY = intermediateY
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) - as.numeric(x[7])
              newX = intermediateX - as.numeric(x[8]) #moves point right
              newY = intermediateY
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 4 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              ##First we need to move the point up

              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              intermediateY = (as.numeric(x[2]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.

              ##Now we need to move the point right

              xySet3 = c(intermediateX,intermediateY,intermediateX,0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet4 = c(intermediateX,intermediateY,0,intermediateY) #This sets the x axis point of contact
              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the right of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2

              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + opp2)
              newY = (intermediateY - adj2)
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              intermediateY = (as.numeric(x[2]) - (adj2 - adj1)) #We know that the animal is moving downward. This reduces as.numeric(x[2]) by the difference between adj2 and adj 1.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the below the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - adj2)
              newY = (intermediateY - opp2)

            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (opp2 - opp1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between opp2 and opp1.
              intermediateY = (as.numeric(x[2]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned above the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + adj2)
              newY = (intermediateY + opp2)
            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (adj2 - adj1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between adj2 and adj1.
              intermediateY = (as.numeric(x[2]) - (opp2 - opp1)) #We know that the individual is moving downwards (towards the negative side of the y axis). This reduces as.numeric(x[2]) by the difference between opp2 and opp1.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the left of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - opp2)
              newY = (intermediateY + adj2)
            }
          }

          if (x[6] == "Up&Left" || x[6] == "UP&LEFT" || x[6] == "up&left" || x[6] == "Up & Left" || x[6] == "UP & LEFT" || x[6] == "up & left" || x[6] == "u&l" || x[6] == "ul" || x[6] == "U&L" || x[6] == "UL"){ #If you want to reposition the reference point at a part of the animal/object ahead of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              intermediateX = as.numeric(x[1]) + as.numeric(x[7]) #moves point up
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY + as.numeric(x[8]) #moves point left
            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              intermediateX = as.numeric(x[1]) - as.numeric(x[7])
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY - as.numeric(x[8]) #moves point left
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) + as.numeric(x[7])
              newX = intermediateX - as.numeric(x[8]) #moves point left
              newY = intermediateY
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) - as.numeric(x[7])
              newX = intermediateX + as.numeric(x[8]) #moves point left
              newY = intermediateY
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 4 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              ##First we need to move the point up

              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              #theta = theta*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              intermediateY = (as.numeric(x[2]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.

              ##Now we need to move the point left

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the right of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2

              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - opp2)
              newY = (intermediateY + adj2)
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (opp2 - opp1))
              intermediateY = (as.numeric(x[2]) - (adj2 - adj1))

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the below the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + adj2)
              newY = (intermediateY + opp2)
            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).
              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (opp2 - opp1))
              intermediateY = (as.numeric(x[2]) + (adj2 - adj1))

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned above the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - adj2)
              newY = (intermediateY - opp2)

            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              #theta = theta*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              hyp2 = hyp1 + as.numeric(x[7]) #The upwards translation essentially amounts to an extension of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (adj2 - adj1))
              intermediateY = (as.numeric(x[2]) - (opp2 - opp1))

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the left of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + opp2)
              newY = (intermediateY - adj2)
            }
          }

          if (x[6] == "Down&Right" || x[6] == "DOWN&RIGHT" || x[6] == "down&right" || x[6] == "Down & Right" || x[6] == "DOWN & RIGHT" || x[6] == "down & right" || x[6] == "d&r" || x[6] == "dr" || x[6] == "D&R" || x[6] == "DR"){ #If you want to reposition the reference point at a part of the animal/object behind of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              intermediateX = as.numeric(x[1]) - as.numeric(x[7])
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY - as.numeric(x[8]) #moves point right

            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              intermediateX = as.numeric(x[1]) + as.numeric(x[7])
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY + as.numeric(x[8]) #moves point right
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) - as.numeric(x[7])
              newX = intermediateX + as.numeric(x[8]) #moves point right
              newY = intermediateY
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) + as.numeric(x[7])
              newX = intermediateX - as.numeric(x[8]) #moves point right
              newY = intermediateY
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 4 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              #theta = theta*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (adj2 - adj1))
              intermediateY = (as.numeric(x[2]) + (opp2 - opp1))

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the right of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2

              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + opp2)
              newY = (intermediateY - adj2)
            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (opp2 - opp1))
              intermediateY = (as.numeric(x[2]) - (adj2 - adj1))

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1
              #Now we need to solve a right triangle that is positioned to the below the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - adj2)
              newY = (intermediateY - opp2)

            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.
              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.
              intermediateX = (as.numeric(x[1]) - (opp2 - opp1))
              intermediateY = (as.numeric(x[2]) + (adj2 - adj1))

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned above the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + adj2)
              newY = (intermediateY + opp2)

            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (adj2 - adj1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between adj2 and adj1.
              intermediateY = (as.numeric(x[2]) - (opp2 - opp1)) #We know that the individual is moving downwards (towards the negative side of the y axis). This reduces as.numeric(x[2]) by the difference between opp2 and opp1.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the left of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - opp2)
              newY = (intermediateY + adj2)

            }

          }

          if (x[6] == "Down&Left" || x[6] == "DOWN&LEFT" || x[6] == "down&left" || x[6] == "Down & Left" || x[6] == "DOWN & LEFT" || x[6] == "down & left" || x[6] == "d&l" || x[6] == "dl" || x[6] == "D&L" || x[6] == "DL"){ #If you want to reposition the reference point at a part of the animal/object behind of the RFID tag/GPS locator

            if(dx > 0 & dy == 0){ #if individuals moved right
              intermediateX = as.numeric(x[1]) - as.numeric(x[7])
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY + as.numeric(x[8]) #moves point right

            }

            if(dx < 0 & dy == 0){ #if individuals moved left
              intermediateX = as.numeric(x[1]) + as.numeric(x[7])
              intermediateY = as.numeric(x[2])
              newX = intermediateX
              newY = intermediateY - as.numeric(x[8]) #moves point right
            }

            if(dx == 0 & dy > 0){ #if animals moved up
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) - as.numeric(x[7])
              newX = intermediateX - as.numeric(x[8]) #moves point right
              newY = intermediateY
            }

            if(dx == 0 & dy < 0){ #if animals moved down
              intermediateX = as.numeric(x[1])
              intermediateY = as.numeric(x[2]) + as.numeric(x[7])
              newX = intermediateX + as.numeric(x[8]) #moves point right
              newY = intermediateY
            }

            ##to ascertain new reference points following diagonal movements, we need to solve 4 right triangles

            if(dx > 0 & dy > 0){ #if animals moved up & right
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.

              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              intermediateY = (as.numeric(x[2]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the right of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2

              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - opp2)
              newY = (intermediateY + adj2)

            }

            if(dx > 0 & dy < 0){ #if animals moved down & right #Note that the 2 differences between this code and the one above (where individuals move up and to the right), are that the relationships between the xy coordinates and the opposite and adjacent sides of the triangles have switched, and y must decrease with adjustment, not increase.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.

              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) + (opp2 - opp1)) #In most cases, this is the same as 0 + opp2. However, in the case that as.numeric(x[1]) is a negative number, the function will still produce accurate results.
              intermediateY = (as.numeric(x[2]) - (adj2 - adj1)) #We know that the animal is moving downward. This reduces as.numeric(x[2]) by the difference between adj2 and adj 1.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the below the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + adj2)
              newY = (intermediateY + opp2)

            }

            if(dx < 0 & dy > 0){ #if animals moved up & left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move down and to the right, but y increases while x decreases.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0)  #This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              adj1 = euc(xySet1) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet2) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.

              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (opp2 - opp1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between opp2 and opp1.
              intermediateY = (as.numeric(x[2]) + (adj2 - adj1)) #In most cases, this is the same as 0 + adj2. However, in the case that as.numeric(x[2]) is a negative number, the function will still produce accurate results.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              adj1 = euc(xySet3) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              opp1 = euc(xySet4) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned above the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (90 - adjAngle)/2
              theta2 = (90 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX - adj2)
              newY = (intermediateY - opp2)


            }

            if(dx < 0 & dy < 0){ #if animals moved down and left #Note that the the relationships between the xy coordinates and the opposite and adjacent sides of the triangles are the same as when individuals move up and to the right, but both y and x decrease.
              xySet1 = c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[1]),0) ##This sets the y axis point of contact as (as.numeric(x[1]), (as.numeric(x[2]) - as.numeric(x[2]))).
              xySet2 = c(as.numeric(x[1]),as.numeric(x[2]),0,as.numeric(x[2])) #This sets the x axis point of contact as ((as.numeric(x[1]) - as.numeric(x[1])), as.numeric(x[2])).

              opp1 = euc(xySet1) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet2) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.

              hyp2 = hyp1 - as.numeric(x[7]) #The downwards translation essentially amounts to an reduction of the triangle's hypotenuse
              opp2 = (hyp2/sin(90*(pi/180)))*sin(theta) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must conert 90 degrees to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              intermediateX = (as.numeric(x[1]) - (adj2 - adj1)) #We know that the individual is moving to the left (towards the negative side of the x axis). This reduces as.numeric(x[1]) by the difference between adj2 and adj1.
              intermediateY = (as.numeric(x[2]) - (opp2 - opp1)) #We know that the individual is moving downwards (towards the negative side of the y axis). This reduces as.numeric(x[2]) by the difference between opp2 and opp1.

              xySet3 = c(intermediateX,intermediateY,intermediateX,0)
              xySet4 = c(intermediateX,intermediateY,0,intermediateY)

              opp1 = euc(xySet3) # To determine the length of the triangle's opposite side, we use the euc function to find the length between the original reference point (xySet1[1:2]) and the x axis (xySet1[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the bottom of the plane (i.e., there are no negative y values), the length of the opposite side will allways equal the reference point's y coordinate.
              adj1 = euc(xySet4) # To determine the length of the triangle's adjacent side, we use the euc function to find the length between the original reference point (xySet2[1:2]) and the y axis (xySet2[3:4]). Note that if the plane's origin point (x = 0, y = 0) is at the leftmost side of the plane (i.e., there are no negative x values), the length of the adjacent side will allways equal the reference point's x coordinate.
              hyp1 =  sqrt((adj1^2) + (opp1^2)) #The triangle's hypotenuse folows the Pythagorean theorum, where a^2 and b^2 = adj1^2 and opp1^2, respectively.

              theta1 = asin(opp1/hyp1) #Sine = opposite/hypotenuse, arcsine is the inverse sine function and can be used to determine the angle across from the opposite side.
              theta1 = theta1*(180/pi) #The asin function returns a value in radians, this code converts that to degrees.
              adjAngle = 90 - theta1

              #Now we need to solve a right triangle that is positioned to the left of the one described above, but shares two vertices with the original. Theta (in degrees) of this new triangle is equal to (180 - adjAngle)/2
              theta2 = (180 - adjAngle)/2
              hyp2 = as.numeric(x[8]) #the length of the hypotenuse of the new triangle is the adjustment distance
              opp2 = (hyp2/sin(90*(pi/180)))*sin((theta2*(pi/180))) #The new opposite side length is calculated according to the law of sines. #Note that because the sin function expects an argument in radian form, we must convert degree values to radians (i.e., multiply it by pi/180)
              adj2 = sqrt((hyp2^2) - (opp2^2)) #The triangle's adjacent side folows the Pythagorean theorum, where c^2 and b^2 = hyp2^2 and opp1^2, respectively.

              #In contrast to translating up or down, opposite and adjacent sides do NOT represent xy coordinates (assuming the plane's origin is in the lower left corner). Rather, they represent necessary adjustments to xy coordinates
              newX = (intermediateX + opp2)
              newY = (intermediateY - adj2)
            }
          }
        }

      CoordAdjusted = data.frame(x.adjusted = newX, y.adjusted = newY)
      return(CoordAdjusted)
    }
    immobAdjustment.point = function(x, locMatrix, immobVec){
      dur = 0
      standDurTest <- FALSE
      while (!standDurTest) { #This means when standDurTest == FALSE #This part of the function determines how long an individual has been immob.
        standDurTest <- immobVec[(x[1] - (1 + dur))] == 0
        standDurTest <- ifelse(is.na(standDurTest) == TRUE, TRUE, standDurTest) #This is for when x[1] -1 relates to the first observation for any individual (because those observations are NA)
        dur = (dur + 1)
      }
      standFrame = data.frame(replaceRow = x[1], replacePoint.x = locMatrix$x.adjusted[(x[1] - dur)], replacePoint.y = locMatrix$y.adjusted[(x[1] - dur)]) #if individuals are immob, their location at a given timestep will be the same as at the previous timestep.
      return(standFrame)
    }
    timeDifference = function(x){
      t1 = x[1]
      t2 = x[2]
      dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
      return(dt)
    }

    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors.
      x <- data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime)
    }

    if(length(x) > 0){ #for some reason using an "else" statement would always result in a table with 0 records...

      if(length(id) > 0){
        if(length(id) == 1 && is.na(match(id[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than id being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "id" values (i.e., if the ids in different list entries are different)
          x$id <- x[,match(id, names(x))]
        }else{ #if length(id) > 1
          x$id = id
        }

      }
      if(length(point.x) > 0){
        if(length(point.x) == 1 && is.na(match(point.x[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$x <- x[,match(point.x, names(x))]
        }else{ #if length(point.x) > 1
          x$x = point.x
        }
      }
      if(length(point.y) > 0){
        if(length(point.y) == 1 && is.na(match(point.y[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$y <- x[,match(point.y, names(x))]
        }else{ #if length(point.y) > 1
          x$y = point.y
        }
      }
      if(length(dateTime) > 0){
        if(length(dateTime) == 1 && is.na(match(dateTime[1], names(x))) == FALSE){ #added 1/14 to accompany the list-processing functionality. If x is a list, rather than point.x being a vector of length(nrow(x)), it may be necessary to designate the colname for intended "point.x" values (i.e., if the x-coordinate values in different list entries are different)
          x$dateTime <- x[,match(dateTime, names(x))]
        }else{ #if length(dateTime) > 1
          x$dateTime = dateTime
        }
      }
    }

    x<-x[order(x$id, x$dateTime),] #data tables should already be ordered like this, but because this function assumes this ordering scheme, I wanted to ensure it was so.
    rownames(x) <- seq(1, nrow(x),1)

    distCoordinates = data.frame(x$x[1:(nrow(x) - 1)], x$y[1:(nrow(x) - 1)], x$x[2:nrow(x)], x$y[2:nrow(x)])
    if(parallel == TRUE){ #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
      cl<-parallel::makeCluster(parallel::detectCores())
      dist = parallel::parApply(cl, distCoordinates,1,euc)
      dist = c(dist, NA) #to make dist the same length as nrow(x)
      parallel::stopCluster(cl)
    }else{
      dist = apply(distCoordinates,1,euc)
      dist = c(dist, NA) #to make dist the same length as nrow(x)
    }

    dx = c((distCoordinates[,1] - distCoordinates[,3])*(-1),NA)
    dy = c((distCoordinates[,2] - distCoordinates[,4])*(-1),NA)

    timesFrame = data.frame(x$dateTime[1:(nrow(x) - 1)], x$dateTime[2:nrow(x)])
    if(parallel == TRUE){
      cl<-parallel::makeCluster(parallel::detectCores())
      dt = parallel::parApply(cl, timesFrame, 1, timeDifference)
      dt = c(dt, NA) #to make dt the same length as nrow(x)
      parallel::stopCluster(cl)
    }else{
      dt = apply(timesFrame, 1, timeDifference)
      dt = c(dt, NA) #to make dt the same length as nrow(x)
    }

    idVec = unique(x$id) ; idSeqVec = x$id; idVecForDtRemoval = idVec[-length(idVec)] #There's no need to remove the max point of the final id point b/c there's already an NA there.
    for(a in idVecForDtRemoval){ #This loop removes the dt value at points describing the last point in each individual's path
      dx[max(which(idSeqVec == a))] = NA
      dy[max(which(idSeqVec == a))] = NA
      dist[max(which(idSeqVec == a))] = NA
      dt[max(which(idSeqVec == a))] = NA
    }
    x$dx = dx
    x$dy = dy
    x$dist = dist
    x$dt = dt

    dataShift <- data.frame(dist = c(NA,x$dist[1:(nrow(x) - 1)]), dx = c(NA,x$dx[1:(nrow(x) - 1)]), dy = c(NA,x$dy[1:(nrow(x) - 1)])) #This is necessary because of the way we calculated/listed these values above. These columns (dist, dx, and dt) refer to the changes a point must make to reach the subsequent point. However, later on in this function, we are not interested in future point alterations. Rather, we need to know how tracked individuals moved during the preceding time step to reach their current point (to determine directionality of movement) #Note that this code shifts values down, but remember that beacuse values are downshifted, the first observation for each id individual will be incorrect. Code below addresses this issue.
    for(b in idVec){ #This code replaces dist, dx, and dy values in the row when a new individual (id) is first observed with NAs to fix the problem noted above. #Note that this loop assumes that the data is sorted by id (i.e., a given id will not repeat in the dataset, after a new id has appeared)
      dataShift[min(which(idSeqVec == b)),] = NA
    }
    xyFrame <- data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, RepositionDir, UpDownRepositionLen, LeftRightRepositionLen)

    if (parallel == TRUE){
      cl<-parallel::makeCluster(parallel::detectCores())
      xyTranslated<-parallel::parApply(cl, xyFrame, 1, translation)
      newCoordinates <- data.frame(data.table::rbindlist(xyTranslated))
      parallel::stopCluster(cl)
    }else{
      xyTranslated <- apply(xyFrame, 1, translation)
      newCoordinates <- data.frame(data.table::rbindlist(xyTranslated))
    }

    RepositionMatrix <- data.frame(matrix(nrow = nrow(x), ncol = (18))) #One row for each observation, Columns listed below
    colnames(RepositionMatrix) <- c("id","x.original","y.original","dist.original", "dx.original", "dy.original","x.adjusted","y.adjusted", "dist.adjusted", "dx.adjusted", "dy.adjusted", "adjustmentDirection","upDownAdjustmentLength","leftRightAdjustmentLength","immob","immobThreshold", "dateTime","dt")
    RepositionMatrix$id = x$id
    RepositionMatrix$x.original = x$x
    RepositionMatrix$y.original = x$y
    RepositionMatrix$dist.original = x$dist
    RepositionMatrix$dx.original = x$dx
    RepositionMatrix$dy.original = x$dy
    RepositionMatrix$adjustmentDirection = RepositionDir
    RepositionMatrix$upDownAdjustmentLength = UpDownRepositionLen
    RepositionMatrix$leftRightAdjustmentLength = LeftRightRepositionLen
    RepositionMatrix$immobThreshold = immobThreshold
    RepositionMatrix$dateTime = x$dateTime
    RepositionMatrix$dt = x$dt
    RepositionMatrix$x.adjusted = newCoordinates$x.adjusted
    RepositionMatrix$y.adjusted = newCoordinates$y.adjusted
    RepositionMatrix$immob = ifelse(xyFrame[,3] > immobThreshold, 0, 1) #if the distance individuals moved was less than / equal to the noted immobThreshold, individuals are said to be "immob," and their position will not change relative to their previous one. (i.e., you assume that any observed movement less than immobThreshold was due to errors or miniscule bodily movements (e.g., head shaking) that are not indicative of actual movement.)

    standlist = which(RepositionMatrix$immob == 1)

    if ((length(standlist) >= 1)){ #To save processing time and reduce chances of errors, this evaluation will not take place if there is only one observation.

      immobVec = RepositionMatrix$immob
      standlistFrame = data.frame(immob = standlist)

      if(parallel == TRUE){
        cl<-parallel::makeCluster(parallel::detectCores())
        immobFrame = data.frame(data.table::rbindlist(parallel::parApply(cl, standlistFrame, 1, immobAdjustment.point, locMatrix = RepositionMatrix, immobVec)))
        parallel::stopCluster(cl)
      }else{
        immobFrame = data.frame(data.table::rbindlist(apply(standlistFrame, 1, immobAdjustment.point, locMatrix = RepositionMatrix, immobVec)))
      }
      if(nrow(immobFrame) > 0){
        RepositionMatrix[immobFrame$replaceRow,7] = immobFrame$replacePoint.x
        RepositionMatrix[immobFrame$replaceRow,8] = immobFrame$replacePoint.y
      }
    }

    newDistCoordinates = data.frame(RepositionMatrix$x.adjusted[1:(nrow(RepositionMatrix) - 1)], RepositionMatrix$y.adjusted[1:(nrow(RepositionMatrix) - 1)], RepositionMatrix$x.adjusted[2:nrow(RepositionMatrix)], RepositionMatrix$y.adjusted[2:nrow(RepositionMatrix)])

    if(parallel == TRUE){ #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
      cl<-parallel::makeCluster(parallel::detectCores())
      newDist = parallel::parApply(cl, newDistCoordinates,1,euc); newDist = c(newDist, NA)
      parallel::stopCluster(cl)
    }else{
      newDist = apply(newDistCoordinates,1,euc); newDist = c(newDist, NA)
    }

    RepositionMatrix$dist.adjusted = newDist

    newdx = c((newDistCoordinates[,1] - newDistCoordinates[,3])*(-1),NA)
    newdy = c((newDistCoordinates[,2] - newDistCoordinates[,4])*(-1),NA)
    RepositionMatrix$dx.adjusted = newdx
    RepositionMatrix$dy.adjusted = newdy

    return(RepositionMatrix)
  }
  list.breaker<-function(x,y,id, dateTime, point.x, point.y, RepositionDir, UpDownRepositionLen, LeftRightRepositionLen, immobThreshold, parallel){
    input<- data.frame(y[unname(unlist(x[1]))])
    reposition<-reposition.generator(input, id, dateTime, point.x, point.y, RepositionDir, UpDownRepositionLen, LeftRightRepositionLen, immobThreshold, parallel)
    return(reposition)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.rep <- apply(breakFrame, 1, list.breaker,y = x,id, dateTime, point.x, point.y, RepositionDir, UpDownRepositionLen, LeftRightRepositionLen, immobThreshold, parallel) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.

    return(list.rep)

  }else{ #if x is a dataFrame

    frame.rep<- reposition.generator(x,id, dateTime, point.x, point.y, RepositionDir, UpDownRepositionLen, LeftRightRepositionLen, immobThreshold, parallel)

    return(frame.rep)
  }

}
