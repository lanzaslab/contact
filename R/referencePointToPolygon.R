#' Create a Rectangular Polygon Using Planar XY Coordinates
#'
#' This function creates a square/rectangular polygon from a single reference point by translating its location multiple times using the same method used in repositionReferencePoint. For example, even though calves in our study were only equiped with RFID tags on their left ear. With this function, we can create polygons that account for the total space used by each individual at each time step.This function is different from similar point-to-polygon functions for two reasons:
#' 1.) It does not assume points lie within the center of the polygon. Rather, the reference point must be a corner of the polygon (Note: "UL" denotes that reference point lies on the upper-left corner of the polygon, "UR" denotes that reference point lies on the upper-right corner of the polygon,"DL" denotes that reference point lies on the down-left corner of the polygon, "DR" denotes that reference point lies on the down-left corner of the polygon). Note that if you want the reference point to be at the center of the polygon, you can first translate the reference point to a central location on tracked individuals using repositionReferencePoint.
#' 2.) Polygon angles/directionality are based on observed movements of tracked individuals.
#'
#' Currently, this function only supports input data with coordinates representing planar ('Euclidean') space (e.g. units of meters).
#'
#' In the output, point1.x and point1.y represent the xy coordinates from the input file. Point2-n coordinates move in a clockwise direction from point1. For example: if point1 is located on the upper left ("UL") corner of the polygon, point2 would be on the upper right corner, point3 on the bottom right, and point 4 on the bottom left.
#'
#' Because this function needs information (dist, dx, dy) from 2 points on an individual's path to work, at least the first point in each individual's path will be removed (the function will report NAs for adjusted locations). Also note that if the distance between an individual's first point in their path and the second one is 0, the function will also report NAs for the second point's adjusted coordinates. The first non-NA values will only be reported for the instance where dist > 0.
#'
#' In the output, if input was previously processed using tempAggregate at a Day or Hour resolutionLevel, dt > secondAgg indicates that tracked individuals were missing in the original dataset on a given day or hour, respectively. In this case, the assumption that individuals are facing a given direction because they moved from the previous timepoint may not be accurate. Consider removing these rows (rows following one with dt > secondAgg; remember that dt indicates the time between recording xy coordinates in row i to row i + 1) from your dataset.
#' @param x Description imminent
#' @param id Vector of length(nrow(data.frame(x))) or singular character data, detailng the relevant colname in x, that denotes what date information will be used. If argument == NULL, datetime.append assumes a column withe colname "id" exists in x. Defaults to NULL.
#' @param point.x Description imminent
#' @param point.y Description imminent
#' @param dateTime Description imminent
#' @param StartLocation Description imminent
#' @param UpDownRepositionLen Description imminent
#' @param LeftRightRepositionLen Description imminent
#' @param CenterPoint Description imminent
#' @param MidPoints Description imminent
#' @param StandingThreshold Description imminent
#' @param parallel Description imminent
#' @keywords data-processing polygon point location planar
#' @export
#' @examples
#' Examples imminent

referencePointToPolygon <-function(x = NULL, id = NULL, dateTime = NULL, point.x = NULL, point.y = NULL, StartLocation = "UL", UpDownRepositionLen = 1, LeftRightRepositionLen = 1, CenterPoint = FALSE, MidPoints = FALSE, StandingThreshold = 0, parallel = TRUE){

  poly.generator<-function(x, id, dateTime, point.x, point.y, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, StandingThreshold, parallel){
    euc=function(x) {
      point1 = x.cor=unlist(c(x[1],x[2]))
      point2 = x.cor=unlist(c(x[3],x[4]))
      euc.dis = raster::pointDistance(p1 = point1, p2 = point2, lonlat = FALSE)
      return(euc.dis)
    }
    timeDifference = function(x){
      t1 = x[1]
      t2 = x[2]
      dt = as.integer(difftime(time1 = t2, time2 = t1, units = "secs"))
      return(dt)
    }
    translation <-function(x){

      dx = as.numeric(x[4])
      dy = as.numeric(x[5])
      dist = as.numeric(x[3])

      if (dist == 0 || is.na(x[1]) == TRUE || is.na(x[2]) == TRUE || is.na(dx) == TRUE || is.na(dy) == TRUE || is.na(dist) == TRUE){ #if individuals did not move or any of the necessary variables are undefined, it reports an NA.
        newX = NA
        newY = NA}else{

          if (x[6] == "Up" || x[6] == "UP" || x[6] == "up"){ #If you want to reposition the reference point at a part of the animal/object ahead of the RFID tag/GPS locator

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

          if (x[6] == "Down" || x[6] == "DOWN" || x[6] == "down"){ #If you want to reposition the reference point at a part of the animal/object behind of the RFID tag/GPS locator

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

          if (x[6] == "Left" || x[6] == "LEFT" || x[6] == "left"){ #If you want to reposition the reference point at a part of the animal/object to the left of the RFID tag/GPS locator

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

          if (x[6] == "Right" || x[6] == "RIGHT" || x[6] == "right"){ #If you want to reposition the reference point at a part of the animal/object to the right of the RFID tag/GPS locator #Note that the calculations are the same as for left translation, but the effects of adjustments on xy coordinates are reversed.

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

          if (x[6] == "Up&Right" || x[6] == "UP&RIGHT" || x[6] == "up&right" || x[6] == "Up & Right" || x[6] == "UP & RIGHT" || x[6] == "up & right"){ #If you want to reposition the reference point at a part of the animal/object ahead of the RFID tag/GPS locator

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

          if (x[6] == "Up&Left" || x[6] == "UP&LEFT" || x[6] == "up&left" || x[6] == "Up & Left" || x[6] == "UP & LEFT" || x[6] == "up & left"){ #If you want to reposition the reference point at a part of the animal/object ahead of the RFID tag/GPS locator

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

          if (x[6] == "Down&Right" || x[6] == "DOWN&RIGHT" || x[6] == "down&right" || x[6] == "Down & Right" || x[6] == "DOWN & RIGHT" || x[6] == "down & right"){ #If you want to reposition the reference point at a part of the animal/object behind of the RFID tag/GPS locator

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
              #theta2 = (180 - adjAngle)/2
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

          if (x[6] == "Down&Left" || x[6] == "DOWN&LEFT" || x[6] == "down&left" || x[6] == "Down & Left" || x[6] == "DOWN & LEFT" || x[6] == "down & left"){ #If you want to reposition the reference point at a part of the animal/object behind of the RFID tag/GPS locator

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
    standingAdjustment.polygon = function(x, locMatrix, standingVec){

      dur = 0
      standDurTest <- FALSE
      while (!standDurTest) { #This means when standDurTest == FALSE #This part of the function determines how long an individual has been standing.
        standDurTest <- standingVec[(x[1] - (1 + dur))] == 0
        standDurTest <- ifelse(is.na(standDurTest) == TRUE, TRUE, standDurTest) #This is for when x[1] -1 relates to the first observation for any individual (because those observations are NA)
        dur = (dur + 1)
      }

      standFrame = data.frame(replaceRow = x[1], replacePoint2.x = locMatrix$cornerPoint2.x[(x[1] - dur)], replacePoint2.y = locMatrix$cornerPoint2.y[(x[1] - dur)], replacePoint3.x = locMatrix$cornerPoint3.x[(x[1] - dur)], replacePoint3.y = locMatrix$cornerPoint3.y[(x[1] - dur)], replacePoint4.x = locMatrix$cornerPoint4.x[(x[1] - dur)], replacePoint4.y = locMatrix$cornerPoint4.y[(x[1] - dur)]) #if individuals are standing, their location at a given timestep will be the same as at the previous timestep. There's no need to replace point1 because it is the pre-translated (a.k.a. original) location of the individual.

      return(standFrame)
    }

    if(length(x) == 0){ #This if statement allows users to input either a series of vectors (id, dateTime, point.x and point.y), a dataframe with columns named the same, or a combination of dataframe and vectors.
      x = data.frame(id = id, x = point.x, y = point.y, dateTime = dateTime)

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
    rownames(x)<-seq(1,nrow(x),1)

    distCoordinates = data.frame(x$x[1:(nrow(x) - 1)], x$y[1:(nrow(x) - 1)], x$x[2:nrow(x)], x$y[2:nrow(x)])
    if(parallel == TRUE){ #This calculates the new distance between adjusted xy coordinates. Reported distances are distances an individual at a given point must travel to reach the subsequent point. Because there is no previous point following the final observation in the dataset, the distance observation at RepositionMatrix[nrow(RepositionMatrix),] is always going to be NA
      cl<-parallel::makeCluster(parallel::detectCores())
      dist = parallel::parApply(cl, distCoordinates,1,euc)
      dist = c(dist, NA) #to make length(dist) == nrow(x)
      parallel::stopCluster(cl)
    }else{
      dist = apply(distCoordinates,1,euc)
      dist = c(dist, NA) #to make length(dist) == nrow(x)
    }

    dx = c((distCoordinates[,1] - distCoordinates[,3])*(-1),NA)
    dy = c((distCoordinates[,2] - distCoordinates[,4])*(-1),NA)

    timesFrame = data.frame(x$dateTime[1:(nrow(x) - 1)], x$dateTime[2:nrow(x)])
    if(parallel == TRUE){
      cl<-parallel::makeCluster(parallel::detectCores())
      dt = parallel::parApply(cl, timesFrame, 1, timeDifference)
      dt = c(dt, NA)  #to make length(dt) == nrow(x)
      parallel::stopCluster(cl)
    }else{
      dt = apply(timesFrame, 1, timeDifference)
      dt = c(dt, NA)  #to make length(dt) == nrow(x)
    }

    idVec = unique(x$id)
    idSeqVec = x$id
    idVecForDtRemoval = idVec[-length(idVec)] #There's no need to remove the max point of the final id point b/c there's already an NA there.
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

    dataShift = data.frame(dist = c(NA,x$dist[1:(nrow(x) - 1)]), dx = c(NA,x$dx[1:(nrow(x) - 1)]), dy = c(NA,x$dy[1:(nrow(x) - 1)])) #This is necessary because of the way we calculated/listed these values above. These columns (dist, dx, and dt) refer to the changes a point must make to reach the subsequent point. However, later on in this function, we are not interested in future point alterations. Rather, we need to know how tracked individuals moved during the preceding time step to reach their current point (to determine directionality of movement) #Note that this code shifts values down, but remember that beacuse values are downshifted, the first observation for each id individual will be incorrect. Code below addresses this issue.

    for(b in idVec){ #This code replaces dist, dx, and dy values in the row when a new individual (id) is first observed with NAs to fix the problem noted above. #Note that this loop assumes that the data is sorted by id (i.e., a given id will not repeat in the dataset, after a new id has appeared)
      dataShift[min(which(idSeqVec == b)),] = NA
    }

    if (StartLocation == "UL"){
      ##xyFrames here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      xyFrame1 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Right", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame2 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Down&Right", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame3 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Down", UpDownRepositionLen, LeftRightRepositionLen)

    }

    if (StartLocation == "UR"){
      ##xyFrames here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      xyFrame1 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Down", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame2 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Down&Left", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame3 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Left", UpDownRepositionLen, LeftRightRepositionLen)

    }

    if (StartLocation == "DR"){
      ##xyFrames here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      xyFrame1 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Left", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame2 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Up&Left", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame3 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Up", UpDownRepositionLen, LeftRightRepositionLen)

    }

    if (StartLocation == "DL"){
      ##xyFrames here calculate the other 3 corners of the polygon (the original xy coordinates are always going to be the first corner point)
      xyFrame1 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Up", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame2 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Up&Right", UpDownRepositionLen, LeftRightRepositionLen)
      xyFrame3 = data.frame(x$x, x$y, dataShift$dist, dataShift$dx, dataShift$dy, "Right", UpDownRepositionLen, LeftRightRepositionLen)

    }

    if (parallel == TRUE){
      cl<-parallel::makeCluster(parallel::detectCores())

      xyTranslated1 = parallel::parApply(cl, xyFrame1, 1, translation)
      newCoordinates1 = data.frame(data.table::rbindlist(xyTranslated1))
      xyTranslated2 = parallel::parApply(cl, xyFrame2, 1, translation)
      newCoordinates2 = data.frame(data.table::rbindlist(xyTranslated2))
      xyTranslated3 = parallel::parApply(cl, xyFrame3, 1, translation)
      newCoordinates3 = data.frame(data.table::rbindlist(xyTranslated3))

      parallel::stopCluster(cl)

    }else{

      xyTranslated1 = apply(xyFrame1, 1, translation)
      newCoordinates1 = data.frame(data.table::rbindlist(xyTranslated1))
      xyTranslated2 = apply(xyFrame2, 1, translation)
      newCoordinates2 = data.frame(data.table::rbindlist(xyTranslated2))
      xyTranslated3 = apply(xyFrame3, 1, translation)
      newCoordinates3 = data.frame(data.table::rbindlist(xyTranslated3))

    }

    polygonMatrix = data.frame(matrix(nrow = nrow(x), ncol = 26)) #One row for each observation, Columns listed below
    colnames(polygonMatrix) = c("id","cornerPoint1.x","cornerPoint1.y","midPoint1.x","midPoint1.y","cornerPoint2.x","cornerPoint2.y","midPoint2.x","midPoint2.y","cornerPoint3.x","cornerPoint3.y","midPoint3.x","midPoint3.y","cornerPoint4.x","cornerPoint4.y","midPoint4.x","midPoint4.y","centroid.x","centroid.y","startLocation","upDownRepositionLength","leftRightRepositionLength","standing","standingThreshold", "dateTime","dt")

    polygonMatrix$id = x$id
    polygonMatrix$cornerPoint1.x = x$x
    polygonMatrix$cornerPoint1.y = x$y
    polygonMatrix$cornerPoint2.x = newCoordinates1$x.adjusted
    polygonMatrix$cornerPoint2.y = newCoordinates1$y.adjusted
    polygonMatrix$cornerPoint3.x = newCoordinates2$x.adjusted
    polygonMatrix$cornerPoint3.y = newCoordinates2$y.adjusted
    polygonMatrix$cornerPoint4.x = newCoordinates3$x.adjusted
    polygonMatrix$cornerPoint4.y = newCoordinates3$y.adjusted
    polygonMatrix$startLocation = StartLocation
    polygonMatrix$upDownRepositionLength = UpDownRepositionLen
    polygonMatrix$leftRightRepositionLength = LeftRightRepositionLen
    polygonMatrix$standingThreshold = StandingThreshold
    polygonMatrix$dateTime = x$dateTime
    polygonMatrix$dt = x$dt

    polygonMatrix$standing = ifelse(xyFrame1[,3] > StandingThreshold, 0, 1) #if the distance individuals moved was less than / equal to the noted StandingThreshold, individuals are said to be "standing," and their position will not change relative to their previous one. (i.e., you assume that any observed movement less than StandingThreshold was due to errors or miniscule bodily movements (e.g., head shaking) that are not indicative of actual movement.)

    standlist = which(polygonMatrix$standing == 1)

    if ((length(standlist) >= 1)){ #To save processing time and reduce chances of errors, this evaluation will not take place if there is only one observation.

      standingVec = polygonMatrix$standing
      standlistFrame = data.frame(standing = standlist)

      if(parallel == TRUE){
        cl<-parallel::makeCluster(parallel::detectCores())
        standingFrame = data.frame(data.table::rbindlist(parallel::parApply(cl, standlistFrame, 1, standingAdjustment.polygon, locMatrix = polygonMatrix, standingVec)))
        parallel::stopCluster(cl)
      }else{
        standingFrame = data.frame(data.table::rbindlist(apply(standlistFrame, 1, standingAdjustment.polygon, locMatrix = polygonMatrix, standingVec)))
      }
      if(nrow(standingFrame) > 0){
        polygonMatrix[standingFrame$replaceRow,match("cornerPoint2.x", names(polygonMatrix))] = standingFrame$replacePoint2.x
        polygonMatrix[standingFrame$replaceRow,match("cornerPoint2.y", names(polygonMatrix))] = standingFrame$replacePoint2.y
        polygonMatrix[standingFrame$replaceRow,match("cornerPoint3.x", names(polygonMatrix))] = standingFrame$replacePoint3.x
        polygonMatrix[standingFrame$replaceRow,match("cornerPoint3.y", names(polygonMatrix))] = standingFrame$replacePoint3.y
        polygonMatrix[standingFrame$replaceRow,match("cornerPoint4.x", names(polygonMatrix))] = standingFrame$replacePoint4.x
        polygonMatrix[standingFrame$replaceRow,match("cornerPoint4.y", names(polygonMatrix))] = standingFrame$replacePoint4.y
      }
    }

    if(CenterPoint == TRUE){
      polygonMatrix$centroid.x = ((polygonMatrix$cornerPoint1.x + polygonMatrix$cornerPoint3.x)/2)
      polygonMatrix$centroid.y = ((polygonMatrix$cornerPoint1.y + polygonMatrix$cornerPoint3.y)/2)
    }else{polygonMatrix <- polygonMatrix[,-c(match("centroid.x", names(polygonMatrix)),match("centroid.y", names(polygonMatrix)))]}

    if(MidPoints == TRUE){
      polygonMatrix$midPoint1.x = ((polygonMatrix$cornerPoint1.x + polygonMatrix$cornerPoint2.x)/2)
      polygonMatrix$midPoint1.y = ((polygonMatrix$cornerPoint1.y + polygonMatrix$cornerPoint2.y)/2)
      polygonMatrix$midPoint2.x = ((polygonMatrix$cornerPoint2.x + polygonMatrix$cornerPoint3.x)/2)
      polygonMatrix$midPoint2.y = ((polygonMatrix$cornerPoint2.y + polygonMatrix$cornerPoint3.y)/2)
      polygonMatrix$midPoint3.x = ((polygonMatrix$cornerPoint3.x + polygonMatrix$cornerPoint4.x)/2)
      polygonMatrix$midPoint3.y = ((polygonMatrix$cornerPoint3.y + polygonMatrix$cornerPoint4.y)/2)
      polygonMatrix$midPoint4.x = ((polygonMatrix$cornerPoint4.x + polygonMatrix$cornerPoint1.x)/2)
      polygonMatrix$midPoint4.y = ((polygonMatrix$cornerPoint4.y + polygonMatrix$cornerPoint1.y)/2)
    }else{polygonMatrix <- polygonMatrix[,-c(match("midPoint1.x", names(polygonMatrix)),match("midPoint1.y", names(polygonMatrix)),match("midPoint2.x", names(polygonMatrix)),match("midPoint2.y", names(polygonMatrix)),match("midPoint3.x", names(polygonMatrix)),match("midPoint3.y", names(polygonMatrix)),match("midPoint4.x", names(polygonMatrix)),match("midPoint4.y", names(polygonMatrix)))]}

    return(polygonMatrix)
  }

  list.breaker<-function(x,y,id, dateTime, point.x, point.y, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, StandingThreshold, parallel){
    input<- data.frame(y[unname(unlist(x[1]))])
    newPolygons<-poly.generator(input, id, dateTime, point.x, point.y, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, StandingThreshold, parallel)
    return(newPolygons)
  }

  if(is.data.frame(x) == FALSE & is.list(x) == TRUE){ #02/02/2019 added the "is.data.frame(x) == FALSE" argument because R apparently treats dataframes as lists.
    breakFrame<- data.frame(seq(1,length(x),1))
    list.poly <- apply(breakFrame, 1, list.breaker,y = x,id, dateTime, point.x, point.y, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, StandingThreshold, parallel) #in the vast majority of cases, parallelizing the subfunctions will result in faster processing than parallelizing the list processing here. As such, since parallelizing this list processing could cause numerous problems due to parallelized subfunctions, this is an apply rather than a parApply or lapply.

    return(list.poly)

  }else{ #if x is a dataFrame
    frame.poly<- poly.generator(x,id, dateTime, point.x, point.y, StartLocation, UpDownRepositionLen, LeftRightRepositionLen, CenterPoint, MidPoints, StandingThreshold, parallel)

    return(frame.poly)
  }
}
