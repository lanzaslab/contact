#' Real-time location data for 20 calves in June 2018
#'
#' A dataset containing planar real-time point locations for 20 calves between
#'    00:00:00 on June 1st, 2018 and 23:59:59 UTC on June 3, 2018. 
#'
#' Calves were approximately 1.5-year-old castrated male cattle (i.e., steer) 
#'    kept in a 30 X 35 m2 pen at the Kansas State University Beef Cattle 
#'    Research Center in Manhattan, KS.
#' 
#' Data collection was supported by U.S. National Institute of Health (NIH) 
#'    grant R01GM117618 as part of the joint National Science 
#'    Foundation-NIH-United States Department of Agriculture Ecology and 
#'    Evolution of Infectious Disease program.
#'
#' @docType data
#' @usage data(calves2018)
#' @format A data frame with 193551 rows and 4 variables:
#' \describe{
#'   \item{calftag}{a unique identifier for each calf}
#'   \item{x}{planar x coordinate}
#'   \item{y}{planar y coordinate}
#'   \item{dateTime}{UTC date and time at which location fix was obtained}
#' }
#' @keywords datasets calves point location planar
#' @references Farthing, T.S., Dawson, D.E., Sanderson, M.W., and Lanzas, 
#'    C. in Review. Accounting for space and uncertainty in real-time-location-
#'    system-derived contact networks. Ecology and Evolution.
#' @examples
#' \donttest{
#' data("calves2018") #alternatively, you may use the command: contact::calves2018
#' head(calves2018)
#' }
"calves2018"