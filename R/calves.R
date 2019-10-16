#' Real-time location data for 10 calves on May 2nd 2016
#'
#' A dataset containing planar real-time point locations for 10 calves between
#'    00:00:00 and 02:00:00 UTC on May 2nd, 2016. These data are a subset of 
#'    the data set published in the supplemental materials of Dawson et al.
#'    2019, and are included here primarily to be used for function-testing 
#'    purposes.
#'
#' Calves were approximately 1.5-year-old beef cattle kept in a 30 X 35 m2 pen 
#'    at the Kansas State University Beef Cattle Research Center in Manhattan, 
#'    KS.
#' 
#' Data collection was supported by U.S. National Institute of Health (NIH) 
#'    grant R01GM117618 as part of the joint National Science 
#'    Foundation-NIH-United States Department of Agriculture Ecology and 
#'    Evolution of Infectious Disease program.
#'
#' @docType data
#' @usage data(calves)
#' @format A data frame with 11118 rows and 5 variables:
#' \describe{
#'   \item{calftag}{a unique identifier for each calf}
#'   \item{x}{planar x coordinate}
#'   \item{y}{planar y coordinate}
#'   \item{time}{UTC time at which location fix was obtained}
#'   \item{date}{date on which fix location occurred}
#' }
#' @keywords datasets calves point location planar
#' @references Dawson, D.E., Farthing, T.S., Sanderson, M.W., and Lanzas, 
#'    C. 2019. Transmission on empirical dynamic contact networks is 
#'    influenced by data processing decisions. Epidemics 26:32-42. 
#' @source \url{https://doi.org/10.1016/j.epidem.2018.08.003}
#' @examples
#' data("calves") #alternatively, you may use the command: contact::calves
#' head(calves)
"calves"