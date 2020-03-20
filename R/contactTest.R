#' Determine if Observed Contacts are More or Less Frequent than in a Random
#'    Distribution (Defunct)
#'
#' This DEFUNCT function was used to determine if tracked individuals in an 
#'    empirical dataset had more or fewer contacts with other tracked 
#'    individuals/specified locations than would be expected at random. The 
#'    function works by comparing an empirically-based contactDur.all or 
#'    contactDur.area function output (emp.input) to the contactDur.all or 
#'    contactDur.area output generated from randomized data (rand.input).
#' 
#' @param ... Any input will return the error message: "'contactTest' is now 
#'    defunct. Please consider using another contact-comparison function 
#'    instead (e.g., \code{\link{contactCompare_chisq}}, 
#'    \code{\link{contactCompare_mantel}}, etc.)."
#' @keywords network-analysis social-network defunct
#' @return Always returns the error message: "'contactTest' is now 
#'    defunct. Please consider using another contact-comparison function 
#'    instead (e.g., \code{\link{contactCompare_chisq}}, 
#'    \code{\link{contactCompare_mantel}}, etc.)."
#'    
#' @export

contactTest<-function(...){
  
  .Defunct(msg = "'contactTest' is now defunct. Please consider using another contact-comparison function instead (e.g., contactCompare_chisq, contactCompare_mantel, etc.)")
  
}
