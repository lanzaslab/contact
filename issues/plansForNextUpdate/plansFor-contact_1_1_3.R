#  1.) Update immobThreshold for-loop in vignette
#    - See if dist2 functions can be made more efficient with similar matrix transformations
#
#  2.) Change the mantel test in contactTest from ape::mantel to ade4::mantel.rtest because it gives more information.
#    - Ensure that all documentation is updated as well (e.g., alternative.hyp can only take the values of "greater," "less," or "two-sided").
#
#  3.) Add a function for creating contact matrices from contactDur output. 
#
#  4.) Update citation info (note: I'm not planning to update the package until the associated paper is accepted for publication).
#
#  5.) Add a fourth shuffle.type to randomizePaths allowing for the Spiegle method to be used while maintaining all days of the empirical set (i.e., shuffle blocks with full replacement). 
#