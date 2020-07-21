
## <!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# contact 1.0.0

  - submitted to CRAN on 10/18/2019

## contact 1.0.1

  - Addressed issues highlighted from CRAN reviewer.
  - Fixed a substantial bug in the “chisq” sub-function of contactTest.

## contact 1.2.0

  - Corrected many bugs in package functions.
  - Improved efficiency (e.g., run time, local memory usage, etc.) of
    package functions
  - Ensured compatability with R 4.0.0.
  - Updated package citation (because our paper was accepted for
    publication\!)
  - Broke contact::contactTest into multiple functions - contactTest is
    now defunct.

## contact 1.2.1

  - Fixed a bug introduced in version 1.2.0, as required by CRAN

## contact 1.2.2

  - Fixed bug associated with temporal sorting of data spanning multiple
    years.
  - Improved efficiency temporal blocking functions.

## contact 1.2.3

  - Fixed bug that forced times in tempAggregate output to begin at
    00:00:00 on the first day in the data set.
  - Added socialEdges function.

## contact 1.2.4

  - Fixed various bugs
  - Realized that the method the findDistThresh function was based on
    was faulty due to a mistake in the paper (Farthing et al. 2020). We
    updated the function to correctly do as intended and updated to
    function documentation to acknowledge the mistake.
  - Drastically improved efficiency of the ntwrkEdges function.

\`\`\`
