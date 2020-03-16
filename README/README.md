
<!-- README.md is generated from README.Rmd. Please edit that file -->

# contact

<!-- badges: start -->

<!-- badges: end -->

This package is designed to process spatiotemporal data into contact and
social networks, and facilitate network analysis by randomizing
individuals’ movement paths and/or related categorical variables. To use
this package, users need only have a dataset containing spatial data
(i.e., latitude/longitude, or planar x & y coordinates), individual IDs
relating spatial data to specific individuals, and date/time information
relating spatial locations to temporal locations. The functionality of
this package ranges from data “cleaning” via multiple filtration
functions, to spatial and temporal data interpolation, and network
creation and summarization. Functions within this package are not
limited to describing interpersonal contacts. Package functions can also
identify and quantify “contacts” between individuals and fixed areas
(e.g., home ranges, water bodies, buildings, etc.). As such, this
package is an incredibly useful resource for facilitating
epidemiological, ecological, ethological and sociological research.

This package was created and is maintained by members of the [Lanzas
Lab](http://www.lanzaslab.org/) at the North Carolina State University
College of Veterinary Medicine.

## Installation

You can install the released version of contact from
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("contact")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("lanzaslab/contact")
```

## Example

*This example is an excerpt from the package vignette.*

Let’s say you want to create a time-aggregated point-location-based
environmental contact networ describing instances when animals,
represented by data points, were in “contact” with a fixed
location/area. As an example of this, here we use *contact* functions to
describe contacts between tracked calves in a single feedlot pen and
their water trough.

``` r
library(contact)

data("calves") #load the calves data set
```

**A.) Ensure all required columns exist in the calves data set (i.e., xy
coordinates, unique individual IDs, dateTime)**

``` r

head(calves)#The calves data set does not have a singular dateTime column. Rather, it has "date" and "time" columns. We must append a dateTime column to the data frame.
#>   calftag     x     y     time       date
#> 1     101 63.38 47.43 00:00:14 05-02-2016
#> 2     101 63.48 46.59 00:00:22 05-02-2016
#> 3     101 62.90 47.07 00:00:26 05-02-2016
#> 4     101 63.27 47.26 00:00:34 05-02-2016
#> 5     101 63.04 46.62 00:00:38 05-02-2016
#> 6     101 63.53 47.24 00:00:42 05-02-2016

calves.dateTime<- contact::datetime.append(x = calves, date = calves$date, time= calves$time, dateTime = NULL, dateFormat = "mdy", dateFake = FALSE, startYear = NULL, tz.in = "UTC", tz.out = NULL, month = FALSE, day = FALSE, year = FALSE, hour = FALSE, minute = FALSE, second = FALSE, daySecond = FALSE, totalSecond = FALSE)
```

**B.) Calculate distances between the water-trough polygon and calves at
each time step.**

First, we must define the location of the water trough. To do this, we
read in point-location data for the water-trough vertices.

``` r

water<- data.frame(x = c(61.43315, 61.89377, 62.37518, 61.82622), y = c(62.44815, 62.73341, 61.93864, 61.67411)) #This is a data frame containing the x and y coordinates of the four trough vertices.
```

As noted in the dist2Area\_df help documention, polygon-vertex
coordinates must be arranged in a particular way. Here we arrange them
accordingly.

``` r

water_poly<-data.frame(matrix(ncol = 8, nrow = 1)) #(ncol = number of vertices)*2
colnum = 0
for(h in 1:nrow(water)){
  water_poly[1,colnum + h] <- water$x[h] #pull the x location for each vertex
  water_poly[1, (colnum + 1 + h)] <- water$y[h] #pull the y location for each vertex
  colnum <- colnum + 1
}
```

Calculate distances between calves and the water polygon at every
timestep.

``` r

water_distance<-contact::dist2Area_df(x = calves.dateTime, y = water_poly, x.id = "calftag", y.id = "water", dateTime = "dateTime", point.x = calves.dateTime$x, point.y = calves.dateTime$y, poly.xy = NULL, parallel = FALSE, dataType = "Point", lonlat = FALSE, numVertices = NULL) #note that the poly.xy and numVertices arguments refer to vertices of polygons in x, not y. Because dataType is "Point," not "Polygon," these arguments are irrelevant here.

head(water_distance)
#>              dateTime totalIndividuals individualsAtTimestep  id dist.to.water
#> 1 2016-05-02 00:00:14               10                     3 101      14.32860
#> 2 2016-05-02 00:00:14               10                     3 104      12.04416
#> 3 2016-05-02 00:00:14               10                     3 108      32.84068
#> 4 2016-05-02 00:00:22               10                     4 101      15.17450
#> 5 2016-05-02 00:00:22               10                     4 104      12.06342
#> 6 2016-05-02 00:00:22               10                     4 108      27.13532
```

**C.) Identify what SpTh value will allow us to capture 99% of contacts,
defined as instances when point-locations were within 0.333 m of the
water trough, given the RTLS accuracy.**

``` r

SpThValues<-contact::findDistThresh(n = 100000, acc.Dist1 = 0.5, acc.Dist2 = NULL, pWithin1 = 90, pWithin2 = NULL, spTh = 0.5) #spTh represents the initially-defined spatial threshold for contact. #spTh represents the initially-defined spatial threshold for contact. Note that we've chosen to use 100,000 in-contact point-location pairs here.

SpThValues #it looks like an adjusted SpTh value of approximately 0.71 m will likely capture 99% of contacts, defined as instances when point-locations were within 0.333 m of the water trough, given the RTLS accuracy. #Note that because these confidence intervals are obtained from distributions generated from random samples, every time this function is run, results will be slightly different. 
#>      mean     5%-CI    10%-CI    15%-CI    20%-CI    25%-CI    30%-CI    35%-CI 
#> 0.7077398 0.7078082 0.7078769 0.7079461 0.7080162 0.7080874 0.7081602 0.7082348 
#>    40%-CI    45%-CI    50%-CI    55%-CI    60%-CI    65%-CI    70%-CI    75%-CI 
#> 0.7083119 0.7083919 0.7084756 0.7085639 0.7086579 0.7087594 0.7088705 0.7089947 
#>    80%-CI    85%-CI    90%-CI    95%-CI    99%-CI       max       TPR 
#> 0.7091379 0.7093102 0.7095342 0.7098779 0.7105498 2.7230331 0.3032000

CI_99<-unname(SpThValues[21]) #we will use this SpTh value moving forward.
```

**D.) Identify time points when calves were within the re-adjusted SpTh
distance from water trough.**

``` r

water_contacts <- contact::contactDur.area(water_distance, dist.threshold=CI_99,sec.threshold=1, blocking = FALSE, equidistant.time = FALSE, parallel = FALSE, reportParameters = TRUE) #Note that because we are not interested in making a time-aggregated network with > 1 temporal levels, we set blocking = FALSE to reduce processing time.

head(water_contacts)
#>   indiv.id area.id contact.id contactDuration    contactStartTime
#> 1      101   water  101-water               1 2016-05-02 00:47:35
#> 2      101   water  101-water               1 2016-05-02 00:47:43
#> 3      101   water  101-water               1 2016-05-02 00:47:47
#> 4      101   water  101-water               1 2016-05-02 00:47:59
#> 5      101   water  101-water               1 2016-05-02 00:48:15
#> 6      101   water  101-water               1 2016-05-02 00:48:36
#>        contactEndTime distThreshold secThreshold equidistant.time
#> 1 2016-05-02 00:47:35     0.7105498            1            FALSE
#> 2 2016-05-02 00:47:43     0.7105498            1            FALSE
#> 3 2016-05-02 00:47:47     0.7105498            1            FALSE
#> 4 2016-05-02 00:47:59     0.7105498            1            FALSE
#> 5 2016-05-02 00:48:15     0.7105498            1            FALSE
#> 6 2016-05-02 00:48:36     0.7105498            1            FALSE
```
