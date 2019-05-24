
#Note that the raw data set provided here has 1000000 fewer observations than the raw data set available in the movebank data repository. This was done to reduce the size below 100mb (the maximum allowable size able to be pushed to Github).

###baboons_raw<- read.csv("Collective movement in wild baboons (data from Strandburg-Peshkin et al. 2015).csv")
###baboons_raw<-baboons_raw[1:(nrow(baboons_raw) - 1000000),]
###save(baboons_raw, file = "baboons_raw.RData")

load("baboons_raw.RData")
#To reduce processing times associated with using our functions on this data set, we will only keep data from the first 2 days of the 14 days presented in the data set.
#so, we first need to quantify days in the data set.
baboon.dateTime<-contact::datetime.append(baboons, dateTime = lubridate::ymd_hms(baboons$timestamp), dateFormat = "ymd", day = TRUE)
#now subset the data set
baboons<-subset(baboon.dateTime, day <= 2)
baboons<-droplevels(baboons)
#now remove the "day" column
baboons<-baboons[,-match("day", names(baboons))]
#finally, the data set reports time points to 1/1000th of the second. This is a bit too precise for our needs. We'll round up each value to the nearest second.
baboons$dateTime <- lubridate::ymd_hms(baboons$timestamp)
#we create the processed data set here
usethis::use_data(baboons)



