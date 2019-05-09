load("calves_raw.RData")
#we're only interested in keeping the first 12 hours-worth of data
data2<-subset(data1, Hour <12)

date <- rep("05-02-2016", nrow(data2)) #create a vector of length nrow(data2) detailing the date data were collected
m <- floor(((data2$Day_Second - 1)/60) - (data2$Hour*60)) #create a vector of length nrow(data2) detailing the minute of each hour fixes were reported. We assumed time started at 00:00:00.
s <- floor((data2$Day_Second - 1) - (m*60) - (data2$Hour*60*60)) #create a vector of length nrow(data2) detailing the second of each minute fixes were reported. We assumed time started at 00:00:00.

#add preceding zeroes to second, minute, and hour values less than 10.
Hour <- ifelse(data2$Hour < 10, paste(0,data2$Hour, sep = ""), data2$Hour)
Minute <- ifelse(m < 10, paste(0,m, sep = ""), m)
Second <- ifelse(s < 10, paste(0,s, sep = ""), s)
time <- (paste(Hour,":",Minute,":",Second, sep = "")) #create the time vector that combines Hour, Minute, and Second Vectors

#Ultimately, we want a data frame showing calftag ids, x location, y location, date, and time columns.
calves <- data2[,1:3] 
calves$time<-time
calves$date<-date
usethis::use_data(calves)
