#load the full 2018 data set
##load("/Volumes/GoogleDrive/My Drive/KSU_RTLS_project/MovementData_2018/Processing_for_Modeling_Paper_01212019/dataset2018.dataframe__Movement2018Processing_tsf_01222019.RData")
##this dataset2018 contains calves from different studies, but we're only interested in keeping the 70 E. coli-pen cattle.
##load the vector of 2018 Ecoli-pen cattle IDs
##load("/Volumes/GoogleDrive/My Drive/KSU_RTLS_project/MovementData_2018/Ecoli2018_RediIDlist.RData")

#dataset2018_Ecoli<-droplevels(dataset2018[dataset2018$ID%in%idlist2018,])

#we're only going to include the June 2018 dates, so we need to get rid of the rest
#dataset2018_Ecoli$month<-lubridate::month(dataset2018_Ecoli$dateTime)
#calves2018<-droplevels(dataset2018_Ecoli[which(dataset2018_Ecoli$month == 6),1:4])
#calves2018<-dataset2018_Ecoli[dataset2018_Ecoli$date%in%dates2018.reduced,1:4]

#save(calves2018, file = "calves2018_raw.RData")

#change "ID," "xloc," and "yloc" columns to "calftag," "x," and "y," respectively.
colnames(calves2018)[1:3]<-c("calftag", "x", "y")
usethis::use_data(calves2018)



