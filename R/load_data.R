#### libraries
library(tigris)
library(rgeos)
library(sp)
library(sf)
library(spdep)
library(spatialreg)
library(MMWRweek)
library(censusapi)
library(plotrix)
library(scales)
library(lubridate)



#### read Covid data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data.path='../data/'
path1=paste(data.path,'FL Case Count Data\\',sep="")
path2=paste(data.path,'CDC_SVI\\',sep="")
path3=paste(data.path,'covid19PVI\\',sep="")


# number of days
N=88
currentDay <- as.Date("07/22/2020",format = "%m/%d/%y")

date2=format(currentDay, "%m%d20%y")
dateVec=c(date2)
dataRead<-read.csv(unzip(paste(path1,"Florida_COVID19_",date2,"_ByZip_CSV.zip",sep=""),
                             "Florida_Cases_Zips_0.csv"))

dataZipCodes=dataRead$ZIP.Code
dataCases=dataRead$Cases_1
dataCases[dataCases=="<5"]=0
dataCases=strtoi(dataCases)
for (i in 2:N){
  # daily data from 08/01/2020 to 10/18/2020. only 08/09/2020 is missing. Make it equal to 08/08 
  day(currentDay)=day(currentDay)+1
  date2=format(currentDay, "%m%d20%y")
  dateVec=c(dateVec,date2)
  
  if(date2!="08092020"){
  dataRead<-read.csv(unzip(paste(path1,"Florida_COVID19_",date2,"_ByZip_CSV.zip",sep=""),
                           "Florida_Cases_Zips_0.csv"))}
  else{
    dataRead<-read.csv(unzip(paste(path1,"Florida_COVID19_","08082020","_ByZip_CSV.zip",sep=""),
                             "Florida_Cases_Zips_0.csv"))
  }
  
  dataZipCodesNew=dataRead$ZIP.Code
  dataCasesNew=dataRead$Cases_1
  dataCasesNew[dataCasesNew=="<5"]=0
  dataCasesNew=strtoi(dataCasesNew)

  dataCases=cbind(dataCases,dataCasesNew)
  dataZipCodes=cbind(dataZipCodes,dataZipCodesNew)  # make sure that all zip code columns (from different data files) are identical
  
}
#### zip level covid data ####
# read "zcta to tract relationship" file 
zctaToTractFile=read.csv("https://www2.census.gov/geo/docs/maps-data/data/rel/zcta_tract_rel_10.txt")
# shape file for FL tracts 
options(tigris_use_cache = TRUE) #ZCTAs can take several minutes to download. To cache the data use this.
fl_tracts=as_Spatial(tracts("Florida",year=2018,cb=TRUE))
# FL tract populations - Retrieve 2010 census tract-level data for any variable within a specific state/county 
Sys.setenv(CENSUS_KEY="0e95aa62e853f18e1d8099dbbc0ca054cdf14f81"); readRenviron("~/.Renviron");Sys.getenv("CENSUS_KEY")
fl2010 <- getCensus(name = "dec/sf1", vintage = 2010,vars = c("H010001"),region = "tract:*", regionin = "state:12")
# Fl county populations
AllCountyNames=fips_codes$county[fips_codes$state_code=='12']  #use tigris package function fips_codes; FL state code is 12
AllCountyCodes=fips_codes$county_code[fips_codes$state_code=='12']
fl2010_cnty <- getCensus(name = "dec/sf1", vintage = 2010,vars = c("H010001"),region = "county:*", regionin = "state:12")
fl2010_cnty$name=AllCountyNames[match(fl2010_cnty$county,AllCountyCodes)]
#zcta to tract for FL
zctaToTractFileFL=zctaToTractFile[which(zctaToTractFile$STATE=='12'),]
(noFLTracts=dim(fl_tracts)[1])
(dim(zctaToTractFileFL))
FLZipCodes=unique(zctaToTractFileFL$ZCTA5)
# how many tracts does each zip code have?
noOfTractsZipCodes=as.matrix(table(factor(zctaToTractFileFL$ZCTA5))) # rownames(noOfTractsZipCodes)[1]
(noFLZipCodes=dim(noOfTractsZipCodes)[1])

#zip codes of case data
uniqueZipCodesC=unique(dataZipCodes[,1])
indCovidCaseofFLzipcodes=match(FLZipCodes,uniqueZipCodesC)  # Many NA's--> These zip codes did not report case counts
# covid daily case count data at zip code level in the order as FLZipCodes
covid_DailyCase=matrix(NA,noFLZipCodes,length(dateVec))
for (i in 1:noFLZipCodes){
    covid_DailyCase[i,]=dataCases[indCovidCaseofFLzipcodes[i],]
}
# replace the NA's resulting due to non-reporting ZC's with 0's
covid_DailyCase[is.na(covid_DailyCase)]=0
rownames(covid_DailyCase)=FLZipCodes

matplot(t(covid_DailyCase),type='l',xlab='date',xaxt='n',ylab='zip code covid counts ',log='y')
seqTime=1:N
axis(1, at=seqTime,labels=dateVec, las=2)

#### tract level covid data ####

# to match the GEOID field in zctaToTractFileFL, create geoID in fl2010
fl2010$geoID=paste(fl2010$state,fl2010$county,fl2010$tract,sep="")
matchTractToPop=match(fl_tracts$GEOID,fl2010$geoID) 
fl_tracts$pop=fl2010$H010001[matchTractToPop]
fl_tracts$validateGEOID=fl2010$geoID[matchTractToPop]
# tracts with non zero populations
indNonZeroTract=which(fl_tracts$pop!=0)
# find zip codes of FL tracts
matchTractToZipCodes=match(fl_tracts$GEOID,as.character(zctaToTractFileFL$GEOID)) 
fl_tracts$ZipCode=zctaToTractFileFL$ZCTA5[matchTractToZipCodes]

noDays=length(dateVec)
daynames=1:noDays
fl_tracts=cbind(fl_tracts, setNames(lapply(daynames, function(x) x=NA), daynames) )

for (i in 1:length(fl_tracts)){
    i1=which(fl_tracts$ZipCode[i]==FLZipCodes)
    if(!is.na(fl_tracts$ZipCode[i])){
    # zip code counts is apportioned to tracts proportional to population
    for (j in 1:noDays){
      fl_tracts[i,12+j]=covid_DailyCase[i1,j]/noOfTractsZipCodes[i1]  
    }
  }
}
covid_dailyCaseTract=fl_tracts@data[,13:(12+noDays)]

matplot(t(covid_dailyCaseTract),type='l',xlab='date',xaxt='n',ylab='tract covid counts ')
seqTime=1:N
axis(1, at=seqTime,labels=dateVec, las=2)


plot(colSums(covid_DailyCase,na.rm = TRUE),type='l',ylim=c(0,1e6),ylab='total florida counts',xlab='day')
lines(colSums(covid_dailyCaseTract,na.rm = TRUE),col='red')
legend("topleft",legend=c('from zip code counts','from tract counts'),col=c('black','red'),lty=c(1,1))

#### read CDC SVI data ####
#read CDC SVI at tract level for FL
CDC_SVI=read.csv(paste(path2,"cdc-sovi-tract 2018-csv-FL.csv",sep =''))
# reorder dataframe to match the shapefile order of the tracts
indToReorderCDC_SVI=match(fl_tracts$GEOID,as.character(CDC_SVI$FIPS))
CDC_SVIr=CDC_SVI[indToReorderCDC_SVI,]
#replace with na'a rows with -999 value in any of the vulnerability variable
CDC_SVIr$SPL_THEME1[which(CDC_SVIr$SPL_THEME1==-999)]=NA
CDC_SVIr$SPL_THEME2[which(CDC_SVIr$SPL_THEME2==-999)]=NA
CDC_SVIr$SPL_THEME3[which(CDC_SVIr$SPL_THEME3==-999)]=NA
CDC_SVIr$SPL_THEME4[which(CDC_SVIr$SPL_THEME4==-999)]=NA
# replace -999 with NA
CDC_SVIr[CDC_SVIr==-999]=0

#### read covid19PVI data - pctbeds, air pollution, pct no insurance,obesity, diabetes, smoking (county level) ####
dataEnvironmental=read.csv(paste(path3,'PVIdata.0722020.csv',sep=""))
dataEnvironmental=dataEnvironmental[,-1] # first column is redundant
# Diabetes prevalence raw value - multiply by 10 
dataEnvironmental$Diabetes=10*dataEnvironmental$Diabetes
# Adult obesity raw value - multiply by 10 
dataEnvironmental$Obesity=10*dataEnvironmental$Obesity
# number of Beds / Population - multiply by 1000 (i.e. number of Beds per 1000 persons)
dataEnvironmental$PctBeds=1000*dataEnvironmental$PctBeds



