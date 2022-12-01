# study the time effect of sally before accounting for demographic and environmental variables 
# Read July to Oct 2020 daily covid data. Mixed models for time effect of hurricane Sally (Sep 11-18, 2020) for all 7 florida regions
#### read data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('load_data.R')
noDays=N
library(stringi)
library(lattice)
library(lme4)
library(nlme)
library(ggplot2)
path1='../data/'
reg.table=read.csv(paste(path1,"FL_REGIONS.csv",sep="")) 
reg.table$FIPS=stri_pad_left(reg.table$FIPS, 3, 0)
#### functions ####
createDataFrame=function(region.no){
  indRegion=which(fl_tracts$COUNTYFP%in%reg.table$FIPS[reg.table$Region==region.no])
  # find tracts of the region with non zero population
  indNonZeroTractRegion=intersect(indRegion,indNonZeroTract)
  # plot study region
  par(mfrow=c(1,1),pty="s",mar=c(3,3,0.5,0)+0.1,mgp=c(1.5,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
  plot(fl_tracts[,1],col='white',border='grey',main="",axes = FALSE)
  plot(fl_tracts[indNonZeroTractRegion,1],border='grey',col='navajowhite2',add=TRUE)
  # matrix for all demographic variables
  X=CDC_SVIr  #(4212  128)
  # create a new data frame where case count for each day is recorded in a column
  y=covid_dailyCaseTract[indNonZeroTractRegion,] # counts
  # demographic covariates
  varnames=c("EP_POV","EP_UNEMP","E_PCI","E_NOHSDP","EP_AGE65","EP_AGE17","EP_DISABL","EP_SNGPNT",
             "EP_MINRTY","EP_LIMENG","EP_MUNIT","EP_MOBILE","EP_CROWD","EP_NOVEH","EP_GROUPQ")
  X[,varnames]=scale(X[,varnames]) # zero mean unit standard deviation scaling
  x=X[indNonZeroTractRegion,]
  noTracts=dim(x)[1]
  # environmental variables for all counties for a given day
  envVarsDay=dataEnvironmental[dataEnvironmental$day==1,]
  # match tracts with counties
  indTractCounty=match(x$COUNTY,envVarsDay$county)
  data.repDays=data.frame(x,count=round(y[,1],),day=1,tractCode=1:noTracts,
                          pop=fl_tracts$pop[indNonZeroTractRegion],envVarsDay[indTractCounty,])
  data.repDays$tractCode=as.factor(data.repDays$tractCode)
  data.repDays$COUNTY=as.factor(data.repDays$COUNTY)
  for (i in 2:noDays){
    envVarsDay=dataEnvironmental[dataEnvironmental$day==i,]
    data.repDays=rbind(data.repDays,
                       data.frame(x,count=round(y[,i],0),day=i,tractCode=1:noTracts,
                                  pop=fl_tracts$pop[indNonZeroTractRegion],envVarsDay[indTractCounty,]))
  }
  # transformed counts
  data.repDays$count_t=log((data.repDays$count+1)/data.repDays$pop) 
  # centering of the week
  data.repDays$day.center=data.repDays$day-mean(data.repDays$day)
  # demographic + environmental covariates
  varnames3=c("EP_POV","EP_UNEMP","E_PCI","E_NOHSDP","EP_AGE65","EP_AGE17","EP_DISABL","EP_SNGPNT",
              "EP_MINRTY","EP_LIMENG","EP_MUNIT","EP_MOBILE","EP_CROWD","EP_NOVEH","EP_GROUPQ",
              "AirPollution","PctNoIns","Diabetes","Obesity","Smoking")
  # create a factor variable for the period of Sally
  beginHurr=which(dateVec=='09112020')
  endHurr=which(dateVec=='09182020')
  periodHurr=seq(beginHurr,endHurr)
  data.repDays$Sally=0
  #data.repDays$Sally[data.repDays$day%in%periodHurr]=1  # shift happens during an interval 
  data.repDays$Sally[data.repDays$day>beginHurr]=1       # sustained shift  
  data.repDays$Sally=as.factor(data.repDays$Sally)
  data.repDays$trendSally=(data.repDays$day-beginHurr)*I(data.repDays$day>beginHurr)     # piecewise linear trend after Sally begin
  data.repDays$trendSallySq=(data.repDays$day-beginHurr)^2*I(data.repDays$day>beginHurr) # quadratic spline after Sally begin
  return(data.repDays)
}

createDataFrameAllFL=function(){
  # matrix for all demographic variables
  X=CDC_SVIr  #(4212  128)
  # create a new data frame where case count for each day is recorded in a column
  y=covid_dailyCaseTract[indNonZeroTract,] # counts
  
  # demographic covariates
  varnames=c("EP_POV","EP_UNEMP","E_PCI","E_NOHSDP","EP_AGE65","EP_AGE17","EP_DISABL","EP_SNGPNT",
             "EP_MINRTY","EP_LIMENG","EP_MUNIT","EP_MOBILE","EP_CROWD","EP_NOVEH","EP_GROUPQ")
  X[,varnames]=scale(X[,varnames]) # zero mean unit standard deviation scaling
  x=X[indNonZeroTract,]
  noTracts=dim(x)[1]
  # environmental variables for all counties for a given day
  envVarsDay=dataEnvironmental[dataEnvironmental$day==1,]
  # regions
  fl_tracts$region=reg.table$Region[match(fl_tracts$COUNTYFP,reg.table$FIPS)]
  # match tracts with counties
  indTractCounty=match(x$COUNTY,envVarsDay$county)
  data.repDays=data.frame(x,count=round(y[,1],),day=1,tractCode=1:noTracts,
                          pop=fl_tracts$pop[indNonZeroTract],envVarsDay[indTractCounty,],
                          region=fl_tracts$region[indNonZeroTract],
                          countyfp=fl_tracts$COUNTYFP[indNonZeroTract])
  data.repDays$tractCode=as.factor(data.repDays$tractCode)
  data.repDays$COUNTY=as.factor(data.repDays$COUNTY)
  for (i in 2:noDays){
    envVarsDay=dataEnvironmental[dataEnvironmental$day==i,]
    data.repDays=rbind(data.repDays,
                       data.frame(x,count=round(y[,i],0),day=i,tractCode=1:noTracts,
                                  pop=fl_tracts$pop[indNonZeroTract],envVarsDay[indTractCounty,],
                                  region=fl_tracts$region[indNonZeroTract],
                                  countyfp=fl_tracts$COUNTYFP[indNonZeroTract]))
  }
  # transformed counts
  data.repDays$count_t=log((data.repDays$count+1)/data.repDays$pop) 
  # centering of the week
  data.repDays$day.center=data.repDays$day-mean(data.repDays$day)
  data.repDays$region=as.factor(data.repDays$region)
  # covariates
  varnames3=c("EP_POV","EP_UNEMP","E_PCI","E_NOHSDP","EP_AGE65","EP_AGE17","EP_DISABL","EP_SNGPNT",
              "EP_MINRTY","EP_LIMENG","EP_MUNIT","EP_MOBILE","EP_CROWD","EP_NOVEH","EP_GROUPQ",
              "AirPollution","PctNoIns","Diabetes","Obesity","Smoking")
  # create a factor variable for the period of Sally
  beginHurr=which(dateVec=='09112020')
  endHurr=which(dateVec=='09182020')
  periodHurr=seq(beginHurr,endHurr)
  data.repDays$Sally=0
  #data.repDays$Sally[data.repDays$day%in%periodHurr]=1  # shift happens during an interval 
  data.repDays$Sally[data.repDays$day>beginHurr]=1       # sustained shift  
  data.repDays$Sally=as.factor(data.repDays$Sally)
  data.repDays$trendSally=(data.repDays$day-beginHurr)*I(data.repDays$day>beginHurr) # piecewise linear trend after Sally begin
  data.repDays$trendSallySq=(data.repDays$day-beginHurr)^2*I(data.repDays$day>beginHurr) # quadratic spline after Sally begin
  #result= list(data.repDays=data.repDays)
  return(data.repDays)
}

createPlots=function (lme_model,covar_name,flag.predictPlot,flag.mapPlot){
  indRegion=which(fl_tracts$COUNTYFP%in%reg.table$FIPS[reg.table$Region==region.no])
  # find tracts of the region with non zero population
  indNonZeroTractRegion=intersect(indRegion,indNonZeroTract)
  #inter.str=paste(covar_name,':trendSally',sep="")
  inter.str=covar_name
  if (flag.predictPlot){
    hat_lme=exp(predict(lme_model))*data.repDays$pop-1
    
    print(xyplot(hat_lme~ day|COUNTY, group=tractCode,
                 type="b", 
                 data=data.repDays.categ,cex=0.5,cex.label=0.5))
    plot(hat_lme,data.repDays.categ$count)
    abline(a=0,b=1,col='grey')
  }
  
  if (flag.mapPlot){
    # map of coefficients: Sally for census tracts of the study region
    ff=ranef(lme_model)
    coeffs=as.numeric(ff$tractCode[[inter.str]])
    #par(mfrow=c(1,1),pty="s",mar=c(.0,.01,0.,0)+0.0001,mgp=c(.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1,cex=1)
    cols <- hcl.colors(4, "heat",alpha=0.8,rev=TRUE)
    brks <- quantile(coeffs, c(0.25,0.50,0.75,0.95))
    plot(fl_tracts[indNonZeroTractRegion,],border='grey')
    plot(fl_tracts[indNonZeroTractRegion,], border=NA,add=TRUE,
         col=alpha(cols[findInterval(coeffs, brks, all.inside=TRUE)],0.5),axes = FALSE)
    title(covar_name)
    leg_str=list()
    for (i in 1:length(brks)-1) {leg_str[i]=paste(as.character(format(brks[i],nsmall=2)),"-",as.character(format(brks[i+1],nsmall=2)))}
    legend("bottomleft",col="white",legend =leg_str, fill = cols, 
           bty = "o", title = "Coefficient",cex=0.75,
           x.intersp = 1, y.intersp = 0.75)
  }
  
}



#### one single level model for each region ####

# loop over all regions
region.no01=1
beginHurr=which(dateVec=='09112020')
for (region.no01 in 1:7){
  data.repDays2=createDataFrame(region.no01)
  lme_mod61=lmer(count_t ~day.center+I(day.center^2)+trendSallySq
                 +(day.center+trendSally|tractCode),data = data.repDays2,REML=FALSE)
  
  lme_mod61q=lmer(count_t ~day.center+I(day.center^2)+trendSallySq
                 +(day.center+I(day.center^2)+trendSallySq|tractCode),data = data.repDays2,REML=FALSE)
  
  lme_mod60=lmer(count_t ~day.center+I(day.center^2)+trendSallySq
                 +(day.center|tractCode),data = data.repDays2,REML=FALSE)
  hat_61lme=exp(predict(lme_mod61))*data.repDays2$pop-1
  print(paste('region=',region.no01,sep=''))
  print(AIC(lme_mod61,lme_mod60,lme_mod61q))
  print(BIC(lme_mod61,lme_mod60,lme_mod61q))
  print(anova(lme_mod60,lme_mod61))
  print(paste('population=',sum(data.repDays2[data.repDays2$day==1,]$pop),sep=''))
  print(paste('number of census tracts=',dim(data.repDays2[data.repDays2$day==1,])[1],sep=''))
  }

# fitted and data for a selected region
region.no01=1  # choose 1 or 4
data.repDays2=createDataFrame(region.no01)
lme_mod61=lmer(count_t ~day.center+I(day.center^2)+trendSallySq
               +(day.center+trendSally|tractCode),data = data.repDays2,REML=FALSE)
lme_mod61_fixed=gls(count_t ~day.center+I(day.center^2)+trendSallySq,data = data.repDays2)
AIC(lme_mod61,lme_mod61_fixed)  # much better- shows the need for random slopes for day and sally

hat_61lme=exp(predict(lme_mod61))*data.repDays2$pop-1

# nlme: Population fitted values, based on the fixed-effects estimates alone, are obtained setting the level argument to 0 (zero).
# lme4: if re.form is ~0 or NA, then the predictions are made at the population level (all random-effect values set to zero)
#hat_61lme_pop=exp(fitted(lme_mod61,re.form=NA))*data.repDays2$pop-1
hat_61lme_pop=exp(fitted(lme_mod61_fixed))*data.repDays2$pop-1


xyplot(count~ day| COUNTY,group=tractCode,data.repDays2, #log((count+1)/pop)
       #layout = c(4,4), 
       #aspect = "xy",
       type = c("p"),pch=1,grid=TRUE,
       xlab = "Day",distribute.type = TRUE,
       ylab = "count",
       abline=c(v=beginHurr),cex.title=0.25
)
xyplot(hat_61lme~ day| COUNTY,group=tractCode,data.repDays2, #log((count+1)/pop)
       #layout = c(4,4), 
       #aspect = "xy",
       type = c("l"),pch=1,grid=TRUE,
       xlab = "Day",distribute.type = TRUE,
       ylab = "count",
       abline=c(v=beginHurr),cex.title=0.25
)

# pred vs observed
plot(hat_61lme,data.repDays2$count,xlab='Predicted',ylab='Observed',col=alpha('grey',0.2),pch=16)
points(hat_61lme_pop,data.repDays2$count,xlab='Predicted',ylab='Observed',col=alpha('darkgrey',0.2),pch=16)

abline(0,1)
rSq=cor(hat_61lme,data.repDays2$count)^2
# pred vs time
plot(count~day.center,data.repDays2)
points(hat_61lme~data.repDays2$day.center,col='red')
points(hat_61lme_pop~data.repDays2$day.center,col='green')




