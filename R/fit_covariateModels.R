# Read July to Oct 2020 daily zip code level covid data. 
# Read census tract level demographic variables and county level health variables
# Choose demographic and health variables for linear mixed models of region 1 (panhandle) 
# and region 4 while including the effect of hurricane Sally (Sep 11-18, 2020)
#### read data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('load_data.R')
noDays=N
library(stringi)
library(lattice)
library(lme4)
library(nlme)
library(caret)
library(corrplot)
library(ggplot2)
library(biscale)
library(cowplot)
library(sf)
#### functions ####
createPlots=function (lme_model,dframe,covar_name,flag.predictPlot,flag.mapPlot){
  inter.str=covar_name
  if (flag.predictPlot){
    hat_lme=exp(predict(lme_model))*dframe$pop-1
    print(xyplot(hat_lme~ day|COUNTY, group=tractCode,
                 layout = c(4,4), 
                 type="l", lwd=2,
                 data=dframe,
                 cex=0.5,cex.label=0.5,
                 ylab='Predicted Count'))  
    plot(hat_lme,dframe$count)
    abline(a=0,b=1,col='grey')
  }
  if (flag.mapPlot){
    # map of coefficients: Sally for census tracts of the study region
    ff=ranef(lme_model)
    if (class(lme_model)[1]=="lmerMod"){
      coeffs=as.numeric(exp(ff$tractCode[[inter.str]]))   # with lme4
    }else if(class(lme_model)[1]=="lme"){
      coeffs=as.numeric(exp(ff[[inter.str]]))            # with nlme
    }
    n_tract=length(indNonZeroTractRegion)
    coeffs_all=matrix(1.0,nrow=n_tract) # coeff=0 means e(coef)=1
    coeffs_all[!(1:n_tract)%in%tractCodeRemove]=coeffs
    #cols <- hcl.colors(4, "heat",alpha=0.8,rev=TRUE)
    cols <- hcl.colors(5, "oranges",alpha=0.8,rev=TRUE)  #see hcl.pals("sequential")
    brks <- quantile(coeffs, c(0.35, 0.50, 0.65, 0.80,0.90))
    par(mfrow=c(1,1),pty="s",mar=c(.0,.01,0.,0)+0.0001,mgp=c(.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1,cex=1)
    plot(fl_tracts[indNonZeroTractRegion,],border='grey')
    #tractCodeRemove are tracts with 0 covid counts therefore coeff's not estimated
    plot(fl_tracts[indNonZeroTractRegion,], border=NA,add=TRUE,
         col=alpha(cols[findInterval(coeffs_all, brks, all.inside=TRUE)],0.5),axes = FALSE)
    title(paste("coefficient of ",covar_name,sep=""),line=-2)
    leg_str=list()
    for (i in 1:length(brks)-1) {
      leg_str[i]=paste(as.character(format(brks[i],digits=3,nsmall=3)),"-",
                       as.character(format(brks[i+1],digits=3,nsmall=3)))}
    legend("bottomleft",col="white",legend =leg_str, fill = cols, 
           bty = "o", title = "Coefficient",cex=0.75,
           x.intersp = 1, y.intersp = 0.75)
  }
}
createBivariate=function (lme_model,dframe,covar_name1,covar_name2,noBreaks){
  inter.str1=covar_name1
  inter.str2=covar_name2
  # map of coefficients: Sally for census tracts of the study region
  ff=ranef(lme_model)
  if (class(lme_model)[1]=="lmerMod"){
    coeffs1=as.numeric(exp(ff$tractCode[[inter.str1]]))
    coeffs2=as.numeric(exp(ff$tractCode[[inter.str2]])) # with lme4
  }else if(class(lme_model)[1]=="lme"){
    coeffs1=as.numeric(exp(ff[[inter.str1]]))  
    coeffs2=as.numeric(exp(ff[[inter.str2]]))  # with nlme
  }
  n_tract=length(indNonZeroTractRegion)
  coeffs_all1=matrix(1.0,nrow=n_tract) 
  coeffs_all1[!(1:n_tract)%in%tractCodeRemove]=coeffs1
  coeffs_all2=matrix(1.0,nrow=n_tract) 
  coeffs_all2[!(1:n_tract)%in%tractCodeRemove]=coeffs2
  fl_tracts_sf=st_as_sf(fl_tracts[indNonZeroTractRegion,], coords = c("x", "y"), crs = 28992, agr = "constant")
  fl_tracts_sf$coef1=as.numeric(coeffs_all1)
  fl_tracts_sf$coef2=as.numeric(coeffs_all2)
  df.map.biclass=bi_class(fl_tracts_sf,x=coef1,y=coef2,style = "quantile", dim = noBreaks)
  
  if (noBreaks==3){mypalette="GrPink"}else{mypalette="GrPink2"}
  map=ggplot() +
    geom_sf(data = df.map.biclass, mapping = aes(fill = bi_class), 
            color = "white", size = 0.1, show.legend = FALSE) +
    bi_scale_fill(pal = mypalette, dim = noBreaks)+  bi_theme()
  
  legend <- bi_legend(pal = mypalette,
                      dim = noBreaks,
                      xlab = paste("Higher ",covar_name1,sep=""),
                      ylab = paste("Higher ",covar_name2,sep=""),
                      size = 8)
  
  # combine map with legend
  finalPlot <- ggdraw() +
    draw_plot(map, 0, 0, 1, 1) +
    draw_plot(legend, 0.1, .1, 0.3, 0.3)
  finalPlot
  return(list(finalPlot=finalPlot,bi_class=df.map.biclass$bi_class))
}


#### create a global data frame ####
path1='../data/'
reg.table=read.csv(paste(path1,"FL_REGIONS.csv",sep="")) 
reg.table$FIPS=stri_pad_left(reg.table$FIPS, 3, 0)
# create a region variable
fl_tracts$region=reg.table$Region[match(fl_tracts$COUNTYFP,reg.table$FIPS)]
#  choose a region between 1,2,.., 7 to consider
region.no=4  # 1 for panhandle.
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
# names of demographic variables of interest
varnames=c("EP_POV","EP_UNEMP","EP_PCI","EP_NOHSDP","EP_AGE65","EP_AGE17","EP_DISABL","EP_SNGPNT",
           "EP_MINRTY","EP_LIMENG","EP_MUNIT","EP_MOBILE","EP_CROWD","EP_NOVEH","EP_GROUPQ")
x=X[indNonZeroTractRegion,]
noTracts=dim(x)[1]
# environmental variables for all counties for the first day
envVarsDay=dataEnvironmental[dataEnvironmental$day==1,]
# replace missing percentBeds with median
envVarsDay$PctBeds[is.na(envVarsDay$PctBeds)]=median(envVarsDay$PctBeds,na.rm = TRUE)
varnames.env=c("AirPollution","PctNoIns","Diabetes","Obesity","PctBeds") # for region 1 and 4
# match tracts with counties
indTractCounty=match(x$COUNTY,envVarsDay$county)
data.repDays=data.frame(x,count=round(y[,1],),day=1,tractCode=1:noTracts,
                        pop=fl_tracts$pop[indNonZeroTractRegion],envVarsDay[indTractCounty,])
data.repDays$tractCode=as.factor(data.repDays$tractCode)
data.repDays$COUNTY=as.factor(data.repDays$COUNTY)
for (i in 2:noDays){
  envVarsDay=dataEnvironmental[dataEnvironmental$day==i,]
  envVarsDay$PctBeds[is.na(envVarsDay$PctBeds)]=median(envVarsDay$PctBeds,na.rm = TRUE)
  data.repDays=rbind(data.repDays,
                          data.frame(x,count=round(y[,i],0),day=i,tractCode=1:noTracts,
                                     pop=fl_tracts$pop[indNonZeroTractRegion],envVarsDay[indTractCounty,]))
}
# find 0 count tracts
aa=data.repDays[c('count','tractCode')]
laa=unique(aa$tractCode)
ss=NA;for(i in 1:length(laa)){ss[i]=sum(aa$count[aa$tractCode==laa[i]])}
tractCodeRemove=aa$tractCode[which(ss==0)]
# transformed counts
data.repDays$count_t=log((data.repDays$count+1)/data.repDays$pop) 
# center the day variable
data.repDays$day.center=data.repDays$day-mean(data.repDays$day)
# covariates
varnames3=c(varnames,varnames.env)
# factor for the time period of Sally
beginHurr=which(dateVec=='09112020')
endHurr=which(dateVec=='09182020')
periodHurr=seq(beginHurr,endHurr)
data.repDays$Sally=0
data.repDays$Sally[data.repDays$day>beginHurr]=1       # sustained shift, for transient shift use  [data.repDays$day%in%periodHurr]=1 
data.repDays$Sally=as.factor(data.repDays$Sally)
data.repDays$trendSally=(data.repDays$day-beginHurr)*I(data.repDays$day>beginHurr) # piecewise linear trend after Sally begin
data.repDays$trendSallySq=(data.repDays$day-beginHurr)^2*I(data.repDays$day>beginHurr) # quadratic spline after Sally begin
#### correlation between covariates ####
correlations=cor(data.repDays[varnames3])
corrplot(correlations,order = "hclust",method='square',type='upper',tl.cex=.5)
highCorr <- findCorrelation(correlations, cutoff = .60) # for region 4, replace POV with NOVEH
if (region.no==4){highCorr[6]=which(varnames3=='EP_NOVEH')}
correlations2=cor(data.repDays[varnames3[-highCorr]])
corrplot(correlations2,order = "hclust",method='square',type='upper',tl.cex=.5)
highCorr2 <- findCorrelation(correlations2, cutoff = .60) # no remaining correlated covariates
print(paste('variables dropped = '));print(varnames3[highCorr])
varnames3drop=varnames3[-highCorr]
#### plot covid trajectories ####
xyplot(count ~ day | COUNTY, group=tractCode,data.repDays[!data.repDays$tractCode%in%tractCodeRemove,], #log((count+1)/pop), # data.repDays[data.repDays$COUNTY!='Calhoun',]
       #layout = c(4,4), 
       type = c("g", "p"),
       xlab = "Day",
       ylab = "count",
       cex=0.5,
       pch=18,
       alpha=0.28,
       abline=c(v=beginHurr))
#### linear mixed effects models - region 1 ####
# 1. as many explanatory variables as possible (beyond optimal model) in the fixed component; no random effect
fmla <- as.formula(paste("count_t ~day.center+I(day.center^2)+trendSallySq+",
                         paste(varnames3drop, collapse= "+",sep="")))
data.repDaysRm=data.repDays[!data.repDays$tractCode%in%tractCodeRemove,]
lme_mod6=gls(fmla,data = data.repDaysRm)
anova(lme_mod6)
## full model in fixed, random coeff for linear spline 
fmla.fixed <- as.formula(paste("count_t ~day.center+I(day.center^2)+trendSallySq+",
                               paste(varnames3drop, collapse= "+",sep="")))
lmc=lmeControl(maxIter=20000,msMaxIter=1000, msMaxEval = 1000)
lme_mod610 <- lme(fmla.fixed,random =list(tractCode=~day.center),
                  data = data.repDaysRm,control=lmc) 
summary(lme_mod610)
anova(lme_mod6,lme_mod610)
lme_mod611 <- lme(fmla.fixed,random =list(tractCode=~day.center+trendSally),
                  data = data.repDaysRm,control=lmc)  
anova(lme_mod610,lme_mod611) #better
(int_mod611=intervals(lme_mod611)) # conf intervals of sigmas
summary(lme_mod611)
createPlots(lme_mod611,data.repDaysRm,"",TRUE,FALSE)
# 2. find optimal structure of the random component
# random slope for EP_GROUPQ         
lme_mod612 <- lme(fmla.fixed,random =list(tractCode=~EP_GROUPQ+day.center+trendSally),
                  data = data.repDaysRm,control=lmc) 
anova(lme_mod611,lme_mod612) # better
(int_mod612=intervals(lme_mod612)) #errors
summary(lme_mod612)
# random slope for EP_GROUPQ and Air Pollution
lme_mod613 <- lme(fmla.fixed,random =list(tractCode=~AirPollution+EP_GROUPQ+day.center+trendSally),
                  data = data.repDaysRm,control=lmc)
anova(lme_mod612,lme_mod613) #not better
# random slope for EP_GROUPQ and EP_LIMENG  
lme_mod614 <- lme(fmla.fixed,random =list(tractCode=~EP_LIMENG+EP_GROUPQ+day.center+trendSally),
                  data = data.repDaysRm,control=lmc) 
anova(lme_mod612,lme_mod614) #not better
# 3. drop insignificant fixed effects based on the optimal random structure (lme_mod612)
lme_mod615 <- lme(count_t ~day.center+I(day.center^2)+trendSallySq+EP_GROUPQ+
                    AirPollution+EP_LIMENG+EP_MUNIT+PctBeds,
                  random =list(tractCode=~EP_GROUPQ+day.center+trendSally),
                  data = data.repDaysRm,control=lmc)
anova(lme_mod612,lme_mod615) #better
intervals(lme_mod615)#errors
summary(lme_mod615)
# lme4 to get intervals
lme_mod617=lmer(count_t ~day.center+I(day.center^2)+trendSallySq+EP_GROUPQ+
                  AirPollution+EP_LIMENG+EP_MUNIT+PctBeds+
                  (EP_GROUPQ+day.center+trendSally|tractCode),data = data.repDaysRm)
summary(lme_mod617)
ci_617=confint(lme_mod617,method='boot',oldNames=FALSE)
AIC(lme_mod617)
# random slope on quadratic spline
lme_mod617q <- lmer(count_t ~day.center+I(day.center^2)+trendSallySq+EP_GROUPQ+
                     AirPollution+EP_LIMENG+EP_MUNIT+PctBeds+
                      (EP_GROUPQ+I(day.center^2)+trendSallySq|tractCode),
                   data = data.repDaysRm)
summary(lme_mod617q)
AIC(lme_mod617q,lme_mod617) #better
ci_617q=confint(lme_mod617q,method='boot',oldNames=FALSE) #todo

# maps with random slopes
#linear spline
createPlots(lme_mod617,data.repDaysRm,"day.center",FALSE,TRUE) #lme_mod615 gives same results
createPlots(lme_mod617,data.repDaysRm,"trendSally",FALSE,TRUE)
createPlots(lme_mod617,data.repDaysRm,"EP_GROUPQ",FALSE,TRUE)
createPlots(lme_mod617,data.repDaysRm,"",TRUE,FALSE)
noBreaks=3
res=createBivariate(lme_mod617,data.repDaysRm,"trendSally","EP_GROUPQ",noBreaks)
png("reg1-TrendGROUPQ.png", width = 8, height =6, units = 'in', res = 600)
print(res$finalPlot)
dev.off()
#quad spline
createPlots(lme_mod617q,data.repDaysRm,"trendSallySq",FALSE,TRUE)
createPlots(lme_mod617q,data.repDaysRm,"EP_GROUPQ",FALSE,TRUE)
# fitted curves and data
lme_mod5.linear=lmer(count_t ~day.center+I(day.center^2)+trendSallySq+
                  (day.center+trendSally|tractCode),data = data.repDaysRm)
lme_mod5.quad=lmer(count_t ~day.center+I(day.center^2)+trendSallySq+
                       (day.center+I(day.center^2)+trendSallySq|tractCode),data = data.repDaysRm)
#hat_lme=exp(predict(lme_mod612))*data.repDaysRm$pop-1        # optimal fixed and random structure, linear temporal
#hat_lme=exp(predict(lme_mod5.linear))*data.repDaysRm$pop-1   # no covar, linear temporal spline
hat_lme=exp(predict(lme_mod5.quad))*data.repDaysRm$pop-1     # no covar, quad temporal spline

df.combined=rbind(data.frame(data.repDaysRm,y=data.repDaysRm$count,type='1data'),
                  data.frame(data.repDaysRm,y=hat_lme,type='0prediction'))
png("reg1-pred and observed.png", width = 8, height =8, units = 'in', res = 600)

hurrTime=rep(beginHurr,2*dim(data.repDaysRm)[1])
ggplot(df.combined, aes(day,y))+ 
  geom_line(aes(linetype=type,colour=tractCode))+ facet_wrap(df.combined$COUNTY,scales='free')+
  theme(legend.position = "none") + geom_vline(xintercept=hurrTime, linetype="dashed", color = "black")+
  xlab("Day") + ylab("Daily Covid case count per tract")
dev.off()
# plot EP_GROUPQ
coeffs=x$EP_GROUPQ            
cols <- hcl.colors(5, "Mint",alpha=0.8,rev=TRUE)
brks <- quantile(coeffs, c(0.25,0.50,0.75,0.80,0.95))
png("reg1-EP_GROUPQ.png", width = 8, height =6, units = 'in', res = 600)
par(mfrow=c(1,1),mar=c(.0,.01,0.,0)+0.0001,mgp=c(.001,.001,0),cex.lab=1,cex.main=1,cex.axis=1,cex=1)
plot(fl_tracts[indNonZeroTractRegion,],border='grey')
plot(fl_tracts[indNonZeroTractRegion,], border=NA,add=TRUE,
     col=alpha(cols[findInterval(coeffs, brks, all.inside=TRUE)],0.5),axes = FALSE)
# show centroids of tracts that fall in the high-high classEP_GROUPQ and TrendSally
points(coordinates(fl_tracts[indNonZeroTractRegion[res$bi_class=='3-3'],]),pch=21,col=NA,bg=alpha('blue',0.5),cex=0.5)
leg_str=list()
for (i in 1:length(brks)-1) {leg_str[i]=paste(as.character(format(brks[i],digits=3,nsmall=2)),
                                              "-",as.character(format(brks[i+1],digits=3,nsmall=2)))}
legend("bottomleft",col="white",legend =leg_str, fill = cols, 
       bty = "o", title = "EP_GROUPQ Value",cex=0.75,
       x.intersp = 1, y.intersp = 0.75)
dev.off()
#fitted vs. residuals, colored according to the tractCode variable:
plot(lme_mod614,col=data.repDaysRm$tractCode)
# Q-Q plot (a straight line indicates normality)
E <- residuals(lme_mod614,type='pearson')
qqnorm(E,col=data.repDaysRm$tractCode,las=1)
qqline(E)
qqnorm(lme_mod614,~resid(.))
plot(lme_mod614, resid(., type = "p")~fitted(.), xlab = "Predicted value", ylab = "Residual",abline = 0)

boxplot(E~tractCode, data = data.repDaysRm, axes = FALSE,ylim = c(-5, 5))
abline(0,0); axis(2)
text(1:length(unique(data.repDays$tractCode)), -2.5, levels(data.repDays$tractCode), cex=0.75, srt=65)
#### linear mixed effects models - region 4 ####
# 1. as many explanatory variables as possible (beyond optimal model) in the fixed component; no random effect
fmla.fixed <- as.formula(paste("count_t ~day.center+I(day.center^2)+trendSallySq+",
                               paste(varnames3drop, collapse= "+",sep="")))
data.repDaysRm=data.repDays[!data.repDays$tractCode%in%tractCodeRemove,]
lme_mod6=gls(fmla.fixed,data = data.repDaysRm)
anova(lme_mod6)
## full model in fixed, random coeff for linear spline 
lme_mod61111=lmer(count_t ~ day.center + I(day.center^2) + trendSallySq + EP_POV + 
                    EP_UNEMP + EP_PCI + EP_AGE17 + EP_DISABL + EP_LIMENG + EP_MUNIT + 
                    EP_MOBILE + EP_CROWD + EP_GROUPQ + PctNoIns + PctBeds+
                    (day.center+trendSally|tractCode),
                  data = data.repDaysRm,REML=FALSE)
summary(lme_mod61111)
createPlots(lme_mod61111,data.repDaysRm,"",TRUE,FALSE)
# 2. optimal random structure. random slope for PctBeds, EP_POV, etc.
lme_mod61211=lmer(count_t ~ day.center + I(day.center^2) + trendSallySq + EP_POV + 
                    EP_UNEMP + EP_PCI + EP_AGE17 + EP_DISABL + EP_LIMENG + EP_MUNIT + 
                    EP_MOBILE + EP_CROWD + EP_GROUPQ + PctNoIns + PctBeds+
                    (PctBeds+day.center+trendSally|tractCode),
                  data = data.repDaysRm,REML=FALSE)
lme_mod61311=lmer(count_t ~ day.center + I(day.center^2) + trendSallySq + EP_POV + 
                    EP_UNEMP + EP_PCI + EP_AGE17 + EP_DISABL + EP_LIMENG + EP_MUNIT + 
                    EP_MOBILE + EP_CROWD + EP_GROUPQ + PctNoIns + PctBeds+
                    (EP_POV+day.center+trendSally|tractCode),
                  data = data.repDaysRm,REML=FALSE)
lme_mod61411=lmer(count_t ~ day.center + I(day.center^2) + trendSallySq + EP_POV + 
                    EP_UNEMP + EP_PCI + EP_AGE17 + EP_DISABL + EP_LIMENG + EP_MUNIT + 
                    EP_MOBILE + EP_CROWD + EP_GROUPQ + PctNoIns + PctBeds+
                    (EP_PCI+day.center+trendSally|tractCode),
                  data = data.repDaysRm,REML=FALSE)

lme_mod61511=lmer(count_t ~ day.center + I(day.center^2) + trendSallySq + EP_POV + 
                    EP_UNEMP + EP_PCI + EP_AGE17 + EP_DISABL + EP_LIMENG + EP_MUNIT + 
                    EP_MOBILE + EP_CROWD + EP_GROUPQ + PctNoIns + PctBeds+
                    (EP_CROWD+day.center+trendSally|tractCode),
                  data = data.repDaysRm,REML=FALSE)
anova(lme_mod61111,lme_mod61211)#not better
anova(lme_mod61111,lme_mod61311)#not better
anova(lme_mod61111,lme_mod61411)#not better
anova(lme_mod61111,lme_mod61511)#better

summary(lme_mod61511)

# 3. drop insignificant fixed effects based on the optimal random structure (lme_mod61511)
lme_mod6161=lmer(count_t ~ day.center + I(day.center^2) + trendSallySq + 
                    EP_DISABL + EP_LIMENG  + EP_CROWD  + PctBeds+ (EP_CROWD+day.center+trendSally|tractCode),
                  data = data.repDaysRm,REML=FALSE)
anova(lme_mod61511,lme_mod6161) #better
summary(lme_mod6161)
ci_6161=confint(lme_mod6161,method='boot',oldNames=FALSE)
# maps with random slopes
createPlots(lme_mod6161,data.repDaysRm,"day.center",FALSE,TRUE)
createPlots(lme_mod6161,data.repDaysRm,"trendSally",FALSE,TRUE)
createPlots(lme_mod6161,data.repDaysRm,"EP_CROWD",FALSE,TRUE)
createPlots(lme_mod6161,data.repDaysRm,"",TRUE,FALSE)
noBreaks=3
res=createBivariate(lme_mod6161,data.repDaysRm,"trendSally","EP_CROWD",noBreaks)
print(res$finalPlot)

# with nlme -omit UNEMP
lme_mod6161n=lme(count_t ~ day.center + I(day.center^2) + trendSallySq + 
                   EP_DISABL + EP_LIMENG  + EP_CROWD  + PctBeds, 
                    random=list(tractCode=~EP_CROWD+day.center+trendSally),
                 data = data.repDaysRm,control=lmc) #converge error
summary(lme_mod6161n) #converge error
# plot EP_CROWD
coeffs=x$EP_CROWD            
cols <- hcl.colors(5, "Red-Purple",alpha=0.8,rev=TRUE)
brks <- quantile(coeffs, c(0.25,0.50,0.75,0.80,0.95))
png("reg4-EP_CROWD.png", width = 5.25, height =7, units = 'in', res = 600)
par(mfrow=c(1,1),mar=c(.0,.01,0.,0)+0.0001,mgp=c(.001,.001,0),cex.lab=1,cex.main=1,cex.axis=1,cex=1)
plot(fl_tracts[indNonZeroTractRegion,],border='grey')
plot(fl_tracts[indNonZeroTractRegion,], border=NA,add=TRUE,
     col=alpha(cols[findInterval(coeffs, brks, all.inside=TRUE)],0.5),axes = FALSE)
# show centroids of tracts that fall in the high-high class for EP_CROWD and TrendSally
points(coordinates(fl_tracts[indNonZeroTractRegion[res$bi_class=='3-3'],]),pch=21,col=NA,bg=alpha('blue',0.5),cex=0.5)
leg_str=list()
for (i in 1:length(brks)-1) {leg_str[i]=paste(as.character(format(brks[i],digits=3,nsmall=2)),
                                              "-",as.character(format(brks[i+1],digits=3,nsmall=2)))}
legend("bottomleft",col="white",legend =leg_str, fill = cols, 
       bty = "o", title = "EP_CROWD Value",cex=0.75,
       x.intersp = 1, y.intersp = 0.75)
# fitted curves and data
lme_mod5.linear=lmer(count_t ~day.center+I(day.center^2)+trendSallySq+
                       (day.center+trendSally|tractCode),data = data.repDaysRm)
hat_lme=exp(predict(lme_mod5.linear))*data.repDaysRm$pop-1   # no covar, linear temporal spline
df.combined=rbind(data.frame(data.repDaysRm,y=data.repDaysRm$count,type='1data'),
                  data.frame(data.repDaysRm,y=hat_lme,type='0prediction'))
hurrTime=rep(beginHurr,2*dim(data.repDaysRm)[1])
ggplot(df.combined, aes(day,y))+ 
  geom_line(aes(linetype=type,colour=tractCode))+ facet_wrap(df.combined$COUNTY,scales='free')+
  theme(legend.position = "none") + geom_vline(xintercept=hurrTime, linetype="dashed", color = "black")+
  xlab("Day") + ylab("Daily Covid case count per tract")

