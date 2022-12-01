# map of regions with Sally track

#track shape file from https://www.nhc.noaa.gov/gis/best_track/al192020_best_track.zip
# cone shape file from https://www.nhc.noaa.gov/gis/
library(tigris)
library(stringi)
library(rgdal)
library(maptools)  #getKMLcoordinates
library(spdep)
library(scales)
#library(plotKML)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path1='../data/'
reg.table=read.csv(paste(path1,"FL_REGIONS.csv",sep="")) 
reg.table$FIPS=stri_pad_left(reg.table$FIPS, 3, 0)

# shape file for FL counties 
options(tigris_use_cache = TRUE) #ZCTAs can take several minutes to download. To cache the data use this.
fl_counties=as_Spatial(counties("Florida",year=2018,cb=TRUE))

# shape files for Sally 
path2='../data/HurricaneSallyTrack/Sally192020_best_track/'
path3='../data/HurricaneSallyTrack/Sally192020_5day_017A/'

# best track
sally.line <- readOGR(dsn = paste(path2, "AL192020_lin.shp",sep=""), stringsAsFactors = F)
sally.points <- readOGR(dsn = paste(path2, "AL192020_pts.shp",sep=""), stringsAsFactors = F)
#advisory 17 a and the uncertainty cone
sally.5day.pgn <- readOGR(dsn = paste(path3, "al192020-017A_5day_pgn.shp",sep=""), stringsAsFactors = F)
sally.5day.ln <- readOGR(dsn = paste(path3, "al192020-017A_5day_lin.shp",sep=""), stringsAsFactors = F)
sally.5day.pts <- readOGR(dsn = paste(path3, "al192020-017A_5day_pts.shp",sep=""), stringsAsFactors = F)


# regions
fl_counties$region=reg.table$Region[match(fl_counties$COUNTYFP,reg.table$FIPS)]
cols=1:7
par(mfrow=c(1,1),mar=c(.0,.01,0.,0)+0.0001,mgp=c(.5,.5,0),cex.lab=1,cex.main=1,cex.axis=1,cex=1)
plot(fl_counties, border='gray',
      col=alpha(cols[fl_counties$region],0.5),axes = FALSE)
leg_str=1:7
legend("bottomleft",col="white",legend =leg_str, fill = alpha(cols,0.5), 
        bty = "o", title = "Region",cex=0.75,
        x.intersp = 1, y.intersp = 0.75,ncol=2)
text(x=coordinates(fl_counties)[,1],y=coordinates(fl_counties)[,2],fl_counties$NAME,cex=0.5) 
plot(sally.5day.ln,add=TRUE)
plot(sally.5day.pts,add=TRUE,col='red')
plot(sally.5day.pgn,add=TRUE)

