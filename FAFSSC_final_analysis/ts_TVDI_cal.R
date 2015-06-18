# TODO: Add comment
# 
# Author: junhui.zhang
###############################################################################
#Packages Declaration
#------------------------------------------------------------------------------
require(doBy)
require(chron)

#User-defined Function
#------------------------------------------------------------------------------

dayofyear <-function(year,month,day){
	year_base<-julian(as.Date(paste(year,1,1,sep='-')))
	julian_date<-julian(as.Date(paste(year,month,day,sep='-')))
	dayofyear1<-(julian_date-year_base+1)
	return(as.numeric(dayofyear1))
}


str2ymd<-function(str_ser,format="%y/%m/%d"){
	#split one string into Date and Time series
#	as.numeric(unlist(split(unlist(strsplit(origin,split='-')),1:3)))->orginc
	unlist(split(unlist(strsplit(as.character(str_ser),split=' ')),1:2)[1])->dates
#	unlist(split(unlist(strsplit(as.character(str_ser),split=' ')),1:2)[2])->times
	unlist(split(unlist(strsplit(as.character(dates),split='/')),1:3)[1])->year
	unlist(split(unlist(strsplit(as.character(dates),split='/')),1:3)[2])->month
	unlist(split(unlist(strsplit(as.character(dates),split='/')),1:3)[3])->day
	doy<-dayofyear(as.numeric(year),as.numeric(month),as.numeric(day))
	results<-data.frame(year=year,month=month,day=day,doy=doy)
	return(results)
}


#------------------------------------------------------------------------------
#								Code Body
#------------------------------------------------------------------------------
read.csv("二道气象数据日值.csv")->ta
ta$doy<-dayofyear(ta$year,ta$month,ta$day)

summaryBy(ta_C~year+doy,data=ta,FUN=mean)->ta.day

#-------
dir(path="./长白山通量数据",pattern='.csv',include.dirs = FALSE,full.names = TRUE)->flst

tmp<-read.csv(flst[1])
tmp<-tmp[,c(1:3,8:9)];names(tmp)<-c("year","month","day","ts5cm_C","sm20cm_pct")
tmp$doy<-dayofyear(tmp$year,tmp$month,tmp$day)
tmp1<-tmp[,c("year","month","day","doy","ts5cm_C","sm20cm_pct")]

for (i in 2:7){
	read.csv(flst[i])->tmp
	tmp<-tmp[,c(1:3,8:9)];names(tmp)<-c("year","month","day","ts5cm_C","sm20cm_pct")
	tmp$doy<-dayofyear(tmp$year,tmp$month,tmp$day)
	tmp<-tmp[,c("year","month","day","doy","ts5cm_C","sm20cm_pct")]
	tmp1<-rbind(tmp1,tmp)
}

for (i in 8:12){
	read.csv(flst[i])->tmp
	names(tmp)<-toupper(names(tmp))
	tmp<-tmp[,c("TIMESTAMP","TS_107_1_AVG","SMOIST_0_1_AVG")]
	mc<-as.character(tmp[,1])
	tmtm<-str2ymd(mc)
	tmp<-data.frame(tmtm$year,tmtm$month,tmtm$day,tmtm$doy,tmp[,2:3])
	names(tmp)<-c("year","month","day","doy","ts5cm_C","sm20cm_pct")
	tmp1<-rbind(tmp1,tmp)	
}

ta_doy<-summaryBy(ts5cm_C+sm20cm_pct~year+doy,data=tmp1)
#ta_doy$id<-paste(ta_doy$year,ta_doy$doy,sep='-')

write.csv(ta_doy,'ts_sm_cbs.csv')

#----------------------------------------------------------------------------------------
#To process METEData 
#-------------------

read.csv("二道气象数据日值.csv")->ta.eda
#names(ta.eda)
ta.eda$doy<-dayofyear(ta.eda$year,ta.eda$month,ta.eda$day)

#ta.eda$id<-paste(ta.eda$year,ta.eda$doy,sep='-')

merge(ta.eda,ta_doy,by=c("year","doy"))->tats_cbs_flx1

write.csv(tats_cbs_flx1,'tats_cbs_flx1.csv')

#---------------------------------------------------------------
#Download the LAI,NDVI,EVI,LSTdata
#
require(MODISTools)
require(NCBI2R)

dat <-data.frame(lat=42.824,long=128.828,end.date=as.Date("2015-06-01"))

## Creating my unique ID for each unique site
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")
products <- "MOD13Q1"     #EVI. NDVI
bands<-GetBands(products)

if(file.exists("MOD13Q1_IndividualPixel") == FALSE) 
		dir.create("MOD13Q1_IndividualPixel")
MODISSubsets(dat, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD13Q1_IndividualPixel",TimeSeriesLength = 15)

products <- "MOD11A2"     #LST
bands<-GetBands(products)
if(file.exists("MOD11A2_IndividualPixel") == FALSE) 
	dir.create("MOD11A2_IndividualPixel")
MODISSubsets(dat, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD11A2_IndividualPixel",TimeSeriesLength = 15)

if(file.exists("MOD15A2_IndividualPixel") == FALSE) 
	dir.create("MOD15A2_IndividualPixel")
products <- "MOD15A2"    #LAI
bands<-GetBands(products)
MODISSubsets(dat, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD15A2_IndividualPixel",TimeSeriesLength = 15)

products <- "MOD09A1"    #RELF
bands<-GetBands(products)
if(file.exists("MOD09A1_IndividualPixel") == FALSE) 
	dir.create("MOD09A1_IndividualPixel")
MODISSubsets(dat, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD09A1_IndividualPixel",TimeSeriesLength = 15)

