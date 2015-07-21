#to search TVDI curves using MODIS MOD13Q1 and MOD11A2


qa<-function(x,...){
  return(c(quantile(x,c(0.01,1),...)))
}

sumfun <- function(x, ...) {
  c(m = mean(x, ...), md = median(x, ...), v = sd(x, ...), l = length(x))
}

library(gregmisc)

q99<-function(x,...){return(c(quantile(x,c(0.99),...)))}
q01<-function(x,...){quantile(x,0.01, ...)}
q50<-function(x, ...){median(x, ...)}


require(doBy)
require(strucchange)
require(breakpoint)
require(robust)
require(Hmisc)
# Step merge LST and NDVI
read.csv("D://东北//data//modis//MOD11A2_Data_QC_Cproj_wll.csv")->MOD11A2.cproj
read.csv("D://东北//data//modis//MOD11A2_Data_QC_ZZZYm.csv")->MOD11A2.zzzy

rbind(MOD11A2.cproj,MOD11A2.zzzy)->MOD11A2
MOD11A2<-MOD11A2[which(MOD11A2$doy %in% c(1,17,33,49,65,81,97,113,129,145,161,177,193,209,225,241,257,273,289,305,321,337,353)),]

read.csv("D://东北//data//modis//MOD13Q1_Data_QC_Cproj_wll.csv")->MOD13Q1.cproj
read.csv("D://东北//data//modis//MOD13Q1_Data_QC_ZZZY.csv")->MOD13Q1.zzzy

rbind(MOD13Q1.cproj,MOD13Q1.zzzy)->MOD13Q1

MOD11A2.a<-MOD11A2[,c("ydID","AYearDoy","year","doy","LST_Day_1km")]
MOD13Q1.a<-MOD13Q1[,c("ydID","AYearDoy","year","doy","X250m_16_days_EVI","X250m_16_days_NDVI")]

merge(MOD11A2.a,MOD13Q1.a)->MOD.LST_VI

month<-matrix(NA,ncol=1,nrow=nrow(MOD.LST_VI))
for (i in 1:nrow(MOD.LST_VI)){month[i]<-ydn2md(MOD.LST_VI$year[i],MOD.LST_VI$doy[i])[[1]]}
MOD.LST_VI<-data.frame(MOD.LST_VI,month)

MOD.LST_VI$TVDI.evi<-matrix(NA,ncol=1,nrow=nrow(MOD.LST_VI))
MOD.LST_VI$TVDI.ndvi<-matrix(NA,ncol=1,nrow=nrow(MOD.LST_VI))
MOD.LST_VI$evi_break<-matrix(NA,ncol=1,nrow=nrow(MOD.LST_VI))
MOD.LST_VI$ndvi_break<-matrix(NA,ncol=1,nrow=nrow(MOD.LST_VI))
lst.min_p<-lst.max_p<-matrix(NA,ncol=1,nrow=nrow(MOD.LST_VI))
#to buils TVDImin TVDImax according to doy using LST-NDVI, LST-EVI regression

doy.lst<-c(97,113,129,145,161,177,193,209,225,241,257,273)
l.dy.lst<-length(doy.lst)
#plot(1,type='n',xlim=c(0,1),ylim=c(200,320))

cut(MOD.LST_VI$X250m_16_days_EVI,20)->MOD.LST_VI$evi_break
cut(MOD.LST_VI$X250m_16_days_NDVI,20)->MOD.LST_VI$ndvi_break


plot(1,xlim=c(0,1),ylim=c(0,1),type='n')
for (i in 1:l.dy.lst){
  
  ii<-which(MOD.LST_VI$doy %in% doy.lst[i])
  
#  cut(MOD.LST_VI$X250m_16_days_EVI,breaks=quantile(MOD.LST_VI$X250m_16_days_EVI,seq(0,1,0.025)))->MOD.LST_VI$evi_break
#  cut(MOD.LST_VI$X250m_16_days_NDVI,breaks=quantile(MOD.LST_VI$X250m_16_days_NDVI,seq(0,1,0.05)))->MOD.LST_VI$ndvi_break[ii]


  summaryBy(X250m_16_days_EVI~evi_break,FUN=median, na.rm=T, data=MOD.LST_VI[ii,])->evi.m
  summaryBy(LST_Day_1km~MOD.LST_VI$evi_break,FUN=range, na.rm=T, data=MOD.LST_VI[ii,])->lst_evi.m
  evi.m<-evi.m[1:20,]
  lst_evi.m<-lst_evi.m[1:20,]
 
  obs.max<-data.frame(x=evi.m[,2],y=lst_evi.m[,3])
  obs.max<-obs.max[!is.na(obs.max$x),]
  
  obs.max1=clowess(x=obs.max$x,y=obs.max$y)
 # ij<-which(obs.max1$y<obs.max$y)
#  obs.max1$y[ij]<-obs.max$y[ij]
#  obs.max1=clowess(x=obs.max1$x,y=obs.max1$y)
#  ij<-which(obs.max1$y<obs.max$y)
#  obs.max1$y[ij]<-obs.max$y[ij]
  obs.max=clowess(x=obs.max1$x,y=obs.max1$y)

  obs.min<-data.frame(x=evi.m[,2],y=lst_evi.m[,2])
  obs.min<-obs.min[!is.na(obs.min$x),]

  obs.min1=clowess(x=obs.min$x,y=obs.min$y)
#  ij<-which(obs.min1$y>obs.min$y)
#  obs.min1$y[ij]<-obs.min$y[ij]
#  obs.min1=clowess(x=obs.min1$x,y=obs.min1$y)
#  ij<-which(obs.min1$y>obs.min$y)
#  obs.min1$y[ij]<-obs.min$y[ij]
  obs.min=clowess(x=obs.min1$x,y=obs.min1$y)
  
#   est_pts<-CE.Normal(as.data.frame(obs.max$y),parallel = TRUE)
#   
#   if(is.na(as.numeric(est_pts[1]))){rlg.max<-rlm(y~x,data=obs.max)}else{
#     loc.idx<-min(est_pts$BP.Loc);
#     rlg.max<-rlm(y~x,data=obs.max[c(loc.idx:nrow(obs.max)),]);
#   }
  
  rlg.max<-lmRob(y~x,data=obs.max,control= lmRob.control(weight=c("Bisquare","Optimal")))

  rlg.min<-lmRob(y~x,data=obs.min,control= lmRob.control(weight=c("Bisquare","Optimal")))
  
  lst.min_p[ii]<-predict(rlg.min,newdata=data.frame(x=MOD.LST_VI$X250m_16_days_EVI[ii]))
  lst.max_p[ii]<-predict(rlg.max,newdata=data.frame(x=MOD.LST_VI$X250m_16_days_EVI[ii]))
  
  MOD.LST_VI$TVDI.evi[ii]<-(MOD.LST_VI$LST_Day_1km[ii]-lst.min_p[ii])/(lst.max_p[ii]-lst.min_p[ii])
  
  points(MOD.LST_VI$X250m_16_days_EVI[ii], MOD.LST_VI$TVDI.evi[ii],pch='.',xlim=c(0,1),ylim=c(0,1),col=i)
#   points(obs.min$x,rlg.min$fitted.values,pch='=',col=2,lwd=2)
#   points(obs.max$x,rlg.max$fitted.values,pch='=',col=2,lwd=2)
#   points(evi.m[,2],lst_evi.m[,2],col=i)
#   points(evi.m[,2],lst_evi.m[,3],pch="。",col=i)
 
  print(summary(rlg.min))
  print(summary(rlg.max))
  
  print(i)
}

plot(1,xlim=c(0,1),ylim=c(0,1),type='n')
for (i in 1:l.dy.lst){
  
  ii<-which(MOD.LST_VI$doy %in% doy.lst[i])
  
  summaryBy(X250m_16_days_NDVI~ndvi_break,FUN=mean, na.rm=T, data=MOD.LST_VI[ii,])->ndvi.m
  summaryBy(LST_Day_1km~MOD.LST_VI$ndvi_break,FUN=range, na.rm=T, data=MOD.LST_VI[ii,])->lst_ndvi.m
  ndvi.m<-ndvi.m[which(!is.na(ndvi.m$X250m_16_days_NDVI.mean)),]
  lst_ndvi.m<-lst_ndvi.m[which(!is.na(ndvi.m$X250m_16_days_NDVI.mean)),]
  
  obs.max<-data.frame(x=ndvi.m[,2],y=lst_ndvi.m[,3])
  obs.max1=clowess(x=obs.max$x,y=obs.max$y)
  ij<-which(obs.max1$y<obs.max$y)
  obs.max1$y[ij]<-obs.max$y[ij]
  obs.max1=clowess(x=obs.max1$x,y=obs.max1$y)
  ij<-which(obs.max1$y<obs.max$y)
  obs.max1$y[ij]<-obs.max$y[ij]
  obs.max=clowess(x=obs.max1$x,y=obs.max1$y)
  
  
  obs.min<-data.frame(x=ndvi.m[,2],y=lst_ndvi.m[,2])
  obs.min<-obs.min[!is.na(obs.min$x),]
  
  obs.min1=clowess(x=obs.min$x,y=obs.min$y)
  ij<-which(obs.min1$y>obs.min$y)
  obs.min1$y[ij]<-obs.min$y[ij]
  obs.min1=clowess(x=obs.min1$x,y=obs.min1$y)
  ij<-which(obs.min1$y>obs.min$y)
  obs.min1$y[ij]<-obs.min$y[ij]
  obs.min=clowess(x=obs.min1$x,y=obs.min1$y)
  
  rlg.max<-lmRob(y~x,data=obs.max,control= lmRob.control(weight=c("Bisquare","Optimal")))
  rlg.min<-lmRob(y~x,data=obs.min,control= lmRob.control(weight=c("Bisquare","Optimal")))
  
  
  lst.min_p[ii]<-predict(rlg.min,newdata=data.frame(x=MOD.LST_VI$X250m_16_days_NDVI[ii]))
  lst.max_p[ii]<-predict(rlg.max,newdata=data.frame(x=MOD.LST_VI$X250m_16_days_NDVI[ii]))
  
  MOD.LST_VI$TVDI.ndvi[ii]<-(MOD.LST_VI$LST_Day_1km[ii]-lst.min_p[ii])/(lst.max_p[ii]-lst.min_p[ii])
  points(MOD.LST_VI$X250m_16_days_NDVI[ii], MOD.LST_VI$TDVI.ndvi[ii],pch='.',xlim=c(0,1),ylim=c(0,1),col=i)
}

write.csv(MOD.LST_VI,"MOD_LST_VI_TVDI.csv",row.names = F)
summaryBy(LST_Day_1km+X250m_16_days_EVI+X250m_16_days_NDVI+TVDI.evi+TVDI.ndvi~ydID+month,data=MOD.LST_VI,FUN=sumfun,na.rm=TRUE)->plot_sum_m
summaryBy(LST_Day_1km+X250m_16_days_EVI+X250m_16_days_NDVI+TVDI.evi+TVDI.ndvi~ydID,data=MOD.LST_VI,FUN=sumfun,na.rm=TRUE)->plot_sum_yr
names(plot_sum_yr)[1]<-'ydID1'

plot_sum_m.w<-reshape(plot_sum_m,direction='wide',idvar=c("ydID"),timevar="month")
for (i in 1:3){
	aa<-expand.grid(c(i),c("TVDI.evi","TVDI.ndvi"),c("m","v","md","l"))
	lst.nm.tf<-paste(aa[,2],aa[,3],aa[,1],sep='.')
	aa1<-expand.grid(c(4),c("TVDI.evi","TVDI.ndvi"),c("m","v","md","l"))
	lst.nm.f<-paste(aa1[,2],aa1[,3],aa1[,1],sep='.')
	for(j in 1:length(lst.nm.tf)){plot_sum_m.w[,lst.nm.tf[j]]<-plot_sum_m.w[,lst.nm.f[j]]}
}

for (i in 10:12){
	aa<-expand.grid(c(i),c("TVDI.evi","TVDI.ndvi"),c("m","v","md","l"))
	lst.nm.tf<-paste(aa[,2],aa[,3],aa[,1],sep='.')
	aa1<-expand.grid(c(9),c("TVDI.evi","TVDI.ndvi"),c("m","v","md","l"))
	lst.nm.f<-paste(aa1[,2],aa1[,3],aa1[,1],sep='.')
	for(j in 1:length(lst.nm.tf)){plot_sum_m.w[,lst.nm.tf[j]]<-plot_sum_m.w[,lst.nm.f[j]]}
}

#names(LmBymth.w)<-c('ydID','lat','long',paste('Lm',1:12,sep=''))
plot_sum<-cbind(plot_sum_m.w,plot_sum_yr)
aa<-expand.grid(c(1:12),c("LST_Day_1km","X250m_16_days_EVI","X250m_16_days_NDVI","TVDI.evi","TVDI.ndvi"),c("m","v","md","l"))
lst.nm<-paste(aa[,2],aa[,3],aa[,1],sep='.')
plot_sum<-plot_sum[,c("ydID",lst.nm,names(plot_sum_yr))]
write.csv(plot_sum,"summLST_VI_TVDI.csv",row.names = F)



read.csv('s_lai.csv')->slai

merge(slai,plot_sum)->flai_tvdi

write.csv(flai_tvdi,"flai_tvdi.csv")

read.csv('bio_soc_cls.csv')->bscr

