#to correcte the influence of canopy on Ts

require(doBy)

factor_Tscorr<-function(lai)
{
  return(1.23*exp(-0.06*lai));
  
}

qa<-function(x,...){
  return(c(range(x,...),mean(x,...)))
}


read.csv("D://东北//data//modis//MOD15A2_Cproj_Data_QC.csv")->lai_cproj
lai_cproj$lai_f<-factor_Tscorr(lai_cproj$Lai_1km)
s.lai_cproj<-summaryBy(Lai_1km+lai_f ~ ydID, FUN= qa, na.rm=T, data=lai_cproj)

read.csv("D://东北//data//modis//MOD15A2_ZZZY_Data_QC.csv")->lai_zzzy
lai_zzzy$lai_f<-factor_Tscorr(lai_zzzy$Lai_1km)
s.lai_zzzy<-summaryBy(Lai_1km+lai_f ~ ydID, FUN= qa, na.rm=T, data=lai_zzzy)

s.lai<-rbind(s.lai_zzzy,s.lai_cproj)
names(s.lai)<-c("ydID","lai0.05","lai0.95","lai_mean","flai0.05","flai0.95","flai_mean")
write.csv(s.lai,"s_lai.csv",row.names=FALSE)



