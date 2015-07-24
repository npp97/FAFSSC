require(doBy)

setwd("D://东北//data//碳专项")
fr.raw<-read.csv("碳专项细根现存量调查20150723.csv")
fr.1<-summaryBy(live_fr_g_m2+dead_fr_g_m2+total_fr_g_m2~pplot+soil_layer,data=fr.raw,FUN=mean,na.rm=T)
fr.2<-summaryBy(live_fr_g_m2.mean+dead_fr_g_m2.mean+total_fr_g_m2.mean~pplot,data=fr.1,FUN=sum,na.rm=T)
write.csv(fr.2,"碳专项细根现存量调查By_plot_20150723.csv",row.names=F)

fr.2<-read.csv("碳专项细根现存量调查By_plot_20150723.csv")
abf<-read.csv('D://东北//data//基础信息//TRMP//afb.csv')

fr.3<-merge(abf,fr.2)

fr.3$ldr<-fr.3$dead_fr_tC_ha/fr.3$live_fr_tC_ha

write.csv()

