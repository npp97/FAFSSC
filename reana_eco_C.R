
require(randomForestSRC)
require(randomForest)
require(zoo)
require(Hmisc)

setwd("D://东北//data//基础信息//RESULT//csv")

read.csv('bio_soc_env_all_gsla.csv')->all.gsl
#----------------------------totc&totc2
rfsrc(TotC_eco~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.totc
plot.variable(gfr.totc,"age_max",partial=T,npts=75)->gr.totc

fr.totc<-data.frame(age=gr.totc$pData[[1]]$x.uniq,totc=gr.totc$pData[[1]]$yhat,totc.se=gr.totc$pData[[1]]$yhat.se)
ylim<-floor(range(fr.totc$age))
totc.a<-data.frame(approx(fr.totc$age,fr.totc$totc,ylim[1]:272),totc.ci=I(2*approx(fr.totc$age,fr.totc$totc.se,ylim[1]:272)$y))
names(totc.a)<-c("age","totc","totc.ci")

##############
rfsrc(TotC_eco2~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.totc2
plot.variable(gfr.totc2,"age_max",partial=T,npts=75)->gr.totc2

fr.totc2<-data.frame(age=gr.totc2$pData[[1]]$x.uniq,totc2=gr.totc2$pData[[1]]$yhat,totc2.se=gr.totc2$pData[[1]]$yhat.se)
ylim<-floor(range(fr.totc2$age))
totc2.a<-data.frame(approx(fr.totc2$age,fr.totc2$totc2,ylim[1]:272),totc2.ci=I(2*approx(fr.totc2$age,fr.totc2$totc2.se,ylim[1]:272)$y))
names(totc2.a)<-c("age","totc","totc.ci")

#---
rfo=randomForest::randomForest(TotC_eco~age_max+lat+long,data=all.gsl,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
partialPlot(rfo,all.gsl,age_max,n.pt=75)->rfr.totc

rfo=randomForest::randomForest(TotC_eco2~age_max+lat+long,data=all.gsl,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
partialPlot(rfo,all.gsl,age_max,n.pt=75)->rfr.totc2

totc2.b<-data.frame(approx(rfr.totc2$x,rfr.totc2$y,ylim[1]:272))
names(totc2.b)<-c("age","rfr.totc2")
####################################################################################
with(fr.totc,errbar(age,totc,totc+2*totc.se,totc-2*totc.se,ylim=c(150,350),pch=21))
with(fr.totc2,errbar(age,totc2,totc2+2*totc2.se,totc2-2*totc2.se,col=3,add=T,pch=19))
grid()
######################################################################################

rfsrc(W_total_t_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.bio
plot.variable(gfr.bio,"age_max",partial=T,npts=75)->gr.bio

fr.bio<-data.frame(age=gr.bio$pData[[1]]$x.uniq,bio=gr.bio$pData[[1]]$yhat,bio.se=gr.bio$pData[[1]]$yhat.se)
ylim<-floor(range(fr.bio$age))
bio.a<-data.frame(approx(fr.bio$age,fr.bio$bio,ylim[1]:272),bio.ci=I(2*approx(fr.bio$age,fr.bio$bio.se,ylim[1]:272)$y))
names(bio.a)<-c("age","bio","bio.ci")

#---
rfo=randomForest::randomForest(W_total_t_ha~age_max+lat+long,data=all.gsl,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
partialPlot(rfo,all.gsl,age_max,n.pt=75)->rfr.bio

##############

rfsrc(W_total_t_ha2~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.bio2
plot.variable(gfr.bio2,"age_max",partial=T,npts=75)->gr.bio2

fr.bio2<-data.frame(age=gr.bio2$pData[[1]]$x.uniq,bio2=gr.bio2$pData[[1]]$yhat,bio2.se=gr.bio2$pData[[1]]$yhat.se)
ylim<-floor(range(fr.bio2$age))
bio2.a<-data.frame(approx(fr.bio2$age,fr.bio2$bio2,ylim[1]:272),bio2.ci=I(2*approx(fr.bio2$age,fr.bio2$bio2.se,ylim[1]:272)$y))
names(bio2.a)<-c("age","bio2","bio2.ci")

#----
rfo=randomForest::randomForest(W_total_t_ha2~age_max+lat+long,data=all.gsl,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
partialPlot(rfo,all.gsl,age_max,n.pt=75)->rfr.bio2

bio2.b<-data.frame(approx(rfr.bio2$x,rfr.bio2$y,ylim[1]:272))
names(bio2.b)<-c("age","rfr.bio2")
####################################################################################
with(fr.bio,errbar(age,bio,bio+2*bio.se,bio-2*bio.se,ylim=c(30,110),pch=21))
with(fr.bio2,errbar(age,bio2,bio2+2*bio2.se,bio2-2*bio2.se,col=2,add=T,pch=19))
with(rfr.bio2,errbar(x,y,y*1.01,y*0.99,col=3,add=T,lty=1,lwd=2))
grid()
######################################################################################

rfsrc(SoilC_1m_t_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.soilc
plot.variable(gfr.soilc,"age_max",partial=T,npts=75)->gr.soilc

fr.soilc<-data.frame(age=gr.soilc$pData[[1]]$x.uniq,soilc=gr.soilc$pData[[1]]$yhat,soilc.se=gr.bio$pData[[1]]$yhat.se)
ylim<-floor(range(fr.soilc$age))
soilc.a<-data.frame(approx(fr.soilc$age,fr.soilc$soilc,ylim[1]:272),soilc.ci=I(2*approx(fr.soilc$age,fr.soilc$soilc.se,ylim[1]:272)$y))
names(soilc.a)<-c("age","soilc","soilc.ci")

#----
rfo=randomForest::randomForest(SoilC_1m_t_ha~age_max+lat+long,data=all.gsl,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
partialPlot(rfo,all.gsl,age_max,n.pt=75)->rfr.soilc

soilc.b<-data.frame(approx(rfr.soilc$x,rfr.soilc$y,ylim[1]:272))
names(soilc.b)<-c("age","rfr.soilc")

#############################################################################################
with(fr.totc2,errbar(age,totc2,totc2+5*totc2.se,totc2-5*totc2.se,col=3,pch=19,ylim=c(30,350)))
with(fr.bio2,errbar(age,bio2,bio2+5*bio2.se,bio2-5*bio2.se,col=2,add=T,pch=21))
with(fr.soilc,errbar(age,soilc,soilc+5*soilc.se,soilc-5*soilc.se,col=4,add=T,pch=22))
with(rfr.totc2,points(y~x,type='l',lwd=2,add=T))
with(rfr.bio2,points(y~x,type='l',lwd=2,add=T))
with(rfr.soilc,points(y~x,type='l',lwd=2,add=T))

grid()

#################################

write.csv(cbind(fr.totc2,fr.bio2,fr.soilc,rfr.soilc,rfr.totc2,rfr.bio2),'ecosystem_carbon_age_2.csv',row.names=F)
write.csv(cbind(totc2.a,bio2.a,soilc.a,totc2.b,bio2.b,soilc.b),'ecosystem_carbon_age_2a.csv',row.names=F)