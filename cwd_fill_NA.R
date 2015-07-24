######################################
#TO fill na of cwd status, 2 kind of output:(1) whole 1450 points using randomforest; (2) subset of 250 points using randomforest 
########################################
#setwd("D://东北//data//碳专项")
#read.csv('ccwd_hlj.csv')->cwd_hlj
#read.csv('yd_hlj.csv')->yd_hlj

#cwd<-summaryBy(ccwd~pplot,data=cwd_hlj,FUN=sum)
#merge(yd_hlj,cwd)->ccwd

#abf<-read.csv('D://东北//data//凋落物-草-灌木生物量碳//all_grs_shrb_lit.csv')
#cwd<-read.csv('cwd.csv')

#abf<-merge(abf,cwd,all.x=T)
#write.csv(abf,'D://东北//data//基础信息//TRMP//cwd_projs.csv')

#---------------------------------------------------------
require(forestFloor)
require(robust)
require(randomForest)
setwd("D://东北//data//基础信息//TRMP")

read.csv('cwd_projs.csv')->cwd.gsl

#TO obtain average curve of CWD~AGE
randomForestSRC::rfsrc(CCWD_tC_ha~age+lat+long,data=cwd.gsl[!is.na(cwd.gsl$CCWD_tC_ha),],importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=3000)->gfr.cwd
randomForestSRC::plot.variable(gfr.cwd,"age",partial=T,npts=70)->gr.cwd

fr.cwd<-data.frame(age=gr.cwd$pData[[1]]$x.uniq,cwd=gr.cwd$pData[[1]]$yhat,cwd.se=gr.cwd$pData[[1]]$yhat.se)
ylim<-floor(range(fr.cwd$age))
cwd.a<-data.frame(approx(fr.cwd$age,fr.cwd$cwd,11:300),ci=I(2*approx(fr.cwd$age,fr.cwd$cwd.se,11:300)$y))
names(cwd.a)<-c("age","cwd","ci")

#TO fill NA of CWD
#
train_fl<-cwd.gsl[!is.na(cwd.gsl$CCWD_tC_ha),]
ii_p<-which(is.na(cwd.gsl$CCWD_tC_ha))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","age")

col.y<-c("CCWD_tC_ha")

Y=train_fl[,col.y]
X=train_fl[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
#ff<- forestFloor(rfo,X)
#print(ff)
#Col=fcol(ff,1)

#plot(ff,col=Col,order_by_importance=TRUE,pch=19)

pp<-predict(rfo)
lmm<-lm(Y~pp)

z<-predict(rfo,newdata=cwd.gsl[,col2prd]);
cwd.gsl$CCWD_tC_ha_prd<-coef(lmm)[1]+coef(lmm)[2]*z

#TO obtain second average curve of CWD~AGE
randomForestSRC::rfsrc(CCWD_tC_ha_prd~age+lat+long,data=cwd.gsl[!is.na(cwd.gsl$CCWD_tC_ha_prd),],importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=3000)->gfr.cwd
randomForestSRC::plot.variable(gfr.cwd,"age",partial=T,npts=70)->gr.cwd

fr2.cwd<-data.frame(age=gr.cwd$pData[[1]]$x.uniq,cwd=gr.cwd$pData[[1]]$yhat,cwd.se=gr.cwd$pData[[1]]$yhat.se)
ylim<-floor(range(fr2.cwd$age))
cwd2.a<-data.frame(approx(fr2.cwd$age,fr2.cwd$cwd,11:300),ci=I(2*approx(fr2.cwd$age,fr2.cwd$cwd.se,11:300)$y))
names(cwd2.a)<-c("age","cwd","ci")

write.csv(cwd.gsl,'all_grs_shrb_litb.csv')

