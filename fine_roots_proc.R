require(doBy)

setwd("D://东北//data//碳专项")
fr.raw<-read.csv("碳专项细根现存量调查20150723.csv")
fr.1<-summaryBy(live_fr_g_m2+dead_fr_g_m2+total_fr_g_m2~pplot+soil_layer,data=fr.raw,FUN=mean,na.rm=T)
fr.2<-summaryBy(live_fr_g_m2.mean+dead_fr_g_m2.mean+total_fr_g_m2.mean~pplot,data=fr.1,FUN=sum,na.rm=T)
write.csv(fr.2,"碳专项细根现存量调查By_plot_20150723.csv",row.names=F)
################################
require(forestFloor)
require(robust)
require(randomForest)
setwd("D://东北//data//基础信息//TRMP")

read.csv('碳专项细根现存量调查By_plot_20150723.csv')->fr.gsl
read.csv('cwd_projs.csv')->cwd.gsl
names(cwd.gsl)[15]<-"age_max"
read.csv('bio_soc_env_all_gsla.csv')->all.gsl

merge(cwd.gsl,fr.gsl)->fr.gsl
#################live fine-roots
train_fl<-fr.gsl[!is.na(fr.gsl$live_fr_tC_ha),]

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","age")

col.y<-c("live_fr_tC_ha")

Y=train_fl[,col.y]
X=train_fl[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
#ff<- forestFloor(rfo,X)
#print(ff)
#Col=fcol(ff,1)

#plot(ff,col=Col,order_by_importance=TRUE,pch=19)
# 
pp<-predict(rfo)
lmm<-lm(Y~pp)

z<-predict(rfo,newdata=all.gsl[,col2prd]);
all.gsl$live_fr_tC_ha<-coef(lmm)[1]+coef(lmm)[2]*z

##############dead fine-roots
train_fl<-fr.gsl[!is.na(fr.gsl$dead_fr_tC_ha),]

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","age")

col.y<-c("dead_fr_tC_ha")

Y=train_fl[,col.y]
X=train_fl[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
#ff<- forestFloor(rfo,X)
#print(ff)
#Col=fcol(ff,1)

#plot(ff,col=Col,order_by_importance=TRUE,pch=19)
# 
pp<-predict(rfo)
lmm<-lm(Y~pp)

z<-predict(rfo,newdata=all.gsl[,col2prd]);
all.gsl$dead_fr_tC_ha<-coef(lmm)[1]+coef(lmm)[2]*z

################total_fr_g_m2
train_fl<-fr.gsl[!is.na(fr.gsl$total_fr_tC_ha),]

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","age")

col.y<-c("total_fr_tC_ha")

Y=train_fl[,col.y]
X=train_fl[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
#ff<- forestFloor(rfo,X)
#print(ff)
#Col=fcol(ff,1)

#plot(ff,col=Col,order_by_importance=TRUE,pch=19)
# 
pp<-predict(rfo)
lmm<-lm(Y~pp)

z<-predict(rfo,newdata=all.gsl[,col2prd]);
all.gsl$total_fr_tC_ha<-coef(lmm)[1]+coef(lmm)[2]*z

write.csv(all.gsl,'bio_soc_env_all_gsla.csv')
