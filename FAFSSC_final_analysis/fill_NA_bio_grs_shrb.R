#------------------------------------------------------------FILL NA　OF GRASS & SHRUB
require(randomForest)
require(forestFloor)
require(robust)

setwd("D://东北//data//基础信息//TRMP")

read.csv('bio_soc_env_all_gsl.csv')->all.gsl

#
#-----------------------BIOMASS
train_bio<-all.gsl[!is.na(all.gsl$W_total_t_ha),]
ii_p<-which(is.na(all.gsl$W_total_t_ha))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","slope","aspect","twi","tpi","W_total_t_ha","lai","age")

col.y<-c("W_total_t_ha")

Y=train_bio[,col.y]
X=train_bio[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
pp<-predict(rfo)
lmm<-rlm(Y~pp)

z<-predict(rfo,newdata=all.gsl[ii_p,col2prd]);
all.gsl$W_total_t_ha[ii_p]<-coef(lmm)[1]+coef(lmm)[2]*z

#----C------------LAI
train_lai<-all.gsl[!is.na(all.gsl$lai),]
ii_p<-which(is.na(all.gsl$lai))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","slope","aspect","twi","tpi","W_total_t_ha","lai","age")

col.y<-c("lai")

Y=train_lai[,col.y]
X=train_lai[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
pp<-predict(rfo)
lmm<-rlm(Y~pp)

z<-predict(rfo,newdata=all.gsl[ii_p,col2prd]);
all.gsl$lai[ii_p]<-coef(lmm)[1]+coef(lmm)[2]*z
#-----------------------SOIL 

train_soilc<-all.gsl[!is.na(all.gsl$SoilC_1m_t_ha),]
ii_p<-which(is.na(all.gsl$SoilC_1m_t_ha))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","slope","aspect","twi","tpi","W_total_t_ha","lai","age")

col.y<-c("SoilC_1m_t_ha")

Y=train_soilc[,col.y]
X=train_soilc[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
pp<-predict(rfo)
lmm<-rlm(Y~pp)

z<-predict(rfo,newdata=all.gsl[ii_p,col2prd]);
all.gsl$SoilC_1m_t_ha[ii_p]<-coef(lmm)[1]+coef(lmm)[2]*z
#-----------------------GRASS C

train_grs<-all.gsl[!is.na(all.gsl$gras_C),]
ii_p<-which(is.na(all.gsl$gras_C))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","slope","aspect","twi","tpi","W_total_t_ha","lai","age")

col.y<-c("gras_C")

Y=train_grs[,col.y]
X=train_grs[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
pp<-predict(rfo)
lmm<-rlm(Y~pp)

z<-predict(rfo,newdata=all.gsl[ii_p,col2prd]);
all.gsl$gras_C[ii_p]<-coef(lmm)[1]+coef(lmm)[2]*z

#
train_shrb<-all.gsl[!is.na(all.gsl$shrb_C),]
ii_p<-which(is.na(all.gsl$shrb_C))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","slope","aspect","twi","tpi","W_total_t_ha","lai","age")

col.y<-c("shrb_C")

Y=train_shrb[,col.y]
X=train_shrb[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
pp<-predict(rfo)
lmm<-rlm(Y~pp)

z<-predict(rfo,newdata=all.gsl[ii_p,col2prd]);
all.gsl$shrb_C[ii_p]<-coef(lmm)[1]+coef(lmm)[2]*z

write.csv(all.gsl,'bio_soc_env_all_gsla.csv')
#


setwd("D://东北//data//基础信息//TRMP")

read.csv('bio_soc_env_all_gsla.csv')->all.gsl
train_fl<-all.gsl[!is.na(all.gsl$litt_Carbon_t_ha),]
ii_p<-which(is.na(all.gsl$litt_Carbon_t_ha))

col2prd<-c("lat","long","elev","slope","aspect","twi","tpi","p_mm","ta_C","W_total_t_ha","lai","age_max")
#col2prd<-c("ta_C","p_mm","slope","aspect","twi","tpi","W_total_t_ha","lai","age")
#col2prd<-c("lat","long","slope","aspect","twi","tpi","W_total_t_ha","lai","age")

col.y<-c("shrb_C")

Y=train_fl[,col.y]
X=train_fl[,col2prd]

rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE,ntree=5000,replace=TRUE,importance=TRUE)
pp<-predict(rfo)
lmm<-rlm(Y~pp)

z<-predict(rfo,newdata=all.gsl[ii_p,col2prd]);
all.gsl$litt_Carbon_t_ha[ii_p]<-coef(lmm)[1]+coef(lmm)[2]*z

write.csv(all.gsl,'bio_soc_env_all_gsla.csv')
