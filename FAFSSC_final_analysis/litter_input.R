lfun <- function(x, ...) {
  c(l = length(x))
}


require(Hmisc)
require(doBy)

read.csv("D://����//data//̼//Biomass_Cproj.csv")->bio.cproj
read.csv("D://����//data//̼//Biomass_zzzy.csv")->bio.zzzy


bio.zzzy<-bio.zzzy[,c("ydID","zwmcbh","xj","Wfoliage_kg","Wstem_kg","Wroots_kg","pft")]
bio.cproj<-bio.cproj[,c("ydID","zwmcbh","xj","Wfoliage_kg","Wstem_kg","Wroots_kg","pft")]
bio.tot<-rbind(bio.zzzy,bio.cproj)

bio.tot$xjbrk<-cut(bio.tot$xj,breaks=c(0,10,20,30,40,50,60,70,100,225))

write.csv(bio.tot,'biomass_tot.csv',row.names=FALSE)

plotlst<-levels(bio.tot$ydID);
nplots<-length(plotlst)
# 
# odnf<-bio.tot$pft %in% c("DNF")
# odbf<-bio.tot$pft %in% c("DBF")
# odbf<-bio.tot$pft %in% c("ENF")



l<- summaryBy(xj ~ ydID+xjbrk, data=bio.tot, FUN=lfun, na.omit=TRUE)
mcarbon<- summaryBy(Wfoliage_kg+Wstem_kg+Wroots_kg ~ ydID + xjbrk, data=bio.tot, FUN=mean, na.omit=T)

#pft<-c("DBF","ENF","DNF")

loss.carbon<-data.frame(matrix(NA,ncol=4,nrow=I(3*nplots)))
names(loss.carbon)<-c("plot","Wfoliage_kg","Wstem_kg","Wroots_kg")
for (iplot in 1:nplots){
    loss.carbon$plot[iplot]<-plotlst[iplot]
    ii_l<-which(l$ydID == plotlst[iplot])
    if(length(l$xj.l[ii_l])< 2){loss.carbon[iplot,2:4]<-0}else
    {
      loss<-(-diff(l$xj.l[ii_l]))
      loss[length(loss)+1]<-loss[length(loss)]
      loss[loss<0]<-0
      loss.carbon[iplot,2:4]<-apply((mcarbon[(mcarbon$ydID %in% plotlst[iplot]),3:5])*loss,2,sum)      
    }
}
#-----------------------------

read.csv('Book8.csv')->b8

rfsrc(LC_tot_tC_ha_yr~age+lat+long,data=b8,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lc
plot.variable(gfr.lc,"age",partial=T,npts=50)->gr.lc

lc.rf<-data.frame(age=gr.lc$pData[[1]]$x.uniq,lc=gr.lc$pData[[1]]$yhat,lc.se=gr.lc$pData[[1]]$yhat.se)
ylim<-floor(range(lc.rf$age))
lc.a<-approx(lc.rf$age,lc.rf$lc,ylim[1]:272)



#-------------------------------------------------
read.csv('afb.csv')->afb

rfsrc(dead_froot_tC_ha~age_max+lat+long,data=afb,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.fr
plot.variable(gfr.fr,"age_max",partial=T,npts=50)->gr.fr

fr.rf<-data.frame(age=gr.fr$pData[[1]]$x.uniq,lc=gr.fr$pData[[1]]$yhat,lc.se=gr.fr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.rf$age))
fr.a<-approx(fr.rf$age,fr.rf$lc,ylim[1]:272)
#-------------------------------------------------------------
require(randomForestSRC)
require(zoo)
setwd("D://东北//data//基础信息//TRMP")

read.csv('bio_soc_env_all_gsl.csv')->all.gsl

rfsrc(grs_tC_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.gr
plot.variable(gfr.gr,"age_max",partial=T,npts=50)->gr.gr

fr.gr<-data.frame(age=gr.gr$pData[[1]]$x.uniq,grs=gr.gr$pData[[1]]$yhat,grs.se=gr.gr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.gr$age))
grs.a<-approx(fr.gr$age,fr.gr$grs,ylim[1]:272)

#------------------------------------------------------
require(randomForestSRC)
require(zoo)
setwd("D://东北//data//基础信息//TRMP")

read.csv('bio_soc_env_all_gsl.csv')->all.gsl

rfsrc(shrb_t_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.shrb
plot.variable(gfr.shrb,"age_max",partial=T,npts=50)->gr.shrb

fr.shrb<-data.frame(age=gr.shrb$pData[[1]]$x.uniq,shrb=gr.shrb$pData[[1]]$yhat,grs.se=gr.shrb$pData[[1]]$yhat.se)
ylim<-floor(range(fr.shrb$age))
shrb.a<-approx(fr.shrb$age,fr.shrb$shrb,ylim[1]:272)

#------------------------------------------------------------FILL NA　OF GRASS & SHRUB
require(randomForest)
require(forestFloor)
require(robust)

setwd("D://东北//data//基础信息//TRMP")

read.csv('bio_soc_env_all_gsl.csv')->all.gsl

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
