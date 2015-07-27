#litters TO merge litters input from 3 files into one list and generate the RF-age 
#relating files :  biobio_soc_env_all_gsla.csv:: leaf_litter; grass_litter,shrub_litter
#				   Book8.csv :: LC_tot_tC_ha_yr
#				   afb.csv:: fine_roots

require(randomForestSRC)
require(zoo)
require(Hmisc)

setwd("D://东北//data//基础信息//RESULT//csv")

read.csv('bio_soc_env_all_gsla.csv')->all.gsl
#----------------------------GRASS
rfsrc(grs_tC_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.gr
plot.variable(gfr.gr,"age_max",partial=T,npts=75)->gr.gr

fr.gr<-data.frame(age=gr.gr$pData[[1]]$x.uniq,grs=gr.gr$pData[[1]]$yhat,grs.se=gr.gr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.gr$age))
grs.a<-data.frame(approx(fr.gr$age,fr.gr$grs,ylim[1]:272),ci=I(2*approx(fr.gr$age,fr.gr$grs.se,ylim[1]:272)$y))
names(grs.a)<-c("age","grs","grs.ci")

#------------------------------------------------------SHRUB
rfsrc(shrb_t_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.shrb
plot.variable(gfr.shrb,"age_max",partial=T,npts=75)->gr.shrb

fr.shrb<-data.frame(age=gr.shrb$pData[[1]]$x.uniq,shrb=gr.shrb$pData[[1]]$yhat,shrb.se=gr.shrb$pData[[1]]$yhat.se)
ylim<-floor(range(fr.shrb$age))
shrb.a<-data.frame(approx(fr.shrb$age,fr.shrb$shrb,ylim[1]:272),ci=I(2*approx(fr.shrb$age,fr.shrb$shrb.se,ylim[1]:272)$y))
names(shrb.a)<-c("age","shrb","shrb.ci")

#------------------leaf_litter
#read.csv('bio_soc_env_all_gsla.csv')->all.gsl
#----------------------------Leaf
rfsrc(W_leaf_kg_tC_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lf
plot.variable(gfr.lf,"age_max",partial=T,npts=75)->gr.lf

fr.lf<-data.frame(age=gr.lf$pData[[1]]$x.uniq,lf=gr.lf$pData[[1]]$yhat,lf.se=gr.lf$pData[[1]]$yhat.se)
ylim<-floor(range(fr.lf$age))
lf.a<-data.frame(approx(fr.lf$age,fr.lf$lf,ylim[1]:272),ci=I(2*approx(fr.lf$age,fr.lf$lf.se,ylim[1]:272)$y))
names(lf.a)<-c("age","lf","lf.ci")
lf.a$lf<-lf.a$lf+shrb.a$shrb*0.15  #assume the ratio of shrub leaf litter is 15%
lf.a$lf.ci<-lf.a$lf.ci+shrb.a$shrb.ci*0.15

#----------------------------------
#read.csv('Book8.csv')->b8 #

rfsrc(LC_tot_tC_ha_yr~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lc
plot.variable(gfr.lc,"age_max",partial=T,npts=75)->gr.lc

lc.rf<-data.frame(age=gr.lc$pData[[1]]$x.uniq,lc=gr.lc$pData[[1]]$yhat,lc.se=gr.lc$pData[[1]]$yhat.se)
ylim<-floor(range(lc.rf$age))
lc.a<-approx(lc.rf$age,lc.rf$lc,ylim[1]:272)
lc.f<-data.frame(age=lc.a$x,lc=lc.a$y,ci=I(2*approx(lc.rf$age,lc.rf$lc.se,ylim[1]:272)$y))
lc.f$lc<-lc.f$lc+lf.a$lf/4
lc.f$ci<-lc.f$ci+lf.a$lf.ci/4   #assume the ratio of foilage:leaf is 1:5
names(lc.f)<-c("age","lc","lc.ci")

########################FINE Roots
#read.csv('bio_soc_env_all_gsla.csv')->all.gsl

rfsrc(live_fr_tC_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.fr
plot.variable(gfr.fr,"age_max",partial=T,npts=75)->gr.fr

fr.rf<-data.frame(age=gr.fr$pData[[1]]$x.uniq,fr=gr.fr$pData[[1]]$yhat,fr.se=gr.fr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.rf$age))
fr.a<-approx(fr.rf$age,fr.rf$fr,ylim[1]:272)
fr.f<-data.frame(age=fr.a$x,fr=fr.a$y,fr.ci=I(2*approx(fr.rf$age,fr.rf$fr.se,ylim[1]:272)$y))

#-------litter
#
#read.csv('bio_soc_env_all_gsla.csv')->all.gsl
rfsrc(litt_Carbon_t_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.sl
plot.variable(gfr.sl,"age_max",partial=T,npts=75)->gr.sl

fr.sl<-data.frame(age=gr.sl$pData[[1]]$x.uniq,sl=gr.gr$pData[[1]]$yhat,sl.se=gr.sl$pData[[1]]$yhat.se)
ylim<-floor(range(fr.sl$age))
sl.a<-data.frame(approx(fr.sl$age,fr.sl$grs,ylim[1]:272),ci=I(2*approx(fr.sl$age,fr.sl$sl.se,ylim[1]:272)$y))
names(sl.a)<-c("age","sl","sl.ci")


litt_inputs<-cbind(fr.gr,fr.shrb,fr.lf,lc.rf,fr.rf,fr.sl)
litt_inputs_i<-cbind(grs.a,shrb.a,lf.a,lc.f,fr.f,sl.a)


#ll<-sl.a
#ll$sl<-sl.a$sl+fr.f$fr
#ll$ci<-sl.a$ci+fr.f$ci

write.csv(litt_inputs,"litt_inputs.csv",row.names=F)
write.csv(litt_inputs_i,"litt_inputs_i.csv",row.names=F)

save.image("litters.RData")
