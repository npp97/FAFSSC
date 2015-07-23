#litters TO merge litters input from 3 files into one list and generate the RF-age 
#relating files :  biobio_soc_env_all_gsla.csv:: leaf_litter; grass_litter,shrub_litter
#				   Book8.csv :: LC_tot_tC_ha_yr
#				   afb.csv:: fine_roots


setwd("D://东北//data//基础信息//TRMP")

read.csv('bio_soc_env_all_gsla.csv')->all.gsl
#----------------------------GRASS
rfsrc(grs_tC_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.gr
plot.variable(gfr.gr,"age_max",partial=T,npts=50)->gr.gr

fr.gr<-data.frame(age=gr.gr$pData[[1]]$x.uniq,grs=gr.gr$pData[[1]]$yhat,grs.se=gr.gr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.gr$age))
grs.a<-data.frame(approx(fr.gr$age,fr.gr$grs,ylim[1]:272),ci=I(2*7*approx(fr.gr$age,fr.gr$grs.se,ylim[1]:272)$y))
names(grs.a)<-c("age","grs","ci")

#------------------------------------------------------SHRUB
rfsrc(shrb_t_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.shrb
plot.variable(gfr.shrb,"age_max",partial=T,npts=50)->gr.shrb

fr.shrb<-data.frame(age=gr.shrb$pData[[1]]$x.uniq,shrb=gr.shrb$pData[[1]]$yhat,shrb.se=gr.shrb$pData[[1]]$yhat.se)
ylim<-floor(range(fr.shrb$age))
shrb.a<-data.frame(approx(fr.shrb$age,fr.shrb$shrb,ylim[1]:272),ci=I(2*7*approx(fr.shrb$age,fr.shrb$shrb.se,ylim[1]:272)$y))
names(shrb.a)<-c("age","shrb","ci")

#------------------leaf_litter
read.csv('bio_soc_env_all_gsla.csv')->all.gsl
#----------------------------Leaf
rfsrc(W_leaf_kg_tC_ha~age_max+lat+long,data=all.gsl,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lf
plot.variable(gfr.lf,"age_max",partial=T,npts=50)->gr.lf

fr.lf<-data.frame(age=gr.lf$pData[[1]]$x.uniq,lf=gr.lf$pData[[1]]$yhat,lf.se=gr.lf$pData[[1]]$yhat.se)
ylim<-floor(range(fr.lf$age))
lf.a<-data.frame(approx(fr.lf$age,fr.lf$lf,ylim[1]:272),ci=I(2*7*approx(fr.lf$age,fr.lf$lf.se,ylim[1]:272)$y))
names(lf.a)<-c("age","lf","ci")
lf.a$lf<-lf.a$lf+shrb.a$shrb*0.15

#----------------------------------
read.csv('Book8.csv')->b8 #

rfsrc(LC_tot_tC_ha_yr~age+lat+long,data=b8,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lc
plot.variable(gfr.lc,"age",partial=T,npts=50)->gr.lc

lc.rf<-data.frame(age=gr.lc$pData[[1]]$x.uniq,lc=gr.lc$pData[[1]]$yhat,lc.se=gr.lc$pData[[1]]$yhat.se)
ylim<-floor(range(lc.rf$age))
lc.a<-approx(lc.rf$age,lc.rf$lc,ylim[1]:272)
lc.f<-data.frame(age=lc.a$x,lc=lc.a$y,ci=I(2*5*approx(lc.rf$age,lc.rf$lc.se,ylim[1]:272)$y))
lc.f$lc<-lc.f$lc+lf.a$lf/5
lc.f$lc.se<-lc.f$ci+lf.a$ci/5

########################FINE Roots
read.csv('afb.csv')->afb

rfsrc(dead_froot_tC_ha~age_max+lat+long,data=afb,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.fr
plot.variable(gfr.fr,"age_max",partial=T,npts=50)->gr.fr

fr.rf<-data.frame(age=gr.fr$pData[[1]]$x.uniq,fr=gr.fr$pData[[1]]$yhat,fr.se=gr.fr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.rf$age))
fr.a<-approx(fr.rf$age,fr.rf$fr,ylim[1]:272)
fr.f<-data.frame(age=fr.a$x,fr=fr.a$y,ci=I(2*5*approx(fr.rf$age,fr.rf$fr,ylim[1]:272)$y))

save.image("litters.RData")
