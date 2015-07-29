
library(MuMIn)
require(randomForestSRC)
require(nlme)

setwd("D://东北//data//基础信息//RESULT//csv")
read.csv("roots_max_min_w.csv")->mmp_l

ii_tot<-which(mmp_l$VARIABLE %in% "TOTC")
rfsrc(root~ta_C+p_mm+veg+bottom,data=mmp_l[ii_tot,],importance="permute.ensemble",keep.inbag=TRUE,replace=TRUE,na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.root.TOTC
plot.variable(gfr.root.TOTC,partial=T)
randomForest(root~ta_C+p_mm+veg+bottom,data=mmp_l[ii_tot,],keep.inbag=TRUE,replace=TRUE,ntree=5000,importance=TRUE)->rf.rtotc


ii_bio<-which(mmp_l$VARIABLE %in% "BIO")
rfsrc(root~ta_C+p_mm+veg+bottom,data=mmp_l[ii_bio,],importance="permute.ensemble",keep.inbag=TRUE,replace=TRUE,na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.root.bio
plot.variable(gfr.root.bio,partial=T)


ii_soilc<-which(mmp_l$VARIABLE %in% "SOILC")
rfsrc(root~ta_C+p_mm+veg+bottom,data=mmp_l[ii_soilc,],importance="permute.ensemble",keep.inbag=TRUE,replace=TRUE,na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.root.soilc
plot.variable(gfr.root.soilc,partial=T)
#------

win.graph()
ii_tot<-which(mmp_l$VARIABLE %in% "TOTC")
rfsrc(top~ta_C+p_mm+veg+bottom,data=mmp_l[ii_tot,],importance="permute.ensemble",keep.inbag=TRUE,replace=TRUE,na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.top.TOTC
plot.variable(gfr.top.TOTC,partial=T)
randomForest(top~ta_C+p_mm+veg+bottom,data=mmp_l[ii_tot,],keep.inbag=TRUE,replace=TRUE,ntree=5000,importance=TRUE)->rf.ttotc

win.graph()
ii_bio<-which(mmp_l$VARIABLE %in% "BIO")
rfsrc(top~ta_C+p_mm+veg+bottom,data=mmp_l[ii_bio,],importance="permute.ensemble",keep.inbag=TRUE,replace=TRUE,na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.top.bio
plot.variable(gfr.top.bio,partial=T)

win.graph()
ii_soilc<-which(mmp_l$VARIABLE %in% "SOILC")
rfsrc(top~ta_C+p_mm+veg+bottom,data=mmp_l[ii_soilc,],importance="permute.ensemble",keep.inbag=TRUE,replace=TRUE,na.action = "na.impute",nimpute = 10, seed= -111,ntree=5000)->gfr.top.soilc
plot.variable(gfr.top.soilc,partial=T)


model.rbio.null<-lme(root~1,random=~1|ZONE,data=mmp_l[mmp_l$VARIABLE=='BIO',])
rbio.model<-gls(log(root)~p_mm+ta_C+veg,data=mmp_l[mmp_l$VARIABLE=='BIO',])
rsoil.model<-gls(log(root)~p_mm+ta_C+veg+bottom,data=mmp_l[mmp_l$VARIABLE=='SOILC',],na.action="na.omit")
rtotc.model<-gls(log(root)~p_mm+ta_C+veg,data=mmp_l[mmp_l$VARIABLE=='TOTC',])
anova(rbio.model)
anova(rsoil.model)
anova(rtotc.model)


tbio2.model<-gls(top~p_mm+ta_C+veg+root,data=mmp_l[mmp_l$VARIABLE=='BIO',])
tsoil2.model<-gls(top~p_mm+ta_C+veg+root,data=mmp_l[mmp_l$VARIABLE=='SOILC',],na.action="na.omit")
ttotc2.model<-gls(top~p_mm+ta_C+veg+root,data=mmp_l[mmp_l$VARIABLE=='TOTC',])
anova(tbio2.model)
anova(tsoil2.model)
anova(ttotc2.model)



All_model.rbio<-lme(ROOT~p_mm.md*ta_C.md+SPEC,random=~LatZ|LonZ/SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
top_model.rbio<-lme(ROOT~p_mm.md+ta_C.md+SPEC,random=~LatZ|LonZ/SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])

All_model.rsoilc<-lme(ROOT~p_mm.md*ta_C.md+SPEC,random=~LatZ|LonZ/SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])
top_model.rsoilc<-lme(ROOT~p_mm.md+ta_C.md+SPEC,random=~LatZ|LonZ/SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])

All_model.rtotc<-lme(ROOT~p_mm.md*ta_C.md+SPEC,random=~LatZ|LatZ/SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])
top_model.rtotc<-lme(ROOT~p_mm.md+ta_C.md+SPEC,random=~LatZ|LatZ/SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])

lm_model.rbio<-gls(ROOT~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
anova(lm_model.rbio)
lm_model.rsoilc<-gls(ROOT~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])
anova(lm_model.rsoilc)
lm_model.rtotc<-gls(ROOT~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])
anova(lm_model.rtotc)
rf.rtotc<-rfsrc(ROOT~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',],importance="permute.ensemble",ntree=5000, seed=-111)
rf.rtot.p<-plot.variable(rf.rtotc,partial=TRUE,ylable='XX')
win.graph()
rf.rsoilc<-rfsrc(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',],importance="permute.ensemble",ntree=5000, seed=-111)
rf.rsoilc.p<-plot.variable(rf.rsoilc,partial=TRUE)
win.graph()
rf.rbio<-rfsrc(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',],importance="permute.ensemble",ntree=5000, seed=-111)
rf.rbio.p<-plot.variable(rf.rbio,partial=TRUE)




AICc_res<-AICc(All_model.rbio,top_model.rbio,All_model.rsoilc,top_model.rsoilc,All_model.rtotc,top_model.rtotc)
AICc_table<-data.frame(Model=row.names(AICc_res),AICc=AICc_res$AICc)
AICc_table$delta<-AICc_table$AICc-AICc_table$AICc[1]
AICc_table$R_squared<-c(r.squaredGLMM(All_model.rbio)[1],r.squaredGLMM(top_model.rbio)[1],r.squaredGLMM(All_model.rsoilc)[1],r.squaredGLMM(top_model.rsoilc)[1],r.squaredGLMM(All_model.rtotc)[1],r.squaredGLMM(top_model.rtotc)[1])
AICc_table

#-----------------------------------------
All_model.tbio<-lme(top~p_mm.md*ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
top_model.tbio<-lme(top~p_mm.md+ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
lm_model<-gls(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
lm_model.2<-gls(top~p_mm.md*ta_C.md,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
lm_model.3<-gls(top~p_mm.md+ta_C.md,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
lm_model.4<-gls(top~p_mm.md+ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])

All_model.tsoilc<-lme(top~p_mm.md*ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])
top_model.tsoilc<-lme(top~p_mm.md+ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])

All_model.ttotc<-lme(top~p_mm.md*ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])
top_model.ttotc<-lme(top~p_mm.md+ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])

AICc_res<-AICc(All_model.tbio,top_model.tbio,All_model.tsoilc,top_model.tsoilc,All_model.ttotc,top_model.ttotc)
AICc_table<-data.frame(Model=row.names(AICc_res),AICc=AICc_res$AICc)
AICc_table$delta<-AICc_table$AICc-AICc_table$AICc[1]
AICc_table$R_squared<-c(r.squaredGLMM(All_model.tbio)[1],r.squaredGLMM(top_model.tbio)[1],r.squaredGLMM(All_model.tsoilc)[1],r.squaredGLMM(top_model.tsoilc)[1],r.squaredGLMM(All_model.ttotc)[1],r.squaredGLMM(top_model.ttotc)[1])
AICc_table


lm_model.tbio<-gls(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
anova(lm_model.tbio)
lm_model.tsoilc<-gls(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])
anova(lm_model.tsoilc)
lm_model.ttotc<-gls(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])
anova(lm_model.ttotc)
rf.ttotc<-rfsrc(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',],importance="permute.ensemble",ntree=10000, seed=-111)
rf.ttot.p<-plot.variable(rf.ttotc,partial=TRUE)
win.graph()
rf.tsoilc<-rfsrc(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',],importance="permute.ensemble",ntree=10000, seed=-111)
rf.tsoilc.p<-plot.variable(rf.tsoilc,partial=TRUE)
win.graph()
rf.tbio<-rfsrc(top~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',],importance="permute.ensemble",ntree=10000, seed=-111)
rf.tbio.p<-plot.variable(rf.tbio,partial=TRUE)



#------------------------------------------------------------------------------------------------------------
vf.bio<-varIdent(form=~1|p_mm.md)
All_model.dbio<-lme(dif~p_mm.md*ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='BIO',],weights=vf.bio)

top_model.dbio<-lme(dif~p_mm.md+ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])

All_model.dsoilc<-lme(dif~p_mm.md*ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])
top_model.dsoilc<-lme(dif~p_mm.md+ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])

All_model.dtotc<-lme(dif~p_mm.md*ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])
top_model.dtotc<-lme(dif~p_mm.md+ta_C.md+SPEC,random=~1|ZONE,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])

AICc_res<-AICc(All_model.dbio,top_model.dbio,All_model.dsoilc,top_model.dsoilc,All_model.dtotc,top_model.dtotc)
AICc_table<-data.frame(Model=row.names(AICc_res),AICc=AICc_res$AICc)
AICc_table$delta<-AICc_table$AICc-AICc_table$AICc[1]
AICc_table$R_squared<-c(r.squaredGLMM(All_model.dbio)[1],r.squaredGLMM(top_model.dbio)[1],r.squaredGLMM(All_model.dsoilc)[1],r.squaredGLMM(top_model.dsoilc)[1],r.squaredGLMM(All_model.dtotc)[1],r.squaredGLMM(top_model.dtotc)[1])
AICc_table


lm_model.dbio<-gls(dif~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])
anova(lm_model.dbio)
lm_model.dsoilc<-gls(dif~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='SOILC',])
anova(lm_model.dsoilc)
lm_model.dtotc<-gam(dif~s(p_mm.md)*s(ta_C.md)+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',])
anova(lm_model.dtotc)

rf.dtotc<-rfsrc(dif~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='TOTC',],importance="permute.ensemble",ntree=10000, seed=-111)
rf.dtot.p<-plot.variable(rf.dtotc,partial=TRUE)

#Precipation 