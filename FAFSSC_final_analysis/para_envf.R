require(nlme)
require(doBy)

sumfun <- function(x, ...) {
	c(tm=mean(x,0.025, ...), m = mean(x, ...), md = median(x, ...), v = sd(x, ...), l = length(x))
}

maxi<-function(x){
	x1<-as.factor(x)
	a<-summaryBy(x~x1,data=data.frame(x,x1),FUN=length)
	return(as.character(a[which.max(a[,2]),1]))
}
zone<-c("C2B","C3B","C3C","B2B","B4C","A2A","A2B","A2D","C4C","A3C","C1B","B4D","B2D","A2C")

env.raw<-read.csv('bio_soc_env_all.csv.csv')

ii<-which(env.raw$cls_max_spec4 %in% zone)

env2ana<-env.raw[ii,]

pt.s<-summaryBy(p_mm+ta_C~cls_max_spec4,data=env2ana,FUN=sumfun)
names(pt.s)[1]<-'ZONE'
sc.s<-summaryBy(soilcode~cls_max_spec4,data=env2ana,FUN=maxi)

read.csv('root_w.csv')->root.w1
root.w<-root.w1[,2:4]
row.names(root.w)<-root.w1[,1]
root.long<-reshape(root.w,idvar="ZONE",ids=row.names(root.w),times=names(root.w),timevar="variable",v.names="ROOT",varying=list(names(root.w)),direction="long")
root.long<-root.long[,c("ZONE","variable","ROOT")]

read.csv('curv_fit.csv')->crv_p
crv_root<-merge(crv_p,root.long,by=c("ZONE","variable"))
#names(crv_root): ZONE,variable,model,a,b,c,d,e,ROOT

ss<-as.data.frame(matrix(NA,ncol=2,nrow=nrow(crv_root)))
names(ss)<-c("bottom","top")
for (i in 1:nrow(crv_p)){
	x=c(15:(floor(crv_root$ROOT[i])+10))
	if(crv_root$model[i] == "Covington_curve"){
		A=crv_root[i,4];B=crv_root[i,5];c0=crv_root[i,6];d=crv_root[i,7];E=crv_root[i,8];
		gf.y=A*x^B*exp(c0*x^d)+E
		ss[i,]<-range(gf.y)
		rm(A,B,c0,d,E,gf.y);
	}
	
	if(crv_root$model[i] == "Logistic"){
		A1=crv_root[i,4];A2=crv_root[i,5];x0=crv_root[i,6];p=crv_root[i,7];
		gf.y<-A2+(A1-A2)/(1+(x/x0)^p)
		ss[i,]<-range(gf.y);
		rm(A1,A2,x0,p,gf.y);
	}
	
	if(crv_root$model[i] == "Slogistic"){
		a=crv_root[i,4];xc0=crv_root[i,5];k=crv_root[i,6]
		gf.y<-a/(1+exp(-k*(x-xc0)))
		ss[i,]<-range(gf.y);
		rm(a,xc0,k,gf.y);
	}
	
	if(crv_root$model[i] == "Srichards"){
		A0=crv_root[i,4];xc1=crv_root[i,5];D0=crv_root[i,6];K0=crv_root[i,7];
		gf.y<-A0*(1+(D0-1)*exp(-K0*(x-xc1)))^(1/(1-D0))
		ss[i,]<-range(gf.y);
		rm(A0,xc1,D0,K0,gf.y);
	}
}

crv_ss<-data.frame(crv_root,ss)

crv_pt<-merge(pt.s,crv_ss)

write.csv(crv_pt,'crv_pt.csv',row.names=F)
#--------------------------------------------------------------

library(MuMIn)

read.csv('crv_pt.csv')->crv_pt
names(crv_pt)
#[1] "ZONE"     "variable" "SPEC"     "p_mm.tm"  "p_mm.md"  "p_mm.v"   "p_mm.l"   "ta_C.tm"  "ta_C.m"   "ta_C.md"  "ta_C.v"   "ta_C.l"   "model"    "a"        "b"        "c"        "d"       
#[18] "e"        "ROOT"     "bottom"   "top"  

crv_pt.ana<-crv_pt[,c("ZONE","variable","LatZ","LonZ","SPEC","p_mm.md","ta_C.md","ROOT","top","dif")]
model.rbio.null<-lme(ROOT~1,random=~1|ZONE,data=crv_pt.ana)
aa<-gls(ROOT~p_mm.md*ta_C.md+SPEC,data=crv_pt.ana[crv_pt.ana$variable=='BIO',])


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

