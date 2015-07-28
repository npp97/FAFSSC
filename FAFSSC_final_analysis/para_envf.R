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


