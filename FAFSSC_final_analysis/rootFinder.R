# TODO: To Search the time of maximum root 0.95
# 
# Author: junhui.zhang
###############################################################################

#
require(rootSolve)
require(grofit)
require(reshape)
#models
rSlogic<-function(x,prob=0.99){return(prob*a-a/(1+exp(-k*(x-xc0))))}
rCoC<-function(x,prob=0.99){return(prob*AA-E-A*x^B*exp(c0*x^d))}
rLogic<-function(x,prob=0.99){return((prob-1)*A2-(A1-A2)/(1+(x/x0)^p))}
rRichd<-function(x,prob=0.99){return(A0*prob-A0*(1+(D0-1)*exp(-K0*(x-xc1)))^(1/(1-D0)))}

read.csv('curv_fit1.csv')->crv_p
#Covington_curve parameter order: A,B,c0,d,E
#Slogic_curce parameter order:a,xc0,k
#Logistic_curve parameter order:A1,A2,x0,p
#Richards_curve parameyer order:A0,xc1,D0,K0

root_rst<-as.data.frame(matrix(NA,ncol=5,nrow=nrow(crv_p)))
names(root_rst)<-c("ZONE","MODEL","VARIABLE","ROOT","Yhat")
root_rst[,1:3]<-crv_p[,1:3]
x=20:400
for (i in 1:nrow(crv_p)){
	if(crv_p$model[i] == "Covington_curve"){
		A=crv_p[i,4];B=crv_p[i,5];c0=crv_p[i,6];d=crv_p[i,7];E=crv_p[i,8];
		grofit.y=A*x^B*exp(c0*x^d)+E
		AA<-quantile(grofit.y,1)
#		AA<-coef(gcFitModel(x,grofit.y)$nls)[1]
		rt<-uniroot.all(rCoC,c(20,400));root_rst[i,4:(4+length(rt)-1)]<-rt
		rm(A,B,c0,d,E,AA);
	}
	
	if(crv_p$model[i] == "Logistic"){
		A1=crv_p[i,4];A2=crv_p[i,5];x0=crv_p[i,6];p=crv_p[i,7];
		root_rst[i,4]<-uniroot(rLogic,c(20,400),extendInt="yes")$root;
		rm(A1,A2,x0,p);
	}
	
	if(crv_p$model[i] == "Slogistic"){
		a=crv_p[i,4];xc0=crv_p[i,5];k=crv_p[i,6]
		root_rst[i,4]<-uniroot(rSlogic,c(20,400),extendInt="yes")$root;
		rm(a,xc0,k);
	}
	
	if(crv_p$model[i] == "Srichards"){
		A0=crv_p[i,4];xc1=crv_p[i,5];D0=crv_p[i,6];K0=crv_p[i,7];
		root_rst[i,4]<-uniroot(rRichd,c(20,400),extendInt="yes")$root;
		rm(A0,xc1,D0,K0);
	}
}

roots<-root_rst[,c("ZONE","VARIABLE","ROOT" )]
roots.w<-reshape(roots,direction='wide',idvar=c("ZONE"),timevar="VARIABLE")
write.csv(roots.w,'root_w1a.csv',row.names=F)