
#by 4p Logistic growth curves
require(lhs)

q75<-function(x,...){return(c(quantile(x,c(0.75),...)))}
q25<-function(x,...){quantile(x,0.25, ...)}
q50<-function(x, ...){median(x, ...)}

parms_gen_lhs<-function (PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 
{
	NUMPARAMS <- length(PARAMETERS)
	ALGORITHM <- tolower(ALGORITHM)
	if (ALGORITHM == "optimum") {
		design <- lhs::optimumLHS(NUMSAMPLES, NUMPARAMS, 
				maxSweeps = 2, eps = 0.1)
	}
	else {
		design <- lhs::randomLHS(NUMSAMPLES, NUMPARAMS)
	}
	for (k in 1:NUMSAMPLES) {
		for (l in 1:NUMPARAMS) {
			lhc_max <- PMAX[l]
			lhc_min <- PMIN[l]
			value <- (design[k, l] * (lhc_max - lhc_min)) + lhc_min
			design[k, l] <- value
		}
	}
	colnames(design) <- c(PARAMETERS)
	return(design)
}

#----------------------------------------------
parms.4log.tot<-c(A1=182.21731,A1.se=3.33417,A2=358.52335,A2.se=2.62196,xc=77.95213,xc.se=1.74168,p=2.8228,p.se=0.18062)
parms.4log.bio<-c(A1=31.10187,A1.se=1.61341,A2=116.77196,A2.se=2.14101,xc=98.87556,xc.se=2.61857,p=2.72947,p.se=0.2088)
parms.4log.soil<-c(A1=153.26201,A1.se=3.13962,A2=243.36767,A2.se=1.6734,xc=66.75651,xc.se=2.50241,p=3.35015,p.se=0.37506)

pr.tot.l<-c(a=(parms.4log.tot["A1"]-2*parms.4log.tot["A1.se"]),xc=(parms.4log.tot["A2"]-2*parms.4log.tot["A2.se"]),d=(parms.4log.tot["xc"]-2*parms.4log.tot["xc.se"]),k=(parms.4log.tot["p"]-2*parms.4log.tot["p.se"]))
pr.tot.u<-c(a=(parms.4log.tot["A1"]+2*parms.4log.tot["A1.se"]),xc=(parms.4log.tot["A2"]+2*parms.4log.tot["A2.se"]),d=(parms.4log.tot["xc"]+2*parms.4log.tot["xc.se"]),k=(parms.4log.tot["p"]+2*parms.4log.tot["p.se"]))
pr.bio.l<-c(a=(parms.4log.bio["A1"]-2*parms.4log.bio["A1.se"]),xc=(parms.4log.bio["A2"]-2*parms.4log.bio["A2.se"]),d=(parms.4log.bio["xc"]-2*parms.4log.bio["xc.se"]),k=(parms.4log.bio["p"]-2*parms.4log.bio["p.se"]))
pr.bio.u<-c(a=(parms.4log.bio["A1"]+2*parms.4log.bio["A1.se"]),xc=(parms.4log.bio["A2"]+2*parms.4log.bio["A2.se"]),d=(parms.4log.bio["xc"]+2*parms.4log.bio["xc.se"]),k=(parms.4log.bio["p"]+2*parms.4log.bio["p.se"]))
pr.soil.l<-c(a=(parms.4log.soil["A1"]-2*parms.4log.soil["A1.se"]),xc=(parms.4log.soil["A2"]-2*parms.4log.soil["A2.se"]),d=(parms.4log.soil["xc"]-2*parms.4log.soil["xc.se"]),k=(parms.4log.soil["p"]-2*parms.4log.soil["p.se"]))
pr.soil.u<-c(a=(parms.4log.soil["A1"]+2*parms.4log.soil["A1.se"]),xc=(parms.4log.soil["A2"]+2*parms.4log.soil["A2.se"]),d=(parms.4log.soil["xc"]+2*parms.4log.soil["xc.se"]),k=(parms.4log.soil["p"]+2*parms.4log.soil["p.se"]))

PARAMETERS = c("A1","A2","xc","p")
NUMSAMPLES=10000
PMIN<-pr.soil.l
PMAX<-pr.soil.u
ALGORITHM <-'normal'
par.soil<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

PARAMETERS = c("A1","A2","xc","p")
NUMSAMPLES=10000
PMIN<-pr.tot.l
PMAX<-pr.tot.u
ALGORITHM <-'normal'
par.tot<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

PARAMETERS = c("A1","A2","xc","p")
NUMSAMPLES=10000
PMIN<-pr.bio.l
PMAX<-pr.bio.u
ALGORITHM <-'normal'
par.bio<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

tot.c<-soil.c<-bio.c<-matrix(NA,nrow=NUMSAMPLES+1,ncol=262)

par.tot<-rbind(par.tot,parms.4log.tot[c("A1","A2","xc","p")])
par.bio<-rbind(par.bio,parms.4log.bio[c("A1","A2","xc","p")])
par.soil<-rbind(par.soil,parms.4log.soil[c("A1","A2","xc","p")])

for (i in 1:10001){
	A1=par.tot[i,1];A2=par.tot[i,2];xc=par.tot[i,3];p=par.tot[i,4]
	fun.tot<-expression(A2+(A1-A2)/(1+(x/xc)^p))
	fx.tot<-D(fun.tot,"x")
	x=11:272	
	tot.c[i,]<-eval(fx.tot)
}

for (i in 1:10001){
	A1=par.bio[i,1];A2=par.bio[i,2];xc=par.bio[i,3];p=par.bio[i,4]
	fun.bio<-expression(A2+(A1-A2)/(1+(x/xc)^p))
	fx.bio<-D(fun.bio,"x")
	x=11:272		
	bio.c[i,]<-eval(fx.bio)	
}
for (i in 1:10001){
	A1=par.soil[i,1];A2=par.soil[i,2];xc=par.soil[i,3];p=par.soil[i,4]
	fun.soil<-expression(A2+(A1-A2)/(1+(x/xc)^p))
	fx.soil<-D(fun.soil,"x")
	x=11:272
	soil.c[i,]<-eval(fx.soil)
}


tot.m<-apply(tot.c,2,FUN=mean,na.rm=T);tot.sd<-apply(tot.c,2,FUN=sd,na.rm=T);
bio.m<-bio.c[10001,];bio.sd<-apply(bio.c,2,FUN=sd,na.rm=T);
soil.m<-soil.c[10001,];soil.sd<-apply(soil.c,2,FUN=sd,na.rm=T);

plot(x,(soil.m/(soil.m+bio.m))*100,type='p',ylim=c(0,100),pch=21,xlab='Stand age(Year)',ylab="Contribution to total carbon stock(%)")
points(x,bio.m/(soil.m+bio.m)*100,type='p',pch=19,col=4)
grid()

pump.cc<-cbind(x,tot.m,tot.sd,bio.m,bio.sd,soil.m,soil.sd)

write.csv(pump.cc,file='d://pump_log_cc.csv',row.names=F)

#
#fun.tot<-expression(358.52+(182.22-358.52)/(1+(x/77.95)^2.8228))
#fun.bio<-expression(116.77+(31.10-116.77)/(1+(x/98.88)^2.72947))
#fun.soil<-expression(243.36767+(153.26201-243.36767)/(1+(x/66.76)^3.35015))
#fx.tot<-D(fun.tot,"x")
#fx.bio<-D(fun.bio,"x")
#fx.soil<-D(fun.soil,"x")
#plot(x,eval(fx.tot),lwd=2,col=1)
#points(x,eval(fx.soil),lwd=2,col=2,lty=2)
#points(x,eval(fx.bio),lwd=2,col=3,lty=3)
#grid()
#a=243.36767
#fun.soil<-expression(a+(153.26201-a)/(1+(x/66.76)^3.35015))
#fx.soil<-D(fun.soil,"x")
#points(x,eval(fx.soil),lwd=2,col=2,lty=2,type='l')
#points(x,eval(fx.soil),lwd=4,col=2,lty=4,type='l')
#points(x,eval(fx.soil),lwd=4,col=4,lty=4,type='l')
#points(x,eval(fx.soil),lwd=4,col=4,lty=1,type='l')
#

