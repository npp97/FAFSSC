# TODO: Add comment
# 
# Author: junhui.zhang
###############################################################################
#by Richards growth curves
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
parms.rich.tot<-c(a=348.81373,a.se=1.07009,xc=73.03939,xc.se=4.93934,d=7.87236,d.se=1.04871,k=0.04568,k.se=0.00491)
parms.rich.bio<-c(a1=109.37542,a1.se=0.72049,xc1=99.32194,xc1.se=4.49261,d1=4.57952,d1.se=0.6059,k1=0.03883,k1.se=0.00441)
parms.rich.soil<-c(a2=240.36798,a2.se=0.77074,xc2=90.3969,xc2.se=10.78625,d2=34.26419,d2.se=19.97034,k2=0.15712,k2.se=0.09012)

pr.tot.l<-c(a=(parms.rich.tot["a"]-2*parms.rich.tot["a.se"]),xc=(parms.rich.tot["xc"]-2*parms.rich.tot["xc.se"]),d=(parms.rich.tot["d"]-2*parms.rich.tot["d.se"]),k=(parms.rich.tot["k"]-2*parms.rich.tot["k.se"]))
pr.tot.u<-c(a=(parms.rich.tot["a"]+2*parms.rich.tot["a.se"]),xc=(parms.rich.tot["xc"]+2*parms.rich.tot["xc.se"]),d=(parms.rich.tot["d"]+2*parms.rich.tot["d.se"]),k=(parms.rich.tot["k"]+2*parms.rich.tot["k.se"]))
pr.bio.l<-c(a=(parms.rich.bio["a1"]-2*parms.rich.bio["a1.se"]),xc=(parms.rich.bio["xc1"]-2*parms.rich.bio["xc1.se"]),d=(parms.rich.bio["d1"]-2*parms.rich.bio["d1.se"]),k=(parms.rich.bio["k1"]-2*parms.rich.bio["k1.se"]))
pr.bio.u<-c(a=(parms.rich.bio["a1"]+2*parms.rich.bio["a1.se"]),xc=(parms.rich.bio["xc1"]+2*parms.rich.bio["xc1.se"]),d=(parms.rich.bio["d1"]+2*parms.rich.bio["d1.se"]),k=(parms.rich.bio["k1"]+2*parms.rich.bio["k1.se"]))
pr.soil.l<-c(a=(parms.rich.soil["a2"]-2*parms.rich.soil["a2.se"]),xc=(parms.rich.soil["xc2"]-2*parms.rich.soil["xc2.se"]),d=(parms.rich.soil["d2"]-2*parms.rich.soil["d2.se"]),k=(parms.rich.soil["k2"]-2*parms.rich.soil["k2.se"]))
pr.soil.u<-c(a=(parms.rich.soil["a2"]+2*parms.rich.soil["a2.se"]),xc=(parms.rich.soil["xc2"]+2*parms.rich.soil["xc2.se"]),d=(parms.rich.soil["d2"]+2*parms.rich.soil["d2.se"]),k=(parms.rich.soil["k2"]+2*parms.rich.soil["k2.se"]))

PARAMETERS = c("a2","xc2","d2","k2")
NUMSAMPLES=10000
PMIN<-pr.soil.l
PMAX<-pr.soil.u
ALGORITHM <-'normal'
par.soil<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

PARAMETERS = c("a","xc","d","k")
NUMSAMPLES=10000
PMIN<-pr.tot.l
PMAX<-pr.tot.u
ALGORITHM <-'normal'
par.tot<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

PARAMETERS = c("a1","xc1","d1","k1")
NUMSAMPLES=10000
PMIN<-pr.bio.l
PMAX<-pr.bio.u
ALGORITHM <-'normal'
par.bio<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

tot.c<-soil.c<-bio.c<-matrix(NA,nrow=NUMSAMPLES+1,ncol=262)

par.tot<-rbind(par.tot,parms.rich.tot[c("a","xc","d","k")])
par.bio<-rbind(par.bio,parms.rich.bio[c("a1","xc1","d1","k1")])
par.soil<-rbind(par.soil,parms.rich.soil[c("a2","xc2","d2","k2")])

for (i in 1:10001){
	a=par.tot[i,1];xc=par.tot[i,2];d=par.tot[i,3];k=par.tot[i,4]
	fun.tot<-expression(a*(1+(d-1)*exp(-1*k*(x-xc)))^(1/(1-d)))
	
	a1=par.bio[i,1];xc1=par.bio[i,2];d1=par.bio[i,3];k1=par.bio[i,4]
	fun.bio<-expression(a1*(1+(d1-1)*exp(-1*k1*(x-xc1)))^(1/(1-d1)))
	
	a2=par.soil[i,1];xc2=par.soil[i,2];d2=par.soil[i,3];k2=par.soil[i,4]
	fun.soil<-expression(a2*(1+(d2-1)*exp(-1*k2*(x-xc2)))^(1/(1-d2)))	
	
	fx.tot<-D(fun.tot,"x")
	fx.bio<-D(fun.bio,"x")
	fx.soil<-D(fun.soil,"x")
	x=11:272
	
	tot.c[i,]<-eval(fx.tot)
	bio.c[i,]<-eval(fx.bio)
	soil.c[i,]<-eval(fx.soil)
}



tot.m<-apply(tot.c,2,FUN=mean,na.rm=T);tot.sd<-apply(tot.c,2,FUN=sd,na.rm=T);
bio.m<-bio.c[10001,];bio.sd<-apply(bio.c,2,FUN=sd,na.rm=T);
soil.m<-soil.c[10001,];soil.sd<-apply(soil.c,2,FUN=sd,na.rm=T);

plot(x,(soil.m/(soil.m+bio.m))*100,type='p',ylim=c(0,100),pch=21,xlab='Stand age(Year)',ylab="Contribution to total carbon stock(%)")
points(x,bio.m/(soil.m+bio.m)*100,type='p',pch=19,col=4)
grid()

pump.cc<-cbind(x,tot.m,tot.sd,bio.m,bio.sd,soil.m,soil.sd)

write.csv(pump.cc,file='d://pump_cc.csv',row.names=F)

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
