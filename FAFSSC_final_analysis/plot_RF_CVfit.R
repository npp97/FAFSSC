#to do RF analysis and Covington's curve fitting on partial AGE response and real data
#Step 1: read data into memmory
require(gdata)
require(FME)
require(sfsmisc)

#-------------------------

Slog_model <- function(parms,x)
{
  with(as.list(parms), 
        return(a/(1+exp(-k*(x-xc)))));
}

## FITTING algorithm 1
Slog_Cost <- function(parms,obs){
  out <- Slog_model(parms, obs$x)
  return((obs$y-out)^2);  # residuals
}


CCurve_model <- function(parms,x)
{
  with(as.list(parms), 
        return(a*x^b*exp(c*x^d)+e));
}

## FITTING algorithm 1
CCurve_Cost <- function(parms,obs) {
  out <- CCurve_model(parms, obs$x)
  return((obs$y-out)^2);  # residuals
}



require(randomForestSRC)

nTree=200;

bio_soc_cls<-read.csv('bio_soc_cls.csv')

plot.cls<-as.character(levels(bio_soc_cls$cls_max_spec4))
ncls<-length(plot.cls)
rst.age_carbon<-data.frame(matrix(NA,nrow = 1, ncol = 4));

names(rst.age_carbon)<-c("VEG.SECTION","VARIABLEs","Plot_Age","CARBON")


for (i in 1:ncls){
  ii<-which(as.character(bio_soc_cls$cls_max_spec4) %in% plot.cls[i])
  data.to.anz<-bio_soc_cls[ii,]
  
  
  gf.veg<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec,data=data.to.anz,importance="permute.ensemble",ntree=nTree,seed=-111)
  gf.s20<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+aspect+slope+tpi+age_max+max_spec,data=data.to.anz,importance="permute.ensemble",ntree=nTree,seed=-111)
  gf.s1m<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+age_max+max_spec,data=data.to.anz,importance="permute.ensemble",ntree=nTree,seed=-111)
  
  gf.veg.p<-plot.variable(gf.veg,partial=TRUE,show.plots = FALSE);
  gf.s20.p<-plot.variable(gf.s20,partial=TRUE,show.plots = FALSE);
  gf.s1m.p<-plot.variable(gf.s1m,partial=TRUE,show.plots = FALSE);
  
  ij<-which(gf.veg.p$xvar.names %in% 'age_max')
  age.veg<-gf.veg.p$pData[[ij]]$x.uniq
  carbon.veg<-gf.veg.p$pData[[ij]]$yhat
  
  ik<-which(gf.s1m.p$xvar.names %in% 'age_max') 
  age.s1m<-gf.s1m.p$pData[[ik]]$x.uniq
  carbon.s1m<-gf.s1m.p$pData[[ik]]$yhat
  
  ih<-which(gf.s20.p$xvar.names %in% 'age_max')   
  age.s20<-gf.s20.p$pData[[ih]]$x.uniq
  carbon.s20<-gf.s20.p$pData[[ih]]$yhat
  
  varb<-rep(c("CBiomass","CSoil1m","CSoil20cm"),c(length(age.veg),length(age.s1m),length(age.s20)));
  vegs<-rep(plot.cls[i],(length(age.veg)+length(age.s1m)+length(age.s20)))
  ages<-c(age.veg,age.s1m,age.s20)  
  carbons<-c(carbon.veg,carbon.s1m,carbon.s20)
  
  rst.age_carbon<-rbind(rst.age_carbon,data.frame(VEG.SECTION=vegs,VARIABLEs=varb,Plot_Age=ages,CARBON=carbons))
  rm(vegs,varb,ages,carbons)
  
  print(plot.cls[i])
}
rst.age_carbon<-rst.age_carbon[-1,]

pdf('./figs/FIG_S6_RF_age_carbon.pdf',paper='A4')
xyplot(CARBON~Plot_Age|VEG.SECTION,groups =VARIABLEs,type=c('l','s','g'),data=rst.age_carbon, prepanel = function(x, y) prepanel.loess(x, y, span = 0.8),cex=0.5,col=1:3,lwd=1.5,lty=1:3,pch=1:3,
       key=list(lines=T,points = TRUE,col=1:3,lwd=1.5,lty=1:3,pch=1:3,cex=1,text=list(levels(as.factor(rst.age_carbon$VEG.SECTION))[c(1,2,3)]),corner=c(0.5,0.975)),
       xlab="Plot Age(Year)" ,ylab='Carbon Density (tC/ha)')
dev.off()

write.csv(rst.age_carbon,'rst.age_carbon.csv')

#Step regression with age using nlstools
require(nlstools)
require(onls)

#the names of output from Slogicals model
#VEG.SECTION, "VARIABLEs", a, k, xc, a.se, k.se, xc.se, rse, R2

fSlogistic.veg <- as.formula(CBiomass ~ a/(1+exp(-k*(age-xc))))
fSlogistic.s1m <- as.formula(CSoil1m ~ a/(1+exp(-k*(age-xc))))
fSlogistic.s20 <- as.formula(CSoil20cm ~  a/(1+exp(-k*(age-xc))))

fSlogistic.veg1 <- as.formula(CBiomass ~ a*age^b*exp(c*age^d)+e)
fSlogistic.s1m1 <- as.formula(CSoil1m ~ a*age^b*exp(c*age^d)+e)
fSlogistic.s201 <- as.formula(CSoil20cm ~ a*age^b*exp(c*age^d)+e)

read.csv("rst.age_carbon_wide.csv")->rst.age_carbon

plot.cls<-as.character(levels(rst.age_carbon$VEG.SECTION))
ncls<-length(plot.cls)
rst.age_carbon.sl<-data.frame(matrix(NA,nrow = ncls*3, ncol = 10));

names(rst.age_carbon.sl)<-c("VEG.SECTION","VARIABLEs", "a", "k", "xc", "a.se", "k.se", "xc.se", "rse", "R2")


for (i in 2:ncls){
  
  ii<-which(as.character(rst.age_carbon$VEG.SECTION) %in% plot.cls[i])
  d2z<-rst.age_carbon[ii,]
  names(d2z)<-c("VEG.SECTION","age","CBiomass","CSoil1m","CSoil20cm")

   s20.onls.slog <- onls(fSlogistic.s20, start = list(a = 50, k=-0.005, xc = 10), data = d2z);
  
  parms0<-c(a=coef(s20.onls.slog)[1],k=coef(s20.onls.slog)[2], xc = coef(s20.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='SANN')
  s20.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='Pseudo')
  

  s201.onls.slog <- onls(fSlogistic.s201, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z);
  
  parms0<-c(a = coef(s201.onls.slog)[1], b=coef(s201.onls.slog)[2], c = coef(s201.onls.slog)[3], d=coef(s201.onls.slog)[4], e=coef(s201.onls.slog)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='SANN')
  s20.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='Pseudo')

  #------------------------------------  
  s1m.onls.slog <- onls(fSlogistic.s20, start = list(a = 50, k=-0.005, xc = 10), data = d2z);
  
  parms0<-c(a=coef(s1m.onls.slog)[1],k=coef(s1m.onls.slog)[2], xc = coef(s1m.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method='SANN')
  s1m.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method='Pseudo')
  
  
  s1m1.onls.slog <- onls(fSlogistic.s1m1, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z);
  
  parms0<-c(a = coef(s1m1.onls.slog)[1], b=coef(s1m1.onls.slog)[2], c = coef(s1m1.onls.slog)[3], d=coef(s1m1.onls.slog)[4], e=coef(s1m1.onls.slog)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method='SANN')
  s1m.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method = 'Pseudo')
  
  s1m1.onls.slog <- onls(fSlogistic.s1m1, start = list(a = coef(st)[1], b=coef(st)[2], c = coef(st)[3], d=coef(st)[4], e=coef(st)[5]), data = d2z);
  #---------------------------------------------------------- 
  veg.onls.slog <- onls(fSlogistic.veg, start = list(a = 50, k=-0.005, xc = 10), data = d2z);
  
  parms0<-c(a=coef(veg.onls.slog)[1],k=coef(veg.onls.slog)[2], xc = coef(veg.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method='SANN')
  veg.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method='Pseudo')
  
  
  veg1.onls.slog <- onls(fSlogistic.veg1, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a = coef(veg1.onls.slog)[1], b=coef(veg1.onls.slog)[2], c = coef(veg1.onls.slog)[3], d=coef(veg1.onls.slog)[4], e=coef(veg1.onls.slog)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method='SANN')
  veg.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method = 'Pseudo',control=list(npop=200,numiter=30000))

   #---------------------------------------------------------- 
  
  
  R2.veg<-(1-veg.boot.slog$rse*(nrow(d2z)-3)/var(d2z$CBiomass));
  R2.s1m<-(1-veg.boot.slog$rse*(nrow(d2z)-3)/var(d2z$CSoil1m));
  R2.s20<-(1-veg.boot.slog$rse*(nrow(d2z)-3)/var(d2z$CSoil20cm));
  
  k=((i-1)*3+1)
  rst.age_carbon.sl[k,3:10]<-c(matrix(veg.boot.slog$estiboot,nrow=1),veg.boot.slog$rse,R2.veg);
  rst.age_carbon.sl[k+1,3:10]<-c(matrix(s1m.boot.slog$estiboot,nrow=1),s1m.boot.slog$rse,R2.s1m);
  rst.age_carbon.sl[k+2,3:10]<-c(matrix(s20.boot.slog$estiboot,nrow=1),s20.boot.slog$rse,R2.s20);
  rst.age_carbon.sl[c(K,K+1,k+2),1]<-plot.cls[i];
  rst.age_carbon.sl[c(K,K+1,k+2),2]<-c("CBiomass","CSoil1m","CSoil20cm");
  
  
  
  veg.onls.slog <- onls(fSlogistic.veg1, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z);
  veg.nls.slog <- onls(fSlogistic.veg1, start = list(a = coef(veg.onls.slog)[1], b=coef(veg.onls.slog)[2], c = coef(veg.onls.slog)[3], d=coef(veg.onls.slog)[4], e=coef(veg.onls.slog)[5]), data = d2z);
  s1m1.onls.slog <- onls(fSlogistic.s1m1, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z);
  s1m1.nls.slog <- nls(fSlogistic.s1m1, start = list(a = coef(s1m.onls.slog)[1], b=coef(s1m.onls.slog)[2], c = coef(s1m.onls.slog)[3], d=coef(s1m.onls.slog)[4], e=coef(s1m.onls.slog)[5]), data = d2z);
  s201.onls.slog <- onls(fSlogistic.s201, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z);
  s201.nls.slog <- nls(fSlogistic.s201, start = list(a = coef(s20.onls.slog)[1], b=coef(s20.onls.slog)[2], c = coef(s20.onls.slog)[3], d=coef(s20.onls.slog)[4], e=coef(s20.onls.slog)[5]), data = d2z);

  veg.boot.slog <- nlsBoot(veg.nls.slog , niter = 1000);
  s1m.boot.slog <- nlsBoot(s1m.nls.slog , niter = 1000);
  s20.boot.slog <- nlsBoot(s20.nls.slog , niter = 1000);

  R2.veg<-(1-veg.boot.slog$rse*(nrow(d2z)-3)/var(d2z$CBiomass));
  R2.s1m<-(1-veg.boot.slog$rse*(nrow(d2z)-3)/var(d2z$CSoil1m));
  R2.s20<-(1-veg.boot.slog$rse*(nrow(d2z)-3)/var(d2z$CSoil20cm));

  k=((i-1)*3+1)
  rst.age_carbon.sl[k,3:10]<-c(matrix(veg.boot.slog$estiboot,nrow=1),veg.boot.slog$rse,R2.veg);
  rst.age_carbon.sl[k+1,3:10]<-c(matrix(s1m.boot.slog$estiboot,nrow=1),s1m.boot.slog$rse,R2.s1m);
  rst.age_carbon.sl[k+2,3:10]<-c(matrix(s20.boot.slog$estiboot,nrow=1),s20.boot.slog$rse,R2.s20);
  rst.age_carbon.sl[c(K,K+1,k+2),1]<-plot.cls[i];
  rst.age_carbon.sl[c(K,K+1,k+2),2]<-c("CBiomass","CSoil1m","CSoil20cm");
}



