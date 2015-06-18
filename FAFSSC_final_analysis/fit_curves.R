#to do RF analysis and Covington's curve fitting on partial AGE response and real data
#Step 1: read data into memmory
require(gdata)
require(FME)
require(Hmisc)
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


require(nlstools)
require(onls)

#the names of output from Slogicals model
#VEG.SECTION, "VARIABLEs", a, k, xc, a.se, k.se, xc.se, rse, R2

fSlogistic.tot <- as.formula(TotC_eco ~ a/(1+exp(-k*(age-xc))))
fSlogistic.veg <- as.formula(CBiomass ~ a/(1+exp(-k*(age-xc))))
fSlogistic.s1m <- as.formula(CSoil1m ~ a/(1+exp(-k*(age-xc))))
fSlogistic.s20 <- as.formula(CSoil20cm ~  a/(1+exp(-k*(age-xc))))

fCCurve.tot <- as.formula(TotC_eco ~ a*age^b*exp(c*age^d)+e)
fCCurve.veg <- as.formula(CBiomass ~ a*age^b*exp(c*age^d)+e)
fCCurve.s1m <- as.formula(CSoil1m ~ a*age^b*exp(c*age^d)+e)
fCCurve.s20 <- as.formula(CSoil20cm ~ a*age^b*exp(c*age^d)+e)

read.csv("bio_soc_cls.csv")->rst.age_carbon

rst.age_carbon<-rst.age_carbon[,c("cls_max_spec4","age_max","TotC_eco","W_total_t_ha" ,"SoilC_1m_t_ha","SoilCden_t_ha_0.20cm")]
names(rst.age_carbon)<-c("VEG.SECTION","age","TotC_eco","CBiomass","CSoil1m","CSoil20cm")
rst.age_carbon$CBiomass<-0.4*rst.age_carbon$CBiomass;

plot.cls<-as.character(levels(rst.age_carbon$VEG.SECTION))
ncls<-length(plot.cls)

print("NPLOS, VEG.SECTION,VARIABLEs, a, k, xc, a.se, k.se, xc.se, rse, R2")
print("MODFIT, VEG.SECTION,VARIABLEs, a, k, xc, a.se, k.se, xc.se, rse, R2")
print("NPLOS, VEG.SECTION,VARIABLEs, a, b, c, d, e, a.se, b.se, c.se, d.se, e.se, rse, R2")
print("MODFIT, VEG.SECTION,VARIABLEs, a, b, c, d, e, a.se, b.se, c.se, d.se, e.se, rse, R2")

sink('fit_log_ccurve.txt')

obiomass<-!is.na(rst.age_carbon$CBiomass)
osoil1m<-!is.na(rst.age_carbon$CSoil1m)
osoil20cm<-!is.na(rst.age_carbon$CSoil20cm)
otot<-!is.na(rst.age_carbon$TotC_eco)
for (i in 1:ncls){
  
  o<-which(as.character(rst.age_carbon$VEG.SECTION) %in% plot.cls[i])
  
  #------------------------------
  ii<-which(o&otot)
  a<-clowess(x=rst.age_carbon$age[ii],y=rst.age_carbon$TotC_eco[ii],f=0.25)
  d2z<-data.frame(age=a$x,TotC_eco=a$y)
  #names(d2z)<-c("VEG.SECTION","age","CBiomass","CSoil1m","CSoil20cm")
  
  tot.onls.slog <- onls(fSlogistic.tot, start = list(a = 50, k=-0.005, xc = 10), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a=coef(tot.onls.slog)[1],k=coef(tot.onls.slog)[2], xc = coef(tot.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$TotC_eco),method='SANN')
  tot.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$TotC_eco),method='Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-LOG","TotC_eco",coef(tot.onls.slog),deviance(tot.onls.slog),logLik(tot.onls.slog)),sep=',')
  print(cat(plot.cls[i],"MODFIT-LOG","TotC_eco",coef(tot.st.slog),deviance(tot.st.slog),tot.st.slog$ssr,sep=','))
  
  tot.onls.ccurve <- onls(fCCurve.tot, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a = coef(tot.onls.ccurve)[1], b=coef(tot.onls.ccurve)[2], c = coef(tot.onls.ccurve)[3], d=coef(tot.onls.ccurve)[4], e=coef(tot.onls.ccurve)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$TotC_eco),method='SANN')
  tot.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$TotC_eco),method='Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-CCU","TotC_eco",coef(tot.onls.ccurve),deviance(tot.onls.ccurve),logLik(tot.onls.ccurve),sep=','))
  print(cat(plot.cls[i],"MODFIT-CCU","TotC_eco",coef(tot.st.ccurve),deviance(tot.st.ccurve),tot.st.ccurve$ssr,sep=','))
  #--------------------------------------
  
  ii<-which(o&osoil20cm)
  a<-clowess(x=rst.age_carbon$age[ii],y=rst.age_carbon$CSoil20cm[ii],f=0.25)
  d2z<-data.frame(age=a$x,CSoil20cm=a$y)
  #names(d2z)<-c("VEG.SECTION","age","CBiomass","CSoil1m","CSoil20cm")
  
  s20.onls.slog <- onls(fSlogistic.s20, start = list(a = 50, k=-0.005, xc = 10), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a=coef(s20.onls.slog)[1],k=coef(s20.onls.slog)[2], xc = coef(s20.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='SANN')
  s20.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-LOG","CSoil20cm",coef(s20.onls.slog),deviance(s20.onls.slog),logLik(s20.onls.slog)),sep=',')
  print(cat(plot.cls[i],"MODFIT-LOG","CSoil20cm",coef(s20.st.slog),deviance(s20.st.slog),s20.st.slog$ssr,sep=','))
  
  s20.onls.ccurve <- onls(fCCurve.s20, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a = coef(s20.onls.ccurve)[1], b=coef(s20.onls.ccurve)[2], c = coef(s20.onls.ccurve)[3], d=coef(s20.onls.ccurve)[4], e=coef(s20.onls.ccurve)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='SANN')
  s20.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil20cm),method='Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-CCU","CSoil20cm",coef(s20.onls.ccurve),deviance(s20.onls.ccurve),logLik(s20.onls.ccurve),sep=','))
  print(cat(plot.cls[i],"MODFIT-CCU","CSoil20cm",coef(s20.st.ccurve),deviance(s20.st.ccurve),s20.st.ccurve$ssr,sep=','))
  
  #------------------------------------  
  ii<-which(o&osoil1m)
  a<-clowess(x=rst.age_carbon$age[ii],y=rst.age_carbon$CSoil1m[ii],f=0.25)
  d2z<-data.frame(age=a$x,CSoil1m=a$y)
  
  d2z<-rst.age_carbon[ii,]
  
  s1m.onls.slog <- onls(fSlogistic.s20, start = list(a = 50, k=-0.005, xc = 10), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a=coef(s1m.onls.slog)[1],k=coef(s1m.onls.slog)[2], xc = coef(s1m.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method='SANN')
  s1m.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method='Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-LOG","CSoil1m",coef(s1m.onls.slog),deviance(s1m.onls.slog),logLik(s1m.onls.slog),sep=','))
  print(cat(plot.cls[i],"MODFIT-LOG","CSoil1m",coef(s1m.st.slog),deviance(s1m.st.slog),s1m.st.slog$ssr,sep=','))
  
  
  s1m.onls.ccurve <- onls(fCCurve.s1m, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z, control = nls.lm.control(maxiter = 10000));
  parms0<-c(a = coef(s1m.onls.ccurve)[1], b=coef(s1m.onls.ccurve)[2], c = coef(s1m.onls.ccurve)[3], d=coef(s1m.onls.ccurve)[4], e=coef(s1m.onls.ccurve)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method='SANN')
  s1m.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CSoil1m),method = 'Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-CCU","CSoil1cm",coef(s1m.onls.ccurve),deviance(s1m.onls.ccurve),logLik(s1m.onls.ccurve),sep=','))
  print(cat(plot.cls[i],"MODFIT-CCU","CSoil1cm",coef(s1m.st.ccurve),deviance(s1m.st.ccurve),s1m.st.ccurve$ssr,sep=','))
  
  #---------------------------------------------------------- 
  ii<-which(o&obiomass)
  a<-clowess(x=rst.age_carbon$age[ii],y=rst.age_carbon$CBiomass[ii],f=0.25)
  d2z<-data.frame(age=a$x,CBiomass=a$y)
  
 # d2z<-rst.age_carbon[ii,]
  #------
  veg.onls.slog <- onls(fSlogistic.veg, start = list(a = 50, k=-0.005, xc = 10), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a=coef(veg.onls.slog)[1],k=coef(veg.onls.slog)[2], xc = coef(veg.onls.slog)[3]);names(parms0)<-c("a","k","xc")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  st<-modFit(f=Slog_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method='SANN')
  veg.st.slog<-modFit(f=Slog_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method='Pseudo')
  
  print(cat(plot.cls[i],"ONLS-LOG","CBiomass",coef(veg.onls.slog),deviance(veg.onls.slog),logLik(veg.onls.slog),sep=','))
  print(cat(plot.cls[i],"MODFIT-LOG","CBiomass",coef(veg.st.slog),deviance(veg.st.slog),veg.st.slog$ssr,sep=','))
  
  
  veg.onls.ccurve <- onls(fCCurve.veg, start = list(a = 50, b=-0.005, c = 0.00001, d=2, e=100), data = d2z, control = nls.lm.control(maxiter = 10000));
  
  parms0<-c(a = coef(veg.onls.ccurve)[1], b=coef(veg.onls.ccurve)[2], c = coef(veg.onls.ccurve)[3], d=coef(veg.onls.ccurve)[4], e=coef(veg.onls.ccurve)[5]);names(parms0)<-c("a","b","c","d","e")
  lower=apply(cbind(0.5*parms0,2.0*parms0),1,min);
  upper=apply(cbind(0.5*parms0,2.0*parms0),1,max);
  
  st<-modFit(f=CCurve_Cost,p=parms0,lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method='SANN')
  veg.st.ccurve<-modFit(f=CCurve_Cost,p=coef(st),lower=lower,upper=upper,obs=data.frame(x=d2z$age,y=d2z$CBiomass),method = 'Pseudo',control=list(numiter=10000))
  
  print(cat(plot.cls[i],"ONLS-CCU","CBiomass",coef(veg.onls.ccurve),deviance(veg.onls.ccurve),logLik(veg.onls.ccurve),sep=','))
  print(cat(plot.cls[i],"MODFIT-CCU","CBiomass",coef(veg.st.ccurve),deviance(veg.st.ccurve),veg.st.ccurve$ssr,sep=','))
  
  #---------------------------------------------------------- 
}

sink()

