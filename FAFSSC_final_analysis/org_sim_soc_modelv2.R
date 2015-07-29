# TODO: to reorganize the simple_soc_model into a function to calibriate the parameters
# 
# Author: junhui.zhang
###############################################################################

require(FME)
require(lhs)

sumfun <- function(x, ...) {
	c(m = mean(x, ...), md = median(x, ...), v = sd(x, ...), l = length(x, ...))
}
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

ssoc<-function(parms,obs){
	
source("D://workDir//git//FAFSSC//FAFSSC_final_analysis//litters_matv2.R")	
#Leaf area index
	
#Parameters defination & Function used
	
	#(1) DOC exported from plot
	#  															dDOC = a_doc*ln(lai)+b_doc; 
	# DOC production, absorb and leaching is the output from complex interactions among hydrology, biochemiclogical, microclimate and biological
	# Up-to-now, there is no aggrement on the formulation and few direct measurement on that, so we formulated a log response function on several reports like Yan et al.,2015, 
	# Peichl et al.,2007�� Close relationship exist between DOCexport and several biological variables,like DBH, above-ground biomass, biomass of fine-roots, LAI etc.
	# For simplicity and sensitive, LAI was selected to describe the dynamics of DOCexpor for the LAI bridges the hydrological, biological and biogeochemicial processes. LAI has 
	# been used as the response variable for precipation interception and indicator of canopy development during restoration, Leaf litter is one of the most important litter that
	# produce DOC.
	a_doc <- parms["a_doc"];
	b_doc <- parms["b_doc"];
	
	#(2) Leaf decomposition
	# assign C/N = 100, TN=1% according to the mean of CarbonPROJECT Data, to calculate the k of leaf litter decomp after
	# Zhang, D., Hui, D., Luo, Y., & Zhou, G. (2008). Rates of litter decomposition 
	# in terrestrial ecosystems: global patterns and controlling factors. Journal of Plant Ecology, 1(2), 85-93.

	# lf_k<-(-2.484 + 0.026*TSoil_C + 0.0287*C/N + 0.461*TN_%)
	# then the lf_k can be shorten as lf_k<-a_lf+b_lf*Tsoil_C 
#	a_lf <- parms$a_lf
#	b_lf <- parms$b_lf
	
	#(3) Fine roots decomposition
	# to calculate the k of fine-root litter decomp after "Silver, W. L., & Miya, R. K. (2001). Global patterns in root decomposition: comparisons of 
	# climate and litter quality effects. Oecologia, 129(3), 407-419."
	# lr_k <- exp(-2.1+0.07*TSoil_C)
	# lr_k <- exp(a_lr + b_lr*TSoil_C)
#	a_lr <- parms$a_lr
#	b_lr <- parms$b_lr
	
	#(4) Soil decomposition
	# to calculate the decomp of soc after Foley et al(1995)
	# Foley, J. A. (1995). An equilibrium model of the terrestrial carbon budget. Tellus B, 47(3), 310-319.
	# and  #Zhou et al.,Zhou, T., Shi, P., Hui, D., & Luo, Y. (2009).#Spatial patterns in temperature sensitivity of soil respiration in China: Estimation
	# with inverse modeling. Science in China Series C: Life Sciences, 52(10), 982-989.
	# q10b=0.54649+0.31165*log(SOC_MgC_ha+3.14537) 
	# or q10c=5.91118*exp(-tsoil$y/11.21539)+1.48429 #Kirschbaum, M. U. (1995). The temperature dependence of soil organic matter decomposition, and the effect of
	#global warming on soil organic C storage. Soil Biology and biochemistry, 27(6), 753-760.
	
	# fts2a<-exp(log(q10a)/10*(TSoil_C-Ts_opt)) Soil temperature correstion
	# fts<-0.8*exp(0.095Ts)
	
	ts_a <- parms["ts_a"]
	ts_b <- parms["ts_b"]
	# fsm<-min(c(1.0,0.25+0.75*(Soil_moisture_[0,1]/AWC)))  Soil moisture correction

	awc <- parms["awc"]

	kbase_dpm<-parms["kbase_dpm"]
	kbase_spm<-parms["kbase_spm"]
	kbase_rpm<-parms["kbase_rpm"]
	
    kbase_som1<-parms["kbase_som1"] 
	kbase_som2<-parms["kbase_som2"]
	kbase_som3<-parms["kbase_som3"]
	kbase_som4<-parms["kbase_som4"]
	
	alf_dpm2som1<-parms["alf_dpm2som1"]
	alf_spm2som2<-parms["alf_spm2som2"]
	alf_rpm2som3<-parms["alf_rpm2som3"]
	
	alf_som1to2<-parms["alf_som1to2"]
	alf_som2to3<-parms["alf_som2to3"]
	alf_som3to4<-parms["alf_som3to4"]

#Driving data tables
#	age <- 11:272; #age
	#soil temperature
		
#ii<-which(age %in% times)
#tsoil_C<-approx(times,tsoil_C[ii],11:272)
#---------------------------------------
#lc=obs$lc
#lai=obs$lai
#tsoil_C=obs$tsoil_C
#sm=obs$sm
#lf=obs$lf
#lr=obs$lr
soc=obs$soc

#to prepare the lcf
x=1:round(parms["turnover.cwd"])
yl<-length(x)
xl=length(lc)

matrix(0,ncol=(length(lc)+yl-1),nrow=length(lc))->deco_lc
for (i in 1:xl){
	k=(i-1)+x	
#	deco_fr[i,k]<-0.3*lr$y[i]*exp(-1*lr_k[i]*x)
	deco_lc[i,k]<-lc[i]/yl
}

lcf<-apply(deco_lc,2,sum)[1:xl]
# to calculate the influence of tsoil and soil moisture
fts<-parms["ts_a"]*exp(parms["ts_b"]*tsoil_C)
fsm<-apply(cbind(rep(1.0,xl),0.25+0.75*sm/awc),1,min)
#fsm[]<-1
fts<-apply(cbind(fts,matrix(1,ncol=1,nrow=length(fts))),1,min)
ef<-fts*fsm
#to part the leaf, root and cwd into dpm, spm and rpm pools
#Calculation code
dpm_lf<-0.39*lf
spm_lf<-0.44*lf
rpm_lf<-0.17*lf

dpm_lr<-0.30*lr
spm_lr<-0.45*lr
rpm_lr<-0.25*lr

spm_lcf<-0.76*lcf
rpm_lcf<-0.24*lcf

som0_1<-0.1*soc[1]*1.6
som0_2<-0.2*soc[1]*1.6
som0_3<-0.3*soc[1]*1.6
som0_4<-0.3*soc[1]*1.6

dpm.leached<-a_doc*log(lf)+b_doc
#dpm.in<-dpm_lf+dpm_lr-dpm.leached
dpm.in<-dpm_lf+dpm_lr
spm.in<-spm_lf+spm_lr+spm_lcf
rpm.in<-rpm_lf+rpm_lr+rpm_lcf

dpm0<-dpm.in[1];spm0<-spm.in[1];rpm0<-rpm.in[1]

spm<-rpm<-dpm<-matrix(0,ncol=1,nrow=xl)
som1<-som2<-som3<-som4<-matrix(0,nrow=length(soc),ncol=1)


for (i in 1:xl){
	dpm[i]<-dpm0+dpm.in[i]
	dpm.out<-dpm[i]*fts[i]*kbase_dpm
	dpm[i]<-dpm[i]-dpm.out
	
	spm[i]<-spm0+spm.in[i]
	spm.out<-spm[i]*fts[i]*kbase_spm
	spm[i]<-spm[i]-spm.out
	
	rpm[i]<-rpm0+rpm.in[i]
	rpm.out<-rpm[i]*fts[i]*kbase_rpm
	rpm[i]<-rpm[i]-rpm.out
	
	dpm0<-dpm[i]
	spm0<-spm[i]
	rpm0<-rpm[i]
	
	som1[i]<-som0_1+dpm.out*alf_dpm2som1
	som1.out<-som1[i]*ef[i]*kbase_som1
	som1[i]<-som1[i]-som1.out

	som2[i]<-som0_2+som1.out*alf_som1to2+spm.out*alf_spm2som2
	som2.out<-som2[i]*ef[i]*kbase_som2
	som2[i]<-som2[i]-som2.out
	
	som3[i]<-som0_3+som2.out*alf_som2to3+rpm.out*alf_rpm2som3
	som3.out<-som3[i]*ef[i]*kbase_som3
	som3[i]<-som3[i]-som3.out
	
	som4[i]<-som0_4+som3.out*alf_som3to4
	som4.out<-som4[i]*ef[i]*kbase_som4
	som4[i]<-som4[i]-som4.out
	
#	som0_1<-som1[i]
#	som0_2<-som2[i]
#	som0_3<-som3[i]
#	som0_4<-som4[i]
}

soc<-som1+som2+som3+som4+rpm+spm+dpm-dpm.leached

return(data.frame(time=11:272,soc=soc))
#Return result
}

ssoc_sim<-function(parms,obs){
	
lc<-obs$lc
lf<-obs$lf
lr<-obs$lr
tsoil_C<-obs$tsoil_C
sm<-obs$sm
lai<-obs$lai
lf2<-obs$lf2

a_doc <- parms["a_doc"];
b_doc <- parms["b_doc"];
	
ts_a <- parms["ts_a"]
ts_b <- parms["ts_b"]

awc <- parms["awc"]
	
kbase_dpm<-parms["kbase_dpm"]
kbase_spm<-parms["kbase_spm"]
kbase_rpm<-parms["kbase_rpm"]
	
kbase_som1<-parms["kbase_som1"] 
kbase_som2<-parms["kbase_som2"]
kbase_som3<-parms["kbase_som3"]
kbase_som4<-parms["kbase_som4"]
	
alf_dpm2som1<-parms["alf_dpm2som1"]
alf_spm2som2<-parms["alf_spm2som2"]
alf_rpm2som3<-parms["alf_rpm2som3"]
	
alf_som1to2<-parms["alf_som1to2"]
alf_som2to3<-parms["alf_som2to3"]
alf_som3to4<-parms["alf_som3to4"]
soc=obs$soc
	
#to prepare the lcf
x=1:round(parms["turnover.cwd"])
yl<-length(x)
xl=length(lc)
	
matrix(0,ncol=(length(lc)+yl-1),nrow=length(lc))->deco_lc
for (i in 1:xl){
	k=(i-1)+x	
#	deco_fr[i,k]<-0.3*lr$y[i]*exp(-1*lr_k[i]*x)
	deco_lc[i,k]<-lc[i]/yl
}
	
lcf<-apply(deco_lc,2,sum)[1:xl]
# to calculate the influence of tsoil and soil moisture
fts<-parms["ts_a"]*exp(parms["ts_b"]*tsoil_C)
fsm<-apply(cbind(rep(1.0,xl),0.25+0.75*sm/awc),1,min)
#fsm[]<-1
fts<-apply(cbind(fts,matrix(1,ncol=1,nrow=length(fts))),1,min)
ef<-fts*fsm
#to part the leaf, root and cwd into dpm, spm and rpm pools
#Calculation code
dpm_lf<-0.39*lf
spm_lf<-0.44*lf
rpm_lf<-0.17*lf
	
dpm_lr<-0.30*lr
spm_lr<-0.45*lr
rpm_lr<-0.25*lr
	
spm_lcf<-0.76*lcf
rpm_lcf<-0.24*lcf
	
som0_1<-0.1*soc[1]*1.6
som0_2<-0.2*soc[1]*1.6
som0_3<-0.3*soc[1]*1.6
som0_4<-0.3*soc[1]*1.6
	
dpm.leached<-a_doc*log(lf2)+b_doc
#dpm.in<-dpm_lf+dpm_lr-dpm.leached
dpm.in<-dpm_lf+dpm_lr
spm.in<-spm_lf+spm_lr+spm_lcf
rpm.in<-rpm_lf+rpm_lr+rpm_lcf
	
dpm0<-dpm.in[1];spm0<-spm.in[1];rpm0<-rpm.in[1]
	
spm<-rpm<-dpm<-matrix(0,ncol=1,nrow=xl)
som1<-som2<-som3<-som4<-matrix(0,nrow=length(soc),ncol=1)
	
	
for (i in 1:xl){
	dpm[i]<-dpm0+dpm.in[i]
	dpm.out<-dpm[i]*fts[i]*kbase_dpm
	dpm[i]<-dpm[i]-dpm.out
	
	spm[i]<-spm0+spm.in[i]
	spm.out<-spm[i]*fts[i]*kbase_spm
	spm[i]<-spm[i]-spm.out
		
	rpm[i]<-rpm0+rpm.in[i]
	rpm.out<-rpm[i]*fts[i]*kbase_rpm
	rpm[i]<-rpm[i]-rpm.out
		
	dpm0<-dpm[i]
	spm0<-spm[i]
	rpm0<-rpm[i]
		
	som1[i]<-som0_1+dpm.out*alf_dpm2som1
	som1.out<-som1[i]*ef[i]*kbase_som1
	som1[i]<-som1[i]-som1.out
		
	som2[i]<-som0_2+som1.out*alf_som1to2+spm.out*alf_spm2som2
	som2.out<-som2[i]*ef[i]*kbase_som2
	som2[i]<-som2[i]-som2.out
		
	som3[i]<-som0_3+som2.out*alf_som2to3+rpm.out*alf_rpm2som3
	som3.out<-som3[i]*ef[i]*kbase_som3
	som3[i]<-som3[i]-som3.out
		
	som4[i]<-som0_4+som3.out*alf_som3to4
	som4.out<-som4[i]*ef[i]*kbase_som4
	som4[i]<-som4[i]-som4.out
		
#	som0_1<-som1[i]
#	som0_2<-som2[i]
#	som0_3<-som3[i]
#	som0_4<-som4[i]
	}
	
soc<-som1+som2+som3+som4+rpm+spm+dpm-dpm.leached
	
return(data.frame(time=11:272,soc=soc))
#Return result
}


ssocCOST<-function(parms,obs){
	out<-ssoc(parms,obs);
#	cost_1<-modCost(model=out,parms,obs=obs);
	return(modCost(model=out,obs=obs))
}

ssocCOST1<-function(parms,obs){
	out<-ssoc(parms,obs);
#	cost_1<-modCost(model=out,parms,obs=obs);
	return(sum((obs$soc-out$soc)^2));  # residuals
}

mdMCMC<-function(ff=ssocCOST,parms,obs,niter=35000,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
	Fit <- modFit(f=ff,p=parms,obs=obs,lower=lower,upper=upper,method='Pseudo');
	parms<-Fit$par
	MCMC1 <- modMCMC(f=ff, p=parms, niter=10000,obs=obs,
			wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
	lower=(summary(MCMC1)[5,])*0.5;
	upper=summary(MCMC1)[7,]*1.5;
	a<-cbind(lower,upper,parms)
	lower<-apply(a,1,min)
	upper<-apply(a,1,max)
	Cov0<-cov(MCMC1$pars)*2.4^2/(length(parms));
	parms=unlist(summary(MCMC1)[6,]);
	MCMC <- modMCMC(f=ff, p=parms, niter=niter,obs=obs,jump=Cov0,ntrydr = 3,
			wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper, burninlength=burninlength);
	return(MCMC);
}
	
pars0 = c(a_doc =-0.4, 
		b_doc= 0.78, 
		ts_a = 0.8 , 
		ts_b = 0.095,
		awc =0.32,
		turnover.cwd=100,
		kbase_dpm = 0.7,
		kbase_spm= 0.07,
		kbase_rpm= 0.014,
		kbase_som1= 0.07,
		kbase_som2= 0.014,
		kbase_som3=0.0014,
		kbase_som4=0.0001,
		alf_dpm2som1=0.61,
		alf_spm2som2=0.45,
		alf_rpm2som3=0.61,
		alf_som1to2=0.72,
		alf_som2to3=0.54,
		alf_som3to4=0.45
);
#
#lower = c(a_doc =-0.5, 
#		b_doc= 0.5, 
#		ts_a = 0.6, 
#		ts_b = 0.05,
#		awc = 0.2,
#		turnover.cwd=100,
#		kbase_dpm = 0.2,
#		kbase_spm= 0.01,
#		kbase_rpm= 0.002,
#		kbase_som1= 0.01,
#		kbase_som2= 0.001,
#		kbase_som3=0.0001,
#		kbase_som4=0.00001,
#		alf_dpm2som1=0.3,
#		alf_spm2som2=0.3,
#		alf_rpm2som3=0.3,
#		alf_som1to2=0.5,
#		alf_som2to3=0.5,
#		alf_som3to4=0.4
#);
#
#upper = c(a_doc =-0., 
#		b_doc= 1.0, 
#		ts_a = 1 , 
#		ts_b = 0.2,
#		awc = 0.5,
#		turnover.cwd=1000,
#		kbase_dpm = 1,
#		kbase_spm= 0.15,
#		kbase_rpm= 0.03,
#		kbase_som1= 0.2,
#		kbase_som2= 0.05,
#		kbase_som3=0.002,
#		kbase_som4=0.0002,
#		alf_dpm2som1=0.7,
#		alf_spm2som2=0.7,
#		alf_rpm2som3=0.85,
#		alf_som1to2=0.9,
#		alf_som2to3=0.8,
#		alf_som3to4=0.7
#);


lower = c(a_doc =-1, 
		b_doc= 0.5, 
		ts_a = 0.6, 
		ts_b = 0.05,
		awc = 0.2,
		turnover.cwd=50,
		kbase_dpm = 0.2,
		kbase_spm= 0.01,
		kbase_rpm= 0.002,
		kbase_som1= 0.01,
		kbase_som2= 0.001,
		kbase_som3=0.0001,
		kbase_som4=0.00001,
		alf_dpm2som1=0.3,
		alf_spm2som2=0.3,
		alf_rpm2som3=0.3,
		alf_som1to2=0.5,
		alf_som2to3=0.5,
		alf_som3to4=0.4
);

upper = c(a_doc =-0., 
		b_doc= 2.0, 
		ts_a = 1 , 
		ts_b = 0.5,
		awc = 0.5,
		turnover.cwd=1000,
		kbase_dpm = 1,
		kbase_spm= 1,
		kbase_rpm= 1,
		kbase_som1= 1,
		kbase_som2= 1,
		kbase_som3=1,
		kbase_som4=1,
		alf_dpm2som1=0.7,
		alf_spm2som2=0.7,
		alf_rpm2som3=0.85,
		alf_som1to2=0.9,
		alf_som2to3=0.8,
		alf_som3to4=0.7
);
source("D://workDir//git//FAFSSC//FAFSSC_final_analysis//litters_matv2.R")	


times<-11:272

obs<-data.frame(time=times,soc=iodata$soc)	

mcfit<-modFit(f=ssocCOST,p=pars0,obs=obs,lower=lower,upper=upper,method="Pseudo")
plot(ssoc(coef(mcfit),obs))
points(obs,col=2)

#y<-(-164.39588*x^(-0.16074)*exp(-1.0748E-06*x^3.01778)+254.96)
#y.sim<-(-191.26625*x^(-0.23265)*exp(-1.300858E-06*x^4.49946)+253.114)
#> points(x,y.sim,pch='+',col=3)
#> y.sim<-(-191.26625*x^(-0.23265)*exp(-1.300858E-09*x^4.49946)+253.114)
#> points(x,y.sim,pch='+',col=3)
#> plot(ssoc(coef(mcfit),obs))
#> points(obs,col=2)
#> points(x,y.sim,pch='+',col=3)
#> points(x,y,pch='+',col=4)

#
#mcfit<-modFit(f=ssocCOST,p=coef(mcfit),obs=obs,lower=lower,upper=upper,method="Pseudo")
#MCMC<-mdMCMC(ff=ssocCOST,parms=coef(mcfit),obs=obs,niter=50000,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=1000);

mcfit<-modFit(f=ssocCOST,p=coef(mcfit),obs=obs,lower=lower,upper=upper,method="Pseudo",control=list(numiter=50000))
#MCMC<-modMCMC(f=ssocCOST,p=coef(mcfit),obs=obs,niter=150000,wvar0=0.1,updatecov=100,jump=0.1,lower=lower,upper=upper,burninlength=5000,ntrydr = 3);
#print(summary(MCMC))
#### MCMC best parameters
##------------------------------------------------------------------------------
#print(MCMC$bestpar)
##mcfit<-modFit(f=ssocCOST,p=pars0,obs=obs, method="Pseudo",lower=lower,upper=upper)

a<-cbind(coef(mcfit)*0.5,lower)
lower=apply(a,1,max)
a<-cbind(coef(mcfit)*1.5,upper)
upper=apply(a,1,min)
a<-cbind(lower,upper)
lower<-apply(a,1,min)
upper<-apply(a,1,max)

PARAMETERS = c("a_doc","b_doc","ts_a","ts_b","awc","turnover.cwd","kbase_dpm","kbase_spm","kbase_rpm","kbase_som1","kbase_som2","kbase_som3","kbase_som4","alf_dpm2som1",
		"alf_spm2som2","alf_rpm2som3","alf_som1to2","alf_som2to3","alf_som3to4");
NUMSAMPLES=10000
PMIN<-lower
PMAX<-upper
ALGORITHM <-'normal'
parms<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

rst<-data.frame(parms,rss=matrix(NA,ncol=1,nrow=nrow(parms)))
for (i in 1:nrow(parms)){
	rst[i,20]<-ssocCOST1(parms[i,1:19],obs)
}

ii<-which(rst[,20]<quantile(rst[,20],0.001))
summary(rst[ii,])

#to narrow the parameters range 
low<-apply(rst[ii,],2,min);upp<-apply(rst[ii,],2,max)
PMIN<-low
PMAX<-upp
NUMSAMPLES=100000
parms<-parms_gen_lhs(PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 

rst2<-data.frame(parms,rss=matrix(NA,ncol=1,nrow=nrow(parms)))
for (i in 1:nrow(parms)){
	rst2[i,20]<-ssocCOST1(parms[i,1:19],obs)
}

ii<-which(rst2[,20]<quantile(rst2[,20],0.0005))
summary(rst2[ii,])

#to check the contribution of each input and soil temperature/moisture influence
iodata.lrf<-iodata
iodata.lcf<-iodata
iodata.lff<-iodata
iodata.soc<-iodata
iodata.tsoil<-iodata
iodata.sm<-iodata
iodata.st<-iodata

iodata.lrf$lr[]<-0; #check contribution of fine roots
iodata.lcf$lc[]<-0; #check contribution of CWD
iodata.lff$lf[]<-0; #check contribution of leaf litter
iodata.soc$soc[]<-0; #check contribution of soc litter
iodata.tsoil$tsoil_C[]<-mean(iodata.tsoil$tsoil_C); #check the influence of soil temperature
iodata.sm$sm[]<-mean(iodata.sm$sm); #check influence of soil moisture
iodata.st$sm[]<-mean(iodata.st$sm);iodata.st$tsoil_C[]<-mean(iodata.st$tsoil_C)

rst.lrf<-cbind(rst2[ii,1:19],matrix(0,ncol=263,nrow=length(ii)))
names(rst.lrf)<-c(names(rst2)[1:19],'rss',paste("year",11:272,sep=''))
rst.tsdi<-rst.smdi<-rst.stdi<-rst.tsoil<-rst.sm<-rst.norm<-rst.lcf<-rst.lff<-rst.soc<-rst.tsoil<-rst.sm<-rst.st<-rst.lrf
rst.leach<-rst.scdi<-rst.lrdi<-rst.lfdi<-rst.lcdi<-rst.lrf
for (i in 1:length(ii)){
	rst.lrf[i,21:282]<-ssoc_sim(unlist(as.list(rst.lrf[i,1:19])),iodata.lrf)$soc
	rst.lcf[i,21:282]<-ssoc_sim(unlist(as.list(rst.lcf[i,1:19])),iodata.lcf)$soc
	rst.lff[i,21:282]<-ssoc_sim(unlist(as.list(rst.lff[i,1:19])),iodata.lff)$soc
	rst.soc[i,21:282]<-ssoc_sim(unlist(as.list(rst.soc[i,1:19])),iodata.soc)$soc	
	rst.tsoil[i,21:282]<-ssoc_sim(unlist(as.list(rst.tsoil[i,1:19])),iodata.tsoil)$soc	
	rst.sm[i,21:282]<-ssoc_sim(unlist(as.list(rst.sm[i,1:19])),iodata.sm)$soc		
	rst.st[i,21:282]<-ssoc_sim(unlist(as.list(rst.st[i,1:19])),iodata.st)$soc	
	rst.norm[i,21:282]<-ssoc_sim(unlist(as.list(rst.norm[i,1:19])),iodata)$soc	
	rst.leach[i,21:282]<-rst.leach[i,"b_doc"]+rst.leach[i,"a_doc"]*log(iodata$lf2)
}

rst.stdi[,21:282]<-rst.norm[,21:282]-rst.st[,21:282]
rst.smdi[,21:282]<-rst.norm[,21:282]-rst.sm[,21:282]
rst.tsdi[,21:282]<-rst.norm[,21:282]-rst.tsoil[,21:282]
rst.lrdi[,21:282]<-rst.norm[,21:282]-rst.lrf[,21:282]
rst.lfdi[,21:282]<-rst.norm[,21:282]-rst.lff[,21:282]
rst.lcdi[,21:282]<-rst.norm[,21:282]-rst.lcf[,21:282]
rst.scdi[,21:282]<-rst.norm[,21:282]-rst.soc[,21:282]


#
#require(randomForestSRC)
#gf<-rfsrc(rss~.,data=rst[,c(1:16,20)],importance="permute.ensemble",nTree=200, seed=-111)
#gfp<-plot.variable(gf,partial=T)

rst.lrf.s<-apply(rst.lrdi,2,sumfun)
rst.lcf.s<-apply(rst.lcdi,2,sumfun)
rst.lff.s<-apply(rst.lfdi,2,sumfun)
rst.soc.s<-apply(rst.scdi,2,sumfun)
rst.norm.s<-apply(rst.norm,2,sumfun)
rst.stdi.s<-apply(rst.stdi,2,sumfun)
rst.smdi.s<-apply(rst.smdi,2,sumfun)
rst.tsdi.s<-apply(rst.tsdi,2,sumfun)
rst.leach.s<-apply(rst.leach,2,sumfun)
rst.s<-rbind(rst.lrf.s[c(1,3,4),21:282],rst.lff.s[c(1,3,4),21:282],rst.lcf.s[c(1,3,4),21:282],rst.stdi.s[c(1,3,4),21:282],rst.smdi.s[c(1,3,4),21:282],rst.tsdi.s[c(1,3,4),21:282],rst.leach.s[c(1,3,4),21:282])

rst.s1<-data.frame(year=11:272,t(rst.s));names(rst.s1)<-c('year','lrf.m','lrf.sd','lrf.l','lff.m','lff.sd','lff.l','lcf.m','lcf.sd','lcf.l','stdi.m','stdi.sd','stdi.l','smdi.m','smdi.sd','smdi.l','tsdi.m','tsdi.sd','tsdi.l','leach.m','leach.sd','leach.l')

write.csv(rst.s1,'rst.s1.csv')
save.image(file='soc_sim_out.RData')

plot(rst.lrf.s[1,21:282],ylim=c(-20,100))
points(rst.lff.s[1,21:282],col=2)
points(rst.lcf.s[1,21:282],col=3)
points(rst.stdi.s[1,21:282],col=4)
points(rst.smdi.s[1,21:282],col=5)
points(rst.tsdi.s[1,21:282],col=6)
grid()

plot(rst.lrf.s[1,21:282],ylim=c(0,120))
points(rst.lff.s[1,21:282]+rst.lrf.s[1,21:282],col=2)
points(rst.lff.s[1,21:282]+rst.lrf.s[1,21:282]+rst.lcf.s[,1,21:282],col=3)
grid()

imin<-which.min(rst2$rss)

cbind(ssoc_sim(unlist(as.list(rst[imin,1:19])),iodata),soc,rst.norm.s)

#Fig S5 scatter age
require(lattice)
Sys.setenv('R_GSCMD'="C://Program Files//gs//gs9.16//bin//gswin64c.exe");
read.csv('carbons_age.csv',colClasses=c('character',"character",'numeric',"character",'numeric'))->cage.long

tiff('./FIG_S5_scatter_age.tif')
ii<-which(cage.long$cls_plot %in% c("C4C","C4D","D1B","C2B","C3B","C3C","B4C","B4D","C1B","B2B","B2D","A2D","A2E","A3B","A3C","A2A","A2B","A2C"))

xyplot(Carbon_t_ha~age_max|cls_plot,groups =Cls_C,type=c('p','smooth','g'),data=cage.long[ii,], prepanel = function(x, y) prepanel.loess(x, y, span = 0.8),cex=0.5,col=1:3,lwd=1.5,lty=1:3,pch=1:3,
		key=list(lines=T,points = TRUE,col=1:3,lty=1:3,lwd=1.5,pch=1:3,cex=1,text=list(levels(as.factor(cage.long$Cls_C))[c(1,2,3)]),corner=c(0.995,0.975)),
		xlab="Stand Age(Year)" ,ylab='Carbon Density (tC/ha)')

dev.off()


