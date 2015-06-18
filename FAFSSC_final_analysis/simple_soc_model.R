#soil carbon budget
#sub1 : litter decomp model
#sub2 : soil decomp model

#
#----------Cal liter decomp
read.csv('soc_litter_input.csv')->llf_lr
approx(x=llf_lr$age,y=llf_lr$lf_tC_ha,11:272)->lf
approx(x=llf_lr$age,y=llf_lr$fr_tC_ha,11:272)->lr
approx(x=llf_lr$age,y=llf_lr$ts,11:272)->tsoil
approx(x=llf_lr$age,y=llf_lr$sm_pct,11:272)->sm
approx(x=llf_lr$age,y=llf_lr$soc_tC_ha,11:272)->soc
approx(x=llf_lr$age,y=llf_lr$fr2_tC_ha,11:272)->lr2
approx(x=llf_lr$age,y=llf_lr$lc_tot_tC_ha,11:272)->lc
approx(x=llf_lr$age,y=llf_lr$lc2_tot_tC_ha,11:272)->lc2

#to calculate the k of fine-root litter decomp after
#Silver, W. L., & Miya, R. K. (2001). Global patterns in root decomposition: comparisons of 
#climate and litter quality effects. Oecologia, 129(3), 407-419.
lr_k<-exp(-2.1+0.07*tsoil$y)


#assign C/N = 100, TN=1% according to the mean of CarbonPROJECT Data,
#to calculate the k of leaf litter decomp after
# Zhang, D., Hui, D., Luo, Y., & Zhou, G. (2008). Rates of litter decomposition 
# in terrestrial ecosystems: global patterns and controlling factors. Journal of Plant Ecology, 1(2), 85-93.

lf_k<-(-2.484 + 0.026*tsoil$y + 0.0287*110 + 0.461*1)

x=1
yl<-length(x)
xl=length(lf$x)

matrix(0,ncol=(length(lf$x)+yl-1),nrow=length(lf$x))->deco_fl
deco_fr<-deco_fl

for (i in 1:xl){
 k=(i-1)+x	
 deco_fl[i,k]<-0.2*lf$y[i]*exp(-1*lf_k[i]*x)
 deco_fl[i,k]<-0.245*lf$y[i]*lf_k[i]
}

d.l<-apply(deco_fl,2,sum)[1:xl]

x=1:2
yl<-length(x)
xl=length(lf$x)

matrix(0,ncol=(length(lf$x)+yl-1),nrow=length(lf$x))->deco_fr
for (i in 1:xl){
	k=(i-1)+x	
#	deco_fr[i,k]<-0.3*lr2$y[i]*exp(-1*lr_k[i]*x)
	deco_fr[i,k]<-0.245*lr2$y[i]/yl
}


d.r<-apply(deco_fr,2,sum)[1:xl]


x=1:600
yl<-length(x)
xl=length(lc$x)

matrix(0,ncol=(length(lf$x)+yl-1),nrow=length(lc$x))->deco_lc
for (i in 1:xl){
	k=(i-1)+x	
#	deco_fr[i,k]<-0.3*lr$y[i]*exp(-1*lr_k[i]*x)
	deco_lc[i,k]<-0.245*lc2$y[i]/yl
}

d.lc<-apply(deco_lc,2,sum)[1:xl]
#d.l<-0.201*lf$y*exp(-1*lf_k)
#d.r<-0.201*lr$y*exp(1*lr_k)

litter_input<-d.l+d.r+d.lc
#litter_input<-d.l+d.r
#--------------------to calculate the decomp of soc after Foley et al(1995)
#Foley, J. A. (1995). An equilibrium model of the terrestrial carbon budget. Tellus B, 47(3), 310-319.

q10a=exp(10*0.204*(1-tsoil$y/36.9))      #Kirschbaum, M. U. (1995). The temperature dependence of soil organic matter decomposition, and the effect of
										#global warming on soil organic C storage. Soil Biology and biochemistry, 27(6), 753-760.
q10b=0.54649+0.31165*log(soc$y+3.14537)  #Zhou et al.,Zhou, T., Shi, P., Hui, D., & Luo, Y. (2009).
										#Spatial patterns in temperature sensitivity of soil respiration in China: Estimation with inverse modeling. 
										#Science in China Series C: Life Sciences, 52(10), 982-989.
q10c=5.91118*exp(-tsoil$y/11.21539)+1.48429 #Kirschbaum, M. U. (1995). The temperature dependence of soil organic matter decomposition, and the effect of
										#global warming on soil organic C storage. Soil Biology and biochemistry, 27(6), 753-760.
fts1<-exp(log(1.9)/10*(tsoil$y-20))
fts2a<-exp(log(q10a)/10*(tsoil$y-25))
fts2b<-exp(log(q10b)/10*(tsoil$y-22.5))
fts2c<-exp(log(q10c)/10*(tsoil$y-25))

fsm<-0.25+0.75*(sm$y/0.3)         #0.35 is parameter from IBIS 6.4, averaged soil water capacity of forest soil
#fsm[]<-1

#
deco1=(1/15+1/750)*fts1*fsm*soc$y
deco2a=(1/15+1/750)*fts2a*fsm*soc$y
deco2b=(1/15+1/750)*fts2b*fsm*soc$y
deco2c=(1/15+1/750)*fts2c*fsm*soc$y

soc_deco<-data.frame(age=sm$x,deco1,deco2a,deco2b,deco2c)

deco<-data.frame(soc_deco,d.l,d.r,d.lc,litter_input,lf=lf$y,lr=lr$y,soc=soc$y)

win.graph()
#par(mfrow=c(1,3))
with(deco,{
#	plot(age,deco1,type='l',lwd=1.5);points(age,deco2a,type='l',lwd=1.5,col=2);points(age,deco2b,type='l',lwd=2,col=3);points(age,deco2c,type='l',col=4)
#	plot(age,I(litter_input-deco1),type='l',lwd=1.5);points(age,I(d.l+d.r-deco2a),type='l',lwd=1.5,col=2);points(age,I(d.l+d.r-deco2b),type='l',lwd=2,col=3);points(age,I(d.l+d.r-deco2c),type='l',col=4)
	plot(age,cumsum(litter_input-deco1)+150,type='l',lwd=1.5,lty=2,col=2,ylim=c(0,400));
#		points(age,cumsum(litter_input-deco2a)+150,type='l',lwd=1.5,col=2);
#		points(age,cumsum(litter_input-deco2b)+150,type='l',lwd=2,col=3);
#		points(age,cumsum(litter_input-deco2c)+150,type='l',col=4);
#		points(age,cumsum(litter_input-deco2c)+150,type='l',col=4);		
		points(age,cumsum(litter_input-0.55*deco2b-0.45*deco1)+150,type='l',col=1,lwd=3);
		points(age,cumsum(d.lc+d.l-0.55*deco2b-0.45*deco1)+150,type='l',lty=2,col=1,lwd=3);
		points(age,cumsum(d.r+d.lc-0.55*deco2b-0.45*deco1)+150,type='l',lty=2,col=2,lwd=3);
		points(age,cumsum(d.r+d.l-0.55*deco2b-0.45*deco1)+150,type='l',lty=2,col=3,lwd=3);
		points(age,cumsum(d.l-0.55*deco2b-0.45*deco1)+150,type='l',lty=3,col=1,lwd=3);	
		points(age,cumsum(d.r-0.55*deco2b-0.45*deco1)+150,type='l',lty=3,col=2,lwd=3);	
		points(age,cumsum(d.lc-0.55*deco2b-0.45*deco1)+150,type='l',lty=3,col=3,lwd=3);			
	points(age,socc)
	grid()
})


