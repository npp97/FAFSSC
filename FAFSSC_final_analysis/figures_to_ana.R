require(randomForestSRC)
require(strucchange)

bio_soc_cls<-read.csv('bio_soc_cls.csv')

bio_soc_cls$FLR=bio_soc_cls$W_foliage_kg_tC_ha+bio_soc_cls$W_leaf_kg_tC_ha+bio_soc_cls$W_root_kg_tC_ha
bio_soc_cls$total_C<-bio_soc_cls$W_total_tC_ha+bio_soc_cls$SoilC_1m_t_ha
bio_soc_cls$total_Cs<-bio_soc_cls$W_total_tC_ha+bio_soc_cls$SoilCden_t_ha_0.20cm


nTree=1000

gf.veg<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+soil_cat+age_max+max_spec,data=bio_soc_cls,importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s20<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+aspect+slope+tpi+soil_cat+age_max+max_spec,data=bio_soc_cls,importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s1m<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+soil_cat+age_max+max_spec,data=bio_soc_cls,importance="permute.ensemble",ntree=nTree,seed=-111)

gf.totC<-rfsrc(total_C~lat+long+aspect+slope+tpi+soil_cat+age_max+max_spec,data=bio_soc_cls,importance="permute.ensemble",ntree=nTree,seed=-111)
gf.totCs<-rfsrc(total_Cs~lat+long+aspect+slope+tpi+soil_cat+age_max+max_spec,data=bio_soc_cls,importance="permute.ensemble",ntree=nTree,seed=-111)

srt_vimp<-function(x)
{
  ximp<-data.frame(xvars=x$xvar.names,vimp=as.numeric(x$importance));
  ximp<-ximp[order(ximp$vimp,decreasing=T),];
  return(ximp);
}

#Fig S1: variable importance

pdf('./figs/FiG_S1_vimps.pdf')
par(mfrow=c(2,3))
ximp.veg<-srt_vimp(gf.veg);dotchart(x=ximp.veg$vimp,labels = ximp.veg$xvar,pch=19);title('Carbon in tree')
ximp.s1m<-srt_vimp(gf.s1m);dotchart(x=ximp.s1m$vimp,labels = ximp.s1m$xvar,pch=19);title('Carbon in top 1 meter soil')
ximp.s20<-srt_vimp(gf.s20);dotchart(x=ximp.s20$vimp,labels = ximp.s20$xvar,pch=19);title("Carbon in top 20cm soil")
ximp.TOC<-srt_vimp(gf.totC);dotchart(x=ximp.TOC$vimp,labels = ximp.TOC$xvar,pch=19);title("Carbon in tree & top 1 meter soil")
ximp.TOCs<-srt_vimp(gf.totCs);dotchart(x=ximp.TOCs$vimp,labels = ximp.TOCs$xvar,pch=19);title("Carbon in tree & top 20cm soil")
dev.off()

win.graph();plot(gf.veg);title(main='gf.veg')
win.graph();plot(gf.s20);title(main='gf.s20')
win.graph();plot(gf.s1m);title(main='gf.s1m')
win.graph();plot(gf.totC);title(main='gf.totC')
win.graph();plot(gf.totCs);title(main='gf.totCs')

gf.veg.v2<-var.select(gf.veg,method='md')
gf.s20.v2<-var.select(gf.s20,method='md')
gf.s1m.v2<-var.select(gf.s1m,method='md')
gf.totC.v2<-var.select(gf.totC,method='md')
gf.totCs.v2<-var.select(gf.totCs,method='md')


win.graph();gf.veg.p<-plot.variable(gf.veg,partial=TRUE);title(main='gf.veg')
win.graph();gf.s20.p<-plot.variable(gf.s20,partial=TRUE);title(main='gf.s20')
win.graph();gf.s1m.p<-plot.variable(gf.s1m,partial=TRUE);title(main='gf.s1m')
win.graph();gf.totC.p<-plot.variable(gf.totC,partial=TRUE);title(main='gf.totC')
win.graph();gf.totCs.p<-plot.variable(gf.totCs,partial=TRUE);title(main='gf.totCs')

save.image('TOC_age-spe-lat-long.RData')

#to Check the influence factors of variables threshold of 100years
ii_l100<-which(bio_soc_cls$age_max<=100)
ii_g100<-which(bio_soc_cls$age_max>=100)

gf.veg.l100<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_l100,],importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s20.l100<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_l100,],importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s1m.l100<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_l100,],importance="permute.ensemble",ntree=nTree,seed=-111)

gf.totC.l100<-rfsrc(total_C~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_l100,],importance="permute.ensemble",ntree=nTree,seed=-111)
gf.totCs.l100<-rfsrc(total_Cs~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_l100,],importance="permute.ensemble",ntree=nTree,seed=-111)

gf.veg.g100<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_g100,],importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s20.g100<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_g100,],importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s1m.g100<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_g100,],importance="permute.ensemble",ntree=nTree,seed=-111)

gf.totC.g100<-rfsrc(total_C~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_g100,],importance="permute.ensemble",ntree=nTree,seed=-111)
gf.totCs.g100<-rfsrc(total_Cs~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_cls[ii_g100,],importance="permute.ensemble",ntree=nTree,seed=-111)


gf.veg.v2.l100<-var.select(gf.veg.l100,method='md')
gf.s20.v2.l100<-var.select(gf.s20.l100,method='md')
gf.s1m.v2.l100<-var.select(gf.s1m.l100,method='md')
gf.totC.v2.l100<-var.select(gf.totC.l100,method='md')
gf.totCs.v2.l100<-var.select(gf.totCs.l100,method='md')


gf.veg.v2.g100<-var.select(gf.veg.g100,method='md')
gf.s20.v2.g100<-var.select(gf.s20.g100,method='md')
gf.s1m.v2.g100<-var.select(gf.s1m.g100,method='md')
gf.totC.v2.g100<-var.select(gf.totC.g100,method='md')
gf.totCs.v2.g100<-var.select(gf.totCs.g100,method='md')

win.graph();gf.veg.p<-plot.variable(gf.veg.g100,partial=TRUE);title(main='gf.veg')
win.graph();gf.s20.p<-plot.variable(gf.s20.g100,partial=TRUE);title(main='gf.s20')
win.graph();gf.s1m.p<-plot.variable(gf.s1m.g100,partial=TRUE);title(main='gf.s1m')



save.image('TOC_age-spe-lat-long.RData')

#To Extract the development of differnt sections


#Figure S2: importance of lat,long and max_spec using TOC

load("veg_soil_rf_envFactors.RData")

require(pracma)

#par(mfrow=c(2,2))
pdf('./figs/FIG_S2_partial_response.pdf')
layout(matrix(c(1,2,3,3),2,2, byrow = FALSE))

par(mar=c(4,4,1,1))
x=gf.totC.p$pData[[2]]$x.uniq
y=gf.totC.p$pData[[2]]$yhat
yerr=gf.totC.p$pData[[2]]$yhat.se
errorbar(x,y,yerr = yerr,xlab='Longitude',cex=1.1,ylab='Partial Variable Importance',pch=19,type='b');
text(118,257,'(A)')

par(mar=c(4,4,1,1))
errorbar(x=gf.totC.p$pData[[3]]$x.uniq,y=gf.totC.p$pData[[3]]$yhat,yerr=gf.totC.p$pData[[3]]$yhat.se,xlab='Latitude',ylab='VIMP',type='b',cex=1.1,pch=19);
text(41.5,257,'(B)')

par(mar=c(4,6,1,1))
ts.m<-data.frame(tsp=apply(matrix(gf.totC.p$pData[[5]]$yhat,ncol=length(levels(gf.totC.p$pData[[5]]$x))),2,FUN=mean))
ts.sd<-data.frame(tsp=apply(matrix(gf.totC.p$pData[[5]]$yhat,ncol=length(levels(gf.totC.p$pData[[5]]$x))),2,FUN=sd))
errorbar(y=1:43,x=ts.m$tsp,xerr=ts.sd$tsp,axes=F,xlab='Partial Variable Importance', ylab='',pch=19,cex=1.1,type='b')
axis(2,at=1:43,labels=as.character(gf.totC.p$pData[[5]]$x.uniq),las=2,cex.axis=0.8)
axis(1)
box()
text(260,43,'(C)')

dev.off()

require(Hmisc)
load('TOC_age-spe-lat-long.RData')
# Figure S3: Importance of AGE using Cbiomass, CSoil1m, Csoil20cm

pdf('./figs/FIG_S3_partial_response_age.pdf')
par(mar=c(4,4,1,1))
x.veg=gf.veg.p$pData[[1]]$x.uniq
y.veg=gf.veg.p$pData[[1]]$yhat
yerr.veg=gf.veg.p$pData[[1]]$yhat.se

x.s1m=gf.s1m.p$pData[[2]]$x.uniq
y.s1m=gf.s1m.p$pData[[2]]$yhat
yerr.s1m=gf.s1m.p$pData[[2]]$yhat.se

x.s20=gf.s20.p$pData[[3]]$x.uniq
y.s20=gf.s20.p$pData[[3]]$yhat
yerr.s20=gf.s20.p$pData[[3]]$yhat.se
ymax=floor(max(c(y.veg,y.s1m,y.s20))/5+1)*5
ymin=floor(min(c(y.veg,y.s1m,y.s20))/5)*5

errbar(x.s1m,y.s1m,yplus = y.s1m+yerr.s1m,yminus=y.s1m-yerr.s1m,cap=0,type='b',pch=19,col=1,xlim=c(0,300),ylim=c(ymin,ymax),xlab='Plot Age(year)',cex=1.1,ylab='Carbon density(tC/ha)');
errbar(x.s20,y.s20,yplus = y.s20+yerr.s20,yminus=y.s20-yerr.s20,cap=0,type='b',pch=20,col=2,xlim=c(0,300),add=T);
errbar(x.veg,y.veg,yplus = y.veg+yerr.veg,yminus=y.veg-yerr.veg,cap=0,type='b',pch=21,col=3,xlim=c(0,300),add=T);
abline(v=c(90,110,140),lwd=1,lty=2,col=1:3)
grid()
legend('bottomright',legend=c("Carbon in 1m Soil","Carbon in 20cm Soil","Biomass Carbon"),pch=19:21,col=1:3)
dev.off()

#Figure S4: Scatterplots between Cbiomass/Csoil1m/Csoil20cm and Age
Sys.setenv('R_GSCMD'="C://Program Files//gs//gs9.16//bin//gswin64c.exe");

#pdf(('./figs/FIG_S4_scatter_age.pdf')
tiff('./figs/FIG_S4_scatter_age.tif',width=2100,height=2100,res=300,compression = "lzw")

plot(SoilC_1m_t_ha~age_max, data = bio_soc_cls, type='p', pch=1,cex=0.85, col=1, xlab='Plot age (year)', ylab= 'Carbon density (tC ha-1)',xlim=c(0,300),ylim=c(0,500))
points(W_total_tC_ha~age_max,data=bio_soc_cls, type='p',pch=2,cex=0.85,col=3)
points(SoilCden_t_ha_0.20cm~age_max,data=bio_soc_cls,type='p',pch=3,cex=0.85,col=2)
legend('topright',legend=c("Csoil1m","Csoil20cm","Cbiomass"),pch=c(1,2,3),col=c(1,2,3))

dev.off()


sumfun <- function(x, ...){
  c(m=mean(x, ...),v=sd(x, ...), l=length(x))
}
rm(dd)
cut(bio_soc_cls$age_max,breaks=quantile(bio_soc_cls$age_max,seq(0,1,0.025)))->cc
#cut(bio_soc_cls$age_max,25)->cc
data.frame(bio_soc_cls,agecls1=cc)->dd
summaryBy(age_max+SoilC_1m_t_ha+W_total_tC_ha+SoilCden_t_ha_0.20cm~agecls1,data=dd,FUN=sumfun,na.rm=T)


#Fig S5 scatter age
read.csv('Carbon_age_long.csv',colClasses=c('character','numeric','numeric',"character","character"))->cage.long

pdf('./figs/FIG_S5_scatter_age.pdf',paper='A4')
xyplot(Carbon_t_ha~age_max|cls_plot,groups =Cls_C,type=c('p','smooth','g'),data=cage.long, prepanel = function(x, y) prepanel.loess(x, y, span = 0.8),cex=0.5,col=1:3,lwd=1.5,lty=1:3,pch=1:3,
      key=list(lines=T,points = TRUE,col=1:3,lwd=1.5,lty=1:3,pch=1:3,cex=1,text=list(levels(as.factor(cage.long$Cls_C))[c(1,2,3)]),corner=c(0.5,0.975)),
      xlab="Plot Age(Year)" ,ylab='Carbon Density (tC/ha)')
dev.off()


