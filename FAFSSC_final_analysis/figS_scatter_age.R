# TODO: Add comment
# 
# Author: junhui.zhang
###############################################################################


require(lattice)
require(reshape)
Sys.setenv('R_GSCMD'="C://Program Files//gs//gs9.16//bin//gswin64c.exe");
read.csv('Carbon_ZONE.csv',colClasses=c('character',rep('numeric',4)))->cage.w
cage.l<-reshape(cage.w,direction='long',v.names=c("Carbon_t_ha"),varying=c("TotC_eco","W_total_tC_ha","SoilC_1m_t_ha"),timevar="Cls_C",time=c("Ecosystem_tC_ha","Biomass_tC_ha","SoilC_tC_ha"))
cage.l<-cage.l[,c("ZONE","Cls_C","age_max","Carbon_t_ha")]

tiff('./FIG_S5_scatter_age.tif',width=7,height=7,res=300,units="in",compression="lzw")
xyplot(Carbon_t_ha~age_max|ZONE,groups =Cls_C,type=c('p','smooth','g'),data=cage.l, prepanel = function(x, y) prepanel.loess(x, y, span = 0.5),cex=0.5,col=1:3,lwd=1.25,lty=1,pch=c(19,2,3),
		key=list(lines=T,points = TRUE,col=1:3,lty=1,lwd=1.5,pch=c(19,2:3),cex=0.9,text=list(levels(as.factor(cage.l$Cls_C))[c(1,2,3)]),corner=c(1,0.975)),
		xlab="Age(Years)" ,ylab=expression('Carbon stock (tC ha'^{-1}*')'))
dev.off()

