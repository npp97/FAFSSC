lfun <- function(x, ...) {
  c(l = length(x))
}


require(Hmisc)
require(doBy)

read.csv("D://东北//data//碳//Biomass_Cproj.csv")->bio.cproj
read.csv("D://东北//data//碳//Biomass_zzzy.csv")->bio.zzzy


bio.zzzy<-bio.zzzy[,c("ydID","zwmcbh","xj","Wfoliage_kg","Wstem_kg","Wroots_kg","pft")]
bio.cproj<-bio.cproj[,c("ydID","zwmcbh","xj","Wfoliage_kg","Wstem_kg","Wroots_kg","pft")]
bio.tot<-rbind(bio.zzzy,bio.cproj)

bio.tot$xjbrk<-cut(bio.tot$xj,breaks=c(0,10,20,30,40,50,60,70,100,225))

write.csv(bio.tot,'biomass_tot.csv',row.names=FALSE)

plotlst<-levels(bio.tot$ydID);
nplots<-length(plotlst)
# 
# odnf<-bio.tot$pft %in% c("DNF")
# odbf<-bio.tot$pft %in% c("DBF")
# odbf<-bio.tot$pft %in% c("ENF")



l<- summaryBy(xj ~ ydID+xjbrk, data=bio.tot, FUN=lfun, na.omit=TRUE)
mcarbon<- summaryBy(Wfoliage_kg+Wstem_kg+Wroots_kg ~ ydID + xjbrk, data=bio.tot, FUN=mean, na.omit=T)

#pft<-c("DBF","ENF","DNF")

loss.carbon<-data.frame(matrix(NA,ncol=4,nrow=I(3*nplots)))
names(loss.carbon)<-c("plot","Wfoliage_kg","Wstem_kg","Wroots_kg")
for (iplot in 1:nplots){
    loss.carbon$plot[iplot]<-plotlst[iplot]
    ii_l<-which(l$ydID == plotlst[iplot])
    if(length(l$xj.l[ii_l])< 2){loss.carbon[iplot,2:4]<-0}else
    {
      loss<-(-diff(l$xj.l[ii_l]))
      loss[length(loss)+1]<-loss[length(loss)]
      loss[loss<0]<-0
      loss.carbon[iplot,2:4]<-apply((mcarbon[(mcarbon$ydID %in% plotlst[iplot]),3:5])*loss,2,sum)      
    }
}
#-----------------------------

read.csv('Book8.csv')->b8

rfsrc(LC_tot_tC_ha_yr~age+lat+long,data=b8,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lc
plot.variable(gfr.lc,"age",partial=T,npts=50)->gr.lc

lc.rf<-data.frame(age=gr.lc$pData[[1]]$x.uniq,lc=gr.lc$pData[[1]]$yhat,lc.se=gr.lc$pData[[1]]$yhat.se)
ylim<-floor(range(lc.rf$age))
lc.a<-approx(lc.rf$age,lc.rf$lc,ylim[1]:272)



#-------------------------------------------------
read.csv('afb.csv')->afb

rfsrc(dead_froot_tC_ha~age_max+lat+long,data=afb,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.fr
plot.variable(gfr.fr,"age_max",partial=T,npts=50)->gr.fr

fr.rf<-data.frame(age=gr.fr$pData[[1]]$x.uniq,lc=gr.fr$pData[[1]]$yhat,lc.se=gr.fr$pData[[1]]$yhat.se)
ylim<-floor(range(fr.rf$age))
fr.a<-approx(fr.rf$age,fr.rf$lc,ylim[1]:272)


