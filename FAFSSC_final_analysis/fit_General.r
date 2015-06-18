 read.csv('bio_soc_env_all.csv.csv')->all
 randomForest(W_total_tC_ha~age_max+lat+long,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T)->gf.b
 randomForest(SoilC_1m_t_ha~age_max+lat+long,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T)->gf.s
 randomForest(TotC_eco~age_max+lat+long,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T)->gf.tot
 partialPlot(gf.tot,all,age_max,col=3,ylim=c(30,350))->g.t
 partialPlot(gf.b,all,age_max,col=3,add=T)->g.b
 partialPlot(gf.s,all,age_max,col=2,add=T)->g.s
 grid(nx=28,ny=30)
 cbind(g.t$x,g.t$y,g.b$y,g.s$y)
 
read.csv('bio_soc_env_all.csv.csv')->all
rfsrc(W_total_tC_ha~age_max+lat+long,data=all,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gf.b
rfsrc(SoilC_1m_t_ha~age_max+lat+long,data=all,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gf.s
rfsrc(TotC_eco~age_max+lat+long,data=all,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gf.t

 plot.variable(gf.t,"age_max",partial=T, npts=50)->g.t
 win.graph()
 plot.variable(gf.b,"age_max",partial=T, npts=50)->g.b
 win.graph()
 plot.variable(gf.s,"age_max",partial=T, npts=50)->g.s
 
 cbind(g.t$pData[[1]]$x.uniq,g.t$pData[[1]]$yhat,g.b$pData[[1]]$yhat,g.s$pData[[1]]$yhat)
  
require(struchange)

 read.csv('bio_soc_env_all.csv.csv')->all
 
randomForest(TVDI.evi.m~age_max+lat+long+slope+aspect+tpi+twi+max_spec,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T)->gf.tvdi
importance(gf.tvdi)
 partialPlot(gf.tvdi,all,age_max,col=3)->g.tvdi


all$ts<-all$ta_C*all$flai_mean


randomForest(ts~age_max+lat+long+slope+aspect+tpi+twi+max_spec,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T)->gf.ts
importance(gf.ts)
 partialPlot(gf.ts,all,age_max,col=3)->g.ts

 cbind(g.ts$x,g.ts$y,g.tvdi$x,g.tvdi$y)
 
 
rfsrc(TVDI.evi.m~age_max+lat+long+slope+aspect+tpi+twi+max_spec,data=all[!is.na(all$SoilC_1m_t_ha),],importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.tvdi
plot.variable(gfr.tvdi,"age_max",partial=T,npts=50)->gr.tdvi

rfsrc(ts~age_max+lat+long+slope+aspect+tpi+twi+max_spec,data=all[!is.na(all$SoilC_1m_t_ha),],importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.ts
plot.variable(gfr.ts,"age_max",partial=T,npts=50)->gr.ts

cbind(gr.ts$pData[[1]]$x.uniq, gr.ts$pData[[1]]$yhat, gr.ts$pData[[1]]$yhat.se, gr.tdvi$pData[[1]]$yhat,gr.tdvi$pData[[1]]$yhat.se)


all$litter_lf<-1*all$W_leaf_kg_tC_ha

randomForest(litter_lf~age_max+lat+long+max_spec,data=all,importance=T)->gf.lf
partialPlot(gf.lf,all,age_max,add=T)->cclf



rfsrc(litter_lf~age_max+lat+long+max_spec,data=all,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gfr.lf
plot.variable(gfr.lf,"age_max",partial=T,npts=50)->gr.lf

data.frame(age=gr.lf$pData[[1]]$x.uniq,lf=gr.lf$pData[[1]]$yhat,lf.se=gr.lf$pData[[1]]$yhat.se,lr=gr.lroot$pData[[1]]$yhat,lr.se=gr.lroot$pData[[1]]$yhat.se)


#---------------------------------Dead fine root input
read.csv('afb.csv')->fine_root_c
randomForest(dead_froot_tC_ha~age_max+lat+long,data=fine_root_c,importance=T)->gf.fr
partialPlot(gf.fr,fine_root_c,age_max,add=T)->cc

rfsrc(dead_froot_tC_ha~age_max+lat+long,data=fine_root_c,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gf.fr
plot.variable(gf.fr,"age_max",partial=T,npts=50)->gr.fr
data.frame(age_fr=gr.fr$x.uniq,fr=gr.fr$yhat)



rfsrc(dead_froot_tC_ha~age_max+lat+long,data=fine_root_c,importance="permute.ensemble",na.action = "na.impute",nimpute = 10, seed= -111,ntree=200)->gf.fr

plot.variable(gf.fr,"age_max",partial=T,npts=50)->gr.fr

data.frame(age_lf=gr.lf$pData[[1]]$x.uniq,lf=gr.lf$pData[[1]]$yhat,age_fr=gr.fr$x.uniq,fr=gr.fr$yhat)

#----------Cal liter decomp
read.csv('deco_lr.csv')->llf_lr
approx(x=llf_lr$age1,y=llf_lr$lf,11:272)->lf
approx(x=llf_lr$age1,y=llf_lr$lr,11:272)->lr
approx(x=llf_lr$age1,y=llf_lr$deco_lf_k,11:272)->lf_k
approx(x=llf_lr$age1,y=llf_lr$deco_lroot_k,11:272)->lr_k


x=1:10
yl<-length(x)
xl=length(lf$x)

# matrix(0,ncol=(length(lf$x)+yl-1),nrow=length(lf$x))->dec_t

# for (i in 1:length(lf$x)){
	# d.l<-0.25*lf$y[i]*exp(-1*lf_k$y[i]*x)
	# d.r<-0.25*lr$y[i]*exp(-1*lr_k$y[i]*x)
	# k=(i-1)+x
	# dec_t[i,c(k)]=d.l+d.r
# }

	d.l<-0.25*lf$y*exp(-1*lf_k$y)
	d.r<-0.25*lr$y*exp(-1*lr_k$y)
	litter_input<-d.l+d.r

#--------------------------------------------------------------
read.csv('deco_lr.csv')->llf_lr
approx(x=llf_lr$age1,y=llf_lr$ts,11:272)->tsoil
approx(x=llf_lr$age1,y=llf_lr$sm,11:272)->sm
approx(x=llf_lr$age1,y=llf_lr$som,11:272)->som

q10=exp(10*0.204*(1-tsoil$y/36.9))
q10=0.54649+0.31165*log(som$y+3.14537)
q10=5.91118*exp(-tsoil$y/11.21539)+1.48429
fts1<-exp(log(2.5)/10*(tsoil$y-20))
fts2<-exp(log(q10)/10*(tsoil$y-20))
#fts2<-exp(-3.3764+0.204*tsoil$y*(1-0.5*tsoil$y/36.9))
fsm<-0.25+0.75*(sm$y/0.35)
fsm[]<-1
deco1=(1/15+1/750)*fts1*fsm*som$y
deco2=(1/15+1/750)*fts2*fsm*som$y

som_deco<-data.frame(age=sm$x,deco1,deco2)

win.graph()
par(mfrow=c(1,2))
plot(som_deco$age,(litter_input-som_deco[,3]))
plot(som_deco$age,cumsum(litter_input-som_deco[,3]))






