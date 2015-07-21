read.csv('pppa.csv')->pppa
require(randomForestSRC)
gf.veg<-rfsrc(LC_tot_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=pppa,importance="permute.ensemble",nTree=200, seed=-111)
gf.veg.p<-plot.variable(gf.veg,partial=TRUE)
plot(gf.veg)
           # Importance   Relative Imp
# max_spec     119.3614         1.0000
# lat           55.3771         0.4639
# long          48.4912         0.4063
# p_mm          42.8308         0.3588
# age_max       36.3392         0.3044
# slope         11.1841         0.0937
# twi            8.3718         0.0701
# ta_C           3.1284         0.0262
# aspect         0.9350         0.0078
# tpi            0.5444         0.0046

gf.veg<-rfsrc(LC_tot_tC_ha~lat+long+aspect+slope+tpi+age_max++twi+p_mm+ta_C,data=pppa,importance="permute.ensemble",nTree=200, seed=-111)
 read.csv('cal_loss_carbon.csv')->closs
 read.csv("bio_soc_env_all.csv.csv")->all
 names(closs)[1]<-'pplot'
 merge(all,closs)->acc
gf.veg<-rfsrc(TVDI.ndvi.m~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=all,importance="permute.ensemble",nTree=200, seed=-111)
plot(gf.veg)
           # Importance   Relative Imp
# ta_C           0.0029         1.0000
# p_mm           0.0010         0.3608
# age_max        0.0008         0.2764
# twi            0.0006         0.2129
# long           0.0004         0.1452
# lat            0.0003         0.1156
# max_spec       0.0001         0.0443
# slope          0.0001         0.0236
# aspect         0.0000         0.0059
# tpi            0.0000         0.0002
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(TVDI.evi.m~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=all,importance="permute.ensemble",nTree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(flai_mean~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=all,importance="permute.ensemble",nTree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(flai_mean~lat+long+aspect+slope+tpi+age_max+max_spec,data=all,importance="permute.ensemble",nTree=200, seed=-111)
plot.variable(gf.veg,partial=T)
plot(gf.veg)

           # Importance   Relative Imp
# long            3e-04         1.0000
# age_max         2e-04         0.7827
# max_spec        2e-04         0.5561
# lat             2e-04         0.5342
# slope           1e-04         0.1793
# aspect          0e+00         0.0545
# tpi             0e+00         0.0169
gf.veg<-rfsrc(flai_mean~lat+long+age_max+max_spec,data=all,importance="permute.ensemble",nTree=200, seed=-111)
gf.veg<-rfsrc(flai_mean~lat+long+age_max+max_spec,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
all$ts_C<-all$ta_C*all$flai_mean
gf.veg<-rfsrc(ts_C~lat+long+age_max+max_spec,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
plot(gf.veg)

           # Importance   Relative Imp
# lat            5.0191         1.0000
# age_max        1.8001         0.3586
# long           1.4363         0.2862
# max_spec       0.4283         0.0853
gf.veg<-rfsrc(ts_C~lat+long+age_max,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(lai_mean~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot(gf.veg)

           # Importance   Relative Imp
# twi            0.0621         1.0000
# p_mm           0.0480         0.7732
# age_max        0.0339         0.5467
# ta_C           0.0299         0.4814
# long           0.0287         0.4617
# lat            0.0234         0.3769
# max_spec       0.0190         0.3067
# slope          0.0030         0.0486
# aspect         0.0020         0.0324
# tpi           -0.0001        -0.0009
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(lai0.95~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(lai0.95~lat+long+aspect+slope+tpi+age_max+max_spec+twi+TVDI.evi.m,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot(gf.veg)

             # Importance   Relative Imp
# TVDI.evi.m       0.8704         1.0000
# long             0.2254         0.2590
# lat              0.1794         0.2061
# twi              0.1481         0.1702
# max_spec         0.1214         0.1395
# slope            0.0347         0.0399
# age_max          0.0299         0.0343
# aspect           0.0258         0.0296
# tpi              0.0099         0.0114
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(lai0.95~lat+long+age_max+twi+TVDI.evi.m+ta_C,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(W_leaf_kg_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi+p_mm+ta_C,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(W_leaf_kg_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(W_leaf_kg_tC_ha~lat+long+age_max+max_spec,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
gf.veg<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance="permute.ensemble",ntree=200, seed=-111)
plot.variable(gf.veg,partial=T)
require(randomForest)
gf.veg<-randomForest(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance=T,ntree=200, seed=-111)
gf.veg<-randomForest(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance=T,ntree=200,nPerm=10)
set.seed(71)
gf.veg<-randomForest(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance=T,ntree=200,nPerm=10)
partialPlot(gf.veg,all,x.var=鈥渁ge_max鈥�)
partialPlot(gf.veg,all,x.var=age_max)
gf.veg<-randomForest(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T,ntree=200,nPerm=10)
partialPlot(gf.veg,all,x.var=age_max)
gf.veg<-randomForest(TotC_eco~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T,ntree=200,nPerm=10)
partialPlot(gf.veg,all,x.var=age_max)
gf.veg<-randomForest(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance=T,ntree=200,nPerm=10)
gf.veg<-randomForest(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T,ntree=200,nPerm=10)
partialPlot(gf.veg,all,x.var=age_max)
win.graph()
gf.veg<-randomForest(TotC_eco~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all[!is.na(all$SoilC_1m_t_ha),],importance=T,ntree=200,nPerm=10)
partialPlot(gf.veg,all,x.var=age_max)
grid()
win.graph()
gf.veg<-randomForest(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+twi,data=all,importance=T,ntree=200,nPerm=10)
partialPlot(gf.veg,all,x.var=age_max)
grid()

#------------------------------------
read.csv('pppa.csv')->pppa
require(randomForestSRC)

win.graph(height=6,width=6)
par(mfrow=c(2,2))
hist(pppa$W_total_tC_ha,20,freq=FALSE,xlab=expression("Biomass Carbon (MgC ha"^{-1}*")"),main='')
hist(pppa$SoilC_1m_t_ha,20,freq=FALSE,xlab=expression("Soil Carbon (MgC ha"^{-1}*")"),main='')
hist(pppa$TotC_eco,20,freq=FALSE,xlab=expression("Total Carbon (MgC ha"^{-1}*")"),main='')
hist(pppa$age_max,20,freq=FALSE,xlab=expression("Plot age (Year)"),main='')

gf.veg<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+age_max+max_spec+soil_gang,data=pppa,importance="permute.ensemble",nTree=200, seed=-111)
#gf.veg.p<-plot.variable(gf.veg,partial=TRUE)
plot(gf.veg)
gf.s1m<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+age_max+max_spec+soil_gang,data=pppa,importance="permute.ensemble",nTree=200, seed=-111)
plot(gf.s1m)
gf.totc<-rfsrc(TotC_eco~lat+long+aspect+slope+tpi+age_max+max_spec+soil_gang,data=pppa,importance="permute.ensemble",nTree=200, seed=-111)
plot(gf.totc)