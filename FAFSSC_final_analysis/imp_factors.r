#require(randomForest)
require(randomForestSRC)
#require(ggRandomForests)
require(FactoMineR)
require(strucchange)
#require(gradientForest)
biomass<-read.csv('biomass_projects.csv')
soc<-read.csv('soc_projects.csv')
geoinfo<-read.csv('age_bkinf_projects_f1.csv')
soil_cls<-read.csv('soil_class_projects.csv')
dom_spec<-read.csv('dom_spec_projects.csv')

# > names(biomass)
 # [1] "pplot"              "projects"            "W_total_t_ha"       "W_leaf_kg_t_ha"     "W_foliage_kg_t_ha" 
 # [6] "W_stem_kg_t_ha"     "W_root_kg_t_ha"     "W_total_tC_ha"      "W_leaf_kg_tC_ha"    "W_foliage_kg_tC_ha"
# [11] "W_stem_kg_tC_ha"    "W_root_kg_tC_ha"   
# > names(dom_spec)
# [1] "pplot"        "projects"     "Tree_no_g0p1" "max.sum"      "max_spec"     "max"         
# > names(geoinfo)
 # [1] "pplot"    "projects" "lat"      "long"     "elev"     "slope"    "aspect"   "twi"      "tpi"      "p_mm"    
# [11] "ta_C"     "xj.m"     "xj.md"    "xj.l"     "age"     
# > names(soc)
# [1] "pplot"                "projects"              "SoilC_1m_t_ha.sum"    "SoilCden_t_ha_0.20cm"
# > names(soil_cls)
 # [1] "pplot"             "projects"           "lat"               "long"              "elev"             
 # [6] "soilcode"          "code_soil_gang"    "code_soil_cat"     "code_soil_sub_cat" "soil_gang"        
# [11] "soil_cat"          "sub_soil_cat"     

soil_cls<-soil_cls[,c(1,2,6:12)]
bkinfo<-merge(geoinfo,soil_cls)
bk_info<-merge(bkinfo,dom_spec)

bio_bk<-merge(bk_info,biomass)
soc_bk<-merge(bk_info,soc)
bio_soc_bk<-merge(bio_bk,soc)

# > names(bk_info)
 # [1] "pplot"             "projects"          "lat"               "long"              "elev"             
 # [6] "slope"             "aspect"            "twi"               "tpi"               "p_mm"             
# [11] "ta_C"              "xj.m"              "xj.md"             "xj.l"              "age"              
# [16] "soilcode"          "code_soil_gang"    "code_soil_cat"     "code_soil_sub_cat" "soil_gang"        
# [21] "soil_cat"          "sub_soil_cat"      "Tree_no_g0p1"      "max.sum"           "max_spec"         
# [26] "max"              
# > 

# > names(bio_bk)
 # [1] "pplot"              "projects"           "lat"                "long"               "elev"              
 # [6] "slope"              "aspect"             "twi"                "tpi"                "p_mm"              
# [11] "ta_C"               "xj.m"               "xj.md"              "xj.l"               "age"               
# [16] "Tree_no_g0p1"       "max.sum"            "max_spec"           "max"                "W_total_t_ha"      
# [21] "W_leaf_kg_t_ha"     "W_foliage_kg_t_ha"  "W_stem_kg_t_ha"     "W_root_kg_t_ha"     "W_total_tC_ha"     
# [26] "W_leaf_kg_tC_ha"    "W_foliage_kg_tC_ha" "W_stem_kg_tC_ha"    "W_root_kg_tC_ha"   
# > names(soc_bk)
 # [1] "pplot"                "projects"             "lat"                  "long"                
 # [5] "elev"                 "slope"                "aspect"               "twi"                 
 # [9] "tpi"                  "p_mm"                 "ta_C"                 "xj.m"                
# [13] "xj.md"                "xj.l"                 "age"                  "Tree_no_g0p1"        
# [17] "max.sum"              "max_spec"             "max"                  "SoilC_1m_t_ha"       
# [21] "SoilCden_t_ha_0.20cm"
# > names(bio_soc_bk)
 # [1] "pplot"                "projects"             "lat"                  "long"                
 # [5] "elev"                 "slope"                "aspect"               "twi"                 
 # [9] "tpi"                  "p_mm"                 "ta_C"                 "xj.m"                
# [13] "xj.md"                "xj.l"                 "age"                  "Tree_no_g0p1"        
# [17] "max.sum"              "max_spec"             "max"                  "W_total_t_ha"        
# [21] "W_leaf_kg_t_ha"       "W_foliage_kg_t_ha"    "W_stem_kg_t_ha"       "W_root_kg_t_ha"      
# [25] "W_total_tC_ha"        "W_leaf_kg_tC_ha"      "W_foliage_kg_tC_ha"   "W_stem_kg_tC_ha"     
# [29] "W_root_kg_tC_ha"      "SoilC_1m_t_ha"        "SoilCden_t_ha_0.20cm"
# > dim(bio_soc_bk)
# [1] 2106   31
# > dim(soc_bk)
# [1] 2112   21
# > dim(bio_bk)
# [1] 2824   29
# > 

bio_soc_bk$FLR=bio_soc_bk$W_foliage_kg_tC_ha+bio_soc_bk$W_leaf_kg_tC_ha+bio_soc_bk$W_root_kg_tC_ha
bio_soc_bk$total_C<-bio_soc_bk$W_total_tC_ha+bio_soc_bk$SoilC_1m_t_ha
bio_soc_bk$total_Cs<-bio_soc_bk$W_total_tC_ha+bio_soc_bk$SoilCden_t_ha_0.20cm

pp<-PCA(bio_soc_cls[,c("lat","long","aspect","slope","tpi","p_mm","ta_C","age_max")])

ii<-which((bio_soc_bk$lat>=40)&(bio_soc_bk$lat<51.5))

bio_soc_bk=bio_soc_bk[ii,]

nTree=200

# gf.totC<-rfsrc(total_C~lat+long+elev+aspect+slope+tpi+p_mm+ta_C+age_max+soil_gang+max_spec+FLR+xj.m,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
# gf.totCs<-rfsrc(total_Cs~lat+long+elev+aspect+slope+tpi+p_mm+ta_C+age_max+soil_gang+max_spec+FLR+xj.m,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
# 
# #gf<-rfsrc(W_total_tC_ha~lat+long+elev+aspect+slope+twi+tpi+p_mm+ta_C+age+soil_cat+soil_gang+Tree_no_g0p1+max_spec,data=bio_bk,importance="permute.ensemble",seed=-111)
# gf.veg<-rfsrc(W_total_tC_ha~lat+long+elev+aspect+slope+tpi+p_mm+ta_C+age_max+soil_cat+soil_gang+max_spec+xj.m,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
# gf.s20<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+elev+aspect+slope+tpi+p_mm+ta_C+age_max+soil_gang+max_spec+W_total_tC_ha+FLR+xj.m,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
# gf.s1m<-rfsrc(SoilC_1m_t_ha~lat+long+elev+aspect+slope+tpi+p_mm+ta_C+age_max+soil_gang+max_spec+W_total_tC_ha+FLR+xj.m,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)

gf.totC<-rfsrc(total_C~lat+long+aspect+slope+tpi+p_mm+ta_C+soil_gang+age_max+max_spec,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
gf.totCs<-rfsrc(total_Cs~lat+long+aspect+slope+tpi+p_mm+ta_C+soil_gang+age_max+max_spec,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)

#gf<-rfsrc(W_total_tC_ha~lat+long+elev+aspect+slope+twi+tpi+p_mm+ta_C+age+soil_cat+soil_gang+Tree_no_g0p1+max_spec,data=bio_bk,importance="permute.ensemble",seed=-111)
gf.veg<-rfsrc(W_total_tC_ha~lat+long+aspect+slope+tpi+soil_cat+soil_gang+age_max+max_spec,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s20<-rfsrc(SoilCden_t_ha_0.20cm~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)
gf.s1m<-rfsrc(SoilC_1m_t_ha~lat+long+aspect+slope+tpi+soil_gang+age_max+max_spec,data=bio_soc_bk,importance="permute.ensemble",ntree=nTree,seed=-111)

gf.veg.v<-var.select(gf.veg,method='vh.vimp')
gf.s20.v<-var.select(gf.s20,method='vh.vimp')
gf.s1m.v<-var.select(gf.s1m,method='vh.vimp')
gf.totC.v<-var.select(gf.totC,method='vh.vimp')
gf.totCs.v<-var.select(gf.totCs,method='vh.vimp')

gf.veg.v2<-var.select(gf.veg,method='md')
gf.s20.v2<-var.select(gf.s20,method='md')
gf.s1m.v2<-var.select(gf.s1m,method='md')
gf.totC.v2<-var.select(gf.totC,method='md')
gf.totCs.v2<-var.select(gf.totCs,method='md')

require(varSelRF)


win.graph();plot(gf.veg);title(main='gf.veg')
win.graph();plot(gf.s20);title(main='gf.s20')
win.graph();plot(gf.s1m);title(main='gf.s1m')
win.graph();plot(gf.totC);title(main='gf.totC')
win.graph();plot(gf.totCs);title(main='gf.totCs')

win.graph();gf.veg.p<-plot.variable(gf.veg,partial=TRUE);title(main='gf.veg')
win.graph();gf.s20.p<-plot.variable(gf.s20,partial=TRUE);title(main='gf.s20')
win.graph();gf.s1m.p<-plot.variable(gf.s1m,partial=TRUE);title(main='gf.s1m')
win.graph();gf.totC.p<-plot.variable(gf.totC,partial=TRUE);title(main='gf.totC')
win.graph();gf.totCs.p<-plot.variable(gf.totCs,partial=TRUE);title(main='gf.totCs')

win.graph();plot.variable(gf.totCs,c("long","lat","max_spec"),partial=TRUE)


plot(bio_soc_bk[,c(3:5,10:16,25,28,38:40)],pch='.') # to check the influence of lat
win.graph(); plot(bio_soc_bk[,c(3:4,6:9,25,28,38:40)],pch='.') # to check the influence terrian parameters

x.long<-(gf.totCs.p$pData[[1]]$x.uniq)
y.long<-(gf.totCs.p$pData[[1]]$yhat)
yse.long<-gf.totCs.p$pData[[1]]$yhat.se

x.lat<-(gf.totCs.p$pData[[4]]$x.uniq)
y.lat<-(gf.totCs.p$pData[[4]]$yhat)
yse.lat<-gf.totCs.p$pData[[4]]$yhat.se

x.max_spec<-(gf.totCs.p$pData[[3]]$x.uniq)
y.max_spec<-apply(matrix(gf.totCs.p$pData[[3]]$yhat,ncol=43),2,FUN=mean)
yse.max_spec<-apply(matrix(gf.totCs.p$pData[[3]]$yhat,ncol=43),2,FUN=sd)
yq2575.max_spec<-apply(matrix(gf.totCs.p$pData[[3]]$yhat,ncol=43),2,FUN=quantile)

save.image('veg_soil_rf_envFactors.RData')

tree.species<-data.frame(tsp=apply(matrix(gf.totCs.p$pData[[3]]$yhat,ncol=length(levels(gf.totCs.p$pData[[3]]$x))),2,FUN=mean))
rownames(tree.species)<-(as.character(gf.totCs.p$pData[[3]]$x.uniq))
lk<-length(levels(gf.totCs.p$pData[[3]]$x))-1
asw<-matrix(0,nrow=lk-1,ncol=1)

require(cluster)
for (i in 2:lk){
	pp<-pam(tree.species,i,diss=FALSE)
	asw[i]<-pp$silinfo$avg.width
}
k.best<-which.max(asw)+1

#K best is 3 or 2 kinds of 
ql<-function(x){
  return(list(qq=(quantile(x)),length(x)))
}
ql<-function(x){
  return(list(min(x),max(x),median(x),length(x)))
}

require(doBy)
read.csv('bio_soc_cls.csv')->bio_soc_cls
# o<-((bio_soc_bk$aspect>45)&(bio_soc_bk$aspect<315))
# ii<-which((bio_soc_cls$p_mm>450)&o)
# bio_soc_cls<-bio_soc_cls[ii,]
summaryBy(age_max~cls_max_spec3,data=bio_soc_cls,FUN=ql)
#summaryBy(age_max~cls_max_spec4,data=bio_soc_cls,FUN=ql)

grp<-as.character(levels(bio_soc_cls$cls_max_spec4))
lgrp=length(grp)

win.graph()
par(mfrow=c(4,5))

for (i in 1:20){
  ii<-which(bio_soc_cls$cls_max_spec4 %in% grp[i])
  #plot(bio_soc_cls$SoilC_1m_t_ha[ii]~bio_soc_cls$age_max[ii],xlable='age',ylable='soc_1m')
  #plot(bio_soc_cls$SoilCden_t_ha_0.20cm[ii]~bio_soc_cls$age_max[ii],xlab='age',ylab='soc_20cm')
  plot(bio_soc_cls$W_total_tC_ha[ii]~bio_soc_cls$age_max[ii],xlab='age',ylab='vegC')
  title(grp[i])
}


