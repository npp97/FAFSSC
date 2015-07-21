# TODO: (1) to merge monthly T&P, LAI and TVDI into sites
#       (2) to generate monthly input of Leaf litter, fr litter and cwd input
# 
# Author: junhui.zhang
###############################################################################

#working dir: "D:\东北\data\基础信息\monthly_inputs"

# required paskages
require(doBy)
require(randomForestSRC)
require(randomForest)

dir<-"D://东北//data//基础信息//monthly_inputs"
setwd(dir)

#[1] "bio_soc_env_all.csv"     "LaiBymth.csv"            "morlity_carbon_loss.csv" "summLST_VI_TVDI.csv"     "TPavgmth.csv"           

#1 TPavgmth.csv
read.csv('TPavgmth.csv')->tp
names(tp)
# "site" "lat"  "long" "Ta1"  "Ta2"  "Ta3"  "Ta4"  "Ta5"  "Ta6"  "Ta7"  "Ta8"  "Ta9"  "Ta10" "Ta11" "Ta12" "P1"   "P2"   "P3"   "P4"   "P5"   "P6"   "P7"   "P8"   "P9"   "P10"  "P11"  "P12" 
read.csv('bio_soc_env_all.csv')->all

site<-all[,c('pplot','lat','long','age_max','aspect','slope','elev','tpi','twi','soil_gang','soil_cat','sub_soil_cat','max_spec')]
names(site)[1]<-'site'
merge(site,tp,by='site')->tp.site
write.csv(tp.site,'tp.site.csv',row.names=F)

tp.site<-read.csv('tp.site.csv')

#TVDI.evi
read.csv('summLST_VI_TVDI.csv')->tvdi
names(tvdi)[1]<-'site'
nm<-c("site",paste("TVDI.evi.m",1:12,sep='.'),paste("TVDI.ndvi.m",1:12,sep='.'))

merge(site,tvdi[,nm],by='site')->tvdi.site
write.csv(tvdi.site,'tvdi.site.csv',row.names=F)

#LAI
read.csv('LaiBymth.csv')->lai
names(lai)[1]<-'site'
merge(site,lai,by='site')->lai.site
lai.site<-lai.site[as.character(unique(lai.site$site)),]
write.csv(lai.site,'lai.site.csv',row.names=F)

#bioC

read.csv('morlity_carbon_loss.csv')->mcl
mcl<-mcl[,c("pplot","lat","long","plot_area_m2")]
all.bio<-all[,c("pplot","lat","long","cls_max_spec4","age_max","xj.l","TotC_eco","SoilC_1m_t_ha","SoilCden_t_ha_0.20cm","W_total_tC_ha","W_leaf_kg_tC_ha","W_foliage_kg_tC_ha","W_stem_kg_tC_ha","W_root_kg_tC_ha" )]
merge(all.bio,mcl,by=c('pplot'),all=T)->all.mcl

merge(all.mcl,all.bio,all=T)->all.bio1

write.csv(all.bio1,'all.bio.csv',row.names=F)

read.csv('all.bio.csv')->bio1.site    #carbon components
names(bio1.site)[1]<-'site'
merge(site,bio1.site,by=c('site','lat','long'))->bio.site
#NOW， WE GOT the full file list scheduled for next calaulation.

read.csv('lai.site.csv')->lai.site  # For ts correction
read.csv('tvdi.site.csv')->tvdi.site  # for soil moisture corection
read.csv('tp.site.csv')->tp.site     #temp and precipation
read.csv('bio.site.csv')->bio.site    #carbon components


#variable to generate : LAI1-LAI12,TVDI.ndvi1-12,TVDI.evi1-12,Ts1-12,Wleaf_kg_tC_ha,W_foliage_kg_tC_ha,W_stem_kg_tC_ha,W_root_kg_tC_ha,ratio_FLR_toBio_pct,stems

#> names(tvdi.site)
#[1] "site"           "lat"            "long"           "age_max"        "aspect"         "slope"          "elev"           "tpi"            "twi"            "soil_gang"      "soil_cat"      
#[12] "sub_soil_cat"   "max_spec"       "TVDI.evi.m.1"   "TVDI.evi.m.2"   "TVDI.evi.m.3"   "TVDI.evi.m.4"   "TVDI.evi.m.5"   "TVDI.evi.m.6"   "TVDI.evi.m.7"   "TVDI.evi.m.8"   "TVDI.evi.m.9"  
#[23] "TVDI.evi.m.10"  "TVDI.evi.m.11"  "TVDI.evi.m.12"  "TVDI.ndvi.m.1"  "TVDI.ndvi.m.2"  "TVDI.ndvi.m.3"  "TVDI.ndvi.m.4"  "TVDI.ndvi.m.5"  "TVDI.ndvi.m.6"  "TVDI.ndvi.m.7"  "TVDI.ndvi.m.8" 
#[34] "TVDI.ndvi.m.9"  "TVDI.ndvi.m.10" "TVDI.ndvi.m.11" "TVDI.ndvi.m.12"
#> 
#> names(lai.site)
#[1] "site"         "lat"          "long"         "age_max"      "aspect"       "slope"        "elev"         "tpi"          "twi"          "soil_gang"    "soil_cat"     "sub_soil_cat"
#[13] "max_spec"     "L95_1"        "L95_2"        "L95_3"        "L95_4"        "L95_5"        "L95_6"        "L95_7"        "L95_8"        "L95_9"        "L95_10"       "L95_11"      
#[25] "L95_12"       "L05_1"        "L05_2"        "L05_3"        "L05_4"        "L05_5"        "L05_6"        "L05_7"        "L05_8"        "L05_9"        "L05_10"       "L05_11"      
#[37] "L05_12"       "L50_1"        "L50_2"        "L50_3"        "L50_4"        "L50_5"        "L50_6"        "L50_7"        "L50_8"        "L50_9"        "L50_10"       "L50_11"      
#[49] "L50_12"       "Lm_1"         "Lm_2"         "Lm_3"         "Lm_4"         "Lm_5"         "Lm_6"         "Lm_7"         "Lm_8"         "Lm_9"         "Lm_10"        "Lm_11"       
#[61] "Lm_12"       
#> 
#> names(bio.site)
#[1] "site"                 "lat"                  "long"                 "aspect"               "slope"                "elev"                 "tpi"                  "twi"                 
#[9] "soil_gang"            "soil_cat"             "sub_soil_cat"         "max_spec"             "plot_area_m2"         "cls_max_spec4"        "age_max"              "xj.l"                
#[17] "stems"                "TotC_eco"             "SoilC_1m_t_ha"        "SoilCden_t_ha_0.20cm" "W_total_tC_ha"        "W_leaf_kg_tC_ha"      "W_foliage_kg_tC_ha"   "W_stem_kg_tC_ha"     
#[25] "W_root_kg_tC_ha"      "ratio_FLR_toBio_pct" 
#>#	
#> names(tp.site)
#[1] "site"         "lat"          "long"         "age_max"      "aspect"       "slope"        "elev"         "tpi"          "twi"          "soil_gang"    "soil_cat"     "sub_soil_cat"
#[13] "max_spec"     "Ta1"          "Ta2"          "Ta3"          "Ta4"          "Ta5"          "Ta6"          "Ta7"          "Ta8"          "Ta9"          "Ta10"         "Ta11"        
#[25] "Ta12"         "P1"           "P2"           "P3"           "P4"           "P5"           "P6"           "P7"           "P8"           "P9"           "P10"          "P11"         
#[37] "P12"  


gf.rflrb<-rfsrc(ratio_FLR_toBio_pct~lat+long+age_max,data=bio.site,importance="permute.ensemble",ntree=200, seed=-111)
gf.rflrb.p<-plot.variable(gf.rflrb,partial=TRUE)
plot(gf.rflrb)

#To check the ratio of rstem

read.csv('rstem.csv')->rstem
names(rstem)
#[1] "age"        "ratio_stem" "ratio_FL"  

r.lo<-loess(ratio_stem~age,span=0.25,data=rstem)
r.lo_p<-predict(r.lo,se=TRUE)
aa<-cbind(rstem,r.lo_p$fit,r.lo_p$se)
write.csv(aa,'d://temp.csv')

#to check the stems

with(bio.site,plot(stems~age_max))
r.lo<-loess(stems~age_max,span=0.25,data=bio.site)
r.lo_p<-predict(r.lo,se=TRUE)
aa<-cbind(bio.site$age_max,bio.site$stems,r.lo_p$fit,r.lo_p$se)
write.csv(aa,'d://temp.csv')

require(HelpersMG)

plot_errbar(x=x, y=y, xlab="axe x", ylab="axe y", bty="n", xlim=c(1,20),
		y.minus=y-1, y.plus=y+1, ylim=c(-3, 3), type="l",lwd=2,
		errbar.y.polygon=TRUE,
		errbar.y.polygon.list=list(border=NA, col=rgb(0.75,0.25,0.1,0.25)))

require(nplr) #for n-parameter logistic curve
require(grofit) #for richards curves


