setwd("D://东北//data//基础信息//RESULT//csv")

read.csv('bio_soc_env_all_gsla.csv')->all.gsl

all.gsl1<-all.gsl

levl<-as.character(levels(as.factor(all.gsl1$cls_max_spec4)))
nlevl<-length(levl)

coff<-as.data.frame(matrix(NA,ncol=4,nrow=nlevl*3))
names(coff)<-c("ZONE","VARIABLE","a","b")
for (i in 1 : nlevl){
	ii<-which(all.gsl1$cls_max_spec4 %in% levl[i])
	coff[(3*(i-1)+1),1:4]<-c(levl[i],"TOTC",coef(lm(TotC_eco2~TotC_eco,data=all.gsl1[ii,])))
	coff[(3*(i-1)+2),1:4]<-c(levl[i],"BIO",coef(lm(W_total_t_ha2~W_total_t_ha,data=all.gsl1[ii,])))
	coff[(3*(i-1)+3),1:4]<-c(levl[i],"SOILC",coef(lm(SoilC_1m_t_ha2~SoilC_1m_t_ha,data=all.gsl1[ii,])))
}

#[1] "A2A" "A2B" "A2C" "A2D" "A3C" "B2B" "B2D" "B2D" "B4C" "B4D" "C1B" "C2B" "C3B" "C3C" "C4C"

#------------------------------------------######################
require(doBy)
clst<-c("A2A","A2B","A2C","A2D","A3C","B2B","B2D","B2D","B4C","B4D","C1B","C2B","C3B","C3C","C4C")

setwd("D://东北//data//基础信息//RESULT//csv")

read.csv('bio_soc_env_all_gsla.csv')->all.gsl

ii<-which(all.gsl$cls_max_spec4 %in% clst)

t.m<-summaryBy(ta_C~cls_max_spec4,data=all.gsl[ii,],FUN=mean,na.rm=T)
p.m<-summaryBy(p_mm~cls_max_spec4,data=all.gsl[ii,],FUN=mean,na.rm=T)

tp<-merge(t.m,p.m)
