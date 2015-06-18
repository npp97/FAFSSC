# TODO: Add comment
# 
# Author: junhui.zhang
###############################################################################
#Packages Declaration
#------------------------------------------------------------------------------
require(doBy)

#User-defined Function
#------------------------------------------------------------------------------

qa<-function(x,...){
	return(c(quantile(x,c(0.01,1),...)))
}

sumfun <- function(x, ...) {
	c(tm=mean(x,0.025, ...), m = mean(x, ...), md = median(x, ...), v = sd(x, ...), l = length(x))
}

#------------------------------------------------------------------------------
#								Code Body
#------------------------------------------------------------------------------
read.csv("碳专项植物成分调查20150610.csv")->cnp

names(cnp)<-c("id","ydID","sample_cat","sample_id","year","month","day","C_con_g_kg","N_con_g_kg","P_con_g_kg","CNratio")

cnp<-cnp[,c("ydID","sample_cat","C_con_g_kg","N_con_g_kg","P_con_g_kg","CNratio")]

cnp.sum<-summaryBy(C_con_g_kg+N_con_g_kg+P_con_g_kg+CNratio~ydID+sample_cat,data=cnp, na.rm=T,FUN=sumfun)

age<-read.csv('bkinfo_for_morility.csv')
names(age)[1]<-'ydID'

merge(cnp.sum,age,by='ydID')->cnp_sum_age
dim(cnp_sum_age)

write.csv(cnp_sum_age,'cnp_sum_age.csv',row.names=F)

