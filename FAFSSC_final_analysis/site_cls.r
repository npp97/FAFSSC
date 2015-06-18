#age estimation with RandomForests
#This verison replace the influence of sp with p_mm 
require(randomForest)
#read data into memory
bk_geoinfo<-read.csv('geoinfo_zzzy.csv')
bk_soil<-read.csv('soil_class_zzzy.csv')
bk_veg<-read.csv("dom_spec_zzzy.csv")

bk<-merge(bk_geoinfo,bk_soil,by=c("pplot","lat","long","elev"))
#bk<-merge(bk,bk_veg,by=c("pplot"))

bkk<-bk[,c("pplot","lat","long","elev","site","slope","aspect","twi","tpi","p_mm","ta_C","code_soil_cat","code_soil_gang","code_soil_sub_cat")]
#data.frame preparation for randomFrest est. of age

#Train the forest using randomForestSRC package
set.seed(111)

cls.rf1 <- randomForest(x=bkk[,c("elev","slope","aspect","twi","tpi","p_mm","ta_C","code_soil_cat","code_soil_sub_cat")],ntree=1500)

b<-predict(cls.rf1,bkk[,c("elev","slope","aspect","twi","tpi","p_mm","ta_C","code_soil_cat","code_soil_sub_cat")])


write.csv(bk_zzzy, file='age_zzzy.csv',row.names=FALSE)
write.csv(bk_cproj, file='age_cproj.csv',row.names=FALSE)