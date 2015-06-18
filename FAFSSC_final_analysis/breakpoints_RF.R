
require(segmented)
require(breakpoint)

bks<-function(obs){
  est_pts<-CE.Normal(as.data.frame(obs$y),Nmax=10)
  if(is.na(as.numeric(est_pts[1]))){a=(matrix(NA,nrow=1,ncol=3));}
  else{
  out.lm<-lm(y~x,data=obs)

  if(est_pts$No.BPs>2){loc.idx<-range(est_pts$BP.Loc);}else{loc.idx=c(est_pts$BP.Loc)}
  o<-segmented(out.lm,seg.Z=~x,psi=list(x=c(obs$x[loc.idx])),
               control=seg.control(display=FALSE))
  a=c(matrix(o$psi[,2:3],nrow=1),coef(o),intercept(o)$x)
  }
  return(a)
}



bks2<-function(obs){
  est_pts<-CE.Normal(as.data.frame(obs$y),Nmax=10)
  if(is.na(as.numeric(est_pts[1]))){a=(matrix(NA,nrow=1,ncol=3));}
  else{
    out.lm<-lm(y~x,data=obs)
    
    if(est_pts$No.BPs>2){loc.idx<-range(est_pts$BP.Loc);}else{loc.idx=c(est_pts$BP.Loc)}
    o<-segmented(out.lm,seg.Z=~x,psi=list(x=c(obs$x[loc.idx])),
                 control=seg.control(display=FALSE))
    a=c(matrix(o$psi[,2:3],nrow=1),coef(o),intercept(o)$x)
  }
  return(a)


bks1<-function(obs){
  est_pts<-CE.Normal(as.data.frame(obs$y),Nmax=3)
  if(is.na(as.numeric(est_pts[1])))
    {a=(matrix(NA,nrow=1,ncol=3));return(a)}
  else{return(obs$x[est_pts$BP.Loc])}
}

#-----------------
read.csv("rst.age_carbon_wide.csv")->rst.age_carbon

plot.cls<-as.character(levels(rst.age_carbon$VEG.SECTION))
ncls<-length(plot.cls)
rst.age_carbon.sl<-data.frame(matrix(NA,nrow = ncls*3, ncol = 10));

names(rst.age_carbon.sl)<-c("VEG.SECTION","VARIABLEs", "a", "k", "xc", "a.se", "k.se", "xc.se", "rse", "R2")

sink('sink_brkpts.txt')

for (i in 1:ncls){
  
  ii<-which(as.character(rst.age_carbon$VEG.SECTION) %in% plot.cls[i])
  d2z<-rst.age_carbon[ii,]
  names(d2z)<-c("VEG.SECTION","age","CBiomass","CSoil1m","CSoil20cm")
  pts.veg<-try(bks(obs=data.frame(approx(x=d2z$age,y=d2z$CBiomass))))
  pts.s1m<-try(bks(obs=data.frame(approx(x=d2z$age,y=d2z$CSoil1m))))
  pts.s20<-try(bks(obs=data.frame(approx(x=d2z$age,y=d2z$CSoil20cm))))

  pts.veg1<-bks1(obs=data.frame(approx(x=d2z$age,y=d2z$CBiomass)))
  pts.s1m1<-bks1(obs=data.frame(approx(x=d2z$age,y=d2z$CSoil1m)))
  pts.s201<-bks1(obs=data.frame(approx(x=d2z$age,y=d2z$CSoil20cm)))
  
  print(c(plot.cls[i],"CBiomass",pts.veg1));print(pts.veg)
  print(c(plot.cls[i],"CSoil1m",pts.s1m1));print(pts.s1m)
  print(c(plot.cls[i],"CSoil20cm",pts.s201));print(pts.s20)

  rm(pts.veg,pts.s20,pts.s1m,pts.veg1,pts.s201,pts.s1m1)
}

sink()
