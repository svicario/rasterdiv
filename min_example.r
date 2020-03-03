#Testing
setwd("~/GitHub/rasterdiv/")
lapply(list.files(pattern="*.R$"), function(x) {source(x)})

library(doSNOW)
library(raster)

a<-matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
b<-raster(a)

##Hill
#Mode single
uno<-Hill(input=a,window=3,alpha=2)
due<-Hill(input=b,window=3,alpha=2)
tre<-Hill(input=a,window=3,alpha=2,nc.cores=2)
quattro<-Hill(input=b,window=3,alpha=2,nc.cores=2)

#Check if output is identical
identical(uno,due,tre,quattro)

#Mode iterative
a<-matrix(c(-0.5,1,0.5,1,0.5,1,1,1,0.5),ncol=3,nrow=3)
b<-raster(a)
uno<-Hill(input=a,window=3,alpha=c(1,5),mode="iterative")
due<-Hill(input=b,window=3,alpha=c(1,5),mode="iterative")
tre<-Hill(input=a,window=3,alpha=c(1,5),nc.cores=2,mode="iterative")
quattro<-Hill(input=b,window=3,alpha=c(1,5),nc.cores=2,mode="iterative")

#Check if outputs are identical
lapply(1:5, function(x) identical(uno[x][[1]],due[x][[1]],tre[x][[1]],quattro[x][[1]]))

#Rao
a<-matrix(c(20,100,10,10,10,100,100,20,20,20,100,10,10,10,100,100,20,20,20,100,10,10,10,100,100,20,20),ncol=3,nrow=9)
b<-raster(a)
rao_n <- Rao(a,distance_m="euclidean",window=3,shannon=FALSE,nc.cores=3,cluster.type="SOCK",na.tolerance=1.0)
rao_o <- spectralrao(a,distance_m="euclidean",window=3,shannon=FALSE,nc.cores=1,cluster.type="SOCK",na.tolerance=1.0)
