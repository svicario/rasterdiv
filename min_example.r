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
tre<-Hill(input=a,window=3,alpha=2,n.process=2)
quattro<-Hill(input=b,window=3,alpha=2,n.process=2)

#Check whether output is identical
identical(uno,due,tre,quattro)

#Mode iterative
uno<-Hill(input=a,window=3,alpha=1:5)
due<-Hill(input=b,window=3,alpha=1:5)
tre<-Hill(input=a,window=3,alpha=1:5,n.process=2)
quattro<-Hill(input=b,window=3,alpha=1:5,n.process=2)

#Check whether outputs are identical
lapply(1:5, function(x) identical(uno[x][[1]],due[x][[1]],tre[x][[1]],quattro[x][[1]]))

##Rao
uno<-Rao(a,distance_m="euclidean",window=3,shannon=FALSE,n.process=1,cluster.type="SOCK",na.tolerance=0)
due<-Rao(b,distance_m="euclidean",window=3,shannon=FALSE,n.process=1,cluster.type="SOCK",na.tolerance=0)

##Renyi
#Mode single
uno<-Renyi(input=a,window=3,alpha=2)
due<-Renyi(input=b,window=3,alpha=2)
tre<-Renyi(input=a,window=3,alpha=2,n.process=2)
quattro<-Renyi(input=b,window=3,alpha=2,n.process=2)
#Check whether output is identical
identical(uno,due,tre,quattro)

#Mode iterative
due<-Hill(input=b,window=3,alpha=0:5)
uno<-Hill(input=a,window=3,alpha=0:5)
tre<-Hill(input=a,window=3,alpha=0:5,n.process=2)
quattro<-Hill(input=b,window=3,alpha=0:5,n.process=2)
#Check whether output is identical
lapply(1:5, function(x) identical(uno[x][[1]],due[x][[1]],tre[x][[1]],quattro[x][[1]]))
