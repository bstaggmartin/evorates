#Note: This file is for development purposes only, and is not run like a traditional test file (it can take
#~20 minutes to run on some machines). This is for making sure models still run correctly under various
#scenarios following updates to the fit.corateBM() function and/or any of the stan scripts it uses.

library(testthat)
set.seed(123)
tree<-phytools::pbtree(n=25)
uni.sim<-gen.corateBM(tree,intra.var=T,n_obs=5)
multi.sim<-gen.corateBM(tree,intra.var=T,n_obs=5,X0=rep(0,3))

#basic models
dat.list<-rep(list(uni.sim,multi.sim),each=4)
Rsig2.list<-rep(c(T,T,F,F),2)
trend.list<-rep(c(F,T,F,T),2)
fit.list<-vector('list',8)
for(i in 1:8){
  fit.list[[i]]<-fit.corateBM(tree,dat.list[[i]]$X,
                              constrain.Rsig2=Rsig2.list[[i]],trend=trend.list[[i]],
                              chains=1,iter=1000)
}
#output expected
for(i in 1:8){
  expect_s3_class(fit.list[[i]],'corateBM_fit')
}
#parameters make sense
##fitted rate values correlated with true ones
for(i in c(3,4,7,8)){
  expect_gt(cor((fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R),0.7)
}

#missing data
mis.uni.inds<-cbind(c(1,13),c(1,1))
mis.uni.dat<-uni.sim$X;mis.uni.dat[mis.uni.inds,]<-NA
mis.multi.inds<-cbind(c(4,7,8,14,14,20),c(2,3,1,1,3,2))
mis.multi.dat<-multi.sim$X;mis.multi.dat[mis.multi.inds]<-NA
mis.dat.list<-rep(list(mis.uni.dat,mis.multi.dat),each=4)
mis.inds.list<-rep(list(mis.uni.inds,mis.multi.inds),each=4)
mis.fit.list<-fit.list
for(i in 1:8){
  mis.fit.list[[i]]<-fit.corateBM(tree,mis.dat.list[[i]],
                                  constrain.Rsig2=Rsig2.list[[i]],trend=trend.list[[i]],
                                  chains=1,iter=1000)
}
#output expected
for(i in 1:8){
  expect_s3_class(mis.fit.list[[i]],'corateBM_fit')
}
#parameters make sense
##fitted rate values correlated with true ones
for(i in c(3,4,7,8)){
  expect_gt(cor((mis.fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R),0.7)
}
##true tip values fall within distribution of fitted ones
for(i in 1:8){
  tmp<-colSums((mis.fit.list[[i]]%chains%'t\\d+$')>
                 rep(dat.list[[i]]$X[mis.inds.list[[i]]],each=500))/500
  expect_false(isTRUE(all.equal(rep(0,length(tmp)),tmp))&isTRUE(all.equal(rep(1,length(tmp)),tmp)))
}

#intraspecific variation
intra.fit.list<-fit.list
for(i in 1:8){
  intra.fit.list[[i]]<-fit.corateBM(tree,dat.list[[i]]$Y,
                                  constrain.Rsig2=Rsig2.list[[i]],trend=trend.list[[i]],
                                  chains=1,iter=1000,
                                  intra.var=T)
}
#output expected
for(i in 1:8){
  expect_s3_class(intra.fit.list[[i]],'corateBM_fit')
}
#parameters make sense
##fitted rate values correlated with true ones
for(i in c(3,4,7,8)){
  expect_gt(cor((intra.fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R),0.7)
}
##fitted trait values correlated with true ones
for(i in 1:8){
  expect_gt(cor(intra.fit.list[[i]]%means%'t\\d+$',as.vector(t(dat.list[[i]]$X))),0.9)
}

#quick check for multiple chains
multchains.fit<-fit.corateBM(tree,dat.list[[i]]$Y,
                             constrain.Rsig2=Rsig2.list[[i]],trend=trend.list[[i]],
                             chains=3,iter=1000,cores=3,
                             intra.var=T)
expect_s3_class(multchains.fit,'corateBM_fit')
#everything looks good, but need to write formal tests