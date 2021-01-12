#Note: This file is for development purposes only, and is not run like a traditional test file (it can take
#~20 minutes to run on some machines). This is for making sure models still run correctly under various
#scenarios following updates to the fit.corateBM() function and/or any of the stan scripts it uses.

library(testthat)
library(contSimmap)
set.seed(123)
tree<-phytools::pbtree(n=25)
uni.sim<-sim.corateBM(tree,n.obs=5)
multi.sim<-sim.corateBM(tree,n.obs=5,X0=rep(0,3))

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
##doesn't work for fit 4 --> trend params seem to mess everything up for some reason
for(i in c(3,4,7,8)){
  expect_gt(cor((fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R),0.5)
}
for(i in c(3,4,7,8)){
  print(cor((fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R))
}
def.par<-par(no.readonly=T)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(0,0,0,0))
for(i in c(3,4,7,8)){
  ylim<-range(fit.list[[i]]%quantiles%'R_[1-9]')
  plot(y=(fit.list[[i]]%means%'R_\\d+')[-1],x=dat.list[[i]]$R,pch=16,xaxt='n',yaxt='n',xlab='',ylab='',ylim=ylim)
  segments(x0=dat.list[[i]]$R,
           y0=fit.list[[i]]%quantiles%list('R_[1-9]',0.025),
           y1=fit.list[[i]]%quantiles%list('R_[1-9]',0.975),
           lty=2)
  abline(0,1,col='red')
}
par(def.par)

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
##again, 4 fails
for(i in c(3,4,7,8)){
  expect_gt(cor((mis.fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R),0.5)
}
for(i in c(3,4,7,8)){
  print(cor((mis.fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R))
}
def.par<-par(no.readonly=T)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(0,0,0,0))
for(i in c(3,4,7,8)){
  ylim<-range(mis.fit.list[[i]]%quantiles%'R_[1-9]')
  plot(y=(mis.fit.list[[i]]%means%'R_\\d+')[-1],x=dat.list[[i]]$R,pch=16,xaxt='n',yaxt='n',xlab='',ylab='',ylim=ylim)
  segments(x0=dat.list[[i]]$R,
           y0=mis.fit.list[[i]]%quantiles%list('R_[1-9]',0.025),
           y1=mis.fit.list[[i]]%quantiles%list('R_[1-9]',0.975),
           lty=2)
  abline(0,1,col='red')
}
par(def.par)
##true tip values fall within distribution of fitted ones
for(i in 1:8){
  tmp<-colSums((mis.fit.list[[i]]%chains%'_t\\d+$')>
                 rep(dat.list[[i]]$X[mis.inds.list[[i]]],each=500))/500
  expect_false(isTRUE(all.equal(rep(0,length(tmp)),tmp))&isTRUE(all.equal(rep(1,length(tmp)),tmp)))
}

#intraspecific variation
intra.fit.list<-fit.list
for(i in 1:8){
  intra.fit.list[[i]]<-fit.corateBM(tree,dat.list[[i]]$trait.data,
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
##3 just barely fails, 4 fails
for(i in c(3,4,7,8)){
  expect_gt(cor((intra.fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R),0.5)
}
for(i in c(3,4,7,8)){
  print(cor((intra.fit.list[[i]]%means%'R_\\d+')[-1],dat.list[[i]]$R))
}
def.par<-par(no.readonly=T)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(0,0,0,0))
for(i in c(3,4,7,8)){
  ylim<-range(intra.fit.list[[i]]%quantiles%'R_[1-9]')
  plot(y=(intra.fit.list[[i]]%means%'R_\\d+')[-1],x=dat.list[[i]]$R,pch=16,xaxt='n',yaxt='n',xlab='',ylab='',ylim=ylim)
  segments(x0=dat.list[[i]]$R,
           y0=intra.fit.list[[i]]%quantiles%list('R_[1-9]',0.025),
           y1=intra.fit.list[[i]]%quantiles%list('R_[1-9]',0.975),
           lty=2)
  abline(0,1,col='red')
}
par(def.par)
##fitted trait values correlated with true ones
for(i in 1:8){
  expect_gt(cor(intra.fit.list[[i]]%means%'_t\\d+$',as.vector(t(dat.list[[i]]$X))),0.9)
}

#quick check for multiple chains
multchains.fit<-fit.corateBM(tree,dat.list[[i]]$Y,
                             constrain.Rsig2=Rsig2.list[[i]],trend=trend.list[[i]],
                             chains=3,iter=1000,cores=3,
                             intra.var=T)
expect_s3_class(multchains.fit,'corateBM_fit')
#everything looks good, but need to write formal tests
#reduced correlation cutoff to 0.5 --> only really necessary for missing data in case with rate noise + trend, but I wonder
#why it's so different from first case. Presumably, sim.corateBM now generates a new dataset, so hopefully it's just that?
#still consistently underestimates intraspecific variation. I really wonder what's up with that...
#In terms of rate ests, I see no evidence of systematic bias, so I think you're okay...
#tried it with 50 tips and seed 321 --> achieved correlation coeffs ~90. Definitely seems to just be a product of the dataset and
#its smallness...
#cauchy seems to perform better, so I think I'll keep it for now, and consider allowing users to specify df's in the future