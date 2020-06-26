post.check<-function(fit,tree,X=NULL,nsim=1000){
  if(length(dim(fit$chains))==3){
    fit<-combine.chains(fit)
  }
  n_e<-nrow(tree$edge)
  samp<-sample(1:dim(fit$chains)[1],nsim,replace=T)
  if(any(grepl('R\\[',dimnames(fit$chains)[[2]]))){
    R_sample<-exp((fit%chains%1:nrow(tree$edge))[samp,])
  }else{
    R_sample<-exp(do.call(cbind,rep(list((fit%chains%'R0')[samp]),n_e)))
  }
  if(all(tree$tip.label%in%dimnames(fit$chains)[[2]])){
    X_sample<-(fit%chains%tree$tip.label)[samp,]
  }else{
    if(is.null(X)){
      stop('need trait vals')
    }
    X_sample<-do.call(rbind,rep(list(X),nsim))
  }
  tmp.tree<-tree
  sim<-vector('numeric',nsim)
  emp<-sim
  for(i in 1:nsim){
    tmp.tree$edge.length<-tree$edge.length*R_sample[i,]
    sim.X<-fastBM(tmp.tree)
    sim.pics2<-pic(sim.X,tree)^2
    sim[i]<-sd(sim.pics2)/mean(sim.pics2)
    emp.pics2<-pic(X_sample[i,],tree)^2
    emp[i]<-sd(emp.pics2)/mean(emp.pics2)
  }
  list(sim=sim,emp=emp)
}

Rsig2.mod<-post.check(test.fit,tree,X=X,nsim=1e4)
bm.mod<-post.check(bm.fit,tree,X=X,nsim=1e4)
trend.mod<-post.check(trend.fit,tree,X=X,nsim=1e4)
full.mod<-post.check(full.fit,tree,X=X,nsim=1e4)

fit.plot<-function(obj){
  dens<-density(obj)
  plot(dens,xlim=range(obj))
  quant.range<-do.call(':',as.list(sapply(quantile(obj,probs=c(0.025,0.975)),
                                          function(ii) which.min(abs(dens$x-ii)))))
  polygon(c(0,dens$y[quant.range],0)~
            rep(dens$x[quant.range],c(2,rep(1,length(quant.range)-2),2)),
          col='black',border=NULL)
}

emp.pics2<-pic(X,tree)^2
emp.stat<-sd(emp.pics2)/mean(emp.pics2)
fit.plot(Rsig2.mod$sim)
abline(v=emp.stat,col='red')
fit.plot(bm.mod$sim)
abline(v=emp.stat,col='red')
fit.plot(trend.mod$sim)
abline(v=emp.stat,col='red')
fit.plot(full.mod$sim)
abline(v=emp.stat,col='red')




sum(Rsig2.mod$sim>emp.stat)/1e4
sum(bm.mod$sim>emp.stat)/1e4
sum(trend.mod$sim>emp.stat)/1e4
sum(full.mod$sim>emp.stat)/1e4

bc<-function(obj){
  tmp<-unlist(obj)
  emp<-density(obj$emp,from=min(tmp),to=max(tmp),bw=diff(range(tmp))/1e6)
  sim<-density(obj$sim,from=min(tmp),to=max(tmp),bw=diff(range(tmp))/1e6)
  xx<-1:length(emp$x)
  dx<-diff(emp$x)[1]
  sum(sapply(xx,function(ii) dx*sqrt(emp$y[xx]*sim$y[xx])))
}

bc<-function(obj){
  tmp<-unlist(obj)
  xx<-seq(min(tmp),max(tmp),length.out=1e3)
  dx<-diff(xx)[1]
  sum<-0
  for(i in 1:(length(xx)-1)){
    sum<-sum+sqrt(sum(obj$emp>=xx[i]&obj$emp<=xx[i+1])*sum(obj$sim>=xx[i]&obj$sim<=xx[i+1]))
  }
  sum
}

bc(bm.mod)

sum(emp.stat<=Rsig2.mod)/1000
sum(emp.stat<=bm.mod)/1000

#between this and looking at Rsig2, one can really quite easily assess evidence for rate heterogeneity!

post.check<-function(fit,tree,X,nsim=1000){
  if(length(dim(fit$chains))==3){
    fit<-combine.chains(fit)
  }
  n_e<-nrow(tree$edge)
  samp<-sample(1:dim(fit$chains)[1],nsim,replace=T)
  if(any(grepl('R\\[',dimnames(fit$chains)[[2]]))){
    R_sample<-exp((fit%chains%1:nrow(tree$edge))[samp,])
  }else{
    R_sample<-exp(do.call(cbind,rep(list((fit%chains%'R0')[samp]),n_e)))
  }
  tmp.tree<-tree
  simpics2<-matrix(NA,length(tree$tip.label)-1,nsim)
  for(i in 1:nsim){
    tmp.tree$edge.length<-tree$edge.length*R_sample[i,]
    simpics2[,i]<-pic(fastBM(tmp.tree),tree)^2
  }
  simpics2
}
Rsig2.mod2<-post.check(test.fit,tree,X)
plot(as.vector(log(Rsig2.mod2))~rep(1:length(emp.pics2),ncol(Rsig2.mod2)),pch=16,col=rgb(0,0,0,0.1))
points(log(emp.pics2),col='red',pch=16,cex=2)
bm.mod2<-post.check(bm.fit,tree,X)
plot(as.vector(log(bm.mod2))~rep(1:length(emp.pics2),ncol(bm.mod2)),pch=16,col=rgb(0,0,0,0.1))
points(log(emp.pics2),col='red',pch=16,cex=2)
fit.plot<-function(obj,emp){
  plot(as.vector(log(obj))~rep(1:length(emp),ncol(obj)),pch=16,col=rgb(0,0,0,0.1))
  points(log(emp),col='red',pch=16,cex=2)
}
trend.mod2<-post.check(trend.fit,tree,X)
fit.plot(trend.mod2,emp.pics2)
full.mod2<-post.check(full.fit,tree,X)
fit.plot(full.mod2,emp.pics2)
