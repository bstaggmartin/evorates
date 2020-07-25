library(contSimmap)
set.seed(123)
tree<-phytools::pbtree(n=100,scale=1)
sim<-gen.corateBM(tree,X0=rnorm(10),Xsig2=rWishart(1,10,diag(10))[,,1],
                  intra.var=T,Ysig2=0.1*rWishart(1,10,diag(10))[,,1],n_obs=rpois(100,2))
plot(sim,tree)
fit<-fit.corateBM(tree,sim$Y,intra.var=T,
                  chains=3,iter=4000,cores=3,refresh=10,control=list(max_treedepth=12))
mat<-matrix(NA,length(sim$X0),length(sim$X0))
mat[lower.tri(mat,diag=T)]<-fit%means%'evocov'
mat2<-mat;mat3<-mat
mat2[lower.tri(mat2,diag=T)]<-fit%quantiles%list('evocov',0.975)
mat3[lower.tri(mat3,diag=T)]<-fit%quantiles%list('evocov',0.025)
plot(c(mat,mat2,mat3)~rep(as.vector(sim$Xsig2),3),pch=rep(c(16,1,1),each=length(mat)))
#to only get post-warmup
for(i in 1:3){
  tmp<-attr(ret@sim$samples[[i]],'sampler_params')
  tmp<-lapply(tmp,function(ii) ii[-(niter-1):0+length(ii)])
  sampler[,,i]<-as.matrix(do.call(cbind,
                                  c(tmp,
                                    list(lp__[,i,]))))
}
#did a 100 tip tree with no intraspecific var with 10 traits in ~10ish minutes and accurately inferred
#covariance measures for the most part
#intra var stress test -- seems to work, but VERY slow. Didn't fully burn in after 1000 generations with
#500 warmup--rates and intracov both mixed very slowly, as you'd expect...
#run for a long time, with thinning...
sim2<-gen.corateBM(tree,X0=rnorm(3),Xsig2=rWishart(1,3,diag(3))[,,1])
fit2<-fit.corateBM(tree,sim2$X,chains=1,iter=1000,refresh=1,intra.var=T)

vec<-rep(NA,100)
vec2<-vec
vec3<-vec
matmat<-matrix(NA,nrow(tree$edge),100)
for(i in 1:100){
  tmp<-opt.fun(tree,sim2$X,trend=F,constrain.Rsig2=F)
  # if(tmp$return_code!=0){
    vec[i]<-tmp$par['Rsig2[1]']
    vec2[i]<-tmp$value
    vec3[i]<-tmp$return_code
    matmat[,i]<-tmp$par[paste('R[',1:nrow(tree$edge),']',sep='')]
  # }
  cat(i,'\n')
}
fit2<-fit.corateBM(tree,sim2$X,iter=1000)
fit3<-fit.corateBM(tree,sim2$X,iter=1000,chains=1,init=list(as.list(tmp$par)))
#optimizing actually working okay...may be worth integrating eventually as an option
#must be find very thin, high-likelihood peaks?
#or does it measure likelihood differently? Seems to settle on similar param ests to Bayesian alg, after
#getting rid of failed attempts

sim<-readRDS('large_dat')
fit<-readRDS('large_fit')
tree<-readRDS('large_tree')
library(contSimmap)
library(ape)

k<-2
plot(0,type='n',xlim=range(sim$Ysig2),ylim=range(fit%chains%'intracov'),
     xlab='true',ylab='posterior distribution')
for(i in 1:k){
  for(j in  1:i){
    dens<-density(fit%chains%paste('X',i,',X',j,'_intracov',sep=''),bw=0.05)
    yy<-c(dens$x,rev(dens$x))
    xx<-0.01*c(dens$y,-1*rev(dens$y))+sim$Ysig2[i,j]
    polygon(yy~xx,border=NA,col=rgb(0.5,0.5,0.5,0.2))
  }
}
points(apply(fit%means%'intracov',1,mean)~sim$Ysig2[lower.tri(sim$Ysig2,diag=T)],pch=16)
abline(0,1,col='red')

plot(0,type='n',xlim=range(sim$Xsig2),ylim=range(fit%chains%'evocov'),
     xlab='true',ylab='posterior distribution')
for(i in 1:k){
  for(j in  1:i){
    dens<-density(fit%chains%paste('X',i,',X',j,'_evocov',sep=''),bw=0.05)
    yy<-c(dens$x,rev(dens$x))
    xx<-0.01*c(dens$y,-1*rev(dens$y))+sim$Xsig2[i,j]
    polygon(yy~xx,border=NA,col=rgb(0.5,0.5,0.5,0.2))
  }
}
points(apply(fit%means%'evocov',1,mean)~sim$Xsig2[lower.tri(sim$Xsig2,diag=T)],pch=16)
abline(0,1,col='red')

plot(0,type='n',xlim=range(sim$R),ylim=range(fit%chains%1:nrow(tree$edge)),
     xlab='true',ylab='posterior distribution')
for(i in 1:nrow(tree$edge)){
  dens<-density(fit%chains%i,bw=0.05)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.05*c(dens$y,-1*rev(dens$y))+sim$R[i]
  polygon(yy~xx,border=NA,col=rgb(0.5,0.5,0.5,0.2))
}
points(apply(fit%means%1:nrow(tree$edge),1,mean)~sim$R,pch=16)
abline(0,1,col='red')

plot(0,type='n',xlim=range(sim$X),ylim=range(fit%chains%tree$tip.label),
     xlab='true',ylab='posterior distribution')
for(i in 1:k){
  for(j in tree$tip.label){
    dens<-density(fit%chains%paste('X',i,'_',j,'$',sep=''),bw=0.05)
    yy<-c(dens$x,rev(dens$x))
    xx<-c(dens$y,-1*rev(dens$y))+sim$X[j,i]
    polygon(yy~xx,border=NA,col=rgb(0.5,0.5,0.5,0.2))
  }
}
points(apply(fit%means%tree$tip.label,1,mean)~as.vector(t(sim$X)),pch=16)
abline(0,1,col='red')

#large_fit took 7.5 hours (4.5 hours warmup, 3 hours sampling, 4000 iterations with first 2000 devoted to
#warmup, with max_treedepth at 12--kinda overkill--had >500 ess for nearly every parameter!)
#no  divergent or max treedepth problems after warmup

tol.el<-sum(tree$edge.length)
wgts<-tree$edge.length/tol.el
plot(0,type='n',xlim=c(1,nrow(tree$edge)),ylim=range(fit%chains%'dev'),
     xlab='branch',ylab='ln(rate) posterior distribution')
segments(x0=rep(1:nrow(tree$edge),each=3),
         y0=as.vector(apply(fit%quantiles%list('dev',c(0,0.025,0.975)),c(1,2),mean)),
         y1=as.vector(apply(fit%quantiles%list('dev',c(0.025,0.975,1)),c(1,2),mean)),
         lty=c(2,1,2),col=c('gray','black','gray'))
t.bg.rate<-sum(sim$R*wgts)
points(sim$R-t.bg.rate,col='blue',pch=16)
abline(h=0,col='red',lty=2)

library(ape)
plot(tree,edge.width=10,cex=0.8)
tmp<-list(R=apply(fit%means%1:nrow(tree$edge),1,mean));class(tmp)<-'corateBM'
plot(tmp,tree,phenogram=F,lwd=10,show.tip.label=F,val.range=range(sim$R))
plot(sim,tree,trait=k,phenogram=F,lwd=10,show.tip.label=F,val.range=range(sim$R))
coords<-get("last_plot.phylo",envir=.PlotPhyloEnv)
trans<-fit$post.probs
colvec<-rep('gray',nrow(tree$edge))
colvec[trans<0.05]<-'blue'
colvec[trans>0.95]<-'red'
colvec<-rgb(trans,trans,trans)
segments(y0=coords$yy[as.vector(t(tree$edge))],
         y1=coords$yy[rep(tree$edge[,2],each=2)],
         x0=coords$xx[rep(tree$edge[,1],each=2)],
         x1=coords$xx[as.vector(t(tree$edge))],
         col=rep(colvec,each=2),
         lwd=3)
par(mfrow=c(1,2))
par(mfrow=c(1,1))
tmp$X<-t(matrix(fit%means%tree$tip.label,k,length(tree$tip.label)))
rownames(tmp$X)<-tree$tip.label
plot(tmp,tree,1,val.range=range(sim$R),lwd=2,alpha=0.5,add=T)
plot(sim,tree,1,val.range=range(sim$R),lwd=4)
XX<-tapply(sim$Y[,1],rownames(sim$Y),mean)[tree$tip.label]
names(XX)<-tree$tip.label
phytools::phenogram(drop.tip(tree,names(XX[is.na(XX)])),XX[!is.na(XX)],add=T,ftype='off')

#got a near false positive for edges 59,60,61--problematic because it recovers pretty rapid edges as
#significantly slow at a 0.1 level! Seems related to lack of edges, odd topology, and uncertainty in species
#means among close relatives

sqrt(diag(sim$Ysig2))/
apply(apply(sim$X,2,range),2,diff)
#the species means for tips 91 and 92 were estimated to be closer than they actually are, and the divergence
#between doesn't seem to reflect their spectacularly high rates of evolution in any case. It seems like a
#fluke due to the hihgly stochastic nature of BM...Note that you could argue you gave the model a bit of
#a curve ball--sd of intraspecific variation is ~2-8% that of the entire observed phenotypic range! Given
#4x that number, the values for a tip could've reasonably spanned a range some 20% of the observed
#phenotypic range in some cases

fit2<-fit.corateBM(tree,sim$X,chains=1,refresh=10)
plot(fit2%means%1:nrow(tree$edge)~sim$R)
abline(0,1)
points(apply(fit%means%1:nrow(tree$edge),1,mean)~sim$R,pch=16)
cor(fit2%means%1:nrow(tree$edge),sim$R)
#eliminating intraspecific variation mitigated, but did not completely solve, the issue
#so the source seems to be a mix of difficulty in accounting for intraspecific variation AND the fact that
#those two tips just are weird! ('we observed what was most likely to have happened' philosophy...)

fit3<-fit.corateBM(tree,sim$Y,chains=1,refresh=10)
plot(fit3%means%1:nrow(tree$edge)~sim$R)
abline(0,1)
points(apply(fit%means%1:nrow(tree$edge),1,mean)~sim$R,pch=16)
cor(fit3%means%1:nrow(tree$edge),sim$R)
#alright, it seems fitting with intraspecific variation when you add noise to the dataset still DRASTICALLY
#improves parameter estimation!
hist(fit3%chains%'R_sig2')

saveRDS(fit3,'shit_fit')
saveRDS(fit2,'good_fit')

sim<-gen.corateBM(tree,Rsig2=0.5,X0=rnorm(2),Xsig2=rWishart(1,2,diag(2))[,,1])
fit<-fit.corateBM(tree,sim$X,iter=1000,refresh=10,chains=1)
