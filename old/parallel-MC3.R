rm(list=ls())
library(doParallel)
library(foreach)
library(phytools)
library(mvnfast)
library(abind)
source("utils.R")
source("gen.corateBM-plot.corateBM.R")
source("exp.relaxed.clock.BM.R")
set.seed(12345)

trees<-pbtree(b=1,d=0.8,t=10,nsim=100,extant.only=T)
no.tips<-sapply(trees,function(ii) length(ii$tip.label))
tree<-trees[which(no.tips==max(no.tips))][[1]]
plot(tree)
corateBM<-gen.corateBM(tree,r0=0,rtrend=0,rvar=0.2,internal=T)
plot(corateBM,tree,log=T,lwd=4,phenogram=T)
recon<-corateBM;recon$traits<-recon$traits[1:length(tree$tip.label)];recon$params$internal=F
plot(recon,tree,log=T,show.tip.label=F,edge.width=4,lwd=4,phenogram=T)
plot(corateBM$edge.rates)

lego<-relaxed.clock.BM(tree,corateBM$traits[1:length(tree$tip.label)],
                              n.iter=1e5,thin=100,report.every=100,win.update=50,block.size=10,tune.period=2e3,
                              adapt.decay=0.5,win=rep(0.5,5),targ.accept=rep(0.4,5),B=0.9,tmp.par.mat="try2",try.swap=60)

##Main problems seem to be from not converging quickly enough: 4 ideas  
###MCMCMC dynamics--allow it to at least initially explore the landscape quickly by raising the ratio to a certain power
###How exactly did clads overcome this issue?
###Is it appropriate to start from ML-esque estimates? It worked mostly for small trees...
###Implement BEAST's adaptive multivariate normal proposal mechanism--there's too high a covariance structure between your params...

no.tips
tree<-trees[which(no.tips==44)][[1]]
plot(tree)
corate<-gen.corateBM(tree,r0=0,rtrend=0,rvar=0.3,internal=T)
plot(corate,tree,log=T,show.tip.label=F,lwd=4,phenogram=F)
cl<-makeCluster(2)
registerDoParallel(cl)
MCMC.list<-foreach(powers=c(1,0.5),
                   .packages=c("phytools","mvnfast","abind")
                   ) %dopar% {relaxed.clock.BM(tree,corate$traits[1:length(tree$tip.label)],
                                               n.iter=2e4,thin=25,report.every=1,
                                               win.update=50,block.size=10,tune.period=1e3,
                                               B=powers)}
stopCluster(cl)

par.mat2<-relaxed.clock.BM(tree,corate$traits[1:length(tree$tip.label)],
                           n.iter=1e5,thin=25,report.every=1,win.update=50,block.size=6,tune.period=1e3,
                           adapt.decay=0.5,win=rep(10,5),B=1)


for(ee in 1:(tree$Nnode*2+1)){
  plot(try2[ee+3,],type='l')
  abline(h=corateBM$ln.node.rates[ee])
}
try<-try2[,100:200]
plot(as.vector(try[4:42,])~rep(corateBM$ln.node.rates,each=101))
