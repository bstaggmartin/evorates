rm(list=ls())
library(phytools)
library(mvnfast)
library(abind)
source("utils.R")
source("gen.corateBM-plot.corateBM.R")
source("relaxed.clock.BM.R")
set.seed(7996)

#Generate tree with autocorrelated rates
tree<-pbtree(n=100)
corateBM<-gen.corateBM(tree,r0=log(10),rtrend=-0.2,rvar=1,internal=T)

#Explore simulated data
plot(corateBM,tree,log=T,edge.width=4)
plot(corateBM,tree,log=T,phenogram=T,lwd=4)
mean(corateBM$edge.rates);var(corateBM$edge.rates)

#Generate tree with one large rate shift
tree<-pbtree(t=10,d=0.5,extant.only=T)
edge.depth<-sapply(1:(tree$Nnode*2+1),function(ii) length(getDescendants(tree,ii)))[tree$edge[,2]]
test.tree<-tree;test.tree$edge.length<-tree$edge.length*edge.depth;plot(test.tree)
wgt.el<-tree$edge.length*edge.depth
shift.pt<-runif(1,0,sum(wgt.el))
shift.edge<-min(which(cumsum(wgt.el)>shift.pt))
if(shift.edge!=1){
  shift.pt<-(shift.pt-cumsum(wgt.el)[shift.edge-1])/edge.depth[shift.edge]
}else{
  shift.pt<-shift.pt/edge.depth[shift.edge]
}
plot(tree)
env<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(x=env$xx[env$edge[shift.edge,1]]+shift.pt,y=env$yy[env$edge[shift.edge,2]],pch=16,col='red')


shifted.edges<-which(tree$edge[,2]%in%getDescendants(tree,tree$edge[shift.edge,2]))
sig<-c(1,10)
rate<-ifelse(1:nrow(tree$edge)%in%shifted.edges,sig[2],sig[1])
rate[shift.edge]<-(shift.pt*sig[1]+(tree$edge.length[shift.edge]-shift.pt)*sig[2])/tree$edge.length[shift.edge]
trans.tree<-tree
trans.tree$edge.length<-trans.tree$edge.length*rate
plot(trans.tree)
x<-fastBM(trans.tree,internal=T)
phenogram(tree,x,ftype='off')
mean(rate);var(rate)

#Run MCMC
par.mat<-relaxed.clock.BM(tree,x[1:length(tree$tip.label)],
                          n.iter=1e4,thin=100,report.every=1,win.update=100,block.size=50,tune.period=1e3)
saveRDS(par.mat,"example_mcmc_output_1lrgshift")
par.mat<-readRDS("example_mcmc_output")

pdf("example_mcmc_output.pdf",width=10,height=10)
edge.rates<-get.edge.rates(tree,par.mat[,seq(100,ncol(par.mat),2)])
par(mfrow=c(2,2))
plot(log(apply(edge.rates,1,median))~log(corateBM$edge.rates),pch=19,xlab="true ln(rate)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=rgb(1,0,0,0.5))
segments(x0=log(corateBM$edge.rates),x1=log(corateBM$edge.rates),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
for(i in sample(1:ncol(edge.rates),5)){
  lines(log(edge.rates[,i])[order(corateBM$edge.rates)]~sort(log(corateBM$edge.rates)),col=rgb(0,0,0,0))
}
abline(0,1,lwd=4)
plot(corateBM,tree,norm.lb=-5,norm.ub=3,log=T,edge.width=3,lwd=3,
     main="true ln(rates)",phenogram=F)
frame()
med.edge.rates<-apply(edge.rates,1,median)
test<-corateBM
test$edge.rates<-med.edge.rates;test$traits<-test$traits[1:length(tree$tip.label)];test$params$internal=F
plot(test,tree,norm.lb=-5,norm.ub=3,log=T,edge.width=3,lwd=3,
     main="median ln(rates) from posterior sample",phenogram=F)
dev.off()


pdf("example_mcmc_output_1lrgshift.pdf",width=10,height=10)
edge.rates<-get.edge.rates(tree,par.mat[,seq(20,ncol(par.mat),1)])
par(mfrow=c(2,2))
plot(log(apply(edge.rates,1,median))~log(rate),pch=19,xlab="true ln(rate)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=rgb(1,0,0,0.5))
segments(x0=log(rate),x1=log(rate),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
for(i in sample(1:ncol(edge.rates),5)){
  lines(log(edge.rates[,i])[order(rate)]~sort(log(rate)),col=rgb(0,0,0,0))
}
abline(0,1,lwd=4)
blah<-list(edge.rates=rate)
class(blah)<-"corateBM"
plot(blah,tree,norm.lb=-4,norm.ub=4,log=T,edge.width=3,lwd=3,
     main="true ln(rates)",phenogram=F)
frame()
med.edge.rates<-apply(edge.rates,1,median)
test<-blah
test$edge.rates<-med.edge.rates;test$traits<-test$traits[1:length(tree$tip.label)];test$params$internal=F
plot(test,tree,norm.lb=-4,norm.ub=4,log=T,edge.width=3,lwd=3,
     main="median ln(rates) from posterior sample",phenogram=F)
dev.off()