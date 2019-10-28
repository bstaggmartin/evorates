rm(list=ls())
library(phytools)
library(mvnfast)
library(abind)
source("utils.R")
source("gen.corateBM-plot.corateBM.R")
source("exp.relaxed.clock.BM.R")
set.seed(7996)

#Generate tree with autocorrelated rates
tree<-pbtree(n=100,scale=1)
corate<-gen.corateBM(tree,rvar=10)
x<-(corate$traits-mean(corate$traits))/sqrt(sum(pic(corate$traits,tree)^2)/19)

par.mat3<-relaxed.clock.BM(tree,x,n.chains=4,try.swap=50,report.chains=1,dT=c(0,0.02,0.04,0.06),report.every=10,
                           inits='informed',win.update=100,block.size=15,
                           win=matrix(3,nrow=5,ncol=4),targ.accept=matrix(0.2,nrow=5,ncol=4))

plot(0,col='white',ylim=c(-1.5,3.5),xlim=c(0,3.1))
for(e in 1:nrow(tree$edge)){
  polygon(c(quantile(burn[tree$edge[e,2]+3,],probs=c(0.025,0.975)),
            quantile(burn[tree$edge[e,1]+3,],probs=c(0.975,0.025)))~
            c(rep(node.depth.edgelength(tree)[tree$edge[e,2]],2),rep(node.depth.edgelength(tree)[tree$edge[e,1]],2)),
          border=NA,col=rgb(0,0,0,0.1))
  segments(x0=node.depth.edgelength(tree)[tree$edge[e,1]],x1=node.depth.edgelength(tree)[tree$edge[e,2]],
          y0=corate$ln.node.rates[tree$edge[e,1]],y1=corate$ln.node.rates[tree$edge[e,2]])
}

#Explore simulated data
plot(corate,tree,log=T,lwd=4)
plot(corate,tree,log=T,phenogram=T,lwd=4)
mean(corate$edge.rates);var(corate$edge.rates)

#Generate tree with one large rate shift
tree<-pbtree(t=6,d=0.5,extant.only=T)
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
                          n.iter=1e4,thin=100,report.every=10,win.update=50,block.size=10,tune.period=1e3)
saveRDS(par.mat,"example_mcmc_output_1lrgshift2")
par.mat<-readRDS("example_mcmc_output")

pdf("example_mcmc_output.pdf",width=10,height=10)
edge.rates<-get.edge.rates(tree,try)
par(mfrow=c(2,2))
plot(log(apply(edge.rates,1,mean))~log(corateBM$edge.rates),pch=19,xlab="true ln(rate)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=rgb(1,0,0,0.5))
segments(x0=log(corateBM$edge.rates),x1=log(corateBM$edge.rates),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
for(i in sample(1:ncol(edge.rates),5)){
  lines(log(edge.rates[,i])[order(corateBM$edge.rates)]~sort(log(corateBM$edge.rates)),col=rgb(0,0,0,0))
}
abline(0,1,lwd=4)
plot(corateBM,tree,norm.lb=-4,norm.ub=4,log=T,lwd=3,
     main="true ln(rates)",phenogram=F)
frame()
med.edge.rates<-apply(edge.rates,1,median)
test<-corateBM
test$edge.rates<-med.edge.rates;test$traits<-test$traits[1:length(tree$tip.label)];test$params$internal=F
plot(test,tree,norm.lb=-4,norm.ub=4,log=T,lwd=3,
     main="median ln(rates) from posterior sample",phenogram=F)
dev.off()


pdf("example_mcmc_output_1lrgshift2.pdf",width=10,height=10)
edge.rates<-get.edge.rates(tree,par.mat[,seq(50,ncol(par.mat),1)])
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
plot(blah,tree,norm.lb=-8,norm.ub=6,log=T,edge.width=3,lwd=3,
     main="true ln(rates)",phenogram=F)
frame()
med.edge.rates<-apply(edge.rates,1,quantile,probs=0.5)
test<-blah
test$edge.rates<-med.edge.rates
plot(test,tree,norm.lb=-8,norm.ub=6,log=T,edge.width=3,lwd=3,
     main="median ln(rates) from posterior sample",phenogram=F)
dev.off()

act.change<-rep(NA,nrow(tree$edge))
for(e in 1:nrow(tree$edge)){
  act.change[e]<-(x[tree$edge[e,2]]-x[tree$edge[e,1]])^2/tree$edge.length[e]
}
plot(log(apply(edge.rates,1,median))~log(act.change),pch=19,xlab="ln(net change)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=ifelse(rate==10,rgb(1,0,0,0.5),rgb(0,0,1,0.5)))
segments(x0=log(act.change),x1=log(act.change),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
abline(0,1,lwd=4);abline(h=0,col="blue",lwd=4);abline(h=log(10),col="red",lwd=4)
text(x=log(act.change),y=log(apply(edge.rates,1,median)),labels=1:nrow(tree$edge),adj=-0.5)
for(i in seq(2.5,3,0.01)){
  points(x=log(act.change)[86:87],y=log(apply(edge.rates,1,median))[86:87],cex=i)
}
plot(tree);add.arrow(tree,tree$edge[86,1])
test2<-blah
test2$edge.rates<-act.change
plot(test2,tree,norm.lb=lb,norm.ub=ub,log=T,edge.width=3,lwd=3,
     "actual ln(net change)",phenogram=F)


pred.x<-c(x[1:length(tree$tip.label)],fastAnc(tree,x[1:length(tree$tip.label)]))
pred.change<-rep(NA,nrow(tree$edge))
for(e in 1:nrow(tree$edge)){
  pred.change[e]<-(pred.x[tree$edge[e,2]]-pred.x[tree$edge[e,1]])^2/tree$edge.length[e]
}
plot(log(apply(edge.rates,1,median))~log(pred.change),pch=19,xlab="ln(predicted net change)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=ifelse(rate==10,rgb(1,0,0,0.5),rgb(0,0,1,0.5)))
segments(x0=log(pred.change),x1=log(pred.change),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
abline(0,1,lwd=4);abline(h=0,col="blue",lwd=4);abline(h=log(10),col="red",lwd=4)
text(x=log(pred.change),y=log(apply(edge.rates,1,median)),labels=1:nrow(tree$edge),adj=-0.5)
for(i in seq(2.5,3,0.02)){
  points(x=log(pred.change)[86:87],y=log(apply(edge.rates,1,median))[86:87],cex=i)
}
test3<-blah
test3$edge.rates<-pred.change
plot(test3,tree,norm.lb=lb,norm.ub=ub,log=T,edge.width=3,lwd=3,
     main="predicted ln(net change)",phenogram=F)

par(mfrow=c(2,2),mar=c(1,1,3,1))
lb<-log(min(pred.change,act.change,med.edge.rates,rate));ub<-log(max(pred.change,act.change,med.edge.rates,rate))
plot(test,tree,norm.lb=lb,norm.ub=ub,log=T,edge.width=3,lwd=3,
     main="median ln(rates) from posterior sample",phenogram=F)
plot(blah,tree,norm.lb=lb,norm.ub=ub,log=T,edge.width=3,lwd=3,
     main="true ln(rates)",phenogram=F)
plot(test3,tree,norm.lb=lb,norm.ub=ub,log=T,edge.width=3,lwd=3,
     main="BM-based ML ln(net change)",phenogram=F)
plot(test2,tree,norm.lb=lb,norm.ub=ub,log=T,edge.width=3,lwd=3,
     main="true ln(net change)",phenogram=F)
#Model output seems (thankfully) relatively insensitive to stochastic variation in actual or predicted rates of change, unless
#these rates of change are a) really anomalous (in the cases of edges 86, 87, and 24) or b) happen to be close to each other on
#the tree (in the cases of edges 94, 95, and 58-61) (since the rate prior is, after all, based on the idea that rates
#phylogenetically covary)