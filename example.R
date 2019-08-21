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

#Run MCMC
par.mat2<-relaxed.clock.BM(tree,corateBM$traits[1:length(tree$tip.label)],
                          n.iter=5e4,thin=200,report.every=10,win.update=100,block.size=51)
saveRDS(par.mat,"example_mcmc_output")


par(mfrow=c(2,2))
plot(log(apply(edge.rates,1,median))~log(corateBM$edge.rates),pch=19,xlab="true ln(rate)",ylab="ln(rate) from posterior sample",
     ylim=c(min(log(apply(edge.rates,1,quantile,probs=0.025))),max(log(apply(edge.rates,1,quantile,probs=0.975)))),cex=1.5,
     col=rgb(1,0,0,0.5))
segments(x0=log(corateBM$edge.rates),x1=log(corateBM$edge.rates),
         y0=log(apply(edge.rates,1,quantile,probs=0.025)),y1=log(apply(edge.rates,1,quantile,probs=0.975)),lwd=2,
         col=rgb(0,0,0,0.5))
for(i in sample(1:ncol(edge.rates),5)){
  lines(log(edge.rates[,i])[order(corateBM$edge.rates)]~sort(log(corateBM$edge.rates)),col=rgb(0,0,0,0.5))
}
abline(0,1,lwd=4)
plot(corateBM,tree,norm.lb=-5,norm.ub=3,log=T,edge.width=3,lwd=3,
     main="true ln(rates)",phenogram=T)
frame()
med.edge.rates<-apply(edge.rates,1,median)
test<-corateBM
test$edge.rates<-med.edge.rates;test$traits<-test$traits[1:length(tree$tip.label)];test$params$internal=F
plot.corateBM(test,tree,norm.lb=-5,norm.ub=3,log=T,edge.width=3,lwd=3,
              main="median ln(rates) from posterior sample",phenogram=F)
