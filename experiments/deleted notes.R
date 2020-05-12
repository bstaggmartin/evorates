library(phytools)
tree<-pbtree(n=100)
mod.tree<-tree
rates.true<-rgamma(nrow(tree$edge),0.4,0.1)
mod.tree$edge.length<-rates.true*mod.tree$edge.length
BM<-fastBM(mod.tree,internal=T)

#Things you figured out (which were unfortunately deleted)
#pics do provide a nice point estimate of the rates along the edges FOR WHICH THEY ARE DIRECTLY ANCESTRAL, but they need to be
#normalized--you think maybe a heirarchical MCMC can do this
#you can also do some nice visualizations by plotting the true history as points at nodes, and your simulation 'credible intervals' via
#segments from a simplified sim object

x<-BM[1:length(tree$tip.label)]
library(contSimmap)
fastpgram(BM,tree)
pics<-pic(x,tree)
edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])^2
plot(edge.pics~rates.true)

#oof, here's the issue...
test<-rgamma(100,1,1)
test2<-rep(1,100)     
sum(dgamma(test,2,2,log=T))
sum(dgamma(test2,2,2,log=T))
#these functions are greedy--dgamma will not be able to discriminate a repeating series of values from points that actually follow a
#distribution
#hmmmmm...
#that seems to be why the discrete gamma is necessary in these cases...
#an mcmc that works directly on the edges might be best, in such a case
#why not drop autocorrelation for now and try to go back to gamma model, but discretize it this time?

#but wait, this should still work, since you're trying to jointly maximize the likelihood of a particular theta being drawn from
#a gamma distribution AND that it produced the observed PIC when plugged into a normal distirbution
#it will be subject to a lot of error, but really 'smooth' things out, like you would expect.
#gamma pars themselves could be estimated in an expectation-maximization framework by iteratively pushing thetas to maximize joint
#lik acoording to gamma they're drawn from and normal deviate they produce, then refitting the gamma distribution to the 'pushed'
#thetas. Theoretically makes sense to me

lnL<-function(ind.cont){
  foo<-function(theta){
    -sum(dnorm(ind.cont,0,theta,log=T),dgamma(theta,1,1,log=T))
  }
  foo
}
normalized<-rep(NA,length(pics))
for(i in 1:length(pics)){
  tmp.fun<-lnL(pics[i])
  normalized[i]<-optimize(tmp.fun,interval=c(0,20))$minimum
}
names(normalized)<-names(pics)
test<-sapply(1:nrow(tree$edge),function(ii) normalized[as.character(tree$edge[ii,1])])
plot(test~rates.true)
abline(0,1)

sim<-simBM(tree,x,n.sim=100,rates=test)
plot(sim,tree)
sim2<-simBM(tree,x,n.sim=100)
plot(sim2,tree)
sim3<-simBM(tree,x,n.sim=100,rates=rates.true)
plot(sim3,tree)
simp<-simplify.sim(sim,probs=c(0.025,0.975))
simp2<-simplify.sim(sim2,probs=c(0.025,0.975))
simp3<-simplify.sim(sim3,probs=c(0.025,0.975))
plot(simp,tree,lines=F,polyg=T,col='black',alph=0.2)
plot(simp2,tree,lines=F,polyg=T,col='black',alph=0.2)
plot(simp3,tree,lines=F,polyg=T,col='black',alph=0.2)
fastpgram(BM,tree,add=T)
colmap<-evolve.colors(tree)
node.colmap<-get.node.colmap(colmap,tree)
plot(BM[-(1:length(tree$tip.label))],col=node.colmap[-(1:length(tree$tip.label))],pch=16,ylim=c(min(simp$nodes),max(simp$nodes)))
draw.errors<-function(sim,tree,offset,colmap){
  segments(x0=1:tree$Nnode+offset,x1=1:tree$Nnode+offset,
           y0=sim$nodes[-(1:length(tree$tip.label)),1],y1=sim$nodes[-(1:length(tree$tip.label)),2],
           col=colmap[-(1:length(tree$tip.label))])
}
draw.errors(simp,tree,0.25,node.colmap)
draw.errors(simp2,tree,0.5,node.colmap)
draw.errors(simp3,tree,0.75,node.colmap)
#in general, probably improves estimation, but increases 'confidence' around the value of particular nodes wayyyy to much for my
#tastes
#one 'hacky' idea would be limit the lower bound to be above the mean pic??? (median--highly skewed distribution)
# test[test<median(pics^2)]<-median(pics^2)
# plot(test~rates.true)
# abline(0,1)
#didn't really work well; how does it compare to pics?
sim4<-simBM(tree,x,n.sim=100,rates=edge.pics)
simp4<-simplify.sim(sim4,probs=c(0.025,0.975))
plot(BM[-(1:length(tree$tip.label))],col=node.colmap[-(1:length(tree$tip.label))],pch=16,ylim=c(min(simp$nodes),max(simp$nodes)))
draw.errors(simp,tree,0.2,node.colmap)
draw.errors(simp2,tree,0.4,node.colmap)
draw.errors(simp3,tree,0.6,node.colmap)
draw.errors(simp4,tree,0.8,node.colmap)
cust.heatmap<-colorRampPalette(c("skyblue2","cyan","chartreuse2","goldenrod","red"))(100)
rate.colmap<-cust.heatmap[round((log(test)-min(log(test),log(rates.true)))/
                                  (max(log(test),log(rates.true))-min(log(test),log(rates.true)))*99+1)]
plot(sim,tree,colmap=rate.colmap,alph=0.2,polyg=T,lines=F)
true.rate.colmap<-cust.heatmap[round((log(rates.true)-min(log(test),log(rates.true)))/
                                       (max(log(test),log(rates.true))-min(log(test),log(rates.true)))*99+1)]
fastpgram(BM,tree,edge.col=alter.cols(true.rate.colmap,0.2),node.col=alter.cols('black',0),add=T,lwd=2)
def.par<-par(no.readonly=T)
par(mfrow=c(1,2))
plot(tree,edge.color=rate.colmap,edge.width=3,show.tip.label=F)
plot(tree,edge.color=true.rate.colmap,edge.width=3,show.tip.label=F)
par(def.par)
#seems to work--the problem (as you thought it would be) is that the rate parameters of sister branches can't be separately, and
#thus low-rate branches adjacent to high-rate branches are 'lumped' to an average rate. This might be less of a problem with the
#autocorrelation (=multivariate normal) model...
#let's try to make an EM alogrithm next
expectation<-function(dat,pars){
  foo<-function(theta){
    -sum(dnorm(dat,0,theta,log=T),dgamma(theta,pars[1],pars[2],log=T))
  }
  foo
}
maximize<-function(theta){
  foo<-function(pars){
    -sum(dgamma(theta,pars[1],pars[2],log=T))
  }
  foo
}
#initialize
dat<-pics
pars<-rexp(2)
theta<-rep(NA,length(dat))
n.iter<-1000
mat<-matrix(NA,ncol=length(pics)+3,nrow=n.iter)
for(i in 1:n.iter){
  #expectation step
  for(pt in 1:length(pics)){
    tmp.fun<-expectation(dat[pt],pars)
    theta[pt]<-optimize(tmp.fun,interval=c(0,1000))$minimum
  }
  #maxmization step
  pars<-optim(pars,maximize(theta))$par
  mat[i,]<-c(theta,pars[1],pars[2],sum(dgamma(theta,pars[1],pars[2],log=T)+dnorm(dat,0,theta)))
}

#damn, it does get greedy and just try to make a thin gamma distribution with all pics at the same value x.x
#maybe try throwing a prior on it? God, this is getting just more and more wonky
#theoretically, though, a multivariate normal MIGHT work, since it will assume that adjacent branches are autocorrelated, thus
#giving us more power to estimate the continuous latent variables underlying the distribution
tree<-pbtree(n=50)
BM<-gen.corateBM(tree,internal=T,rvar=3)
x<-BM$traits[1:length(tree$tip.label)]
opt.corate<-function(x,tree){
  pics<-pic(x,tree)
  node.heights<-node.depth.edgelength(tree)
  node.heights<-node.heights[-(1:(length(tree$tip.label)+1))]
  phylo.vcv<-vcvPhylo(tree,anc.nodes=T)
  phylo.vcv<-phylo.vcv[-(1:length(tree$tip.label)),-(1:length(tree$tip.label))]
  expectation<-function(dat,pars,nhgts=node.heights,pvcv=phylo.vcv){
    foo<-function(theta){
      -sum(dnorm(dat,0,exp(theta),log=T),dmvn(theta[-1],theta[1]+nhgts*pars[1],pars[2]*pvcv,log=T))
    }
    foo
  }
  maximize<-function(theta,nhgts=node.heights,pvcv=phylo.vcv){
    foo<-function(pars){
      -sum(dmvn(theta[-1],theta[1]+nhgts*pars[1],pars[2]*pvcv,log=T))
    }
    foo
  }
  pars<-c(rnorm(1),rexp(1))
  pars.old<-rep(0,2)
  dat<-pics
  theta<-rnorm(49)
  theta.old<-rep(0,49)
  while(any(round(c(pars-pars.old,theta-theta.old),6)!=0)){
    tmp.fun<-expectation(dat,pars)
    theta.old<-theta
    theta<-optim(theta,tmp.fun)$par
    tmp.fun<-maximize(theta)
    pars.old<-pars
    pars<-optim(pars,tmp.fun)$par #I think the problem is here, because it only considers the thetas, not their relationship to
                                  #observed data (i.e., pics)
    cat(pars,mean(theta),var(theta),'\n')
  }
  list(pars=pars,theta=theta)
}

indices<-sapply(1:tree$Nnode+length(tree$tip.label),function(ii) which(tree$edge[,1]==ii))
comp.true.rates<-sapply(1:ncol(indices),function(ii) mean(BM$edge.rates[indices[,ii]]))

edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])^2
theta<-rowMeans(mat2)
names(theta)<-names(pics)
test<-exp(sapply(1:nrow(tree$edge),function(ii)theta[as.character(tree$edge[ii,1])]))
# plot(edge.pics~BM$edge.rates)
# plot(test~BM$edge.rates)
# abline(0,1)

rate.colmap<-cust.heatmap[round((log(test)-min(log(test),log(BM$edge.rates)))/
                                  (max(log(test),log(BM$edge.rates))-min(log(test),log(BM$edge.rates)))*99+1)]
true.rate.colmap<-cust.heatmap[round((log(BM$edge.rates)-min(log(test),log(BM$edge.rates)))/
                                  (max(log(test),log(BM$edge.rates))-min(log(test),log(BM$edge.rates)))*99+1)]
par(mfrow=c(1,2))
plot(tree,edge.color=rate.colmap,edge.width=4)
plot(tree,edge.color=true.rate.colmap,edge.width=4)
median(exp(theta))
median(pics^2)

#it definitely works!!! It just chooses a local optima each time...but maybe you could use this in a full Bayesian algorithm?
#oh, odd--you get the same behavior where it always says 'early burst' dynamic...maybe you should eliminate rtrend from the model...
#fuck, fixed error, turned out it suffers the same issues with pulling the mvn closer and closer until you get the infinited lik
#effect
par(def.par)
sim<-simBM(tree,x,n.sim=100,rates=test)
plot(sim,tree,col='gray',alph=0.1)
fastpgram(BM$traits,tree,add=T)


expectation<-function(dat,pars,nhgts=node.heights,pvcv=phylo.vcv){
  foo<-function(theta){
    -sum(dnorm(dat,0,exp(theta),log=T),dmvn(theta[-1],theta[1]+nhgts*pars[2],exp(pars[1])*pvcv,log=T))
  }
  foo
}
maximize<-function(dat,nhgts=node.heights,pvcv=phylo.vcv){
  foo<-function(pars){
    -sum(dnorm(dat,0,exp(pars[-(1:2)]),log=T),dmvn(pars[-(1:3)],pars[3]+nhgts*pars[2],exp(pars[1])*pvcv,log=T))
  }
}
# pars[tt$edge[order(tt$edge[,2][tt$edge[,2]>length(tt$tip.label)]),1]-length(tt$tip.label)+3]+
#   tt$edge.length[order(tt$edge[,2][tt$edge[,2]>length(tt$tip.label)])]*pars[2],
# tt$edge.length[order(tt$edge[,2][tt$edge[,2]>length(tt$tip.label)])]*pars[1]))

tmp.fun<-maximize(dat)
mat<-sapply(1:100,function(ii) optim(rnorm(51,0,10),tmp.fun,method='SANN')$par[1:2])
# mat<-matrix(NA,ncol=100,nrow=2)
# for(i in 1:100){
#   root<-rnorm(1,0,3)
#   var<-rlnorm(1,0,3)
#   trend<-rnorm(1,0,3)
#   rand.theta<-rmvn(1,root+node.heights*trend,var*phylo.vcv)
#   rand.theta[rand.theta>10]<-10;rand.theta[rand.theta< -10]<- -10
#   pars<-c(var,trend,root,rand.theta)
#   mat[,i]<-optim(pars,tmp.fun)$par[1:2]
# }

#new idea! Do what you did above, then optimize expectations of thetas under all of them?
start<-proc.time()
mat2<-sapply(1:100,function(ii) optim(rnorm(49,0,5),expectation(dat,mat[,ii]),method='SANN')$par)
proc.time()-start

plot(log(pics^2)~log(comp.true.rates),pch=16,col='red',xlab='true',ylab='predicted',ylim=c(min(log(pics^2),mat2),max(log(pics^2),mat2)))
points(as.vector(mat2)~log(rep(comp.true.rates,100)),pch=16,col=rgb(0,0,0,0.1));abline(0,1)
points(log(pics^2)~log(comp.true.rates),pch=16,col='red')
points(apply(mat2,1,median)~log(comp.true.rates),pch=16,col='blue')


#doesn't seem to help :(
#YES SIMULATED ANNEALING FTW
#okay, you're close. This works, but it's contingent upon getting good pars values, which ISN'T a trivial problem. With your
#method, it tended to always return parameters that are simply a function of the rate params you randomly generated, and it's VERY
#easy to break the dmvn function by giving it a wonky multivariate normal datapoint...
#Okay, the problem was simply drawing a negative rvar value, by guessing the log of rvar rather than rvar itself, the function works
#just fine!!!
logLik<-maximize(dat)
plot(apply(rbind(mat,mat2),2,logLik))
abline(v=46)
plot(mat2[,45]~log(comp.true.rates))
mat[,46]
abline(0,1)
#might it be as simple as eliminating rtrend?

plot(log(pics^2)~log(comp.true.rates),pch=16,col='red',xlab='true',ylab='predicted')





#dump from .Rhistory
plot(lnL(edge.pics,df=seq(0.1,100,0.1))~seq(0.1,100,0.1))
plot(lnL(edge.pics,df=seq(0.1,100,0.1)))
plot(lnL(edge.pics,df=seq(0.1,10,0.1)))
lnL(edge.pics,df=seq(0.1,10,0.1))
#pics can be thought of as realization of pics~chisq(df=rate) #why has no one done anything with this before?
#pics can be thought of as realization of pics~chisq(df=rate) #why has no one done anything with this before?
testy<-pic(test,tree2)
testy<-pic(test,tree2)
testy<-pic(test,tree)
testy<-pic(trest,tree)
plot(edge.pics~rates.true)
# edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
# names(pics)<-NULL
# pics<-c(tip_pics,pics)
# tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])
pics<-pic(x,tree)^2
sim<-simBM(tree,x,n.sim=100,rates=edge.pics)
plot(edge.pics~rates.true)
# edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
# names(pics)<-NULL
# pics<-c(tip_pics,pics)
# tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])
pics<-pic(x,tree)^2
sim<-simBM(tree,x,n.sim=100,rates=edge.pics)
plot(edge.pics~rates.true)
# edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
# names(pics)<-NULL
# pics<-c(tip_pics,pics)
# tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])
pics<-pic(x,tree)^2
plot(edge.pics~rates.true)
# edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
# names(pics)<-NULL
# pics<-c(tip_pics,pics)
# tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])
pics<-pic(x,tree)^2
plot(edge.pics~rates.true)
# edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
# names(pics)<-NULL
# pics<-c(tip_pics,pics)
# tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
edge.pics<-sapply(1:nrow(tree$edge),function(ii) pics[as.character(tree$edge[ii,1])])
pics<-pic(x,tree)^2
plot(x=edge.pics,y=apply(predicted,1,median))
plot(x=edge.pics,y=apply(predicted,1,mean))
plot(edge.pics~apply(predicted,1,mean))
sim<-simBM(tree,x,n.sim=100,rates=edge.pics)
plot(edge.pics~rates.true)
edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
names(pics)<-NULL
pics<-c(tip_pics,pics)
tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
pics<-pic(x,tree)^2
plot(edge.pics~rates.true)
sim<-simBM(tree,x,n.sim=100,rates=edge.pics)
plot(edge.pics~rates.true)
edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
names(pics)<-NULL
pics<-c(tip_pics,pics)
tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
pics<-pic(x,tree)^2
sim2<-simBM(tree,x,n.sim=100,rates=edge.pics)
sim<-simBM(tree,x,n.sim=100,rates=edge.pics)
plot(edge.pics~rates.true)
edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
names(pics)<-NULL
pics<-c(tip_pics,pics)
tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
pics<-pic(x,tree)^2
plot(edge.pics~rates.true)
edge.pics
edge.pics<-apply(cbind(pics[tree$edge[,1]],pics[tree$edge[,2]]),1,mean)
cbind(pics[tree$edge[,1]],pics[tree$edge[,2]])
names(pics)<-NULL
pics<-c(tip_pics,pics)
tip_pics
tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[as.character(tree$edge[which(tree$edge[,2]==ii),1])])
pics[tree$edge[which(tree$edge[,2]==1),1]]
pics[tree$edge[which(tree$edge[,2]==1),1]]
tip_pics
tip_pics<-sapply(1:length(tree$tip.label),function(ii) pics[tree$edge[which(tree$edge[,2]==ii),1]])
pics
pics<-pic(x,tree)^2
plot(BM[-(1:length(tree$tip.label))],col=node.colmap[-(1:length(tree$tip.label))],pch=16,ylim=c(min(simp$nodes),max(simp$nodes)))
segments(x0=1:tree$Nnode,x1=1:tree$Nnode,y0=simp2$nodes[-(1:length(tree$tip.label)),1],y1=simp2$nodes[-(1:length(tree$tip.label)),3],
         col=node.colmap[-(1:length(tree$tip.label))])