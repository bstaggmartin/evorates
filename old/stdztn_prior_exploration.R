tree<-pbtree(n=100,scale=1)
n.sims<-100
rates<-c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)
sig<-matrix(0,ncol=length(rates),nrow=n.sims)
for(i in 1:length(rates)){
  x<-fastBM(tree,sig2=rates[i],nsim=n.sims)
  x<-apply(x,2,function(xx) (xx-mean(xx))/sd(xx))
  #sig[,i]<-sapply(1:n.sims,function(ii) phylosig(tree,x[,ii],method='K'))
  sig[,i]<-sapply(1:n.sims,function(ii) sum(pic(x[,ii],tree,scaled=F)^2)/(length(tree$tip.label)-1))
}

grad<-rainbow(length(rates),alpha=0.25)
plot(density(sig[,1]),xlim=c(0,1),ylim=c(0,2))
for(i in 1:length(rates)){
  polygon(density(sig[,i]),col=grad[i])
}
legend('topright',col=grad,pch=16,legend=as.character(rates))

#when standardizing values and tree height, rate param seems fairly constant for a tree of given size (as you expected); larger trees
#associated with lower rates than high ones, with a 25-tip tree having ~0.6-0.7 avg. and 500-tip tree having a more consistent 0.3.
#there should be some way to derive an empirical expectation, right? Anyways, assuming that the average rate falls somewhere around 0.5
#seems like a safe assumption for vague priors--more informed than 1, anyways!

#it SHOULD be 1, given the behavior of brownian motion--the problem is clustering of data points (due to phylogenetic autocorrelation)
#infalting sd estimates. An alternative strategy is to simply standardize trait values by a 'phylogenetically-corrected' standard dev.
#(e.g., the estimated brownian motion rate!). Then the prior for r0 may simply be centered on 0 with a high level of vagueness say like
#1e6 to 1e-6...

#something odd here about fastBM--accumulated variance higher than expected at low values and lower than expected at high values
#you forgot about the squaring! Of course, duh
ntax<-100
tree<-pbtree(n=ntax,scale=1)
n.sims<-100
rates<-c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)
x<-array(NA,dim=c(ntax,n.sims,length(rates)))
est.rates<-matrix(NA,ncol=length(rates),nrow=n.sims)
for(i in 1:length(rates)){
  x[,,i]<-fastBM(tree,sig2=rates[i],nsim=n.sims)
  x[,,i]<-apply(x[,,i],2,function(xx) (xx-mean(xx))/
                  sqrt(sum(pic(xx,tree)^2)/(ntax-1)))
  est.rates[,i]<-apply(x[,,i],2,function(xx) sum(pic(xx,tree)^2)/(ntax-1))
}
#now it's all ones. Beautiful

iter<-1
plot(density(x[,1,iter]),xlim=c(min(x[,,iter]),max(x[,,iter])),col=rgb(0,0,0,0))
for(i in 1:n.sims){
  polygon(density(x[,i,iter]),border=rgb(0,0,0,0.1))
}

grad<-rainbow(length(rates),alpha=0.25)
plot(density(est.rates[,1]),xlim=c(0,2),ylim=c(0,10))
for(i in 1:length(rates)){
  polygon(density(est.rates[,i]),col=grad[i])
}
legend('topright',col=grad,pch=16,legend=as.character(rates))

#null expectations
##r0: log(0), not sure yet (maybe like sig2=log(100))
##rtrend: log(0), not sure yet (maybe like sig2=log(1000))
##rvar: exp(1)? (see below), not sure yet--log(1000 ought to do it, methinks)
##root: 0, not sure yet (maybe like sig2=6?)

#I'll just look at the distribution of PICs (as estimates of node rates to get the estimated rate of PIC variation accumulation
#and PIC trending under BM...)

ntax<-100
tree<-pbtree(n=ntax,scale=1)
n.sims<-100
rates<-c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3)
x<-array(NA,dim=c(ntax,n.sims,length(rates)))
pics<-array(NA,dim=c(2*ntax-1,n.sims,length(rates)))
get.nodewise.pics<-function(x,tree){
  tmp.pics<-log(pic(x,tree)^2)
  c(tmp.pics[as.character(sapply(1:ntax,function(tip) tree$edge[which(tree$edge[,2]==tip),1]))],tmp.pics)
}
est.rorates<-matrix(NA,ncol=length(rates),nrow=n.sims)
for(i in 1:length(rates)){
  x[,,i]<-fastBM(tree,sig2=rates[i],nsim=n.sims)
  x[,,i]<-apply(x[,,i],2,function(xx) (xx-mean(xx))/
                  sqrt(sum(pic(xx,tree)^2)/(ntax-1)))
  pics[,,i]<-apply(x[,,i],2,get.nodewise.pics,tree=tree)
  est.rorates[,i]<-apply(pics[1:ntax,,i],2,function(xx) sum(pic(xx,tree)^2)/(ntax-1))
}

grad<-rainbow(length(rates),alpha=0.25)
plot(density(log(est.rorates[,1])),xlim=c(0,20),ylim=c(0,1))
for(i in 1:length(rates)){
  polygon(density(log(est.rorates[,i])),col=grad[i])
}
legend('topright',col=grad,pch=16,legend=as.character(rates))
mean(log(est.rorates));sd(log(est.rorates))
#on a log scale, the expectation seems to be around e...I might be making a mistake, but that also could be a mathematical thing...

#playing with proposed prior defaults
#rvar
plot(dnorm(seq(-log(1e6),log(1e6),0.01),exp(1),log(100),log=T)~exp(seq(-log(1e6),log(1e6),0.01)),type='l',xlim=c(0,1e4))
#basically seems to translate into 'we would be surprised by 1800-fold changes in rate; maybe bump up to 1000 to be safe?
     