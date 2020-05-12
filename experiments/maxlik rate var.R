#function for maximizing the likelihood of pars (a trend and variance associated with the evolution of rates, as well as latent,
#unobserved rate parameters associated with each observed phylogenetically independent contrast [PIC]) given the observed PICs and the
#structure of the tree (described via the phylogenetic variance covariance matrix for internal nodes ONLY [excluding root])
maximize<-function(pics,pvcv=phylo.vcv){
  foo<-function(pars){
    -sum(dnorm(pics,0,exp(pars[-(1:2)]),log=T),dmvn(pars[-(1:3)],pars[3]+diag(pvcv)*pars[2],exp(pars[1])*pvcv,log=T))
  }
}
#function for maximizing the likelihood of theta (latent, unobserved rate parameters associated with each observed PIC) given the
#observed PICs, trend and variance associated with the evolution of rates, and the structure of the tree (described via the phylogenetic
#variance covariance matrix for internal nodes ONLY [excluding root])
expectation<-function(pics,rvar,rtrend,pvcv=phylo.vcv){
  foo<-function(theta){
    -sum(dnorm(pics,0,exp(theta),log=T),dmvn(theta[-1],theta[1]+diag(pvcv)*rtrend,exp(rvar)*pvcv,log=T))
  }
  foo
}
#in general, this is a difficult optimization problem due to the number of parameter, but we make more identifiable by declaring that
#the latent rate parameters are associated with some 
revo<-rnorm(2,0,10)
theta<-rnorm(49,0,10)
while(any(round((revo-revo.old||theta-theta.old),5)!=0)){
  tmp.fun<-maximize(pics,phylo.vcv)
  revo<-optim(c(revo,theta),tmp.fun,method='SANN')$par[1:2]
  tmp.fun<-expectation(pics,revo[1],revo[2],phylo.vcv)
  theta<-optim(theta,tmp.fun,method='SANN')$par
  theta.old<-theta
  revo.old<-revo
}

x<-BM$traits[1:length(tree$tip.label)]

#experiment with standardization
phy<-tree
phy$edge.length<-phy$edge.length/max(node.depth.edgelength(tree))
phy$edge.length
phylo.vcv<-vcvPhylo(phy,anc.nodes=T);phylo.vcv<-phylo.vcv[-(1:length(tree$tip.label)),-(1:length(tree$tip.label))]
phylo.vcv
pics.std<-pic(x,phy)
comp.true.rates.std<-comp.true.rates/max(node.depth.edgelength(tree))/sqrt(mean(pics.std^2))
x.std<-(x-mean(x))/sqrt(mean(pics.std^2))
pics.std<-pic(x.std,phy)
comp.true.rates.std
revo<-rnorm(51,0,10)
while(any(round(revo-revo.old,5)!=0)){
  tmp.fun<-maximize(pics.std,phylo.vcv)
  revo<-optim(revo,tmp.fun,method='SANN')$par
  revo.old<-revo
}
#didn't improve anything, it inflates rvar to be higher since the tree gets shorter but the log relationships between rates/pics
#doesn't really change



revo1<-revo
theta1<-theta
#seems to work out...
#maybe it would be better if we standardized trait values and phylogeny height?
