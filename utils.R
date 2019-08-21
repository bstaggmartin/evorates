#Draw from a bactrian distribution (see Yang et al. 2013 for details)
##works by drawing from a binomial distribution, then using that draw as an indicator to draw from either a "left" or "right" normal
##distribution
rbac<-function(n,x=0,sd=1,m=0.9){
  mu<-ifelse(rbinom(n,1,0.5)>0.5,x-m*sd,x+m*sd)
  rnorm(n,mu,sqrt(1-m^2)*sd)
}

#Extract edge-wise rates from the output of a relaxed.clock.BM mcmc
##helpful since the autocorrelated rate model output gives (log) node-wise rate parameters, rather than the edge-wise rates (which
##can be calculated by averaging across adjacent nodes)
get.edge.rates<-function(tree,par.mat){
  exp(sapply(1:ncol(par.mat),function(ii) {
    apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),1,mean)
  }))
}
