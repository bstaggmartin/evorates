#Draw from a bactrian distribution (see Yang et al. 2013 for details)
##works by drawing from a binomial distribution, then using that draw as an indicator to draw from either a "left" or "right" normal
##distribution
rbac<-function(n,x=0,sd=1,m=0.9){
  mu<-ifelse(rbinom(n,1,0.5)>0.5,x-m*sd,x+m*sd)
  rnorm(n,mu,sqrt(1-m^2)*sd)
}

#Extract edge-wise rates from the output of a relaxed.clock.BM mcmc
##helpful since the autocorrelated rate model output gives (log) node-wise rate parameters, rather than the edge-wise rates (which
##can be approximated by averaging across adjacent nodes)
get.edge.rates<-function(tree,par.mat){
  exp(sapply(1:ncol(par.mat),function(ii) {
    apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),1,mean)
  }))
}

#Function to sample the indices of vector as "evenly" as possible by halving the vector into finer and finer increments
##I'm hoping this will make simulated BM bridges more noisy, a la true BM, but I still need to test this; note that this
##function is, unfortunately, much slower than sample, but I have yet to think of a more efficient way to do this...
sample.evenly<-function(vec){
  samp<-c(0,round(length(vec)/2),length(vec)+1)
  while(!all(1:length(vec)%in%samp)){
    space<-mean(sort(samp,decreasing=T)[-length(samp)]-sort(samp,decreasing=T)[-1])
    samp<-unique(c(samp,round(sort(samp)[-length(samp)]+space/2)))
  }
  vec[samp[-c(1,3)]]
}
##Update: did not seem to make the bridges more noisy...but keeping the function anyways; I commented out a line of code in the
##'SimBM' function so I can easily switch between random sampling and this type of sampling when simulating BM bridges