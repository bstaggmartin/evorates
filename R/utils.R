#Draw from a bactrian distribution (see Yang et al. 2013 for details)
##works by drawing from a binomial distribution, then using that draw as an indicator to draw from either a "left" or "right" normal
##distribution
#' @export
rbac<-function(n,x=0,sd=1,m=0.9){
  mu<-ifelse(rbinom(n,1,0.5)>0.5,x-m*sd,x+m*sd)
  rnorm(n,mu,sqrt(1-m^2)*sd)
}

#Extract edge-wise rates from the output of a relaxed.clock.BM mcmc
##helpful since the autocorrelated rate model output gives (log) node-wise rate parameters, rather than the edge-wise rates (which
##can be approximated by averaging across adjacent nodes)
#' @export
get.edge.rates<-function(tree,par.mat){
  exp(sapply(1:ncol(par.mat),function(ii) {
    apply(cbind(par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,1]],par.mat[4:(nrow(par.mat)-1),ii][tree$edge[,2]]),1,mean)
  }))
}

#Function to sample the indices of vector as "evenly" as possible by halving the vector into finer and finer increments
##I'm hoping this will make simulated BM bridges more noisy, a la true BM, but I still need to test this; note that this
##function is, unfortunately, much slower than sample, but I have yet to think of a more efficient way to do this...
#' @export
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

#fast phenogram function--I just find myself using this over and over again. It's not as 'safe' as Liam's function and doesn't do
#tip labels, but it's a nice in how simple and flexible it is. Also works with the coloring system I created.
#' @export
fastpgram<-function(x,tree,add=F,node.col='black',edge.col='black',...){
  if(length(x)==length(tree$tip.label)){
    x<-c(x,fastAnc(tree,x))
  }
  if(add){
    points(x~node.depth.edgelength(tree),col=node.col,...)
  }else{
    plot(x~node.depth.edgelength(tree),col=node.col,...)
  }
  segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
           y0=x[tree$edge[,1]],y1=x[tree$edge[,2]],col=edge.col,...)
}

#get coordinates of nodes from phylogenetic plotting environment
#' @export
get.node.coords<-function(){
  tmp<-get("last_plot.phylo",envir=.PlotPhyloEnv)[c('xx','yy')]
  cbind(tmp[[1]],tmp[[2]])
}

#simplify simBM object to quantiles of simulations
#' @export
simplify.sim<-function(sim,probs=c(0.025,0.5,0.975)){
  sim$nodes<-t(apply(sim$nodes,1,quantile,probs=probs))
  if(!is.null(sim$edges)){
    sim$edges<-aperm(apply(sim$edges,c(1,2),quantile,probs=probs,na.rm=T),c(2,3,1))
  }
  sim
}
