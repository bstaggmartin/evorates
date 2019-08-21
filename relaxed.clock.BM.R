relaxed.clock.BM<-function(tree,x,n.iter=1e5,thin=100,report.every=100,
                           r0.prior.var=log(1e6),rtrend.prior.var=log(1e6),rvar.prior.var=log(1e6),root.prior.var=1e6,
                           block.size=5,win=rep(1,5),prop.mix=0.9,
                           tune.period=1e4,win.update=50,adapt.decay=0.5,targ.accept=rep(0.234,5)){
  
  
  ##CREATING SOME HELPFUL FUNCTIONS##
  #Sums edge-wise var-cov matrices, stored in the C array, scaling each matrix according to vector of scalars, rates.
  trans.C<-function(C,rates){
    apply(array(sapply(1:length(rates),function(e) C[,,e]*rates[e]),dim=dim(C)),c(1,2),sum)
  }
  
  #Propose a new set of values by adding a bactrian distributed variate to currently accepted set of values, pars. Use
  #update.indices to indicate which parameters are having new values proposed. Note that the rvar parameter, pars[2], is handled
  #specially since it has new values proposed via sliding scale, rather than sliding window, moves.
  prop.vals<-function(pars,update.indices,win,tt=tree){
    indicator<-ifelse((1:length(pars))%in%update.indices,1,0)
    indicator[2]<-0
    pars.prime<-rbac(length(pars),pars,
                     indicator*c(win[1:3],rep(win[5],length(tt$tip.label)),win[4],rep(win[5],tt$Nnode-1))/2,prop.mix)
    if(2%in%update.indices){
      pars.prime[2]<-exp(rbac(1,0,win[2]/2,prop.mix))*pars[2]
    }
    pars.prime
  }
  
  #Get the (log) likelihood & prior ratio between a proposed set of parameters, pars.prime, and the currently accepted set of
  #parameters, pars. Use update.indices to indicate which parameters are having new values proposed, and the function will
  #automatically simplify likelihood ratio calculations.
  get.R<-function(pars,pars.prime,update.indices,
                  xx=x,CC=C,tt=tree,nn=adj.nodes,ee=adj.edges,
                  r0.pri=r0.prior.var,rtrend.pri=rtrend.prior.var,rvar.pri=rvar.prior.var,root.pri=root.prior.var){
    #update "simple priors"
    if(1%in%update.indices){
      rtrend.prior<-dnorm(pars.prime[1],0,rtrend.pri,log=T)-dnorm(pars[1],0,rtrend.pri,log=T)
    }else{rtrend.prior<-0}
    ##note that a scale factor is calculated and added to the prior ratio for rvar; this is because this parameter is updated
    ##according to a scaling, rather than sliding, move, so the proposal distribution is asymmetric and the hastings ratio is
    ##not 1
    if(2%in%update.indices){
      sf<-pars.prime[2]/pars[2]
      rvar.prior<-dexp(pars.prime[2],1/sqrt(rvar.pri),log=T)-dexp(pars[2],1/sqrt(rvar.pri),log=T)+log(sf)
    }else{rvar.prior<-0}
    if(3%in%update.indices){
      root.prior<-dnorm(pars.prime[3],0,r0.pri,log=T)-dnorm(pars[3],0,r0.pri,log=T)
    }else{root.prior<-0}
    if((length(tt$tip.label)+4)%in%update.indices){
      r0.prior<-dnorm(pars.prime[length(tt$tip.label)+4],0,r0.pri,log=T)-dnorm(pars[length(tt$tip.label)+4],0,r0.pri,log=T)
    }else{r0.prior<-0}
    #update "complex" rate prior
    ##if either rtrend or rvar is updated, the likelihood of each node rate parameter is affected and must be re-calculated
    if(any((1:2)%in%update.indices)){
      rate.prior<-dmvn(pars.prime[4:length(pars.prime)][-(length(tt$tip.label)+1)],
                       mu=pars.prime[length(tt$tip.label)+4]+
                         pars.prime[1]*node.depth.edgelength(tt)[-(length(tt$tip.label)+1)],
                       sigma=trans.C(CC[-(length(tt$tip.label)+1),-(length(tt$tip.label)+1),],
                                     rep(pars.prime[2],dim(CC)[3])),log=T)-
        dmvn(pars[4:length(pars.prime)][-(length(tt$tip.label)+1)],
             mu=pars[length(tt$tip.label)+4]+
               pars[1]*node.depth.edgelength(tt)[-(length(tt$tip.label)+1)],
             sigma=trans.C(CC[-(length(tt$tip.label)+1),-(length(tt$tip.label)+1),],
                           rep(pars[2],dim(CC)[3])),log=T)
    ##if neither rtrend nor rvar is updated, but some node rate parameters are updated, then the likelihood of the rate
    ##parameters for these nodes and adjacent nodes must be re-calculated
    }else if(any((4:length(pars))%in%update.indices)){
      nodes<-unique(as.vector(nn[(update.indices-3)[update.indices>3],]))
      nodes<-nodes[(nodes!=length(tt$tip.label)+1)&(!(is.na(nodes)))]
      rate.prior<-dmvn(pars.prime[nodes+3],
                       mu=pars.prime[length(tt$tip.label)+4]+
                         pars.prime[1]*node.depth.edgelength(tt)[nodes],
                       sigma=trans.C(CC[nodes,nodes,],
                                     rep(pars.prime[2],dim(CC)[3])),log=T)-
        dmvn(pars[nodes+3],
             mu=pars[length(tt$tip.label)+4]+
               pars[1]*node.depth.edgelength(tt)[nodes],
             sigma=trans.C(CC[nodes,nodes,],
                           rep(pars[2],dim(CC)[3])),log=T)
    }else{rate.prior<-0}
    #update likelihood
    ##note that, since ancestral trait values are not being explicitly estimated, these multivariate normal likelihoods cannot
    ##be simplified in the same way they are when updating the rate prior
    if(any((3:length(pars))%in%update.indices)){
      ##create a vector to store edge-wise rates by averaging the node-wise rate parameters at adjacent nodes; if any node rate
      ##parameters are being updated, update the edge-wise rates for each edge adjacent to the effected nodes
      edge.rates.prime<-edge.rates<-exp(apply(cbind(pars[4:length(pars)][tt$edge[,1]],pars[4:length(pars)][tt$edge[,2]]),1,mean))
      if(any((4:length(pars))%in%update.indices)){
        edges<-as.vector(ee[(update.indices-3)[update.indices>3],])
        edges<-edges[!(is.na(edges))]
        edge.rates.prime[edges]<-exp(sapply(edges,function(e) mean(c(pars.prime[tt$edge[e,1]+3],pars.prime[tt$edge[e,2]+3]))))
      }
      log.lik<-dmvn(xx,
                    mu=rep(pars.prime[3],length(xx)),
                    sigma=trans.C(CC[1:length(xx),1:length(xx),],edge.rates.prime),log=T)-
        dmvn(xx,
             mu=rep(pars[3],length(xx)),
             sigma=trans.C(CC[1:length(xx),1:length(xx),],edge.rates),log=T)
    }else{log.lik<-0}
    #return sum of all (log) prior/likelihood ratios
    sum(rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior,log.lik)
  }
  
  #Get the (log) posterior probability of a currently accepted set of parameters, pars
  get.post<-function(pars,xx=x,CC=C,tt=tree,
                     r0.pri=r0.prior.var,rtrend.pri=rtrend.prior.var,rvar.pri=rvar.prior.var,root.pri=root.prior.var){
    #"simple" priors
    rtrend.prior<-dnorm(pars[1],0,rtrend.pri,log=T)
    rvar.prior<-dexp(pars[2],1/sqrt(rvar.pri),log=T)
    root.prior<-dnorm(pars[3],0,root.pri,log=T)
    r0.prior<-dnorm(pars[length(tt$tip.label)+4],0,r0.pri,log=T)
    #"complex" prior
    rate.prior<-dmvn(pars[4:length(pars)][-(length(tt$tip.label)+1)],
                     mu=pars[length(tt$tip.label)+4]+pars[1]*node.depth.edgelength(tt)[-(length(tt$tip.label)+1)],
                     sigma=trans.C(CC[-(length(tt$tip.label)+1),-(length(tt$tip.label)+1),],rep(pars[2],dim(CC)[3])),log=T)
    #getting edge-wise rates and calculating likelihood
    edge.rates<-exp(apply(cbind(pars[4:length(pars)][tt$edge[,1]],pars[4:length(pars)][tt$edge[,2]]),1,mean))
    log.lik<-dmvn(xx,
                  mu=rep(pars[3],length(xx)),
                  sigma=trans.C(CC[1:length(xx),1:length(xx),],edge.rates),log=T)
    #return the sum
    sum(rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior,log.lik)
  }
  
  
  ##EXTRACTING RELEVANT INFO FROM TREE##
  #Creating array of edge-wise var-cov matrices for given tree
  C<-array(0,dim=c(rep(tree$Nnode*2+1,2),nrow(tree$edge)))
  for(e in 1:nrow(tree$edge)){
    D<-unique(c(tree$edge[e,2],getDescendants(tree,tree$edge[e,2])))
    C[as.matrix(cbind(expand.grid(D,D),rep(e,length(D)^2)))]<-tree$edge.length[e]
  }
  #Creating look-up tables for adjacent nodes/edges for each node in given tree
  P.edges<-sapply(1:(tree$Nnode*2+1),function(n) which(tree$edge[,2]==n))
  P.edges[lengths(P.edges)<1]<-NA;P.edges<-unlist(P.edges)
  P.nodes<-tree$edge[P.edges,1]
  D.edges<-sapply(1:(tree$Nnode*2+1),function(n) which(tree$edge[,1]==n))
  D.edges[lengths(D.edges)<1]<-list(c(NA,NA));D.edges<-matrix(unlist(D.edges),ncol=2,byrow=T)
  D.nodes<-matrix(tree$edge[as.vector(D.edges),2],ncol=2)
  adj.nodes<-cbind(1:(tree$Nnode*2+1),P.nodes,D.nodes);adj.edges<-cbind(P.edges,D.edges)
  #Using the phylogenetic-independent contrasts to choose somewhat informed starting parameter values for the MCMC chain, with
  #rtrend being set to 0, rvar being set to the variance in (log) pic's standardized by total height of the tree, root being set
  #to the mean trait value, and node rate parameters being set to their associated (log) pic's. Tip rate parameters are just set
  #to the (log) pic's of their ancestral nodes. In future versions of this function, there will be a way to specify user-defined
  #starting values, or just randomly sample from prior distributions to allow for better assessment of convergence.
  pics<-pic(x,tree)^2
  pars<-c(0,var(log(pics))/node.depth.edgelength(tree)[1],mean(x),
          c(log(pics)[as.character(tree$edge[which(tree$edge[,2]%in%(1:length(tree$tip.label))),1])],log(pics)))
  names(pars)<-NULL
  
  
  ##INITIALIZING MCMC ALGORITHM##
  n.iter<-n.iter-n.iter%%thin
  par.mat<-matrix(NA,nrow=5+tree$Nnode*2,ncol=n.iter/thin) 
  n.blocks<-floor(length(pars)/block.size)
  if(length(pars)%%block.size!=0){
    block.sizes<-c(rep(block.size,n.blocks),length(pars)%%block.size)
    n.blocks<-n.blocks+1
  }else{
    block.sizes<-rep(block.size,n.blocks)
  }
  accept<-rep(0,5)
  prop<-rep(0,5)
  for(i in 2:n.iter){
    blocks<-split(sample(1:length(pars)),rep(1:n.blocks,block.sizes))
    for(j in 1:length(blocks)){
      pars.prime<-prop.vals(pars,blocks[[j]],win)
      R<-get.R(pars,pars.prime,blocks[[j]])
      n.par.updates<-c(ifelse(c(1:3,length(tree$tip.label)+4)%in%blocks[[j]],1,0),
                       sum(blocks[[j]]%in%c(4:(length(tree$tip.label)+3),(length(tree$tip.label)+5):length(pars))))
      if(exp(R)>runif(1,0,1)){
        pars<-pars.prime
        accept<-accept+n.par.updates
      }
      prop<-prop+n.par.updates
    }
    if(i%%report.every==0){
      tmp.edge.rates<-exp(apply(cbind(pars[4:length(pars)][tree$edge[,1]],pars[4:length(pars)][tree$edge[,2]]),1,mean))
      cat(i,pars[length(tree$tip.label)+4],pars[1:3],mean(tmp.edge.rates),var(tmp.edge.rates),get.post(pars),accept/prop,"\n")
    }
    if(i%%thin==0){
      par.mat[,i/thin]<-c(pars,get.post(pars))
    }
    if(i<=tune.period&i%%win.update==0){
      n.update<-i/win.update
      cat("updating proposal windows...\n")
      old.win<-win
      win<-exp(log(win)+ifelse(accept/prop-targ.accept>0,1,-1)*min(0.1,n.update^-adapt.decay))
      cat("r0: ",old.win[4]," -> ",win[4],"\n","rtrend: ",old.win[1]," -> ",win[1],"\n","rvar: ",old.win[2]," -> ",win[2],"\n",
          "root: ",old.win[3]," -> ",win[3],"\n","rate: ",old.win[5]," -> ",win[5],"\n")
      accept<-rep(0,5)
      prop<-rep(0,5)
    }
  }
  
  
  return(par.mat)
}
