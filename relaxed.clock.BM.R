relaxed.clock.BM<-function(tree,x,n.iter=1e5,
                           r0.prior.var=1e6,rtrend.prior.var=1e6,rvar.prior.var=1e6,root.prior.var=1e6,
                           block.size=5,win=rep(1,5),tune.period=1e4,thin=100,win.update=50,adapt.decay=0.5,report.every=100){
  
  
  ##CREATING SOME HELPFUL FUNCTIONS##
  #sums edge-wise var-cov matrices, scaling each matrix according to vector of scalars
  trans.C<-function(C,rates){
    return(apply(array(sapply(1:length(rates),function(e) C[,,e]*rates[e]),dim=dim(C)),c(1,2),sum))
  }
  
  get.post<-function(pars,xx=x,CC=C,tt=tree,
                     r0.pri=r0.prior.var,rtrend.pri=rtrend.prior.var,rvar.pri=rvar.prior.var,root.pri=root.prior.var){
    edge.rates<-exp(apply(cbind(pars[4:length(pars)][tt$edge[,1]],pars[4:length(pars)][tt$edge[,2]]),1,mean))
    log.lik<-dmvnorm(xx,
                     mean=rep(pars[3],length(xx)),
                     sigma=trans.C(CC[1:length(xx),1:length(xx),],edge.rates),log=T)
    rtrend.prior<-dnorm(pars[1],0,rtrend.pri,log=T)
    rvar.prior<-dexp(pars[2],1/sqrt(rvar.pri),log=T)
    root.prior<-dnorm(pars[3],0,root.pri,log=T)
    r0.prior<-dnorm(pars[length(tt$tip.label)+4],0,r0.pri,log=T)
    rate.prior<-dmvnorm(pars[4:length(pars)][-(length(tt$tip.label)+1)],
                        mean=pars[length(tt$tip.label)+4]+pars[1]*node.depth.edgelength(tt)[-(length(tt$tip.label)+1)],
                        sigma=trans.C(CC[-(length(tt$tip.label)+1),-(length(tt$tip.label)+1),],rep(pars[2],dim(CC)[3])),log=T)
    return(sum(log.lik,rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior))
  }
  
  
  ##EXTRACTING RELEVANT INFO FROM TREE##
  C<-array(0,dim=c(rep(tree$Nnode*2+1,2),nrow(tree$edge)))
  for(e in 1:nrow(tree$edge)){
    D<-unique(c(tree$edge[e,2],getDescendants(tree,tree$edge[e,2])))
    C[as.matrix(cbind(expand.grid(D,D),rep(e,length(D)^2)))]<-tree$edge.length[e]
  }
  
  pics<-pic(x,tree)^2
  
  pars<-c(0,var(log(pics))/node.depth.edgelength(tree)[1],mean(x),
          fastBM(tree,a=mean(log(pics)),sig2=var(log(pics))/node.depth.edgelength(tree)[1],internal=T))
  names(pars)<-NULL
  edge.rates<-exp(apply(cbind(pars[4:length(pars)][tree$edge[,1]],pars[4:length(pars)][tree$edge[,2]]),1,mean))
  
  
  ##INITIALIZING MCMC ALGORITHM##
  n.iter<-n.iter-n.iter%%thin
  par.mat<-matrix(NA,nrow=5+tree$Nnode*2,ncol=n.iter/thin) 
  n.blocks<-floor(tree$Nnode*2/block.size)
  if(nrow(tree$edge)%%block.size!=0){
    block.sizes<-c(rep(block.size,n.blocks),(tree$Nnode*2)%%block.size)
    n.blocks<-n.blocks+1
  }else{
    block.sizes<-rep(block.size,n.blocks)
  }
  accept<-rep(0,5)
  prop<-rep(0,5)
  
  for(i in 2:n.iter){
    
    #rtrend update#
    theta.prime<-runif(1,pars[1]-win[1]/2,pars[1]+win[1]/2)
    R<-dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
               mean=pars[length(tree$tip.label)+4]+theta.prime*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)-
      dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
              mean=pars[length(tree$tip.label)+4]+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
              sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)+
      dnorm(theta.prime,0,rtrend.prior.var,log=T)-dnorm(pars[1],0,rtrend.prior.var,log=T)
    if(exp(R)>runif(1,0,1)){
      pars[1]<-theta.prime
      accept[1]<-accept[1]+1
    }
    prop[1]<-prop[1]+1
    
    #rvarupdate#
    sf<-exp(runif(1,-win[2]/2,win[2]/2));theta.prime<-sf*pars[2]
    R<-dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
               mean=pars[length(tree$tip.label)+4]+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(theta.prime,dim(C)[3])),log=T)-
      dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
              mean=pars[length(tree$tip.label)+4]+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
              sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)+
      dexp(theta.prime,1/sqrt(rvar.prior.var),log=T)-dexp(pars[2],1/sqrt(rvar.prior.var),log=T)+
      log(sf)
    if(exp(R)>runif(1,0,1)){
      pars[2]<-theta.prime
      accept[2]<-accept[2]+1
    }
    prop[2]<-prop[2]+1
    
    #root update#
    theta.prime<-runif(1,pars[3]-win[3]/2,pars[3]+win[3]/2)
    R<-dmvnorm(x,
               mean=rep(theta.prime,length(x)),
               sigma=trans.C(C[1:length(x),1:length(x),],edge.rates),log=T)-
      dmvnorm(x,
              mean=rep(pars[3],length(x)),
              sigma=trans.C(C[1:length(x),1:length(x),],edge.rates),log=T)+
      dnorm(theta.prime,0,root.prior.var,log=T)-dnorm(pars[3],0,root.prior.var,log=T)
    if(exp(R)>runif(1,0,1)){
      pars[3]<-theta.prime
      accept[3]<-accept[3]+1
    }
    prop[3]<-prop[3]+1
    
    #r0 update#
    theta.prime<-runif(1,pars[length(tree$tip.label)+4]-win[4]/2,pars[length(tree$tip.label)+4]+win[4]/2)
    R<-dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
               mean=theta.prime+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
               sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)-
      dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
              mean=pars[length(tree$tip.label)+4]+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
              sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)+
      dnorm(theta.prime,0,r0.prior.var,log=T)-dnorm(pars[length(tree$tip.label)+4],0,r0.prior.var,log=T)
    if(exp(R)>runif(1,0,1)){
      pars[length(tree$tip.label)+4]<-theta.prime
      edge.rates<-exp(apply(cbind(pars[4:length(pars)][tree$edge[,1]],pars[4:length(pars)][tree$edge[,2]]),1,mean))
      accept[4]<-accept[4]+1
    }
    prop[4]<-prop[4]+1
    
    blocks<-split(sample((1:(tree$Nnode*2+1))[-(length(tree$tip.label)+1)]),rep(1:n.blocks,block.sizes))
    theta.prime<-pars[4:length(pars)]
    for(j in 1:length(blocks)){
      sf<-exp(runif(block.sizes[j],-win[5],win[5]));theta.prime[blocks[[j]]]<-sf*pars[blocks[[j]]+3]
      edge.rates.prime<-exp(apply(cbind(theta.prime[tree$edge[,1]],theta.prime[tree$edge[,2]]),1,mean))
      R<-dmvnorm(x,
                 mean=rep(pars[3],length(x)),
                 sigma=trans.C(C[1:length(x),1:length(x),],edge.rates.prime),log=T)-
        dmvnorm(x,
                mean=rep(pars[3],length(x)),
                sigma=trans.C(C[1:length(x),1:length(x),],edge.rates),log=T)+
        dmvnorm(theta.prime[-(length(tree$tip.label)+1)],
                mean=pars[length(tree$tip.label)+4]+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
                sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)-
        dmvnorm(pars[4:length(pars)][-(length(tree$tip.label)+1)],
                mean=pars[length(tree$tip.label)+4]+pars[1]*node.depth.edgelength(tree)[-(length(tree$tip.label)+1)],
                sigma=trans.C(C[-(length(tree$tip.label)+1),-(length(tree$tip.label)+1),],rep(pars[2],dim(C)[3])),log=T)+
        sum(log(sf))
      
      
      if(exp(R)>runif(1,0,1)){
        pars[4:length(pars)]<-theta.prime
        edge.rates<-edge.rates.prime
        accept[5]<-accept[5]+block.sizes[j]
      }
      prop[5]<-prop[5]+block.sizes[j]
    }
    
    
    if(i%%report.every==0){
      cat(i,pars[length(tree$tip.label)+4],pars[1:3],mean(edge.rates),var(edge.rates),get.post(pars),accept/prop,"\n")
    }
    
    if(i%%thin==0){
      par.mat[,i/thin]<-c(pars,get.post(pars))
    }
    
    if(i<=tune.period&i%%win.update==0){
      n.update<-i/win.update
      cat("updating proposal windows...\n")
      
      old.win<-win
      win<-exp(log(win)+ifelse(accept/prop-c(0.44,0.44,0.44,0.44,0.44)>0,1,-1)*min(0.1,n.update^-adapt.decay))
      
      cat("r0: ",old.win[1]," -> ",win[1],"\n","rtrend: ",old.win[2]," -> ",win[2],"\n","rvar: ",old.win[3]," -> ",win[3],"\n",
          "root: ",old.win[4]," -> ",win[4],"\n","rate: ",old.win[5]," -> ",win[5],"\n")
      
      accept<-rep(0,5)
      prop<-rep(0,5)
    }
    
    
  }
  return(par.mat)
}