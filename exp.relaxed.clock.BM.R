#I've started playing around with priors more such that it's more feasible to sample from prior distribution for starting points and
#really dig into whether things have 'converged'
##One issue I'm realizing is that, while I think some parameters (like r0, rtrend, and rvar) should be heavily penalized to prevent
##exponential 'blow-up' of rates, another promising approach might be to simply constrain the upper and lower limits of the values
##node-wise rates may attain...but having such a 'compound prior' is currently beyond my expertise in how to implement. Maybe I
##should talk to Gideon if he knows anything about it.

relaxed.clock.BM<-function(tree,x,n.iter=1e5,thin=100,inits='random',report.every=100, #basic MCMC pars
                           n.chains=1,sample.chains=1,report.chains=1,try.swap=Inf,dT=0, #multi-chain pars
                           r0.prior.var=log(1e2),rtrend.prior.var=log(1e1),rvar.prior.var=log(1e1),root.prior.var=1e3, #prior pars
                           block.size=10,win=matrix(1,ncol=n.chains,nrow=5),prop.mix=0.9, #prop. pars
                           tune.period=1e4,win.update=50,adapt.decay=0.5,targ.accept=matrix(0.3,ncol=n.chains,nrow=5),report.win=T, #adaptive prop, pars
                           tmp.par.mat="tmp.par.mat"){
  
  
  ##10/3: trans.C and get.post have been updated to be able to 'vectorize' over multiple sets of parameters/rates given a matrix##
  ##where each column corresponds to a set of parameters/rates##
  ##10/7: get.R is on its way, and prop.vals is done now
  ##10/14: get.R done and seems to be working and faster than two calls to get.post() by a 0.5-2 milliseconds!
  ##10/14: Unfortunately, the new get.R seems to be slower than 3 calls to get.R by ~0.5 milliseconds...but the syntax is much
  ##nicer--it will be interesting to see how that scales up with larger trees (all this was tested with 30 tip tree)
  ##10/14 SUMMARY: get.R IS consistently faster than using two get.post's to get the likelihood ratio, and the multi-chain
  ##get.post seems just slightly faster than the single-chain get.post. However, the multi-chain get.R is slightly slower than the
  ##single-chain get.R...though, it does allow for the 'ugly bits' of code to be handled internally and 'wraps up' everything
  ##quite nicely for the main code of the MCMC chain. I should see if there's any additional tricks I might use to further speed up
  ##the multi-chain get.R
  ###This will allow users to easily specify multiple independent chains or an MCMCMC approach with heated chains, though it###
  ###won't allow for parallelization...###
  ##CREATING SOME HELPFUL FUNCTIONS##
  #Some helper functions that are later called by the get.R function--they help in processing indices taken from the adjacent nodes
  #and edges lookup tables
  ex.nas<-function(x){
    x[!is.na(x)]
  }
  ex.root<-function(x,root=num.tip+1){
    x[x!=root]
  }
  
  #Sums edge-wise var-cov matrices, stored in the C array, scaling each matrix according to vector of scalars.
  ##You could gain some speed here by forgoing the sapply?--checked w/ microbenchmark, actually slowed it down!##
  ##Might still want to check into making trans.C output results as lists rather than arrays such that it works more smoothly##
  ##with the later mapply functions that call it...##10/7: did it!
  trans.C<-function(C,rates){
    #quick fix--may want to recode function to better deal with accidental coercion of C array to vector later
    #happens in cases where only one node is updated at a time--C[n,n,] turns into a vector...
    if(length(dim(C))<3){
      if(is.vector(rates)){
        sum(C*rates)
      }else{
        lapply(1:ncol(rates),function(ii) sum(C*rates[,ii]))
      }
    }else{
      if(is.vector(rates)){
        apply(sweep(C,3,rates,'*'),c(1,2),sum)
      }else{
        lapply(1:ncol(rates),function(ii) apply(sweep(C,3,rates[,ii],'*'),c(1,2),sum))
      }
    }
  }
  
  #Propose a new set of values by adding a bactrian distributed variate to currently accepted set of values, pars. Use
  #update.indices to indicate which parameters are having new values proposed. Note that the rvar parameter, pars[2], is handled
  #specially since it has new values proposed via sliding scale, rather than sliding window, moves.
  prop.vals<-function(pars,update.indices,win,nchain=n.chains,ntip=num.tip,nnod=num.node,mix=prop.mix){
    if(nchain==1){
      indicator<-ifelse((1:(nnod+3))%in%update.indices,1,0)
      indicator[2]<-0
      pars.prime<-rbac(length(pars),pars,
                       indicator*c(win[1:3],rep(win[5],ntip),win[4],rep(win[5],ntip-2))/2,mix)
      if(any(update.indices==2)){
        pars.prime[2]<-exp(rbac(1,0,win[2]/2,mix))*pars[2]
      }
      pars.prime
    }else{
      indicator<-sapply(1:nchain,function(ii) ifelse((1:(nnod+3))%in%as.vector(update.indices[,ii]),1,0))
      var.ud<-which(indicator[2,]==1)
      indicator[2,]<-0
      pars.prime<-matrix(rbac(length(pars),pars,
                              indicator*c(win[1:3,],rep(win[5,],ntip),win[4,],rep(win[5,],ntip-2))/2,mix),ncol=nchain)
      if(length(var.ud)>0){
        pars.prime[2,var.ud]<-exp(rbac(length(var.ud),0,win[2,]/2,mix))*pars[2,var.ud]
      }
      pars.prime
    }
  }
  
  #Get the (log) likelihood & prior ratio between a proposed set of parameters, pars.prime, and the currently accepted set of
  #parameters, pars. Use update.indices to indicate which parameters are having new values proposed, and the function will
  #automatically simplify likelihood ratio calculations.
  get.R<-function(pars,pars.prime,update.indices,nchain=n.chains,
                  xx=x,CC=C.array,cc=c.mat,em=edge.mat,ntip=num.tip,nnod=num.node,nn=adj.nodes,ee=adj.edges,
                  r0.pri=r0.prior.var,rtrend.pri=rtrend.prior.var,rvar.pri=rvar.prior.var,root.pri=root.prior.var){
    if(nchain==1){
      #update "simple priors"
      if(any(update.indices==1)){
        rtrend.prior<-dnorm(pars.prime[1],0,rtrend.pri,log=T)-dnorm(pars[1],0,rtrend.pri,log=T)
      }else{rtrend.prior<-0}
      ##note that a scale factor is calculated and added to the prior ratio for rvar; this is because this parameter is updated
      ##according to a scaling, rather than sliding, move, so the proposal distribution is asymmetric and the hastings ratio is
      ##not 1
      if(any(update.indices==2)){
        sf<-pars.prime[2]/pars[2]
        rvar.prior<-dexp(pars.prime[2],1/sqrt(rvar.pri),log=T)-dexp(pars[2],1/sqrt(rvar.pri),log=T)+log(sf)
      }else{rvar.prior<-0}
      if(any(update.indices==3)){
        root.prior<-dnorm(pars.prime[3],0,root.pri,log=T)-dnorm(pars[3],0,root.pri,log=T)
      }else{root.prior<-0}
      if(any(update.indices==(ntip+4))){
        r0.prior<-dnorm(pars.prime[ntip+4],0,r0.pri,log=T)-dnorm(pars[ntip+4],0,r0.pri,log=T)
      }else{r0.prior<-0}
      #update "complex" rate prior
      ##if either rtrend or rvar is updated, the likelihood of each node rate parameter is affected and must be re-calculated
      if(any(update.indices==1)|any(update.indices==2)){
        rate.prior<-dmvn(pars.prime[1:nnod+3][-(ntip+1)],
                         mu=pars.prime[ntip+4]+pars.prime[1]*diag(cc)[-(ntip+1)],
                         sigma=cc[-(ntip+1),-(ntip+1)]*pars.prime[2],log=T)-
          dmvn(pars[1:nnod+3][-(ntip+1)],
               mu=pars[ntip+4]+pars[1]*diag(cc)[-(ntip+1)],
               sigma=cc[-(ntip+1),-(ntip+1)]*pars[2],log=T)
        ##if neither rtrend nor rvar is updated, but some node rate parameters are updated, then the likelihood of the rate
        ##parameters for these nodes and adjacent nodes must be re-calculated
      }else if(any((4:length(pars))%in%update.indices)){
        nodes<-unique(as.vector(nn[(update.indices-3)[update.indices>3],]))
        nodes<-nodes[(nodes!=ntip+1)&(!(is.na(nodes)))]
        rate.prior<-dmvn(X=pars.prime[nodes+3],
                         mu=pars.prime[ntip+4]+pars.prime[1]*diag(cc)[nodes],
                         sigma=cc[nodes,nodes]*pars.prime[2],log=T)-
          dmvn(X=pars[nodes+3],
               mu=pars[ntip+4]+pars[1]*diag(cc)[nodes],
               sigma=cc[nodes,nodes]*pars[2],log=T)
      }else{rate.prior<-0}
      #update likelihood
      ##note that, since ancestral trait values are not being explicitly estimated, these multivariate normal likelihoods cannot
      ##be simplified in the same way they are when updating the rate prior
      if(any((3:length(pars))%in%update.indices)){
        ##create a vector to store edge-wise rates by averaging the node-wise rate parameters at adjacent nodes; if any node rate
        ##parameters are being updated, update the edge-wise rates for each edge adjacent to the effected nodes
        edge.rates.prime<-edge.rates<-exp(apply(cbind(pars[1:nnod+3][em[,1]],pars[1:nnod+3][em[,2]]),1,mean))
        edges<-as.vector(ee[(update.indices-3)[update.indices>3],])
        edges<-edges[!(is.na(edges))]
        edge.rates.prime[edges]<-exp(sapply(edges,function(e) mean(c(pars.prime[em[e,1]+3],pars.prime[em[e,2]+3]))))
        log.lik<-dmvn(X=xx,
                      mu=rep(pars.prime[3],ntip),
                      sigma=trans.C(CC,edge.rates.prime),log=T)-
          dmvn(X=xx,
               mu=rep(pars[3],ntip),
               sigma=trans.C(CC,edge.rates),log=T)
      }else{log.lik<-0}
      #return sum of all (log) prior/likelihood ratios
      sum(rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior,log.lik)
    }else{
      if(any(update.indices==1)){
        rtrend.prior<-dnorm(pars.prime[1,],0,rtrend.pri,log=T)-dnorm(pars[1,],0,rtrend.pri,log=T)
      }else{rtrend.prior<-rep(0,nchain)}
      if(any(update.indices==2)){
        sf<-pars.prime[2,]/pars[2,]
        rvar.prior<-dexp(pars.prime[2,],1/sqrt(rvar.pri),log=T)-dexp(pars[2,],1/sqrt(rvar.pri),log=T)+log(sf)
      }else{rvar.prior<-rep(0,nchain)}
      if(any(update.indices==3)){
        root.prior<-dnorm(pars.prime[3,],0,root.pri,log=T)-dnorm(pars[3,],0,root.pri,log=T)
      }else{root.prior<-rep(0,nchain)}
      if(any(update.indices==(ntip+4))){
        r0.prior<-dnorm(pars.prime[ntip+4,],0,r0.pri,log=T)-dnorm(pars[ntip+4,],0,r0.pri,log=T)
      }else{r0.prior<-rep(0,nchain)}
      all.ud<-sapply(1:nchain,function(ii) any(update.indices[,ii]==1)|any(update.indices[,ii]==2))
      n.all<-sum(all.ud);alls<-which(all.ud)
      prt.ud<-sapply(1:nchain,function(ii) any((1:nnod+3)%in%update.indices[,ii]))&!all.ud
      n.prt<-sum(prt.ud);prts<-which(prt.ud)
      rate.prior<-rep(0,nchain)
      if(n.all>0){
        rate.prior[all.ud]<-mapply(dmvn,
                                   X=split(pars.prime[c(1:ntip,(ntip+2):nnod)+3,all.ud],
                                           rep(1:n.all,each=nnod-1)),
                                   mu=split(pars.prime[ntip+4,all.ud]+pars.prime[1,all.ud]*rep(diag(cc)[-(ntip+1)],each=n.all),
                                            rep(1:n.all,nnod-1)),
                                   sigma=lapply(1:n.all,function(ii) cc[-(ntip+1),-(ntip+1)]*pars[2,alls[ii]]),
                                   log=T)-
          mapply(dmvn,
                 X=split(pars[c(1:ntip,(ntip+2):nnod)+3,all.ud],
                         rep(1:n.all,each=nnod-1)),
                 mu=split(pars[ntip+4,all.ud]+pars[1,all.ud]*rep(diag(cc)[-(ntip+1)],each=n.all),
                          rep(1:n.all,nnod-1)),
                 sigma=lapply(1:n.all,function(ii) cc[-(ntip+1),-(ntip+1)]*pars[2,alls[ii]]),
                 log=T)
      }
      if(n.prt>0){
        nodes<-lapply(split(update.indices[,prt.ud],rep(1:n.prt,each=nrow(update.indices))),
                             function(ii) ex.root(ex.nas(unique(as.vector(nn[(ii-3)[ii>3],])))))
        rate.prior[prt.ud]<-mapply(dmvn,
                                   X=lapply(1:n.prt,function(ii) 
                                     pars.prime[nodes[[ii]]+3,prts[ii]]),
                                   mu=lapply(1:n.prt,function(ii) 
                                     pars.prime[ntip+4,prts[ii]]+pars.prime[1,prts[ii]]*diag(cc)[nodes[[ii]]]),
                                   sigma=lapply(1:n.prt,function(ii) 
                                     cc[nodes[[ii]],nodes[[ii]]]*pars.prime[2,prts[ii]]),
                                   log=T)-
          mapply(dmvn,
                 X=lapply(1:n.prt,function(ii) 
                   pars[nodes[[ii]]+3,prts[ii]]),
                 mu=lapply(1:n.prt,function(ii) 
                   pars[ntip+4,prts[ii]]+pars[1,prts[ii]]*diag(cc)[nodes[[ii]]]),
                 sigma=lapply(1:n.prt,function(ii) 
                   cc[nodes[[ii]],nodes[[ii]]]*pars[2,prts[ii]]),
                 log=T)
      }
      trt.ud<-sapply(1:nchain,function(ii) any((0:nnod+3)%in%update.indices[,ii]))
      n.trt<-sum(trt.ud);trts<-which(trt.ud)
      log.lik<-rep(0,nchain)
      if(n.trt>0){
        edge.rates.prime<-edge.rates<-sapply(trts,function(ii) 
          exp(apply(cbind(pars[1:nnod+3,ii][em[,1]],pars[1:nnod+3,ii][em[,2]]),1,mean)))
        edges<-lapply(split(update.indices[,trt.ud],rep(1:n.trt,each=nrow(update.indices))),
                      function(ii) ex.nas(unique(as.vector(ee[(ii-3)[ii>3],]))))
        columns<-rep(1:n.trt,lengths(edges));par.columns<-trts[columns];edges<-unlist(edges)
        edge.rates.prime[cbind(edges,columns)]<-
          exp(mapply(function(e,c) mean(c(pars.prime[em[e,1]+3,c],pars.prime[em[e,2]+3,c])),e=edges,c=par.columns))
        log.lik[trt.ud]<-mapply(dmvn,
                                X=rep(list(xx),n.trt),
                                mu=split(rep(pars.prime[3,trt.ud],each=ntip),rep(1:n.trt,each=ntip)),
                                sigma=trans.C(CC,edge.rates.prime),
                                log=T)-
          mapply(dmvn,
                 X=rep(list(xx),n.trt),
                 mu=split(rep(pars[3,trt.ud],each=ntip),rep(1:n.trt,each=ntip)),
                 sigma=trans.C(CC,edge.rates),
                 log=T)
      }
      apply(cbind(rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior,log.lik),1,sum)
    }
  }
  
  #Get the (log) posterior probability of a currently accepted set of parameters, pars
  get.post<-function(pars,nchain=n.chains,xx=x,CC=C.array,cc=c.mat,em=edge.mat,ntip=num.tip,nnod=num.node,
                     r0.pri=r0.prior.var,rtrend.pri=rtrend.prior.var,rvar.pri=rvar.prior.var,root.pri=root.prior.var){
    #"simple" priors
    if(nchain==1){
      rtrend.prior<-dnorm(pars[1],0,rtrend.pri,log=T)
      rvar.prior<-dexp(pars[2],1/sqrt(rvar.pri),log=T)
      root.prior<-dnorm(pars[3],0,root.pri,log=T)
      r0.prior<-dnorm(pars[ntip+4],0,r0.pri,log=T)
      #"complex" prior
      rate.prior<-dmvn(pars[1:nnod+3][-(ntip+1)],
                       mu=pars[ntip+4]+pars[1]*diag(cc)[-(ntip+1)],
                       sigma=cc[-(ntip+1),-(ntip+1)]*pars[2],log=T)
      #getting edge-wise rates and calculating likelihood
      edge.rates<-exp(apply(cbind(pars[1:nnod+3][em[,1]],pars[1:nnod+3][em[,2]]),1,mean))
      log.lik<-dmvn(xx,
                    mu=rep(pars[3],ntip),
                    sigma=trans.C(CC,edge.rates),log=T)
      #return the sum
      sum(rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior,log.lik)
    }else{
      rtrend.prior<-dnorm(pars[1,],0,rtrend.pri,log=T)
      rvar.prior<-dexp(pars[2,],1/sqrt(rvar.pri),log=T)
      root.prior<-dnorm(pars[3,],0,root.pri,log=T)
      r0.prior<-dnorm(pars[ntip+4,],0,r0.pri,log=T)
      #"complex" prior
      rate.prior<-mapply(dmvn,
                         X=split(pars[1:nnod+3,][-(ntip+1),],
                                 rep(1:nchain,each=nnod-1)),
                         mu=split(pars[ntip+4,]+pars[1,]*rep(diag(cc)[-(ntip+1)],each=nchain),
                                  rep(1:nchain,nnod-1)),
                         sigma=lapply(1:nchain,function(ii) cc[-(ntip+1),-(ntip+1)]*pars[2,ii]),
                         log=T)
      #getting edge-wise rates and calculating likelihood
      edge.rates<-sapply(1:nchain,function(ii) 
        exp(apply(cbind(pars[1:nnod+3,ii][em[,1]],pars[1:nnod+3,ii][em[,2]]),1,mean)))
      log.lik<-mapply(dmvn,
                      X=rep(list(xx),nchain),
                      mu=split(rep(pars[3,],each=ntip),rep(1:ncol(pars),each=ntip)),
                      sigma=trans.C(CC,edge.rates),
                      log=T)
      #return the sums
      apply(cbind(rtrend.prior,rvar.prior,root.prior,r0.prior,rate.prior,log.lik),1,sum)
    }
  }
  
  
  ##EXTRACTING RELEVANT INFO FROM TREE##
  num.tip<-length(tree$tip.label);num.node<-num.tip+tree$Nnode
  edge.mat<-tree$edge
  #Creating array of edge-wise var-cov matrices for given tree
  ##Simplify this--you only need the *whole* vcv matrix for each node, you don't need an array w/ slice for each edge
  ##You only need an array for the trait values, which only involves the *extant* tips!
  ##This will help you cut down on array size and memory usage!!!
  #Tip vcv array, sliced up by edge
  C.array<-array(0,dim=c(rep(num.tip,2),nrow(edge.mat)))
  for(e in 1:nrow(edge.mat)){
    D<-getDescendants(tree,edge.mat[e,2]);D<-D[D<=num.tip]
    C.array[as.matrix(cbind(expand.grid(D,D),rep(e,length(D)^2)))]<-tree$edge.length[e]
  }
  #All node vcv matrix
  node.hgts<-node.depth.edgelength(tree)
  c.mat<-matrix(node.hgts[mrca(tree,full=T)],ncol=num.node)
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
  if(inits=='informed'){
    pics<-pic(x,tree)^2
    pars<-c(0,var(log(pics))/diag(c.mat)[1],mean(x),
            c(log(pics)[as.character(edge.mat[which(edge.mat[,2]%in%(1:num.tip)),1])],log(pics)))
    names(pars)<-NULL
    pars<-matrix(rep(pars,n.chains),ncol=n.chains,nrow=length(pars))
  }else if(inits=='random'){
    pars<-matrix(NA,ncol=n.chains,nrow=num.node+3)
    pars[c(1:3,num.tip+4),]<-rbind(rnorm(n.chains,0,rtrend.prior.var),
                                   rexp(n.chains,1/sqrt(rvar.prior.var)),
                                   rnorm(n.chains,0,root.prior.var),
                                   rnorm(n.chains,0,r0.prior.var))
    for(e in 1:nrow(edge.mat)){
      pars[edge.mat[e,2]+3,]<-rnorm(n.chains,pars[edge.mat[e,1]+3,]+pars[1,]*tree$edge.length[e],pars[2,]*tree$edge.length[e])
    }
  }
  #make sure users can specify their own init values too!
  
  
  ##INITIALIZING MCMC ALGORITHM##
  n.iter<-n.iter-n.iter%%thin
  par.mat<-array(NA,dim=c(num.node+4,n.iter/thin,length(sample.chains)))
  n.blocks<-floor((num.node+3)/block.size)
  if((num.node+3)%%block.size!=0){
    block.sizes<-c(rep(block.size,n.blocks),(num.node+3)%%block.size)
    n.blocks<-n.blocks+1
  }else{
    block.sizes<-rep(block.size,n.blocks)
  }
  block.dif<-block.sizes[1]-block.sizes[length(block.sizes)]
  accept<-matrix(0,nrow=5,ncol=n.chains)
  prop<-accept
  if(length(dT)==n.chains){
    betas<-1/(1+dT)
  }else{
    betas<-1/(1+dT*(1:n.chains-1))
  }
  if(n.chains>1&!is.infinite(try.swap)){
    Ts<-sort(unique(dT))
    which.Ts<-lapply(Ts,function(ii) which(dT==ii))
    pair.Ts<-lapply(1:(length(which.Ts)-1),function(ii) expand.grid(which.Ts[[ii]],which.Ts[[ii+1]]))
    pairs<-cbind(unlist(lapply(pair.Ts,function(ii) ii[,1])),unlist(lapply(pair.Ts,function(ii) ii[,2])))
  }
  for(i in 2:n.iter){
    blocks<-asplit(array(sapply(1:n.chains,function(ii) 
      c(sample(1:(num.node+3)),rep(NA,block.dif))),dim=c(block.size,n.blocks,n.chains)),2)
    blocks[[length(blocks)]]<-blocks[[length(blocks)]][1:block.sizes[length(block.sizes)],]
    for(j in 1:n.blocks){
      pars.prime<-prop.vals(pars,blocks[[j]],win)
      R<-get.R(pars,pars.prime,blocks[[j]])
      n.par.updates<-sapply(1:n.chains,function(ii)
        c(ifelse(c(1:3,num.tip+4)%in%blocks[[j]][,ii],1,0),sum(blocks[[j]][,ii]%in%c(1:num.tip+3,(num.tip+2):num.node+3))))
      prop<-prop+n.par.updates
      accept.cond<-exp(betas*R)>runif(n.chains,0,1)
      if(any(accept.cond)){
        pars[,accept.cond]<-pars.prime[,accept.cond]
        accept[,accept.cond]<-accept[,accept.cond]+n.par.updates[,accept.cond]
      }
    }
    if(i%%try.swap==0){
      pair<-pairs[sample(1:nrow(pairs),1),]
      R<-betas[pair[1]]*get.post(pars[,pair[2]],nchain=1)+betas[pair[2]]*get.post(pars[,pair[1]],nchain=1)-
        (betas[pair[1]]*get.post(pars[,pair[1]],nchain=1)+betas[pair[2]]*get.post(pars[,pair[2]],nchain=1))
      if(exp(R)>runif(1,0,1)){
        tmp.pars<-pars[,pair[1]]
        pars[,pair[1]]<-pars[,pair[2]]
        pars[,pair[2]]<-tmp.pars
        cat("chains ",pair[1]," & ",pair[2],"swapped!\n")
      }
    }
    if(i%%report.every==0){
      tmp.edge.rates<-exp(sapply(report.chains,function(ii) 
        apply(cbind(pars[1:num.node+3,ii][edge.mat[,1]],pars[1:num.node+3,ii][edge.mat[,2]]),1,mean)))
      for(chain in report.chains){
        cat(sprintf('%9.0f',i),' ',sprintf('%3.0f',chain),' ',sprintf('%10.3f',pars[c(1:3,num.tip+4),chain]),' ',
            sprintf('%12.3f',c(mean(tmp.edge.rates[,chain]),var(tmp.edge.rates[,chain]))),' ',
            sprintf('%12.3f',get.post(pars[,chain],nchain=1)),
            sprintf('%4.3f',accept[,chain]/prop[,chain]),'\n')
      }
    }
    if(i%%thin==0){
      par.mat[,i/thin,]<-c(pars[,sample.chains],get.post(pars[,sample.chains],nchain=length(sample.chains)))
      assign(tmp.par.mat,par.mat,envir=.GlobalEnv)
    }
    if(i<=tune.period&i%%win.update==0){
      n.update<-i/win.update
      cat("updating proposal windows...\n")
      old.win<-win
      win<-exp(log(win)+ifelse(accept/prop-targ.accept>0,1,-1)*min(0.1,n.update^-adapt.decay))
      for(chain in report.chains){
        cat(sprintf('%3.0f',chain),' ',
            sprintf('%8.8s','rtrend: '),sprintf('%6.3f',old.win[1,chain]),' -> ',sprintf('%6.3f',win[1,chain]),'\n',
            sprintf('%13.8s','rvar: '),sprintf('%6.3f',old.win[2,chain]),' -> ',sprintf('%6.3f',win[2,chain]),'\n',
            sprintf('%13.8s','root: '),sprintf('%6.3f',old.win[3,chain]),' -> ',sprintf('%6.3f',win[3,chain]),'\n',
            sprintf('%13.8s','r0: '),sprintf('%6.3f',old.win[4,chain]),' -> ',sprintf('%6.3f',win[4,chain]),'\n',
            sprintf('%13.8s','rate: '),sprintf('%6.3f',old.win[5,chain]),' -> ',sprintf('%6.3f',win[5,chain]),'\n')
      }
      accept<-matrix(0,nrow=5,ncol=n.chains)
      prop<-accept
    }
  }
  
  
  par.mat
}
