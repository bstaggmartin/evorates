#' @importFrom truncnorm rtruncnorm
#' @export
simBM<-function(tree,x,n.sim=1000,along.branch=T,res=100,return.MLE=T,rates=NULL,
                lb=-Inf,ub=Inf,lb.reflective=F,ub.reflective=F,zero.brlen=1e-8){
  ##BASIC ERROR CHECKS AND FIXES##
  if(!inherits(tree,'phylo')){
    stop('tree should be object of class \"phylo\"')
  }
  if(!inherits(x,'numeric')){
    stop('x should be of class \"numeric\"')
  }
  if(!is.vector(x)){
    stop('x looks like a matrix; simBM does not yet support multivariate trait evolution')
  }
  if(ub<lb){
    stop('badly-defined boundaries; upper boundary must be greater than lower boundary')
  }
  if(lb>min(x)|ub<max(x)){
    warning('range of x exceeds either upper or lower boundary;
            simBM still worked, but the results are likely wonky (and most definitely unrealistic)')
  }
  n.sim<-round(n.sim);res<-round(res)
  if(n.sim<1|res<1){
    stop('n.sim or res is less than 1')
  }
  ####
  
  ##EXTRACTING USEFUL PARAMETERS##
  tree<-ape::reorder.phylo(tree)
  ntax<-length(tree$tip.label)
  edges<-tree$edge
  edge.lens<-tree$edge.length
  if(any(edge.lens==0)){
    warning('tree has branch lengths of 0, which pose mathematical difficulties;
            simBM set branch lengths of 0 to',zero.brlen,'.')
    edge.lens[edge.lens==0]<-zero.brlen
  }
  if(length(x)!=ntax){
    stop('differing numbers of elements in x and tips in tree;
         simBM does not yet support intraspecific variation or missing values--
         please prune your tree and/or average trait values accordingly')
  }
  if(is.null(names(x))){
    warning('x has no names;
            simBM assumed elements of x are given in the same order as tip labels in tree')
    names(x)<-tree$tip.label
  }
  tip.states<-x[tree$tip.label]
  sigma<-mean(ape::pic(x[1:ntax],phy=multi2di(tree))^2)
  if(!lb.reflective){
    lb.trunc<-lb
  }else{
    lb.trunc<- -Inf
  }
  if(!ub.reflective){
    ub.trunc<-ub
  }else{
    ub.trunc<-Inf
  }
  ####
  
  ##IF THERE ARE REFLECTIVE BOUNDARIES, DEFINING FUNCTIONS TO TAKE CARE OF THAT##
  if(lb.reflective|ub.reflective){
    refl.b.calc<-function(x,x.base,sig,
                          app.lb=lb,app.ub=ub,hard.lb=lb.trunc,hard.ub=ub.trunc,refl.lb=lb.reflective,refl.ub=ub.reflective){
      if(refl.lb&!(refl.ub)){
        while(!(all(x>app.lb))){
          x<-ifelse(x<app.lb,app.lb-(x-app.lb),x)
          if(!(all(x<app.ub))){
            problem.indices<-which(x>app.ub)
            x[problem.indices]<-rtruncnorm(length(problem.indices),hard.lb,hard.ub,x.base,sqrt(sig))
          }
        }
      }
      if(refl.ub&!(refl.lb)){
        while(!(all(x<app.ub))){
          x<-ifelse(x>app.ub,app.ub-(x-app.ub),x)
          if(!(all(x>app.lb))){
            problem.indices<-which(x<app.lb)
            x[problem.indices]<-rtruncnorm(length(problem.indices),hard.lb,hard.ub,x.base,sqrt(sig))
          }
        }
      }
      if(refl.lb&refl.ub){
        while(!(all(x>app.lb&x<app.ub))){
          x<-ifelse(x<app.lb,app.lb-(x-app.lb),x)
          x<-ifelse(x>app.ub,app.ub-(x-app.ub),x)
        }
      }
      x
    }
    refl.b.calc.interp<-function(x,x.base,sig,time.pts,time.pt,
                                 app.lb=lb,app.ub=ub,hard.lb=lb.trunc,hard.ub=ub.trunc,refl.lb=lb.reflective,refl.ub=ub.reflective){
      if(refl.lb&!(refl.ub)){
        while(!(all(x>app.lb))){
          x<-ifelse(x<app.lb,app.lb-(x-app.lb),x)
          if(!(all(x<app.ub))){
            problem.indices<-which(x>app.ub)
            x[problem.indices]<-interp(length(problem.indices),x.base[,problem.indices],time.pts,time.pt[problem.indices],
                                       sqrt(sig),hard.lb,hard.ub)
          }
        }
      }
      if(refl.ub&!(refl.lb)){
        while(!(all(x<app.ub))){
          x<-ifelse(x>app.ub,app.ub-(x-app.ub),x)
          if(!(all(x>app.lb))){
            problem.indices<-which(x<app.lb)
            x[problem.indices]<-interp(length(problem.indices),x.base[,problem.indices],time.pts,time.pt[problem.indices],
                                       sqrt(sig),hard.lb,hard.ub)
          }
        }
      }
      if(refl.lb&refl.ub){
        while(!(all(x>app.lb&x<app.ub))){
          x<-ifelse(x<app.lb,app.lb-(x-app.lb),x)
          x<-ifelse(x>app.ub,app.ub-(x-app.ub),x)
        }
      }
      x
    }
  }
  ####
  
  #modify edge lengths
  if(!is.null(rates)){
    edge.lens.true<-edge.lens
    edge.lens<-edge.lens*rates/sigma
  }
  
  ##TIPs-TO-ROOT TREE TRAVERSAL##
  mu<-rep(NA,nrow(edges)+1)
  p<-rep(NA,nrow(edges)+1)
  for(e in nrow(edges):0){
    if(e==0){
      n<-ntax+1
      d<-edges[which(edges[,1]==n),2]
      p.a<-sum(p[d])
      mu[n]<-sum(mu[d]*p[d]/p.a)
      t<-0
      p[n]<-p.a/(1+t*p.a)
      break
    }
    n<-edges[e,2]
    if(length(which(edges[,1]==n))==0){
      mu[n]<-tip.states[n]
      t<-edge.lens[e]
      p[n]<-1/t
    }
    if(length(which(edges[,1]==n))>0){
      d<-edges[which(edges[,1]==n),2]
      p.a<-sum(p[d])
      mu[n]<-sum(mu[d]*p[d]/p.a)
      t<-edge.lens[e]
      p[n]<-p.a/(1+t*p.a)
    }
  }
  ####
  
  ##ROOT-TO-TIPS TREE TRAVERSAL WITH NOISE##
  mu.sim<-matrix(NA,nrow=nrow(edges)+1,ncol=n.sim)
  mu.sim[1:ntax,]<-rep(mu[1:ntax],n.sim)
  mu.sim[ntax+1,]<-rtruncnorm(n.sim,lb.trunc,ub.trunc,mu[ntax+1],sqrt(sigma*1/p[ntax+1]))
  p.sim<-(p[edges[,2]]/(1-edge.lens*p[edges[,2]])+1/edge.lens)[order(edges[,2])]
  p.sim<-append(p.sim,p[ntax+1],ntax)
  if(lb.reflective|ub.reflective){
    mu.sim[ntax+1,]<-refl.b.calc(mu.sim[ntax+1,],mu[ntax+1],sqrt(sigma*1/p[ntax+1]))
  }
  for(e in 1:nrow(edges)){
    n<-edges[e,2]
    if(length(which(edges[,1]==n))==0){
      next
    }
    a.n<-edges[e,1]
    t<-edge.lens[e]
    mu.sim[n,]<-mu[n]*p[n]*t+mu.sim[a.n,]-mu.sim[a.n,]*p[n]*t
    mu.sim[n,]<-rtruncnorm(n.sim,lb.trunc,ub.trunc,mu.sim[n,],sqrt(sigma*1/p.sim[n]))
    if(lb.reflective|ub.reflective){
      mu.sim[n,]<-refl.b.calc(mu.sim[n,],mu[n],sqrt(sigma*1/p.sim[n]))
    }
  }
  ####
  
  #unmodify edge lengths
  if(!is.null(rates)){
    edge.lens<-edge.lens.true
  }
  #make rates vector all sigma if no rate heterogeneity
  if(is.null(rates)){
    rates<-rep(sigma,nrow(edges))
  }
  
  ##SIMULATING ANAGENETIC CHANGE ALONG BRANCHES##
  if(along.branch){
    time.vec<-seq(0,max(node.depth.edgelength(tree)),length.out=res)
    ###Defining functions for interpolating Brownian Bridges between two points
    get.last.pt<-function(ii,tt,tts){
      max(which(!is.na(ii)&tts<tt))
    }
    get.next.pt<-function(ii,tt,tts){
      min(which(!is.na(ii)&tts>tt))
    }
    interp<-function(sims,mus,time.pts,time.pt,sig,int.lb.trunc=lb.trunc,int.ub.trunc=ub.trunc){
      mus<-as.matrix(mus)
      last.pt<-mapply(get.last.pt,
                      ii=split(mus,rep(1:sims,each=nrow(mus))),
                      tt=split(time.pts[time.pt],1:length(time.pt)),
                      tts=rep(list(time.pts),sims))
      next.pt<-mapply(get.next.pt,
                      ii=split(mus,rep(1:sims,each=nrow(mus))),
                      tt=split(time.pts[time.pt],1:length(time.pt)),
                      tts=rep(list(time.pts),sims))
      rtruncnorm(sims,int.lb.trunc,int.ub.trunc,
                            (time.pts[time.pt]-time.pts[last.pt])/(time.pts[next.pt]-time.pts[last.pt])*
                              (mus[cbind(next.pt,1:sims)]-mus[cbind(last.pt,1:sims)])+mus[cbind(last.pt,1:sims)],
                            sqrt(sig*(time.pts[time.pt]-time.pts[last.pt])*
                              (time.pts[next.pt]-time.pts[time.pt])/
                              (time.pts[next.pt]-time.pts[last.pt])))
    }
    ###
    mu.mat<-array(NA,dim=c(nrow(edges),length(time.vec),n.sim))
    for(e in 1:nrow(edges)){
      n1<-edges[e,1];n2<-edges[e,2]
      t1<-node.depth.edgelength(tree)[n1];t2<-node.depth.edgelength(tree)[n2]
      int.ts<-which(time.vec>t1&time.vec<t2)
      if(length(int.ts>0)){
        start.mus<-mu.sim[n1,];end.mus<-mu.sim[n2,]
        int.mus<-mu.mat[e,int.ts,]
        #have to make sure vector is treated 'vertically' when only simulating 1 contSimmap
        if(is.vector(int.mus)&n.sim==1){
          int.mus<-as.matrix(int.mus)
        }
        mu.tmp<-rbind(start.mus,int.mus,end.mus)
        t.tmp<-c(t1,time.vec[int.ts],t2)
        ###special case when there is only 1 intervening time point between two nodes
        if(length(int.ts)==1){
          mu.tmp[2,]<-interp(sims=n.sim,mus=mu.tmp,time.pts=t.tmp,time.pt=rep(2,n.sim),sig=rates[e])
          if(lb.reflective|ub.reflective){
            mu.tmp[2,]<-refl.b.calc.interp(mu.tmp[2,],mu.tmp,rates[e],t.tmp,2)
          }
        ###normal case when there are multiple intervening time points
        }else{
          ord.tmp<-sapply(1:n.sim,function(ii) sample(x=2:(length(t.tmp)-1),size=length(t.tmp)-2))
          for(i in 1:nrow(ord.tmp)){
            tmp.indices<-cbind(ord.tmp[i,],1:n.sim)
            mu.tmp[tmp.indices]<-interp(sims=n.sim,mus=mu.tmp,time.pts=t.tmp,time.pt=ord.tmp[i,],sig=rates[e])
            if(lb.reflective|ub.reflective){
              mu.tmp[tmp.indices]<-refl.b.calc.interp(mu.tmp[tmp.indices],mu.tmp,rates[e],t.tmp,ord.tmp[i,])
            }
          }
        }
        mu.tmp<-mu.tmp[-c(1,nrow(mu.tmp)),]
        mu.mat[e,int.ts,]<-mu.tmp
      }
    }
  }
  ####
  
  #re-modify edge lengths
  if(!is.null(rates)){
    edge.lens.true<-edge.lens
    edge.lens<-edge.lens*rates/sigma
  }
  
  ##ROOT-TO-TIPS TREE TRAVERSAL TO GET OVERALL ML ESTIMATES##
  if(return.MLE){
    mu.true<-mu
    p.true<-p
    for(e in 1:nrow(edges)){
      n<-edges[e,2]
      if(length(which(edges[,1]==n))==0){
        next
      }
      a.n<-edges[e,1]
      t<-edge.lens[e]
      mu.true[n]<-mu.true[n]*p.true[n]*t+mu.true[a.n]-mu.true[a.n]*p.true[n]*t
      p.true[n]<-p.true[n]/(1-t*p.true[n])+(p.true[a.n]-p.true[n])/(1+t*(p.true[a.n]-p.true[n]))
    }
  }
  ####
  
  ##OUTPUT##
  if(!return.MLE&!along.branch){
    out<-list(nodes=mu.sim)
  }else if(return.MLE&!along.branch){
    out<-list(nodes=mu.sim,MLE=mu.true,MLE.se=c(rep(0,ntax),sigma*1/p.true[ntax+(1:tree$Nnode)]))
  }else if(!return.MLE&along.branch){
    out<-list(nodes=mu.sim,edges=mu.mat,ts=time.vec)
  }else{
    out<-list(nodes=mu.sim,edges=mu.mat,ts=time.vec,MLE=mu.true,MLE.se=c(rep(0,ntax),sigma*1/p.true[ntax+(1:tree$Nnode)]))
  }
  class(out)<-"simBM"
  out
  ####
}
