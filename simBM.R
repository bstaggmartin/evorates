simBM<-function(tree,x,n.sim=1000,along.branch=T,res=0.01,return.MLE=T,
                lb=-Inf,ub=Inf,lb.reflective=F,ub.reflective=F){
  
  ##EXTRACTING USEFUL PARAMETERS##
  ntax<-length(tree$tip.label)
  edges<-tree$edge
  edge.lens<-tree$edge.length
  tip.states<-x[tree$tip.label];names(tip.states)<-1:ntax
  get.no.D<-function(node,tree){
    tmp<-getDescendants(tree,node)
    length(tmp[tmp<=length(tree$tip.label)])
  }
  no.D<-sapply(1:(tree$Nnode*2+1),get.no.D,tree=tree)
  sigma<-mean(pic(x[1:ntax],phy=tree)^2)
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
  
  ##ROOT-TO-TIPS TREE TRAVERSAL WITH NOISE##
  mu.sim<-matrix(rep(mu,n.sim),nrow=nrow(edges)+1,ncol=n.sim)
  p.sim<-p
  for(i in 1:n.sim){
    mu.sim[ntax+1,i]<-rtruncnorm(1,lb.trunc,ub.trunc,mu.sim[ntax+1,i],sigma*1/p[ntax+1])
    
    ##REFLECTIVE BOUNDS CALCULATIONS##
    if(lb.reflective&!(ub.reflective)){
      while(mu.sim[ntax+1,i]<lb){
        if(mu.sim[ntax+1,i]<lb){
          mu.sim[ntax+1,i]<-lb-(mu.sim[ntax+1,i]-lb)
        }
        if(mu.sim[ntax+1,i]>ub){
          mu.sim[ntax+1,i]<-rtruncnorm(1,lb.trunc,ub.trunc,mu.sim[ntax+1,i],sigma*1/p[ntax+1])
        }
      }
    }
    if(ub.reflective&!(lb.reflective)){
      while(mu.sim[ntax+1,i]>ub){
        if(mu.sim[ntax+1,i]>ub){
          mu.sim[ntax+1,i]<-ub-(mu.sim[ntax+1,i]-ub)
        }
        if(mu.sim[ntax+1,i]<lb){
          mu.sim[ntax+1,i]<-rtruncnorm(1,lb.trunc,ub.trunc,mu.sim[ntax+1,i],sigma*1/p[ntax+1])
        }
      }
    }
    if(lb.reflective&ub.reflective){
      while(mu.sim[ntax+1,i]<lb|mu.sim[ntax+1,i]>ub){
        if(mu.sim[ntax+1,i]<lb){
          mu.sim[ntax+1,i]<-lb-(mu.sim[ntax+1,i]-lb)
        }
        if(mu.sim[ntax+1,i]>ub){
          mu.sim[ntax+1,i]<-ub-(mu.sim[ntax+1,i]-ub)
        }
      }
    }
    ####
    
    for(e in 1:nrow(edges)){
      n<-edges[e,2]
      if(length(which(edges[,1]==n))==0){
        next
      }
      a.n<-edges[e,1]
      t<-edge.lens[e]
      mu.sim[n,i]<-mu.sim[n,i]*p[n]*t+mu.sim[a.n,i]-mu.sim[a.n,i]*p[n]*t
      if(i==1){
        p.sim[n]<-p.sim[n]/(1-t*p.sim[n])+1/t #+1/t since the ancestral node trait value is now fixed
      }
      mu.sim[n,i]<-rtruncnorm(1,lb.trunc,ub.trunc,mu.sim[n,i],sigma*1/p.sim[n])
      #sigma*1/p.sim[n]*sqrt(t/(t+t.d))*log(length(tree$tip.label)/(no.D[n]+1)))
      
      ##REFLECTIVE BOUNDS CALCULATIONS##
      if(lb.reflective&!(ub.reflective)){
        while(mu.sim[n,i]<lb){
          if(mu.sim[n,i]<lb){
            mu.sim[n,i]<-lb-(mu.sim[n,i]-lb)
          }
          if(mu.sim[n,i]>ub){
            mu.sim[n,i]<-rtruncnorm(1,lb.trunc,ub.trunc,mu.sim[n,i],sigma*1/p[n])
          }
        }
      }
      if(ub.reflective&!(lb.reflective)){
        while(mu.sim[n,i]>ub){
          if(mu.sim[n,i]>ub){
            mu.sim[n,i]<-ub-(mu.sim[n,i]-ub)
          }
          if(mu.sim[n,i]<lb){
            mu.sim[n,i]<-rtruncnorm(1,lb.trunc,ub.trunc,mu.sim[n,i],sigma*1/p[n])
          }
        }
      }
      if(lb.reflective&ub.reflective){
        while(mu.sim[n,i]<lb|mu.sim[n,i]>ub){
          if(mu.sim[n,i]<lb){
            mu.sim[n,i]<-lb-(mu.sim[n,i]-lb)
          }
          if(mu.sim[n,i]>ub){
            mu.sim[n,i]<-ub-(mu.sim[n,i]-ub)
          }
        }
      }
      ####
      
    }
  }
  
  ##SIMULATING ANAGENETIC CHANGE ALONG BRANCHES##
  if(along.branch){
    time.vec<-seq(0,max(node.depth.edgelength(tree)),length.out=1/res)
    interp<-function(sims,mus,time.pts,time.pt,sig,int.lb.trunc=lb.trunc,int.ub.trunc=ub.trunc){
      mus<-as.matrix(mus)
      last.pt<-max(which(!is.na(mus[,1])&time.pts<time.pts[time.pt]))
      next.pt<-min(which(!is.na(mus[,1])&time.pts>time.pts[time.pt]))
      return(rtruncnorm(sims,int.lb.trunc,int.ub.trunc,
                        mus[last.pt,]+(time.pts[time.pt]-time.pts[last.pt])/
                          (time.pts[next.pt]-time.pts[last.pt])*
                          (mus[next.pt,]-mus[last.pt,]),
                        sig*(time.pts[time.pt]-time.pts[last.pt])*
                          (time.pts[next.pt]-time.pts[time.pt])/
                          (time.pts[next.pt]-time.pts[last.pt])))
    }
    mu.mat<-array(NA,dim=c(nrow(edges),length(time.vec),n.sim))
    for(e in 1:nrow(edges)){
      n1<-edges[e,1];n2<-edges[e,2]
      t1<-node.depth.edgelength(tree)[n1];t2<-node.depth.edgelength(tree)[n2]
      int.ts<-which(time.vec>t1&time.vec<t2)
      if(length(int.ts>0)){
        start.mus<-mu.sim[n1,];end.mus<-mu.sim[n2,]
        mu.tmp<-rbind(start.mus,mu.mat[e,int.ts,],end.mus)
        t.tmp<-c(t1,time.vec[int.ts],t2)
        ord.tmp<-sample(x=2:(length(t.tmp)-1),size=length(t.tmp)-2)
        #ord.tmp<-sample.evenly(2:(length(t.tmp)-1))
        if(length(ord.tmp)==1){
          ord.tmp<-2
        }
        for(t in ord.tmp){
          mu.tmp[t,]<-interp(sims=n.sim,mus=mu.tmp,time.pts=t.tmp,time.pt=t,sig=sigma)
          
          ##REFLECTIVE BOUNDS CALCULATIONS##
          if(lb.reflective&!(ub.reflective)){
            while(!(all(mu.tmp[t,]>lb))){
              mu.tmp[t,]<-ifelse(mu.tmp[t,]<lb,lb-(mu.tmp[t,]-lb),mu.tmp[t,])
              if(!(all(mu.tmp[t,]<ub))){
                problem.indices<-which(mu.tmp[t,]>ub)
                mu.tmp[t,problem.indices]<-interp(sims=length(problem.indices),mus=mu.tmp[,problem.indices],
                                                  time.pts=t.tmp,time.pt=t,sig=sigma)
              }
            }
          }
          if(ub.reflective&!(lb.reflective)){
            while(!(all(mu.tmp[t,]<ub))){
              mu.tmp[t,]<-ifelse(mu.tmp[t,]>ub,ub-(mu.tmp[t,]-ub),mu.tmp[t,])
              if(!(all(mu.tmp[t,]>lb))){
                problem.indices<-which(mu.tmp[t,]<lb)
                mu.tmp[t,problem.indices]<-interp(sims=length(problem.indices),mus=mu.tmp[,problem.indices],
                                                  time.pts=t.tmp,time.pt=t,sig=sigma)
              }
            }
          }
          if(lb.reflective&ub.reflective){
            while(!(all(mu.tmp[t,]>lb&mu.tmp[t,]<ub))){
              mu.tmp[t,]<-ifelse(mu.tmp[t,]<lb,lb-(mu.tmp[t,]-lb),mu.tmp[t,])
              mu.tmp[t,]<-ifelse(mu.tmp[t,]>ub,ub-(mu.tmp[t,]-ub),mu.tmp[t,])
            }
          }
          ####
          
        }
        mu.tmp<-mu.tmp[-c(1,nrow(mu.tmp)),]
        mu.mat[e,int.ts,]<-mu.tmp
      }
    }
    if(!return.MLE){
      out<-list(nodes=mu.sim,edges=mu.mat,ts=time.vec)
      class(out)<-"simBM"
      return(out)
    }
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
    if(!along.branch){
      out<-list(nodes=mu.sim,MLE=mu.true,MLE.se=c(rep(0,ntax),sigma*1/p.true[ntax+(1:tree$Nnode)]))
      class(out)<-"simBM"
      return(out)
    }else{
      out<-list(nodes=mu.sim,edges=mu.mat,ts=time.vec,MLE=mu.true,MLE.se=c(rep(0,ntax),sigma*1/p.true[ntax+(1:tree$Nnode)]))
      class(out)<-"simBM"
      return(out)
    }
  }
  out<-list(nodes=mu.sim)
  class(out)<-"simBM"
  return(out)
}