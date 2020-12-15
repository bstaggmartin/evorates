#simulate trait and rate data under an autocorrelated Brownian motion model
#' @export
sim.corateBM<-function(tree,R0=0,Rsig2=1,X0=0,Rmu=0,evocov=1,trait.names=NULL,
                       n.obs=rep(1,length(tree$tip.label)),intra.var=F,intracov=0.2^2,
                       anc.states=F){
  tree<-.format.tree(tree)$tree
  k<-max(sapply(list(X0,evocov,intracov),NROW))
  if(length(X0)<k){
    warning('X0 implies lower number of traits than evocov and/or intracov: recycling X0 to match other inputs')
    X0<-rep(X0,length.out=k)
  }
  if(NROW(evocov)<k){
    warning('evocov implies lower number of traits than X0 and/or intracov: recycled evocov (assuming no covariance if undeclared) to match other inputs')
  }else if(length(dim(evocov))==0&k>1){
    warning('evocov is a vector: turned evocov into a matrix with no covariance')
  }
  evocov<-.coerce.to.cov.mat(evocov,k)
  chol.evocov<-t(chol(evocov))
  if(is.null(trait.names)){
    trait.names<-rep(NA,length.out=k)
  }
  def.trait.names<-paste('X',1:k,sep='_')
  trait.names<-ifelse(is.na(trait.names),def.trait.names,trait.names)
  
  n<-length(tree$tip.label)
  edge<-tree$edge
  n_e<-nrow(edge)
  ndepths<-node.depth.edgelength(tree)
  elen<-tree$edge.length
  X<-matrix(NA,n+tree$Nnode,k)
  rownames(X)<-c(tree$tip.label,n+1:tree$Nnode)
  colnames(X)<-trait.names
  X[n+1,]<-X0
  n_R<-vector('numeric',max(edge))
  n_R[n+1]<-R0
  R<-vector('numeric',n_e)
  seed<-matrix(rnorm(n_e*(k+2)),nrow=n_e)
  for(i in 1:n_e){
    t<-elen[i]
    n_R[edge[i,2]]<-n_R[edge[i,1]]+seed[i,1]*sqrt(Rsig2*t)
    R[i]<-mean(n_R[edge[i,]])+seed[i,2]*sqrt(Rsig2*t/12)
    if(Rmu!=0&t!=0){
      R[i]<-R[i]-log(abs(Rmu))-log(t)+log(abs(diff(exp(Rmu*ndepths[edge[i,]]))))
    }
    X[edge[i,2],]<-X[edge[i,1],]+sqrt(t*exp(R[i]))*chol.evocov%*%seed[i,-(1:2)]
  }
  scalar<-1/mean(diag(evocov))
  evocov<-evocov*scalar
  R<-R-log(scalar);R0<-R0-log(scalar)
  out<-list('R_0'=R0,'R_sig2'=Rsig2,'X_0'=X0,'R'=R,'R_mu'=Rmu,'tree'=tree)
  if(anc.states){
    out$X<-X
  }else{
    out$X<-as.matrix(X[1:n,])
  }
  if(k>1){
    out$evocov<-evocov
    out$k<-k
  }
  if(any(n.obs>1)){
    intra.var<-T
  }
  if(intra.var){
    if(NROW(intracov)<k){
      warning('intracov implies lower number of traits than X0 and/or evocov: recycled intracov (assuming no covariance if undeclared) to match other inputs')
    }else if(length(dim(intracov))==0&k>1){
      warning('intracov is a vector: turned intracov into a matrix with no covariance')
    }
    intracov<-.coerce.to.cov.mat(intracov,k)
    chol.intracov<-t(chol(intracov))
    if(is.null(names(n.obs))){
      n.obs<-rep(n.obs,length.out=n)
    }else{
      n.obs<-n.obs[tree$tip.label]
      n.obs[is.na(n.obs)]<-1
      names(n.obs)<-tree$tip.label
    }
    seed<-matrix(rnorm(k*sum(n.obs)),k,sum(n.obs))
    Y<-X[rep(1:n,n.obs),]+t(chol.intracov%*%seed)
    rownames(Y)<-rep(tree$tip.label,n.obs)
    out$trait.data<-Y
    out$intracov<-intracov
    out$n.obs<-n.obs
  }else{
    out$trait.data<-as.matrix(X[1:n,])
  }
  class(out)<-'corateBM'
  out
}


#simulate trait and rate data under an autocorrelated Brownian motion model
#REMOVE oprate and wnrate stuff!
#' @export
sim.corateBM.old<-function(tree,R0=0,Rsig2=1,X0=0,Rmu=0,evocov=1,Ralpha=0,Rtheta=0,WN=F,
                           n.obs=rep(1,length(tree$tip.label)),intra.var=F,intracov=0.2^2,
                           anc.states=F){
  tree<-.format.tree(tree)$tree
  k<-max(sapply(list(X0,evocov,intracov),NROW))
  if(length(X0)<k){
    warning('X0 implies lower number of traits than evocov and/or intracov: recycling X0 to match other inputs')
    X0<-rep(X0,length.out=k)
  }
  if(NROW(evocov)<k){
    warning('evocov implies lower number of traits than X0 and/or intracov: recycled evocov (assuming no covariance if undeclared) to match other inputs')
  }else if(length(dim(evocov))==0&k>1){
    warning('evocov is a vector: turned evocov into a matrix with no covariance')
  }
  if(Ralpha>0&Rmu!=0){
    warning('trended rate evolution (Rmu != 1) is current unsupported under OU rate evolution models (alpha > 0): Rmu set to 0')
    Rmu<-0
  }
  evocov<-.coerce.to.cov.mat(evocov,k)
  chol.evocov<-t(chol(evocov))
  
  n<-length(tree$tip.label)
  edge<-tree$edge
  edge.t.ranges<-matrix(node.depth.edgelength(tree)[edge],ncol=2)
  n_e<-nrow(edge)
  elen<-tree$edge.length
  X<-matrix(NA,n+tree$Nnode,k)
  if(WN){
    R<-rnorm(n_e,R0+rowMeans(matrix(node.depth.edgelength(tree)[edge],ncol=2))*Rmu,sqrt(Rsig2/(elen+1)))
  }else{
    n_R<-X[,1]
    n_R[n+1]<-R0
    R<-vector('numeric',n_e)
  }
  rownames(X)<-c(tree$tip.label,n+1:tree$Nnode)
  X[n+1,]<-X0
  #could speed things up here by generating only a single seed...
  for(i in 1:n_e){
    t<-elen[i]
    if(!WN){
      if(t==0){
        mu<-n_R[edge[i,1]]
        var<-0
      }else if(Ralpha>0){
        mu<-n_R[edge[i,1]]*exp(-Ralpha*t)+Rtheta*(1-exp(-Ralpha*t))
        var<-Rsig2/(2*Ralpha)*(1-exp(-2*Ralpha*t))
      }else{
        mu<-n_R[edge[i,1]]
        var<-Rsig2*t
      }
      n_R[edge[i,2]]<-rnorm(1,mu,sqrt(var))
      if(elen[i]!=0){
        if(Ralpha>0){
          mu<-Rtheta+(sum(n_R[edge[i,]])-2*Rtheta)*(cosh(Ralpha*t)-1)/(Ralpha*t*sinh(Ralpha*t))
          var<-Rsig2/(Ralpha^2*t)*((2-2*cosh(Ralpha*t))/(Ralpha*t*sinh(Ralpha*t))+1)
          if(is.nan(var)){
            var<-0
          }else if(var>Rsig2*t/12|var<0){
            var<-Rsig2*t/12
          }
        }else{
          mu<-mean(n_R[edge[i,]])+Rmu*t/2
          var<-Rsig2*t/12
        }
      }
      R[i]<-rnorm(1,mu,sqrt(var))
    }
    seed<-rnorm(k)
    X[edge[i,2],]<-X[edge[i,1],]+sqrt(t*exp(R[i]))*chol.evocov%*%seed
  }
  scalar<-1/mean(diag(evocov))
  evocov<-evocov*scalar
  R<-R-log(scalar);R0<-R0-log(scalar)
  out<-list('R_0'=R0,'R_sig2'=Rsig2,'R_alpha'=Ralpha,'R_theta'=Rtheta,'X_0'=X0,'R'=R,'R_mu'=Rmu,'tree'=tree,'WN'=WN)
  if(anc.states){
    out$X<-X
  }else{
    out$X<-as.matrix(X[1:n,])
  }
  if(k>1){
    out$evocov<-evocov
    out$k<-k
  }
  if(any(n.obs>1)){
    intra.var<-T
  }
  if(intra.var){
    if(NROW(intracov)<k){
      warning('intracov implies lower number of traits than X0 and/or evocov: recycled intracov (assuming no covariance if undeclared) to match other inputs')
    }else if(length(dim(intracov))==0&k>1){
      warning('intracov is a vector: turned intracov into a matrix with no covariance')
    }
    intracov<-.coerce.to.cov.mat(intracov,k)
    chol.intracov<-t(chol(intracov))
    if(is.null(names(n.obs))){
      n.obs<-rep(n.obs,length.out=n)
    }else{
      n.obs<-n.obs[tree$tip.label]
      n.obs[is.na(n.obs)]<-1
      names(n.obs)<-tree$tip.label
    }
    seed<-matrix(rnorm(k*sum(n.obs)),k,sum(n.obs))
    Y<-X[rep(1:n,n.obs),]+t(chol.intracov%*%seed)
    rownames(Y)<-rep(tree$tip.label,n.obs)
    out$trait.data<-Y
    out$intracov<-intracov
    out$n.obs<-n.obs
  }else{
    out$trait.data<-as.matrix(X[1:n,])
  }
  class(out)<-'corateBM'
  out
}
