#simulate trait and rate data under an autocorrelated Brownian motion model
#' @export
gen.corateBM<-function(tree,R0=0,Rsig2=1,X0=0,Rmu=0,Xsig2=1,
                       intra.var=F,n_obs=rep(1,length(tree$tip.label)),Ysig2=0.2^2,
                       anc.states=F){
  k<-max(sapply(list(X0,Xsig2,Ysig2),NROW))
  if(length(X0)<k){
    warning('X0 implies lower number of traits than Xsig2 and/or Ysig2: recycling X0 to match other inputs')
    X0<-rep(X0,length.out=k)
  }
  if(NROW(Xsig2)<k){
    warning('Xsig2 implies lower number of traits than X0 and/or Ysig2: recycled Xsig2 (assuming no covariance if undeclared) to match other inputs')
  }else if(length(dim(Xsig2))==0&k>1){
    warning('Xsig2 is a vector: turned Xsig2 into a matrix with no covariance')
  }
  Xsig2<-coerce.to.cov.mat(Xsig2,k)
  chol.Xsig2<-t(chol(Xsig2))
  n<-length(tree$tip.label)
  n_e<-nrow(tree$edge)
  X<-matrix(NA,n+tree$Nnode,k)
  n_R<-X[,1]
  R<-rep(NA,n_e)
  rownames(X)<-c(tree$tip.label,n+1:tree$Nnode)
  X[n+1,]<-X0;n_R[n+1]<-R0
  for(i in 1:n_e){
    n_R[tree$edge[i,2]]<-rnorm(1,n_R[tree$edge[i,1]]+Rmu*tree$edge.length[i],sqrt(tree$edge.length[i]*Rsig2))
    R[i]<-rnorm(1,mean(n_R[tree$edge[i,]]),sqrt(Rsig2*tree$edge.length[i]/12))
    seed<-rnorm(k)
    X[tree$edge[i,2],]<-X[tree$edge[i,1],]+sqrt(tree$edge.length[i]*exp(R[i]))*chol.Xsig2%*%seed
  }
  scalar<-1/mean(diag(Xsig2))
  Xsig2<-Xsig2*scalar
  R<-R-log(scalar);R0<-R0-log(scalar)
  if(anc.states){
    out<-list('R0'=R0,'Rsig2'=Rsig2,'X0'=X0,'R'=R,'X'=X,'Rmu'=Rmu)
  }else{
    out<-list('R0'=R0,'Rsig2'=Rsig2,'X0'=X0,'R'=R,'X'=as.matrix(X[1:n,]),'Rmu'=Rmu)
  }
  if(k>1){
    out$Xsig2<-Xsig2
    out$k<-k
  }
  if(intra.var){
    if(NROW(Ysig2)<k){
      warning('Ysig2 implies lower number of traits than X0 and/or Xsig2: recycled Ysig2 (assuming no covariance if undeclared) to match other inputs')
    }else if(length(dim(Ysig2))==0&k>1){
      warning('Ysig2 is a vector: turned Ysig2 into a matrix with no covariance')
    }
    Ysig2<-coerce.to.cov.mat(Ysig2,k)
    chol.Ysig2<-t(chol(Ysig2))
    n_obs<-rep(n_obs,length.out=n)
    seed<-matrix(rnorm(k*sum(n_obs)),k,sum(n_obs))
    Y<-X[rep(1:n,n_obs),]+t(chol.Ysig2%*%seed)
    rownames(Y)<-rep(tree$tip.label,n_obs)
    out$Y<-Y
    out$Ysig2<-Ysig2
    out$n_obs<-n_obs
  }
  class(out)<-'corateBM'
  out
}
