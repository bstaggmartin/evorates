#simulate trait and rate data under an autocorrelated Brownian motion model
#' @export
gen.corateBM<-function(tree,R0=0,Rsig2=1,X0=0,Rmu=0,evocov=1,
                       intra.var=F,n_obs=rep(1,length(tree$tip.label)),intracov=0.2^2,
                       anc.states=F){
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
    X[tree$edge[i,2],]<-X[tree$edge[i,1],]+sqrt(tree$edge.length[i]*exp(R[i]))*chol.evocov%*%seed
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
  if(intra.var){
    if(NROW(intracov)<k){
      warning('intracov implies lower number of traits than X0 and/or evocov: recycled intracov (assuming no covariance if undeclared) to match other inputs')
    }else if(length(dim(intracov))==0&k>1){
      warning('intracov is a vector: turned intracov into a matrix with no covariance')
    }
    intracov<-.coerce.to.cov.mat(intracov,k)
    chol.intracov<-t(chol(intracov))
    n_obs<-rep(n_obs,length.out=n)
    seed<-matrix(rnorm(k*sum(n_obs)),k,sum(n_obs))
    Y<-X[rep(1:n,n_obs),]+t(chol.intracov%*%seed)
    rownames(Y)<-rep(tree$tip.label,n_obs)
    out$trait.data<-Y
    out$intracov<-intracov
    out$n_obs<-n_obs
  }else{
    out$trait.data<-as.matrix(X[1:n,])
  }
  class(out)<-'corateBM'
  out
}
