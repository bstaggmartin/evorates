#returns error if correlation matrices are not square
#returns 0 if matrix is square but not valid correlation matrix (i.e., diagonal elements not equal to 1, non-diagonal
#elements outside of [-1,1] interval, non positive-definite)
#' @export
dlkj<-function(cor,nu=1,log=F){
  K<-dim(cor)[1]
  if(K!=dim(cor)[2]){
    stop('correlation matrices must be square')
  }
  if(is.matrix(cor)){
    tmp<-array(NA,c(K,K,1))
    tmp[,,]<-cor
  }
  comp<-rep(1,K)
  inds<-1:(K-1)
  norm.const<-sum(lgamma(nu+(K-1)/2)-(inds/2)*log(pi)-lgamma(nu+(K-inds-1)/2))
  foo<-function(x,nu,comp){
    tmp<-x[lower.tri(x)]
    if(any(tmp< -1)|any(tmp>1)|isFALSE(all.equal(diag(x),comp))){
      -Inf
    }else{
      test<-try(chol(x),silent=T)
      if(inherits(test,'try-error')){
        -Inf
      }else{
        determinant(x)$modulus*(nu-1)
      }
    }
  }
  out<-apply(cor,3,foo,nu=nu,comp=comp)+norm.const
  if(!log){
    out<-exp(out)
  }
  out
}
sig<-aperm(rethinking::rlkjcorr(1e3,5,2),c(2,3,1))
dlkj(sig,nu=2)
plot(dlkj(sig,nu=2,log=T)~apply(sig,3,rethinking::dlkjcorr,eta=2))
#yep, just differs by a constant

#returns 0 if any rows are not unit simplices
#' @export
ddiri<-function(simplex,alpha=1,log=F){
  if(is.null(dim(simplex))){
    simplex<-t(simplex)
  }
  K<-ncol(simplex)
  n<-nrow(simplex)
  alpha<-rep(alpha,length.out=K)
  norm.const<-lgamma(sum(alpha))-sum(lgamma(alpha))
  tmp<-rowSums(log(simplex)*(rep(alpha,each=n)-1))
  tmp[isFALSE(all.equal(rowSums(simplex),rep(1,n)))]<- -Inf
  out<-tmp+norm.const
  if(!log){
    out<-exp(out)
  }
  out
}
# simplex<-gtools::rdirichlet(5,rep(1,3))
# ddiri(simplex,rep(2,3))
# gtools::ddirichlet(simplex,rep(2,3))
# #ha! I think your function is the exact same, just about