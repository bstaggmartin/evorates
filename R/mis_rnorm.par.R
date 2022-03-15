#' @export
rnorm.par<-function(niter=1,n=1,nchains=1,simplify=TRUE,fit=NULL,chain.names=paste('chain',seq_len(nchains))){
  if(!is.null(fit)){
    if(.is.par(fit)){
      fit<-.make.par.3D(fit)
      dims<-dim(fit)
      niter<-dims[1]
      nchains<-dims[length(dims)]
      tmp<-dimnames(fit)
      chain.names<-tmp[[length(tmp)]]
    }else if(inherits(fit,'evorates_fit')){
      niter<-fit$sampler.control$iter-fit$sampler.control$warmup
      nchains<-fit$sampler.control$chains
      tmp<-attr(fit$chains,'chains')
      if(is.null(tmp)){
        tmp<-dimnames(fit$chains)
        chain.names<-tmp[[length(tmp)]]
      }else{
        chain.names<-tmp
      }
    }
  }
  out<-array(rnorm(niter*n*nchains),
             c(niter,n,nchains),
             c(iterations=list(NULL),
               parameters=list(paste('N(0;1)',if(n==1) '' else 1:n,sep='')),
               chains=list(chain.names)))
  attr(out,'param_type')<-'chains'
  out<-.add.par.class(out)
  if(simplify){
    out<-.simplify.par(out)
  }
  out
}