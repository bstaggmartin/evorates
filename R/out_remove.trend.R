#remove trends from rates in evorates_fit object
#could make more efficient by pre-selecting rather than post-selecting
#also check out get.bg.rate --> there might be some ways to speed this up
#maybe add support for inits? Probably would be rarely used...

#' @export
remove.trend<-function(fit,
                       element=c('chains','quantiles','means','MAPs','diagnostics'),
                       select=NULL,select.extra=NULL,
                       simplify=TRUE){
  if(!inherits(fit,'evorates_fit')){
    stop("fit must be a fitted evolving rates model fit (class 'evorates_fit')")
  }
  tree<-fit$call$tree
  try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics')),silent=TRUE)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a evolving rates model fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', or 'diagnostics'")
  }
  element<-try.element
  if(element=='quantiles'&!is.null(fit$quantiles)&is.null(select.extra)){
    select.extra<-as.numeric(gsub('%','',dimnames(fit$quantiles)[[1]]))/100
  }
  e<-nrow(tree$edge)
  try.Rsig2<-try(fit%chains%'^R_sig2$',silent=TRUE)
  if(inherits(try.Rsig2,'try-error')){
    warning('fit has constrained R_sig2: residual rates will be identical to rate at root of tree')
    tmp<-.expand.element(chains)
    dims<-dim(tmp)
    tmp.dimnames<-dimnames(tmp)
    tmp.dimnames[[2]]<-paste('R',1:e,sep='_')
    tmp<-array(.int.chains(fit,0),c(dims[1],e,dims[3]),tmp.dimnames)
    tmp[,tree$edge.length==0,]<-NA
  }else{
    tmp<-.int.chains(fit,1:e)
    try.Rmu<-try(fit%chains%'^R_mu$',silent=TRUE)
    if(inherits(try.Rmu,'try-error')){
      warning('fit has constrained R_mu: residual rates will be identical to estimated rates')
    }else{
      edgerans<-tree$edge
      edgerans[,]<-ape::node.depth.edgelength(tree)[edgerans]
      trend.contrib<-aperm(tmp,c(1,3,2))
      trend.contrib[,,]<-aperm(.int.chains(fit,'^R_mu$'),c(1,3,2))
      trend.contrib<-aperm(trend.contrib,c(3,1,2))
      elen<-rep(NA,e)
      elen[tree$edge.length>0]<-tree$edge.length[tree$edge.length>0]
      trend.contrib<-log(abs(trend.contrib))+log(elen)-
        log(abs(exp(trend.contrib*edgerans[,2])-exp(trend.contrib*edgerans[,1])))
      trend.contrib<-aperm(trend.contrib,c(2,1,3))
      tmp<-tmp+trend.contrib
    }
  }
  
  tmp.fit<-fit
  tmp.fit$chains<-tmp
  if(element!='chains'){
    if(element=='diagnostics'){
      tmp.dimnames<-dimnames(tmp.fit$param.diags)
      tmp.fit$param.diags<-array(NA,c(4,dim(tmp)[2],nchains),
                                 c(tmp.dimnames[1],dimnames(tmp)[2],tmp.dimnames[3]))
      tmp.fit$param.diags[2,,]<-apply(tmp.fit$chains,c(2,3),rstan::ess_bulk)
      tmp.fit$param.diags[3,,]<-apply(tmp.fit$chains,c(2,3),rstan::ess_tail)
      tmp.fit$param.diags[4,,]<-apply(tmp.fit$chains,c(2,3),rstan::Rhat)
    }else{
      if(element=='quantiles'&!is.null(fit$quantiles)){
        report.quantiles<-as.numeric(gsub('%','',dimnames(fit$quantiles)[[1]]))/100
      }
      tmp.fit[[element]]<-NULL
    }
  }
  
  if(is.null(select)){
    select<-dimnames(tmp)[[2]]
  }
  if(element=='chains'|element=='quantiles'|element=='diagnostics'&!is.null(select.extra)){
    select<-list(select,select.extra)
  }
  out<-do.call(paste('.int.',element,sep=''),list(fit=tmp.fit,select=select))
  if(simplify){
    out<-.simplify.element(out)
  }
  out
}
