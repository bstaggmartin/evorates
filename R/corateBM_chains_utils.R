#need to fix stop condition because of diagnostics element
#' @export
select.chains<-function(fit,chains){
  if(length(dim(fit$chains))==2){
    stop('correlated BM fit appears to have already been simplified to a single chain')
  }
  if(is.character(chains)){
    chains.exist<-chains%in%dimnames(fit$chains)[[3]]
    if(all(!chains.exist)){
      stop('none of the specified chains found')
    }
    if(any(!chains.exist)){
      warning(paste(chains[which(!chains.exist)],collapse=', '),' not found')
      chains<-chains[which(chains.exist)]
    }
  }else if(is.numeric(chains)){
    chains.exist<-chains<=dim(fit$chains)[3]
    if(all(chains>dim(fit$chains)[3])){
      stop('none of the specified chains found')
    }
    if(any(chains>dim(fit$chains)[3])){
      warning(paste(paste('chain',chains[which(!chains.exist)]),collapse=', '),' not found')
      chains<-chains[which(chains.exist)]
    }
  }else{
    stop('the chains argument must be a character vector specifying chain names or a numeric vector specifying chain indices')
  }
  fit$chains<-fit$chains[,,chains]
  fit$diagnostics$sampler<-fit$diagnostics$sampler[,,chains]
  fit$diagnostics$params<-fit$diagnostics$params[,,chains]
  if(!is.null(fit$quantiles)){
    report.quantiles<-as.numeric(substr(dimnames(fit$quantiles)[[1]],0,
                                        nchar(dimnames(fit$quantiles)[[1]])-1))/100
    fit$quantiles<-fit$quantiles[,,chains]
    if(length(dim(fit$quantiles))==0){
      attr(fit$quantiles,'quantiles')<-paste(report.quantiles*100,'%',sep='')
    }
  }
  if(!is.null(fit$means)){
    fit$means<-fit$means[,chains]
  }
  if(!is.null(fit$MAPs)){
    fit$MAPs<-fit$MAPs[,chains]
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-fit$post.probs[,chains]
  }
  if(length(chains)==1){
    if(is.character(chains)){
      attr(fit,'chains')<-chains
    }else{
      attr(fit,'chains')<-paste('chain',chains)
    }
  }
  fit
}

#automatically excludes warmup not included in chains and inits--doesn't make sense anymore
#' @export
combine.chains<-function(fit,exclude.all.warmup=T){
  if(length(dim(fit$chains))==2){
    warning("correlated rate BM fit appears to have already been simplified to a single chain: nothing done")
  }else{
    if(exclude.all.warmup){
      chains.len<-dim(fit$chains)[1]
      fit$chains<-fit$chains[((chains.len-fit$diagnostics$iter+fit$diagnostics$warmup)+1):chains.len,,]
    }
    diags.len<-dim(fit$diagnostics$sampler)[1]
    chains.len<-dim(fit$chains)[1]
    if(diags.len-chains.len!=0){
      fit$diagnostics$sampler<-reduce.array(fit$diagnostics$sampler,1:(diags.len-chains.len))
    }
    fit$chains<-do.call(rbind,asplit(fit$chains,3))
    names(dimnames(fit$chains))<-c('iterations','parameters')
    fit$diagnostics$sampler<-do.call(rbind,asplit(fit$diagnostics$sampler,3))
    names(dimnames(fit$diagnostics$sampler))<-c('iterations','parameters')
    fit$diagnostics$params<-fit$diagnostics$params[,,1]
    fit$diagnostics$params[1,]<-NA
    fit$diagnostics$params[2,]<-apply(fit$chains,2,rstan::ess_bulk)
    fit$diagnostics$params[3,]<-apply(fit$chains,2,rstan::ess_tail)
    fit$diagnostics$params[4,]<-apply(fit$chains,2,rstan::Rhat)
    if(!is.null(fit$quantiles)){
      report.quantiles<-as.numeric(substr(dimnames(fit$quantiles)[[1]],0,
                                          nchar(dimnames(fit$quantiles)[[1]])-1))/100
      fit$quantiles<-apply(fit$chains,2,quantile,probs=report.quantiles)
      if(length(dim(fit$quantiles))==0){
        attr(fit$quantiles,'quantiles')<-paste(report.quantiles*100,'%',sep='')
      }
    }
    if(!is.null(fit$means)){
      fit$means<-apply(fit$chains,2,mean)
    }
    if(!is.null(fit$MAPs)){
      fit$MAPs<-fit$chains[which.max(fit$diagnostics$sampler[,'lp__']),]
    }
    if(!is.null(fit$post.probs)){
      fit$post.probs<-apply(fit%chains%'_dev',2,function(ii) sum(ii>0)/length(ii))
    }
    fit
  }
}

#for getting rid of warmup in chains and related elements
#making it only reduce sampler
#' @export
exclude.warmup<-function(fit,sampler=T){
  if(length(dim(fit$chains))==2){
    simplified<-T
  }else{
    simplified<-F
  }
  chains.len<-dim(fit$chains)[1]
  incl.inds<-((chains.len-fit$diagnostics$iter+fit$diagnostics$warmup)+1):chains.len
  if(all(incl.inds==(1:chains.len))){
    if(sampler){
      diags.len<-dim(fit$diagnostics$sampler)[1]
      if(diags.len<fit$diagnostics$iter){
        warning('warmup is already excluded from chains and sampler: nothing done')
      }else{
        fit$diagnostics$sampler<-coerce.to.array(fit$diagnostics$sampler)
        fit$diagnostics$sampler<-reduce.array(fit$diagnostics$sampler,1:fit$diagnostics$warmup)
        if(simplified){
          fit$diagnostics$sampler<-fit$diagnostics$sampler[,,1]
        }
      }
    }else{
      warning('warmup is already excluded from chains: nothing done')
    }
  }else{
    fit$chains<-coerce.to.array(fit$chains)
    fit$chains<-fit$chains[incl.inds,,]
    chains.len<-dim(fit$chains)[1]
    fit$diagnostics$params<-coerce.to.array(fit$diagnostics$params)
    fit$diagnostics$params[2,,]<-apply(fit$chains,c(2,3),rstan::ess_bulk)
    fit$diagnostics$params[3,,]<-apply(fit$chains,c(2,3),rstan::ess_tail)
    fit$diagnostics$params[4,,]<-apply(fit$chains,c(2,3),rstan::Rhat)
    if(sampler){
      fit$diagnostics$sampler<-coerce.to.array(fit$diagnostics$sampler)
      fit$diagnostics$sampler<-reduce.array(fit$diagnostics$sampler,1:fit$diagnostics$warmup)
    }
    if(!is.null(fit$quantiles)){
      report.quantiles<-as.numeric(substr(dimnames(fit$quantiles)[[1]],0,
                                          nchar(dimnames(fit$quantiles)[[1]])-1))/100
      fit$quantiles<-apply(fit$chains,c(2,3),quantile,probs=report.quantiles)
      if(length(dim(fit$quantiles))==0){
        attr(fit$quantiles,'quantiles')<-paste(report.quantiles*100,'%',sep='')
      }
    }
    if(!is.null(fit$means)){
      fit$means<-apply(fit$chains,c(2,3),mean)
    }
    if(!is.null(fit$MAPs)){
      nchain<-dim(fit$chains)[3]
      diags.len<-dim(fit$diagnostics$sampler)[1]
      if(diags.len-chains.len>0){
        trimmed.sampler.params<-reduce.array(fit$diagnostics$sampler,1:(diags.len-chains.len))
      }else{
        trimmed.sampler.params<-fit$diagnostics$sampler
      }
      MAP.inds<-sapply(1:nchain, function(ii) which.max(trimmed.sampler.params[,'lp__',ii]))
      MAPs<-sapply(1:nchain,function(ii) fit$chains[MAP.inds[ii],,ii])
      dimnames(MAPs)<-dimnames(fit$chains)[-1]
      fit$MAPs<-MAPs
    }
    if(!is.null(fit$post.probs)){
      fit$post.probs<-apply(fit%chains%'_dev',c(2,3),function(ii) sum(ii>0)/length(ii))
    }
    if(simplified){
      fit$chains<-fit$chains[,,1]
      fit$diagnostics$sampler<-fit$diagnostics$sampler[,,1]
      fit$diagnostics$params<-fit$diagnostics$params[,,1]
    }
  }
  fit
}

#it works, but you have to think of the implications of the sampler excluding warmup--is that okay?
#still gotta work on it--returns error in this case.
#I think it works now--gotta double-check when I'm less tired...