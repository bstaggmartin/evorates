#' @export
select.chains<-function(fit,chains,simplify=T){
  #input processing
  if(length(dim(fit$chains))==2){
    if(simplify){
      return(fit)
    }else{
      for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
        fit[[i]]<-.expand.element(fit[[i]],simplify=T)
      }
    }
  }
  if(is.numeric(chains)){
    chains<-paste('chain',chains)
  }
  chains.exist<-chains%in%dimnames(fit$chains)[[3]]
  if(all(!chains.exist)){
    stop('none of the specified chains found')
  }
  if(any(!chains.exist)){
    warning(paste(chains[which(!chains.exist)],collapse=', '),' not found')
    chains<-chains[which(chains.exist)]
  }
  fit$sampler.control$chains<-length(chains)
  inds<-match(chains,dimnames(fit$chains)[[3]])
  #subset chains
  for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
    fit[[i]]<-.index.element(fit[[i]],inds,length(dim(fit[[i]])))
    if(simplify){
      fit[[i]]<-.simplify.element(fit[[i]])
    }
  }
  fit
}

#automatically excludes iterations not included in chains and inits--doesn't make sense anymore
#' @export
combine.chains<-function(fit,simplify=T){
  #process input
  if(length(dim(fit$chains))==2){
    if(simplify){
      return(fit)
    }else{
      for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
        fit[[i]]<-.expand.element(fit[[i]],simplify=T)
      }
      chains.len<-dim(fit$chains)[1]
      diags.len<-dim(fit$sampler.params)[1]
      if(diags.len-chains.len>0){
        fit$sampler.params<-.index.element(fit$sampler.params,1:(diags.len-chains.len),1,T)
      }
      return(fit)
    }
  }else{
    fit$sampler.control$chains<-1
    #trim sampler.params to iterations included in chains
    chains.len<-dim(fit$chains)[1]
    diags.len<-dim(fit$sampler.params)[1]
    if(diags.len-chains.len>0){
      fit$sampler.params<-.index.element(fit$sampler.params,1:(diags.len-chains.len),1,T)
    }
    #combine chains, sampler.params, and param.diags
    out.chains<-paste(dimnames(fit$chains)[[3]],collapse=', ')
    for(i in c('chains','sampler.params')){
      tmp<-do.call(rbind,asplit(fit[[i]],3))
      fit[[i]]<-array(tmp,c(dim(tmp),1),c(dimnames(fit[[i]])[-3],chains=out.chains))
    }
    fit$param.diags<-.index.element(fit$chains,1:4,1)
    dimnames(fit$param.diags)<-c(diagnostics=list(c('inits','bulk_ess','tail_ess','Rhat')),
                                 dimnames(fit$chains)[-1])
    fit$param.diags[1,,]<-NA
    fit$param.diags[2,,]<-apply(fit$chains,c(2,3),rstan::ess_bulk)
    fit$param.diags[3,,]<-apply(fit$chains,c(2,3),rstan::ess_tail)
    fit$param.diags[4,,]<-apply(fit$chains,c(2,3),rstan::Rhat)
    #redo any other elements
    for(i in names(fit)[names(fit)%in%c('quantiles','means','MAPs')]){
      fit[[i]]<-NULL
      fit[[i]]<-do.call(paste('.int.',i,sep=''),list(fit=fit,select='.|dev'))
      if(i=='means'|i=='MAPs'){
        fit[[i]]<-array(fit[[i]],dim(fit[[i]])[-1],dimnames(fit[[i]])[-1])
      }
    }
    if(!is.null(fit$post.probs)){
      fit$post.probs<-apply(.int.chains(fit,'R_\\d+_dev'),c(2,3),function(ii) sum(ii>0)/length(ii))
    }
    #simplify as needed
    if(simplify){
      for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
        fit[[i]]<-.simplify.element(fit[[i]])
      }
    }
    fit
  }
}

#for getting rid of warmup in chains and related elements
#' @export
exclude.warmup<-function(fit,warmup=fit$sampler.control$warmup,sampler=T){
  #process input
  if(length(dim(fit$chains))==2){
    simplified<-T
  }else{
    simplified<-F
  }
  if(simplified){
    for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
      fit[[i]]<-.expand.element(fit[[i]],simplify=T)
    }
  }
  target.len<-fit$sampler.control$iter-warmup
  diags.len<-dim(fit$sampler.params)[1]
  #trim sampler as needed if desired
  if(sampler){
    if(diags.len-target.len>0){
      fit$sampler.params<-.index.element(fit$sampler.params,1:(diags.len-target.len),1,T)
    }else if(diags.len-target.len<0){
      warning('desired number of iterations (',target.len,') is longer than number of available iterations (',diags.len,') in sampler.params element: no iterations excluded from sampler.params')
    }
  }
  #trim chains as needed if desired
  chains.len<-dim(fit$chains)[1]
  if(chains.len-target.len<=0){
    if(chains.len-target.len<0){
      warning('desired number of iterations (',target.len,') is longer than number of available iterations (',chains.len,') in chains element: no iterations excluded from chains')
    }
    if(!simplified){
      return(fit)
    }else{
      for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
        fit[[i]]<-.simplify.element(fit[[i]])
      }
      return(fit)
    }
  }else{
    fit$chains<-.index.element(fit$chains,1:(chains.len-target.len),1,T)
  }
  #redo parameter diagnostics (other than inits)
  fit$param.diags[2,,]<-apply(fit$chains,c(2,3),rstan::ess_bulk)
  fit$param.diags[3,,]<-apply(fit$chains,c(2,3),rstan::ess_tail)
  fit$param.diags[4,,]<-apply(fit$chains,c(2,3),rstan::Rhat)
  #redo any other elements
  for(i in names(fit)[names(fit)%in%c('quantiles','means','MAPs')]){
    fit[[i]]<-NULL
    fit[[i]]<-do.call(paste('.int.',i,sep=''),list(fit=fit,select='.|dev'))
    if(i=='means'|i=='MAPs'){
      fit[[i]]<-array(fit[[i]],dim(fit[[i]])[-1],dimnames(fit[[i]])[-1])
    }
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-apply(.int.chains(fit,'R_\\d+_dev'),c(2,3),function(ii) sum(ii>0)/length(ii))
  }
  #simplify as needed
  if(simplified){
    for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
      fit[[i]]<-.simplify.element(fit[[i]])
    }
  }
  fit
}
#seems to work for various situations now--things will get funky as you include thinning, though...

#' @export
thin.chains<-function(fit,thin=2){
  #process input
  if(thin<=1){
    return(fit)
  }
  if(length(dim(fit$chains))==2){
    simplified<-T
  }else{
    simplified<-F
  }
  if(simplified){
    for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
      fit[[i]]<-.expand.element(fit[[i]],simplify=T)
    }
  }
  #get indices of iterations to be included, modify sampler.control element accordingly
  niter<-fit$sampler.control$iter
  incl.inds<-setNames(rep(list(seq(1,niter,thin)),2),c('chains','sampler.params'))
  chains.len<-dim(fit$chains)[1]
  diags.len<-dim(fit$sampler.params)[1]
  incl.inds$chains<-incl.inds$chains-niter+chains.len
  incl.inds$sampler.params<-incl.inds$sampler.params-niter+diags.len
  incl.inds<-lapply(incl.inds,function(ii) ii[ii>0])
  fit$sampler.control$iter<-floor(fit$sampler.control$iter/thin)
  fit$sampler.control$warmup<-floor(fit$sampler.control$warmup/thin)
  fit$sampler.control$thin<-fit$sampler.control$thin*thin
  #thin the chains
  for(i in c('chains','sampler.params')){
    fit[[i]]<-.index.element(fit[[i]],incl.inds[[i]],1)
  }
  #redo parameter diagnostics (other than inits)
  fit$param.diags[2,,]<-apply(fit$chains,c(2,3),rstan::ess_bulk)
  fit$param.diags[3,,]<-apply(fit$chains,c(2,3),rstan::ess_tail)
  fit$param.diags[4,,]<-apply(fit$chains,c(2,3),rstan::Rhat)
  #redo any other elements
  for(i in names(fit)[names(fit)%in%c('quantiles','means','MAPs')]){
    fit[[i]]<-NULL
    fit[[i]]<-do.call(paste('.int.',i,sep=''),list(fit=fit,select='.|dev'))
    if(i=='means'|i=='MAPs'){
      fit[[i]]<-array(fit[[i]],dim(fit[[i]])[-1],dimnames(fit[[i]])[-1])
    }
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-apply(.int.chains(fit,'R_\\d+_dev'),c(2,3),function(ii) sum(ii>0)/length(ii))
  }
  #simplify as needed
  if(simplified){
    for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
      fit[[i]]<-.simplify.element(fit[[i]])
    }
  }
  fit
}
