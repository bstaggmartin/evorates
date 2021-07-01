#' Select chains in a fitted evolving rates model
#' 
#' 
#' This function extracts specific chains from the output of the function \code{fit.evorates}.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
#' @param chains A character or numeric vector. If of class character, the given text is matched to names of
#' the available chains in \code{fit}. If of class numeric, the function assumes \code{chains} refers to the
#' name "chain \code{chains}".
#' @param simplify TRUE or FALSE value: should the resulting elements of \code{fit} have dimensions of length
#' 1 collapsed and stored as attributes?
#' 
#' 
#' @return an object of class "\code{evorates_fit}". All previously-existing elements in \code{fit} will be
#' included.
#' 
#' 
#' @family chain management
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #get chain 1 and 2
#' select.chains(example.fit,c('chain 1','chain 2'))
#' #or
#' select.chains(example.fit,1:2)
#' 
#' 
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
#' Combine chains in a fitted evolving rates model
#' 
#' 
#' This function combines all chains from the output of the function \code{fit.evorates} in sequential order.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
#' @param simplify TRUE or FALSE value: should the resulting elements of \code{fit} have dimensions of length
#' 1 collapsed and stored as attributes?
#' 
#' 
#' @return an object of class "\code{evorates_fit}". All previously-existing elements in \code{fit} will be
#' included.
#' 
#' 
#' @details By sequential, I mean that the chains beginning of the second chain will follow the end of the
#' first chain, and so on. No permutation business is attempted. It is obviously unwise to run this function
#' with a \code{evorates_fit} including chains that have not yet reached a stationary distribution or include
#' obvious burn-in/warmup iterations.
#' 
#' The resulting chain name will simply be all the previous chain names separated by commas (e.g.,
#' "chain 1, chain 2, chain 3, ...).
#' 
#' Any iterations included in the sampler.params element but not in the chains element are
#' automatically excluded. Otherwise, issues would result from other functions which assuming any
#' difference between the number of iterations in sampler.params and chains are due to iterations missing from
#' the beginning of chains (i.e., burn-in/warmup). Any warmup iterations remaining in the chains element,
#' however, will be kept.
#' 
#' After combining chains, the resulting composite 'markov chain' no longer corresponds to any one set of
#' initialization values. As such, the 'inits' diagnostic from the param.diags element is set to NA.
#' 
#' 
#' @family chain management
#' 
#' 
#' @examples
#' #requires example fitted model object
#' combine.chains(example.fit)
#' 
#' 
#' @export
combine.chains<-function(fit,simplify=T){
  #process input
  if(length(dim(fit$chains))==2){
    if(simplify){
      return(fit)
    }else{
      for(i in names(fit)[!(names(fit)%in%c('sampler.control','call'))]){
        fit[[i]]<-.expand.element(fit[[i]],simplify=simplify)
      }
      chains.len<-dim(fit$chains)[1]
      diags.len<-dim(fit$sampler.params)[1]
      if(diags.len-chains.len>0){
        fit$sampler.params<-.index.element(fit$sampler.params,1:(diags.len-chains.len),1,T)
      }
      fit$sampler.control$warmup<-0
      fit$sampler.control$iter<-chains.len
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
    fit$sampler.control$warmup<-0
    fit$sampler.control$iter<-chains.len
    fit
  }
}

#for getting rid of warmup in chains and related elements
#' Exclude warmup in a fitted evolving rates model
#' 
#' 
#' This function gets rid of an arbitrary number of iterations at the beginning of all chains from the output
#' of the function \code{fit.evorates}.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
#' @param warmup the number of iterations to be excluded, counting from the 1st iteration, even if it isn't
#' stored in \code{fit}.
#' @param sampler TRUE or FALSE value: should the iterations be excluded from the \code{sampler.params} element
#' too?
#' 
#' 
#' @details Keep in \code{warmup} counts form the first iteration. For example, if \code{fit} has 2000
#' iterations with the first 1000 as warmup, you would set \code{warmup} to 1500 to keep the last 500
#' iterations.
#' 
#' 
#' @return an object of class "\code{evorates_fit}". All previously-existing elements in \code{fit} will be
#' included.
#' 
#' 
#' @family chain management
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #exclude warmup iterations from sampler
#' exclude.warmup(example.fit)
#' #exclude some non-warmup iterations
#' exclude.warmup(example.fit,1100)
#' #exclude some non-warmup iterations, but keep in sampler.params
#' exclude.warmup(example.fit,1100,F)
#' 
#' 
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
    fit$chains<-.index.element(fit$chains,1:(chains.len-target.len),1,TRUE)
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

#' Thin chains in a fitted evolving rates model
#' 
#' 
#' This function thins all chains from the output of the function \code{fit.evorates} by keeping only every
#' \code{thin}th iteration.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
#' @param thin 
#' 
#' 
#' @details Keep in \code{warmup} counts form the first iteration. For example, if \code{fit} has 2000
#' iterations with the first 1000 as warmup, you would set \code{warmup} to 1500 to keep the last 500
#' iterations.
#' 
#' 
#' @return an object of class "\code{evorates_fit}". All previously-existing elements in \code{fit} will be
#' included.
#' 
#' 
#' @family chain management
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #keep every 3rd iteration
#' thin.chains(example.fit,3)
#' 
#' 
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
