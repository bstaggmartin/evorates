#almost certainly needs cleaning up

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
select.chains<-function(fit,chains,simplify=TRUE){
  #input processing
  par.inds<-which(names(fit)!='call'&names(fit)!='sampler.control')
  for(i in par.inds){
    fit[[i]]<-.expand.par(fit[[i]])
  }
  chain.nms<-dimnames(fit$chains)[[3]]
  if(is.numeric(chains)){
    chains<-chain.nms[chains]
  }
  inds<-pmatch(chains,chain.nms)
  probs<-is.na(inds)
  if(any(probs)){
    inds<-inds[!probs]
    warning('Some chain selections were out of bounds and ignored')
  }
  fit$sampler.control$chains<-length(chains)
  #subset chains
  for(i in par.inds){
    type<-.get.par.type(fit[[i]])
    fit[[i]]<-.add.par.class(fit[[i]][,,inds,drop=FALSE])
    attr(fit[[i]],'param_type')<-type
    if(simplify){
      fit[[i]]<-.simplify.par(fit[[i]])
    }
  }
  fit
}

#automatically excludes iterations not included in chains and inits--doesn't make sense anymore
#would break for a BM fit with 1 parameter--need to work on this
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
combine.chains<-function(fit,simplify=TRUE){
  #input processing
  par.inds<-which(names(fit)!='call'&names(fit)!='sampler.control')
  nchains<-fit$sampler.control$chains
  fit$sampler.control$chains<-1
  fit$sampler.control$warmup<-0
  #warmup won'ts always be excluded...below is safer
  fit$chains<-.expand.par(fit$chains)
  fit$sampler.control$iter<-dim(fit$chains)[1]*nchains
  
  if(nchains==1){ #handle cases where there's only 1 chain already
    #trim sampler.params if sampler.params still includes warmup iterations not in chains
    sampler.params<-.make.par.3D(fit[['sampler.params']])
    diag.len<-dim(sampler.params)[1]
    targ.len<-dim(.make.par.3D(fit[['chains']]))[1]
    diff.len<-diag.len-targ.len
    if(diff.len){
      sampler.params<-.add.par.class(sampler.params[-seq_len(diff.len),,,drop=FALSE])
      attr(sampler.params,'param_type')<-'chains'
      if(simplify){
        sampler.params<-.simplify.par(sampler.params)
      }
      fit[['sampler.params']]<-sampler.params
    }
    #expand as necessary
    if(!simplify){
      for(i in par.inds){
        fit[[i]]<-.expand.par(fit[[i]])
      }
    }
  }else{ #handle cases where chains do need to be combined
    full<-.call.op('chains',fit,c('.',names(fit$sampler.params)))
    outdims<-dim(full)
    outdims[1]<-outdims[1]*outdims[3]
    outdims[3]<-1
    outnames<-dimnames(full)
    outnames[[3]]<-paste(outnames[[3]],collapse=', ')
    outdims<-outdims[c(1,3,2)]
    outnames<-outnames[c(1,3,2)]
    full<-aperm(array((aperm(full,c(1,3,2))),outdims,outnames),c(1,3,2))
    break.ind<-outdims[3]-9
    fit$chains<-.add.par.class(full[,seq_len(break.ind),,drop=FALSE])
    fit$sampler.params<-.add.par.class(full[,seq_len(9)+break.ind,,drop=FALSE])
    attr(fit$chains,'param_type')<-attr(fit$sampler.params,'param_type')<-'chains'
    for(i in c('quantiles','diagnostics')){
      if(!is.null(fit[[i]])){
        extra.select<-dimnames(.make.par.3D(fit[[i]]))[[1]]
        fit[[i]]<-NULL
        fit[[i]]<-.call.op(i,fit,list('.',extra.select),FALSE)
        if(i=='diagnostics'){
          fit[[i]][extra.select=='inits',,]<-NA
        }
      }
    }
    #just use averaging for mean-based param_blocks...
    outnames<-outnames[c(1,3,2)]
    for(i in c('means','post.probs')){
      if(!is.null(fit[[i]])){
        nms<-names(fit[[i]])
        outdims<-c(1,length(nms),1)
        outnames[[2]]<-nms
        fit[[i]]<-.add.par.class(array(.rowMeans(fit[[i]],outdims[2],nchains),outdims,outnames))
        attr(fit[[i]],'param_type')<-'means'
      }
    }
  }
  
  #simplify as needed
  if(simplify){
    for(i in par.inds){
      fit[[i]]<-.simplify.par(fit[[i]])
    }
  }
  
  #output
  fit
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
#' @details Keep in mind \code{warmup} counts form the first iteration. For example, if \code{fit} has 2000
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
exclude.warmup<-function(fit,warmup=fit$sampler.control$warmup,sampler=TRUE,simplify=TRUE){
  #process input
  fit$sampler.control$warmup<-warmup
  niter<-fit$sampler.control$iter
  targ.len<-niter-warmup
  par.inds<-match('chains',names(fit))
  if(sampler){
    fit$sampler.control$warmup<-0
    fit$sampler.control$iter<-targ.len
    par.inds<-c(par.inds,match('sampler.params',names(fit)))
  }
  
  fit<-.proc.dim1.mods(fit,par.inds,targ.len,NULL,NULL)
  
  #simplify as needed
  if(simplify){
    par.inds<-which(names(fit)!='call'&names(fit)!='sampler.control')
    for(i in par.inds){
      fit[[i]]<-.simplify.par(fit[[i]])
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
thin.chains<-function(fit,thin=2,simplify=TRUE){
  #process input
  niter<-fit$sampler.control$iter
  thin<-round(thin[1])
  thin.inds<-seq.int(1,niter,thin)
  fit$sampler.control$iter<-floor(fit$sampler.control$iter/thin)
  fit$sampler.control$warmup<-floor(fit$sampler.control$warmup/thin)
  fit$sampler.control$thin<-fit$sampler.control$thin*thin
  par.inds<-match(c('chains','sampler.params'),names(fit))
  
  fit<-.proc.dim1.mods(fit,par.inds,NULL,thin.inds,niter)
  
  #simplify as needed
  if(simplify){
    par.inds<-which(names(fit)!='call'&names(fit)!='sampler.control')
    for(i in par.inds){
      fit[[i]]<-.simplify.par(fit[[i]])
    }
  }
  fit
}

.proc.dim1.mods<-function(fit,par.inds,targ.len,thin.inds,niter){
  for(i in par.inds){
    fit[[i]]<-.make.par.3D(fit[[i]])
    diag.len<-dim(fit[[i]])[1]
    if(is.null(thin.inds)){
      diff.len<-diag.len-targ.len
      if(diff.len<0){
        stop('Number of iterations to exclude exceeds the total number of iterations in ',
             if(i==1) 'chains' else 'sampling parameters')
      }
      inds<-if(diff.len) -seq_len(diff.len) else substitute()
    }else{
      tmp<-niter-diag.len
      inds<-thin.inds[thin.inds>tmp]-tmp
    }
    fit[[i]]<-.add.par.class(fit[[i]][inds,,,drop=FALSE])
    attr(fit[[i]],'param_type')<-'chains'
  }
  for(i in c('quantiles','means','diagnostics')){
    if(!is.null(fit[[i]])){
      extra.select<-dimnames(.make.par.3D(fit[[i]]))[[1]]
      if(i=='diagnostics'){
        tmp.inds<-extra.select=='inits'
        inits.save<-fit[[i]][tmp.inds,,,drop=FALSE]
      }
      fit[[i]]<-NULL
      fit[[i]]<-.call.op(i,fit,list('.',extra.select),FALSE)
      if(i=='diagnostics'){
        fit[[i]][tmp.inds,,]<-inits.save
      }
    }
  }
  if(!is.null(fit$post.probs)){
    Rdevs<-.call.op('chains',fit,'^Rdev_[1-9][0-9]*$|Rdev_[1-9][0-9]*$')
    Rdevs[Rdevs==0]<-NA
    Rdevs<-Rdevs>0
    #add rate deviation posterior probabilities
    #should never get instances where deviations are perfectly 0...but just in case
    fit$post.probs<-.call.op('means',list(chains=Rdevs),'.',FALSE)
    fit$post.probs[is.infinite(fit$post.probs)|is.nan(fit$post.probs)]<-0.5
  }
  fit
}