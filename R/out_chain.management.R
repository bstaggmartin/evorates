#almost certainly needs cleaning up

#' Select chains in a fitted evorates model
#' 
#' 
#' This function subsets fitted evorates model to specific chains.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param chains A character or numeric vector specifying which chains to select. Notably, \code{NA}
#' (or out of bound) selections are not allowed and automatically ignored with a warning.
#' @param simplify \code{TRUE} or \code{FALSE}: should the resulting \code{param_block} array(s) in
#' \code{fit} be  simplified? If \code{TRUE} (the default), dimensions of length 1 in the result are
#' automatically collapsed, with corresponding information stored as attributes (this is the default
#' behavior of param_block operators).
#' 
#' @return An object of class "\code{evorates_fit}". All previously-existing param_block arrays stored
#' in \code{fit} will be included.
#' 
#' 
#' @family chain management
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #three ways to get the second and first chain
#' fit <- select.chains(cet_fit, c('chain 1','chain 2'))
#' fit <- select.chains(cet_fit, c(1,2))
#' fit <- select.chains(cet_fit, -c(3,4))
#' 
#' #note handling of NAs
#' fit <- select.chains(cet_fit, c('chain 1', NA))
#' fit <- select.chains(cet_fit, c('chain 1', 'chain 5'))
#' fit <- select.chains(cet_fit, c(1,5))
#' 
#' #here's what happens when you do vs. don't simplify the result
#' select.chains(cet_fit, 1)$means; select.chains(cet_fit, 1, simplify = FALSE)$means
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

#' Combine chains in a fitted evorates model
#' 
#' 
#' This function combines all chains in fitted evorates model sequentially into a single, larger chain.
#' Generally, this function should only be run after confirming that (via Rhat diagnostics) that each
#' chain adequately converged!
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param simplify \code{TRUE} or \code{FALSE}: should the resulting \code{param_block} array(s) in
#' \code{fit} be  simplified? If \code{TRUE} (the default), dimensions of length 1 in the result are
#' automatically collapsed, with corresponding information stored as attributes (this is the default
#' behavior of param_block operators).
#' 
#' 
#' @return An object of class "\code{evorates_fit}". All previously-existing \code{param_block} arrays stored
#' in \code{fit} will be included.
#' 
#' 
#' @details Chains are combined sequentially, meaning that the beginning of the second chain will follow
#' the end of the first chain, the third will follow the second, and so on--no permutation business is
#' attempted. The resulting chain name will simply be all the previous chain names separated by commas
#' (e.g., "\code{chain 1, chain 2, chain 3, <...>}"). Because of this, initial values in the parameter
#' diagnostics param_block (if present) are set to \code{NA} and any remaining "warmup iterations"
#' (defined as iterations present in \code{fit$sampler.params} but not in \code{fit$chains}) are
#' removed, with all retained iterations reclassified (perhaps misleadingly in some cases)
#' as non-warmup iterations.
#' 
#' 
#' @family chain management
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #make sure all chains co
#' 
#' 
#' @export
combine.chains<-function(fit,simplify=TRUE){
  #input processing
  par.inds<-which(names(fit)!='call'&names(fit)!='sampler.control')
  nchains<-fit$sampler.control$chains
  fit$sampler.control$chains<-1
  fit$sampler.control$warmup<-0
  #warmup won't always be excluded...below is safer
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
      fit[[i]]<-.make.par.3D(fit[[i]])
      extra.select<-dimnames(fit[[i]])[[1]]
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