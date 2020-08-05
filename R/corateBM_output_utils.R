####OPERATORS####


#need to update %MAPs%, %chains%, %means%, %sampler% with updated functions to be more efficient
#I think I'm done
#no need for %params% or %post.probs%--have comparison and diagnostic functions to take care of that

#' Extract posterior distributions from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting sampled parameter values for particular parameters from
#' output of the function \code{fit.corateBM}.
#'
#'
#' @param fit An object of class "\code{corateBM_fit}", the output of a call to \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to parameter
#' names using regular expressions. If of class numeric, the numbers are matched to indices of rate parameters
#' ("\code{R_i}"), including the root ("\code{R_0}"). Exhibits special behavior with regards to rate deviation
#' parameters ("\code{R_i_dev}") and commas ("\code{,}"); see details below.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of iterations,
#' then parameters, then chains, collapsing any dimensions of length 1 and storing associated name
#' information in attributes.
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma ("\code{,}"), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore 
#' ("\code{_}"). In practice, this helps extract parameters related to covariance matrices, as
#' "\code{trait1,trait2_evocov}" is automatically converted to "\code{trait2,trait1_evocov}". This is helpful
#' because corateBM model outputs only store values for the lower triangle of the covariance matrix, since
#' they are symmetric. Can return unpredictable results in cases where "\code{,}" or "\code{_}" are part of
#' trait names or tip labels. Use double-backslashes ("\code{\\\\,}" or "\code{\\\\_}") to indicate "escaped"
#' commas and underscores that the function should ignore during this swapping behavior. Note that you can
#' also use "\code{\\\\,}" to do encode regular expression searches of the form "\code{{m,n}}".
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' "\code{_dev}" are by default excluded unless text with "\code{dev}" is explicitly included at the end of a
#' text element the select argument.
#' 
#' 
#' @family corateBM_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object and tree
#' #get rate posterior distributions
#' example.fit%chains%1:nrow(example.tree$edge)
#' #',' behavior
#' example.fit%chains%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%chains%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%chains%c('.','dev') #all parameters
#' 
#' 
#' @export
`%chains%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %chains% operator only accepts fitted correlated rate BM fits (class 'corateBM_fit') on left hand side")
  }
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  out<-.int.chains(fit,select)
  .simplify.element(out)
}

#takes ~760 milliseconds! (holy shit)
#check if def.report.quantiles are in an already existing quantiles object
#if not, select parameters, then extract quantiles!
#improved speed 760-fold (als created helpful new function, .simplify.array...)
#' Extract posterior distribution quantiles from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting quantiles of sampled parameter values for particular
#' parameters from output of the function \code{fit.corateBM}. Doesn't require a pre-stored \code{quantiles}
#' element.
#'
#'
#' @param fit An object of class "\code{corateBM_fit}", the output of a call to \code{fit.corateBM}.
#' @param select A list with two elements: 1) A character or numeric vector. If of class character, the given
#' text is matched to parameter names using regular expressions. If of class numeric, the numbers are matched
#' to indices of rate parameters ("\code{R_i}"), including the root ("\code{R_0}"). Exhibits special behavior
#' with regards to rate deviation parameters ("\code{R_i_dev}") and commas ("\code{,}"); see details below.
#' 2) A numeric vector of values between 0 and 1 specifying the quantiles to extract. If unsupplied, quantiles
#' are taken from \code{fit$quantiles}, if it exists. If \code{fit$quantiles} doesn't exist, the operator
#' defaults to 2.5, 50, and 97.5% quantiles. If only a vector is supplied, the operator assumes that vector
#' corresponds to the first element of the list.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of quantiles,
#' then parameters, then chains, collapsing any dimensions of length 1 and storing associated name
#' information in attributes. When the output is of total length 1, defaults to being named according to
#' parameters and stores quantiles and chains information as attributes.
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma ("\code{,}"), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore 
#' ("\code{_}"). In practice, this helps extract parameters related to covariance matrices, as
#' "\code{trait1,trait2_evocov}" is automatically converted to "\code{trait2,trait1_evocov}". This is helpful
#' because corateBM model outputs only store values for the lower triangle of the covariance matrix, since
#' they are symmetric. Can return unpredictable results in cases where "\code{,}" or "\code{_}" are part of
#' trait names or tip labels. Use double-backslashes ("\code{\\\\,}" or "\code{\\\\_}") to indicate "escaped"
#' commas and underscores that the function should ignore during this swapping behavior. Note that you can
#' also use "\code{\\\\,}" to do encode regular expression searches of the form "\code{{m,n}}".
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' "\code{_dev}" are by default excluded unless text with "\code{dev}" is explicitly included at the end of
#' a text element the select argument.
#' 
#' 
#' @family corateBM_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object and tree
#' #get rate posterior distribution quantiles
#' example.fit%quantiles%1:nrow(example.tree$edge)
#' #',' behavior
#' example.fit%quantiles%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%quantiles%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%quantiles%c('.','dev') #all parameters
#' #non-default quantile extraction
#' example.fit%quantiles%list('R_0',0.432)
#' 
#' 
#' @export
`%quantiles%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %quantiles% operator only accepts fitted correlated rate BM fits (class 'corateBM_fit') on left hand side")
  }
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  out<-.int.quantiles(fit,select)
  .simplify.element(out)
}

#' Extract posterior distribution means from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting means of sampled parameter values for particular parameters
#' from output of the function \code{fit.corateBM}. Doesn't require a pre-stored \code{means} element.
#'
#'
#' @param fit An object of class "\code{corateBM_fit}", the output of a call to \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to parameter
#' names using regular expressions. If of class numeric, the numbers are matched to indices of rate parameters
#' ("\code{R_i}"), including the root ("\code{R_0}"). Exhibits special behavior with regards to rate deviation
#' parameters ("\code{R_i_dev}") and commas ("\code{,}"); see details below.
#' 
#' 
#' @return A numeric vector or matrix. The dimensions will always go in the order of parameters then chains,
#' collapsing any dimensions of length 1 and storing associated name information in attributes.
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma ("\code{,}"), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore 
#' ("\code{_}"). In practice, this helps extract parameters related to covariance matrices, as
#' "\code{trait1,trait2_evocov}" is automatically converted to "\code{trait2,trait1_evocov}". This is helpful
#' because corateBM model outputs only store values for the lower triangle of the covariance matrix, since
#' they are symmetric. Can return unpredictable results in cases where "\code{,}" or "\code{_}" are part of
#' trait names or tip labels. Use double-backslashes ("\code{\\\\,}" or "\code{\\\\_}") to indicate "escaped"
#' commas and underscores that the function should ignore during this swapping behavior. Note that you can
#' also use "\code{\\\\,}" to do encode regular expression searches of the form "\code{{m,n}}".
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' "\code{_dev}" are by default excluded unless text with "\code{dev}" is explicitly included at the end of a
#' text element the select argument.
#' 
#' 
#' @family corateBM_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object and tree
#' #get rate posterior distributions
#' example.fit%means%1:nrow(example.tree$edge)
#' #',' behavior
#' example.fit%means%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%means%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%means%c('.','dev') #all parameters
#' 
#' 
#' @export
`%means%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %means% operator only accepts fitted correlated rate BM fits (class 'corateBM_fit') on left hand side")
  }
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  out<-.int.means(fit,select)
  .simplify.element(out)
}

#' Extract max a posteriori parameter estimates from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting sampled parameter values from iterations that exhibited the
#' highest posterior probability for each chain in the output of the function \code{fit.corateBM}. Doesn't
#' require a pre-stored \code{MAPs} element.
#'
#'
#' @param fit An object of class "\code{corateBM_fit}", the output of a call to \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to parameter
#' names using regular expressions. If of class numeric, the numbers are matched to indices of rate parameters
#' ("\code{R_i}"), including the root ("\code{R_0}"). Exhibits special behavior with regards to rate deviation
#' parameters ("\code{R_i_dev}") and commas ("\code{,}"); see details below.
#' 
#' 
#' @return A numeric vector or matrix. The dimensions will always go in the order of parameters then chains,
#' collapsing any dimensions of length 1 and storing associated name information in attributes.
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma ("\code{,}"), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore 
#' ("\code{_}"). In practice, this helps extract parameters related to covariance matrices, as
#' "\code{trait1,trait2_evocov}" is automatically converted to "\code{trait2,trait1_evocov}". This is helpful
#' because corateBM model outputs only store values for the lower triangle of the covariance matrix, since
#' they are symmetric. Can return unpredictable results in cases where "\code{,}" or "\code{_}" are part of
#' trait names or tip labels. Use double-backslashes ("\code{\\\\,}" or "\code{\\\\_}") to indicate "escaped"
#' commas and underscores that the function should ignore during this swapping behavior. Note that you can
#' also use "\code{\\\\,}" to do encode regular expression searches of the form "\code{{m,n}}".
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' "\code{_dev}" are by default excluded unless text with "\code{dev}" is explicitly included at the end of a
#' text element the select argument.
#' 
#' Note that finding the maximum a posteriori (MAP) parameter estimates is not the goal of Bayesian inference
#' per se, and that there is no guaruntee that MAP parameter estimates for a given chain will even be close to
#' the "true" MAP, even if the chain otherwise seems "healthy". In the package author's experience, passing a
#' correlated rate Brownian Motion model to the \code{rstan} function \code{optimizing} indeed appears to
#' frequently find narrow posterior probability 'peaks' that substantially exceed the highest posterior
#' probabilities found via \code{fit.corateBM}. These peaks conincide with similar parameter estimates to the
#' output of \code{fit.corateBM} but with inflated \code{R_sig2} estimates (often substantially exceeding
#' simulated values). Nonetheless, some users may find that MAP parameter estimates provide interesting
#' information about the posterior distribution surface and/or the model fit, particularly in cases where
#' they drastically differ from posterior distribution means.
#' 
#' 
#' @family corateBM_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object and tree
#' #get rate posterior distributions
#' example.fit%MAPs%1:nrow(example.tree$edge)
#' #',' behavior
#' example.fit%MAPs%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%MAPs%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%MAPs%c('.','dev') #all parameters
#' 
#' 
#' @export
`%MAPs%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %MAPs% operator only accepts fitted correlated rate BM fits (class 'corateBM_fit') on left hand side")
  }
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  out<-.int.MAPs(fit,select)
  .simplify.element(out)
}

#' Extract MCMC sampler parameters from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting parameters used to tune the behavior of the stan-based
#' Markov chain Monte Carlo (MCMC) sampler from output of the function \code{fit.corateBM}. Helpful for
#' assessing problematic MCMC behavior.
#'
#'
#' @param fit An object of class "\code{corateBM_fit}", the output of a call to \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to parameter
#' names using regular expressions. If of class numeric, the numbers are taken simply taken as indices. See
#' details for names/order of available parameters.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of quantiles,
#' then parameters, then chains, collapsing any dimensions of length 1 and storing associated name
#' information in attributes. Note that sampler parameters include warmup iterations. If output from the
#' \code{\%chains\%} operator includes \code{n} iterations, these will correspond to the last \code{n}
#' iterations included in output from the \code{\%sampler\%} operator.
#' 
#' 
#' @details The names of the available parameters are "\code{accept_stat__}", "\code{stepsize__}",
#' "\code{treedepth__}", "\code{n_leapfrog__}", "\code{divergent__}", "\code{energy__}", and "\code{lp__}",
#' in that order. See Stan documentation for more information on what these parameters mean.
#' 
#' 
#' @family corateBM_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #get log posterior probability (+ some arbitrary constant) of iterations in mcmc sampler
#' example.fit%sampler%'lp'
#' #same result using numeric select
#' example.fit%sampler%7
#' #get all sampler parameters
#' example.fit%sampler%1:7
#' 
#' 
#' @seealso \code{\link[rstan]{check_hmc_diagnostics}}
#' 
#' 
#' @export
`%sampler%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %sampler% operator only accepts fitted correlated rate BM fits (class 'corateBM_fit') on left hand side")
  }
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  out<-.int.sampler(fit,select)
  .simplify.element(out)
}

#' Extract posterior distribution diagnostics from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting diagnostic summary statistics of posterior distributions for
#' particular parameters from output of the function \code{fit.corateBM}.
#'
#'
#' @param fit An object of class "\code{corateBM_fit}", the output of a call to \code{fit.corateBM}.
#' @param select A list with two elements: 1) A character or numeric vector. If of class character, the given
#' text is matched to parameter names using regular expressions. If of class numeric, the numbers are matched
#' to indices of rate parameters ("\code{R_i}"), including the root ("\code{R_0}"). Exhibits special behavior
#' with regards to rate deviation parameters ("\code{R_i_dev}") and commas ("\code{,}"); see details below.
#' 2) A character or numeric vector. If of class character, the given text is matched to diagnostic names
#' using regular expressions (with no special behavior). If of class numeric, the numbers are taken simply
#' taken as indices. The names of the available diagnostics are "\code{inits}", "\code{bulk_ess}", 
#' "\code{tail_ess}", and "\code{Rhat}", in that order. See details for definitions of what these diagnostics
#' mean. If only a vector (or list of length 1) is supplied, the operator assumes that vector corresponds to
#' the first element of the list and defaults to extracting all available diagnostics.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of diagnostics,
#' then parameters, then chains, collapsing any dimensions of length 1 and storing associated name
#' information in attributes. When the output is of total length 1, defaults to being named according to
#' parameters and stores diagnostics and chains information as attributes.
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma ("\code{,}"), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore 
#' ("\code{_}"). In practice, this helps extract parameters related to covariance matrices, as
#' "\code{trait1,trait2_evocov}" is automatically converted to "\code{trait2,trait1_evocov}". This is helpful
#' because corateBM model outputs only store values for the lower triangle of the covariance matrix, since
#' they are symmetric. Can return unpredictable results in cases where "\code{,}" or "\code{_}" are part of
#' trait names or tip labels. Use double-backslashes ("\code{\\\\,}" or "\code{\\\\_}") to indicate "escaped"
#' commas and underscores that the function should ignore during this swapping behavior. Note that you can
#' also use "\code{\\\\,}" to do encode regular expression searches of the form "\code{{m,n}}".
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' "\code{_dev}" are by default excluded unless text with "\code{dev}" is explicitly included at the end of
#' a text element the select argument.
#' 
#' Insert explanation of diagnostics here
#' 
#' 
#' @family corateBM_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object and tree
#' #get all rate posterior distribution diagnostics
#' example.fit%diagnostics%1:nrow(example.tree$edge)
#' #',' behavior
#' example.fit%diagnostics%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%diagnostics%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%diagnostics%c('.','dev') #all parameters
#' #specific diagnostics extraction
#' example.fit%diagnostics%list('R_0','0.432,'inits') #initial values
#' example.fit%diagnostics%list('R_0','0.432,'ess') #bulk and tail effective sample sizes
#' 
#' 
#' @export
`%diagnostics%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %diagnostics% operator only accepts fitted correlated rate BM fits (class 'corateBM_fit') on left hand side")
  }
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  out<-.int.diagnostics(fit,select)
  .simplify.element(out)
}


####CHAIN MANAGEMENT####


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
      warning('desired number of iterations (',target.len,') is longer than number of available iterations (',diags.len,') in sampler.params element: no iterations excluded')
    }
  }
  #trim chains as needed if desired
  chains.len<-dim(fit$chains)[1]
  if(chains.len-target.len<=0){
    if(chains.len-target.len<0){
      warning('desired number of iterations (',target.len,') is longer than number of available iterations (',chains.len,') in chains element: no iterations excluded')
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
  chains.len<-dim(fit$chains)[1]
  diags.len<-dim(fit$sampler.params)[1]
  incl.inds<-rep(list(seq(1,chains.len,thin)),2)
  names(incl.inds)<-c('chains','sampler.params')
  diff.iter<-chains.len-length(incl.inds$chains)
  true.warmup<-chains.len-fit$sampler.control$iter+fit$sampler.control$warmup
  if(true.warmup>0){
    diff.warmup<-diff.iter-length(incl.inds$chains[incl.inds$chains>true.warmup])
  }else{
    diff.warmup<-0
  }
  fit$sampler.control$iter<-fit$sampler.control$iter-diff.iter
  fit$sampler.control$warmup<-fit$sampler.control$warmup-diff.warmup
  #account for any iterations in parameter diagnostics not inlcuded in chains
  if(diags.len-chains.len>0){
    incl.inds$sampler.params<-c(1:(diags.len-chains.len),incl.inds$sampler.params+diags.len-chains.len)
  }
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

####OTHER####


#' @export
get.cov.mat<-function(fit,traits=colnames(fit$call$trait.data),
                      element=c('chains','quantiles','means','MAPs','diagnostics'),
                      type=c('evocov','intracov'),
                      output.list=F,select.extra=NULL,
                      simplify=T){
  if(!inherits(fit,'corateBM_fit')){
    stop("fit must be a fitted correlated rates BM fit (class 'corateBM_fit')")
  }
  if(ncol(fit$call$trait.data)==1){
    stop('fit appears to be for univariate trait dataset: there are no trait covariance parameters')
  }
  try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics')),silent=T)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', or 'diagnostics'")
  }
  element<-try.element
  try.type<-try(match.arg(type,c('evocov','intracov')),silent=T)
  if(inherits(try.type,'try-error')){
    stop(type," is not an available option for covariance type: please specifiy either 'evocov' for evolutionary covariance or 'intracov' for intraspecific covariance")
  }
  type<-try.type
  if(type=='intracov'&sum(grepl('intracov',dimnames(fit$chains)[['parameters']]))==0){
    warning("covariance type 'intracov' selected, but no intraspecific variance modeled in corateBM_fit: defaulted to 'evocov'")
    type<-'evocov'
  }
  if(is.numeric(traits)){
    traits<-colnames(fit$call$trait.data)[traits]
  }
  #check if any trait names not available
  traits.exist<-traits%in%colnames(fit$call$trait.data)
  if(all(!traits.exist)){
    stop('none of the specified traits found')
  }
  if(any(!traits.exist)){
    warning(paste(traits[which(!traits.exist)],collapse=', '),' not found')
    traits<-traits[which(traits.exist)]
  }
  k<-length(traits)
  new.traits<-gsub(',','\\\\,',traits)
  new.traits<-gsub('_','\\\\_',new.traits)
  select<-paste(rep(new.traits,each=k),',',
                rep(new.traits,k),paste('_',type,sep=''),
                sep='')[lower.tri(matrix(NA,k,k),diag=T)]
  if(element=='quantiles'|element=='diagnostics'&!is.null(select.extra)){
    select<-list(select,select.extra)
  }
  tmp<-do.call(paste('.int.',element,sep=''),list(fit=fit,select=select))
  out<-array(NA,c(k,k,dim(tmp)[1],dim(tmp)[3]),c(rep(list(parameters=traits),2),dimnames(tmp)[c(1,3)]))
  for(i in 1:k){
    for(j in 1:i){
      tmp.name<-grep(paste(paste(traits[c(i,j)],',',traits[c(j,i)],'_',type,sep=''),collapse='|'),
                     dimnames(tmp)[[2]])
      out[i,j,,]<-tmp[,tmp.name,]
      if(i!=j){
        out[j,i,,]<-tmp[,tmp.name,]
      }
    }
  }
  if(output.list){
    out<-asplit(out,4)
    out<-lapply(out,function(ii) if(dim(ii)[3]==1) ii else asplit(ii,3))
    if(length(out)==1){
      attr(out[[1]],'chains')<-names(out)
      out<-out[[1]]
    }
  }else if(simplify){
    new.dims<-c(dim(out)[1:2],ifelse(dim(out)[3:4]==1,NA,dim(out)[3:4]))
    new.dimnames<-dimnames(out)[!is.na(new.dims)]
    new.out<-array(out,new.dims[!is.na(new.dims)],new.dimnames)
    for(i in which(is.na(new.dims))){
      attr(new.out,names(dimnames(out))[i])<-dimnames(out)[[i]]
    }
    out<-new.out
  }
  out
}
#7/27: %quantiles% for whole large_fit array is only longer than a simple apply implementation by
#~200 milliseconds...