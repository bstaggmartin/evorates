#need to update %MAPs%, %chains%, %means%, %sampler% with updated functions to be more efficient

#' Extract posterior distributions from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting sampled parameter values for particular parameters
#' from output of the function \code{fit.corateBM}.
#'
#'
#' @param fit An object of class 'corateBM_fit', which is the typical output of a call to
#' \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to
#' parameter names using regular expressions. If of class numeric, the numbers are matched to indices
#' of rate parameters (\code{R_i}), include root (\code{R_0}). Exhibits special behavior with regards to
#' rate deviation (\code{R_i_dev}) and commas (','); see details below.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of iterations,
#' then parameters, then chains, collapsing any dimensions of length 1. When \code{select} matches to a
#' single parameter, the parameter name is stored as an attribute labelled "parameters".
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma (','), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore
#' ('_'). In practice, this helps extract parameters related to covariance matrices, as
#' 'trait1,trait2_evocov' is automatically converted to 'trait2,trait1_evocov'. This is helpful because
#' corateBM model outputs only store values for the lower triangle of the covariance matrix, since they are
#' symmetric. Can return unpredictable results in cases where ',' or '_' are part of trait names or
#' tip labels. Use double-backslashes ('\\,'/'\\_') to indicate 'escaped' commas and underscores that the
#' function should ignore during this swapping behavior. Note that you can also use '\\,' to do encode
#' expression searches of the form '{m,n}'.
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' '_dev' are by default excluded unless text with 'dev' is explicitly included at the end of a text
#' element the select argument.
#' 
#' 
#' @seealso
#' \code{\link{\%quantiles\%}}, \code{\link{\%means\%}}, \code{\link{\%MAPs\%}}, \code{\link{\%sampler\%}}
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
  fit[['chains']]<-expand.element(fit[['chains']])
  out<-int.op(fit,'chains',select)
  simplify.element(out)
}

#takes ~760 milliseconds! (holy shit)
#check if def.report.quantiles are in an already existing quantiles object
#if not, select parameters, then extract quantiles!
#improved speed 760-fold (als created helpful new function, simplify.array...)
#' Extract posterior distribution quantiles from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting quantiles of sampled parameter values for
#' particular parameters from output of the function \code{fit.corateBM}. Doesn't require a pre-stored
#' quantiles element.
#'
#'
#' @param fit An object of class 'corateBM_fit', which is the typical output of a call to
#' \code{fit.corateBM}.
#' @param select A list with two elements: 1) character or numeric vector. If of class character, the given
#' text is matched to parameter names using regular expressions. If of class numeric, the numbers are
#' matched to indices of rate parameters (\code{R_i}), including root (\code{R_0}). Exhibits special behavior
#' with regards to rate deviation parameters (\code{R_i_dev}) and commas (','); see details below. 2) A
#' numeric vector of values between 0 and 1 specifying the quantiles to extract. If unsupplied, quantiles
#' are taken from \code{fit$quantiles}, if it exists. If \code{fit$quantiles} doesn't exist, the operator
#' defaults to 2.5, 50, and 97.5% quantiles. If only a vector is supplied, the operator assumes that vector
#' corresponds to the first element of the list.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of quantiles,
#' then parameters, then chains, collapsing any dimensions of length 1. When the first element of
#' \code{select} matches to a single parameter, the parameter name is stored as an attribute labelled 
#' "parameters". When there is a single quantile being extracted, it is stored as an attribute labelled
#' "quantiles". When a single value is extracted, the result is named by its quantile with a 'parameters'
#' attribute.
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma (','), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore
#' ('_'). In practice, this helps extract parameters related to covariance matrices, as
#' 'trait1,trait2_evocov' is automatically converted to 'trait2,trait1_evocov'. This is helpful because
#' corateBM model outputs only store values for the lower triangle of the covariance matrix, since they are
#' symmetric. Can return unpredictable results in cases where ',' or '_' are part of trait names or
#' tip labels. Use double-backslashes ('\\,'/'\\_') to indicate 'escaped' commas and underscores that the
#' function should ignore during this swapping behavior. Note that you can also use '\\,' to do encode
#' expression searches of the form '{m,n}'.
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' '_dev' are by default excluded unless text with 'dev' is explicitly included at the end of a text
#' element the select argument.
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
#' @seealso
#' \code{\link{\%chains\%}}, \code{\link{\%means\%}}, \code{\link{\%MAPs\%}}, \code{\link{\%sampler\%}}
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
  if(is.null(fit[['quantiles']])){
    def.report.quantiles<-c(0.025,0.5,0.975)
  }else{
    fit[['quantiles']]<-expand.element(fit[['quantiles']])
    def.report.quantiles<-as.numeric(substr(dimnames(fit[['quantiles']])[[1]],0,
                                            nchar(dimnames(fit[['quantiles']])[[1]])-1))/100
  }
  if(is.list(select)){
    if(length(select)<2){
      select[[2]]<-def.report.quantiles
    }
  }else{
    select<-list(select,def.report.quantiles)
  }
  fit[['chains']]<-expand.element(fit[['chains']])
  tmp<-int.op(fit,'chains',select[[1]])
  out<-array(NA,dim=c(length(select[[2]]),dim(tmp)[-1]))
  dimnames(out)<-c('quantiles'=list(paste(select[[2]]*100,'%',sep='')),
                   dimnames(tmp)[-1])
  if(!is.null(fit[['quantiles']])){
    matches<-match(select[[2]],def.report.quantiles)
    out[(1:dim(out)[1])[!is.na(matches)],,]<-
      fit[['quantiles']][matches[!is.na(matches)],dimnames(tmp)[[2]],]
    select[[2]]<-select[[2]][is.na(matches)]
  }
  if(length(select[[2]])>0){
    out[is.na(out[,1,1]),,]<-apply(tmp,c(2,3),quantile,probs=select[[2]])
  }
  simplify.element(out)
}

#' Extract posterior distribution means from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting means of sampled parameter values for particular
#' parameters from output of the function \code{fit.corateBM}. Doesn't require a pre-stored means element.
#'
#'
#' @param fit An object of class 'corateBM_fit', which is the typical output of a call to
#' \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to
#' parameter names using regular expressions. If of class numeric, the numbers are matched to indices
#' of rate parameters (\code{R_i}), including root (\code{R_0}). Exhibits special behavior with regards to
#' rate deviation parameters (\code{R_i_dev}) and commas (','); see details below.
#' 
#' 
#' @return A numeric vector or matrix. The dimensions will always go in the order of parameters then chains,
#' collapsing any dimensions of length 1. When \code{select} matches to a single parameter, the parameter
#' name is stored as an attribute labelled "parameters".
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma (','), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore
#' ('_'). In practice, this helps extract parameters related to covariance matrices, as
#' 'trait1,trait2_evocov' is automatically converted to 'trait2,trait1_evocov'. This is helpful because
#' corateBM model outputs only store values for the lower triangle of the covariance matrix, since they are
#' symmetric. Can return unpredictable results in cases where ',' or '_' are part of trait names or
#' tip labels. Use double-backslashes ('\\,'/'\\_') to indicate 'escaped' commas and underscores that the
#' function should ignore during this swapping behavior. Note that you can also use '\\,' to do encode
#' expression searches of the form '{m,n}'.
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' '_dev' are by default excluded unless text with 'dev' is explicitly included at the end of a text
#' element the select argument.
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
#' @seealso
#' \code{\link{\%chains\%}}, \code{\link{\%quantiles\%}}, \code{\link{\%MAPs\%}}, \code{\link{\%sampler\%}}
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
  if(is.null(fit[['means']])){
    fit[['chains']]<-expand.element(fit[['chains']])
    tmp<-int.op(fit,'chains',select)
    out<-array(apply(tmp,c(2,3),mean),c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['means']]<-expand.element(fit[['means']])
    out<-int.op(fit,'means',select)
  }
  simplify.element(out)
}

#' Extract max a posteriori parameter estimates from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting sampled parameter values from iterations that exhibited
#' the highest posterior probability for each chain in the output of the function \code{fit.corateBM}.
#' Doesn't require a pre-stored MAPs element.
#'
#'
#' @param fit An object of class 'corateBM_fit', which is the typical output of a call to
#' \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to
#' parameter names using regular expressions. If of class numeric, the numbers are matched to indices
#' of rate parameters (\code{R_i}), including root (\code{R_0}). Exhibits special behavior with regards to
#' rate deviation parameters (\code{R_i_dev}) and commas (','); see details below.
#' 
#' 
#' @return A numeric vector or matrix. The dimensions will always go in the order of parameters then chains,
#' collapsing any dimensions of length 1. When \code{select} matches to a single parameter, the parameter
#' name is stored as an attribute labelled "parameters".
#' 
#' 
#' @details These operators become really useful with some knowledge on how to use regular expressions.
#' Please consult Rstudio's
#' \href{https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf}{Regular Expressions Cheat Sheet}
#' to learn more about regular expressions.
#' 
#' For convenience, when faced with text including a single comma (','), this operator will attempt to
#' switch the text preceding the comma with any text comeing after the comma and before an underscore
#' ('_'). In practice, this helps extract parameters related to covariance matrices, as
#' 'trait1,trait2_evocov' is automatically converted to 'trait2,trait1_evocov'. This is helpful because
#' corateBM model outputs only store values for the lower triangle of the covariance matrix, since they are
#' symmetric. Can return unpredictable results in cases where ',' or '_' are part of trait names or
#' tip labels. Use double-backslashes ('\\,'/'\\_') to indicate 'escaped' commas and underscores that the
#' function should ignore during this swapping behavior. Note that you can also use '\\,' to do encode
#' expression searches of the form '{m,n}'.
#' 
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' '_dev' are by default excluded unless text with 'dev' is explicitly included at the end of a text
#' element the select argument.
#' 
#' Note that finding the maximum a posteriori (MAP) parameter estimates is not the goal of Bayesian
#' inference per se, and that there is no guaruntee that MAP parameter estimates for a given chain will
#' even be close to the "true" MAP, even if the chain otherwise seems "healthy". In the package author's
#' experience, passing a correlated rate Brownian Motion model to the \code{rstan} function
#' \code{optimizing} indeed appears to frequently find narrow posterior probability 'peaks' that
#' substantially exceed the highest posterior probabilities found via \code{fit.corateBM}. These peaks
#' conincide with similar parameter estimates to the output of \code{fit.corateBM} but with inflated
#' \code{R_sig2} estimates (often substantially exceeding simulated values). Nonetheless, some users may
#' find that MAP parameter estimates provide interesting information about the posterior distribution
#' surface and/or the model fit, particularly in cases where they drastically differ from posterior
#' distribution means.
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
#' @seealso
#' \code{\link{\%chains\%}}, \code{\link{\%quantiles\%}}, \code{\link{\%means\%}}, \code{\link{\%sampler\%}}
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
  if(is.null(fit[['MAPs']])){
    fit[['chains']]<-expand.element(fit[['chains']])
    tmp<-int.op(fit,'chains',select)
    chains.len<-dim(tmp)[1]
    nchain<-dim(tmp)[3]
    lp<-as.matrix('%sampler%'(fit,'lp'))
    diags.len<-nrow(lp)
    if(diags.len-chains.len>0){
      lp<-as.matrix(lp[-(1:(diags.len-chains.len)),])
    }
    MAP.inds<-sapply(1:nchain, function(ii) which.max(lp[,ii]))
    out<-array(sapply(1:nchain,function(ii) tmp[MAP.inds[ii],,ii]),
               c(1,dim(tmp)[-1]),dimnames(tmp))
  }else{
    fit[['MAPs']]<-expand.element(fit[['MAPs']])
    out<-int.op(fit,'MAPs',select)
  }
  simplify.element(out)
}

#' Extract MCMC sampler parameters from a fitted correlated rates Brownian Motion model
#'
#'
#' This is an operator for efficiently extracting parameters used to tune the behavior of the stan-based
#' Markov chain Monte Carlo (MCMC) sampler from output of the function \code{fit.corateBM}. Helpful for
#' assessing problematic MCMC behavior.
#'
#'
#' @param fit An object of class 'corateBM_fit', which is the typical output of a call to
#' \code{fit.corateBM}.
#' @param select A character or numeric vector. If of class character, the given text is matched to
#' parameter names using regular expressions. If of class numeric, the numbers are taken simply taken
#' as indices. See details for names/order of available parameters.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array. The dimensions will always go in the order of iterations,
#' then parameters, then chains, collapsing any dimensions of length 1. When \code{select} matches to a
#' single parameter, the parameter name is stored as an attribute labelled "parameters". Note that
#' sampler parameters include warmup iterations. If output from the \code{\%chains\%} operator includes n
#' iterations, these will correspond to the last n iterations included in output from the \code{\%sampler\%}
#' opertator.
#' 
#' 
#' @details The names of the available parameters are 'accept_stat__', 'stepsize__', treedepth__',
#' 'n_leapfrog__', 'divergent__', 'energy__', and 'lp__', in that order. See Stan documentation for more
#' information on what these parameters mean.
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
#' @seealso
#' \code{\link{\%chains\%}}, \code{\link{\%quantiles\%}}, \code{\link{\%means\%}}, \code{\link{\%MAPs\%}},
#' \code{\link[rstan]{check_hmc_diagnostics}}
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
  fit[['sampler']]<-expand.element(fit[['diagnostics']][['sampler']])
  if(is.numeric(select)){
    select<-paste(c('accept_stat','stepsize','treedepth','n_leapfrog','diverget','energy','lp'),
                  '__',sep='')[select]
  }
  out<-int.op(fit,'sampler',select)
  simplify.element(out)
}

#internal operator -- select parameters out of fit$element based on edge numbers or regular expressions
#' @keywords internal
#' @export
int.op<-function(fit,element,select){
  select<-as.vector(select)
  param.names<-dimnames(fit[[element]])[[2]]
  if(is.numeric(select)){
    if(element=='sampler'){
      select<-param.names[select]
    }else{
      n.R<-suppressWarnings(
        max(as.numeric(substr(param.names,regexpr('R_',param.names)+2,nchar(param.names))),na.rm=T)
      )
      if(all(select>n.R|select<0)){
        stop("all numbers out of range for edge rate (R) parameters: are edge indices based on right tree?")
      }
      if(any(select>n.R|select<0)){
        warning(paste(select[select>n.R|select<0],collapse=', '),
                " out of range for edge rate (R) parameters: are edge indices based on right tree?")
        select<-select[select<=n.R&select>=0]
      }
      select<-paste('R_',sort(unique(select)),'$',sep='')
    }
  }else if(!is.character(select)){
    stop("the %",element,"% operator only accepts numeric or character vectors on right hand side")
  }
  select<-gsub('\\\\,',paste(rep('~',17),collapse=''),select)
  select<-gsub('\\\\_',paste(rep('@',17),collapse=''),select)
  if(any(grepl(',',select))){
    tmp<-strsplit(select,split=',')
    if(any(lengths(tmp)==3)){
      warning('text before and after comma not swapped for ',
              paste(select[which(lengths(tmp)==3)],collapse=', '),
              " since more than 1 comma was found: for any commas that are part of trait/tip names, please use '\\,' instead")
    }
    if(any(lengths(tmp)==2)){
      connect.inds<-which(lengths(tmp)==2)
      connect.inds<-cbind(connect.inds,length(select)+connect.inds)
      new.select<-tmp[connect.inds[,1]]
      new.select.l<-sapply(new.select,'[',1)
      new.select.r<-sapply(new.select,'[',2)
      new.select.r<-strsplit(new.select.r,'_')
      new.select<-sapply(1:length(new.select),
                         function(ii) paste(paste(new.select.r[[ii]][1],new.select.l[[ii]][1],sep=','),
                                            paste(new.select.r[[ii]][-1],collapse='_'),
                                            sep=if(length(new.select.r[[ii]])>1) '_' else ''))
      select<-c(select,new.select)
    }else{
      connect.inds<-NULL
    }
  }else{
    connect.inds<-NULL
  }
  select<-gsub('~{17}',',',select)
  select<-gsub('@{17}','_',select)
  tmp<-lapply(select,function(ii) grep(ii,param.names))
  forbidden.inds<-grep('_dev$',param.names)
  for(i in grep('dev$',select,invert=T)){
    tmp[[i]]<-tmp[[i]][!(tmp[[i]]%in%forbidden.inds)]
  }
  if(all(lengths(tmp)==0)){
    stop("couldn't find any corresponding parameters")
  }
  if(any(lengths(tmp)==0)){
    problem.select<-which(lengths(tmp)==0)
    if(!is.null(connect.inds)){
      for(i in 1:nrow(connect.inds)){
        if(any(lengths(tmp[connect.inds[i,]])!=0)){
          problem.select<-problem.select[-which(problem.select%in%connect.inds[i,])]
        }
      }
    }
    if(length(problem.select)>0){
      warning("couldn't find parameters corresponding to: ",
              unique(paste(select[problem.select],collapse=', ')))
    }
  }
  inds<-sort(unique(unlist(tmp)))
  index.array(fit[[element]],inds,2)
}
