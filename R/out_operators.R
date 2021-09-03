#' Extract posterior distributions from a fitted Evolving Rates model
#'
#'
#' This is an operator for extracting sampled values for particular parameters from the output of \code{fit.evorates}
#' or \code{output.evorates}
#'
#'
#' @param fit An object of class "\code{evorates_fit}"
#' @param select A list with two elements: 1) A character or numeric vector. If of class character, the given
#' text is matched to parameter
#' names using regular expressions. If of class numeric, the numbers are matched to indices of rate parameters
#' ("\code{R_i}"), including the root ("\code{R_0}"). Exhibits special behavior with regards to rate deviation
#' parameters ("\code{R_i_dev}") and commas ("\code{,}"); see details below. 2) A vector of integers specifying
#' particular iterations to extract from the chain. If unsupplied, all iterations will be extracted.
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
#' Since rate deviation parameters are oftentimes redundant with rate parameters, parameter names ending in
#' "\code{_dev}" are by default excluded unless text with "\code{dev}" is explicitly included at the end of a
#' text element the select argument.
#' 
#' (NOTE: THIS BEHAVIOR STILL EXISTS, BUT MULTIVARIATE EVORATES MODELS ARE DEPRECIATED FOR THE TIME BEING, SO NO
#' EVORATES MODELS ACTUALLY INCLUDE COVARIANCE MATRIX PARAMETERS...YET) For convenience, when faced with text
#' including a single comma ("\code{,}"), this operator will attempt to
#' switch the text preceding the comma with any text coming after the comma and before an underscore 
#' ("\code{_}"). In practice, this helps extract parameters related to covariance matrices, as
#' "\code{trait1,trait2_evocov}" is automatically converted to "\code{trait2,trait1_evocov}". This is helpful
#' because evorates model outputs only store values for the lower triangle of the covariance matrix, since
#' they are symmetric. Can return unpredictable results in cases where "\code{,}" or "\code{_}" are part of
#' trait names or tip labels. Use double-backslashes ("\code{\\\\,}" or "\code{\\\\_}") to indicate "escaped"
#' commas and underscores that the function should ignore during this swapping behavior. Note that you can
#' also use "\code{\\\\,}" to do encode regular expression searches of the form "\code{{m,n}}".
#' 
#' 
#' @family evorates_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #get rate posterior distributions
#' example.fit%chains%1:nrow(example.fit$call$tree$edge)
#' #',' behavior
#' example.fit%chains%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%chains%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%chains%c('.','dev') #all parameters
#' #non-default iteration extraction
#' example.fit%xchains%list('R_0',c(23,45,52))
#' 
#' 
#' @export
`%chains%`<-function(fit,select){
  is.char<-try(is.character(select),silent=TRUE)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  flag<-FALSE
  if(.is.evorates.element(fit)){
    type<-.get.element.type(fit)
    if(type!='chains'){
      flag<-TRUE
    }else{
      fit<-list(chains=fit)
    }
  }else if(!inherits(fit,'evorates_fit')){
    flag<-TRUE
  }
  if(flag){
    stop("the %chains% operator only accepts fitted evolving rates model fits (class 'evorates_fit') or loose chains elements on left hand side")
  }
  out<-.int.chains(fit,select)
  .simplify.element(out)
}

#takes ~760 milliseconds! (holy shit)
#check if def.report.quantiles are in an already existing quantiles object
#if not, select parameters, then extract quantiles!
#improved speed 760-fold (als created helpful new function, .simplify.array...)

#ALTER THIS TO ALLOW STORAGE OF LOOSE QUANTILE ELEMENTS SO PARTICULAR QUANTILES CAN BE EXTRACTED

#' Extract posterior distribution quantiles from a fitted evolving rates model
#'
#'
#' This is an operator for efficiently extracting quantiles of sampled parameter values for particular
#' parameters from output of the function \code{fit.evorates}. Doesn't require a pre-stored \code{quantiles}
#' element.
#'
#'
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
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
#' because evorates model outputs only store values for the lower triangle of the covariance matrix, since
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
#' @family evorates_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #get rate posterior distribution quantiles
#' example.fit%quantiles%1:nrow(example.fit$call$tree$edge)
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
  is.char<-try(is.character(select),silent=TRUE)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  flag<-FALSE
  if(.is.evorates.element(fit)){
    type<-.get.element.type(fit)
    if(type!='chains'&type!='quantiles'){
      flag<-TRUE
    }else{
      if(type=='chains'){
        fit<-list(chains=fit)
      }else{
        fit<-list(quantiles=fit)
      }
    }
  }else if(!inherits(fit,'evorates_fit')){
    flag<-TRUE
  }
  if(flag){
    stop("the %chains% operator only accepts fitted evolving rates model fits (class 'evorates_fit') or loose chains/quantiles elements on left hand side")
  }
  out<-.int.quantiles(fit,select)
  .simplify.element(out)
}

#' Extract posterior distribution means from a fitted evolving rates model
#'
#'
#' This is an operator for efficiently extracting means of sampled parameter values for particular parameters
#' from output of the function \code{fit.evorates}. Doesn't require a pre-stored \code{means} element.
#'
#'
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
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
#' because evorates model outputs only store values for the lower triangle of the covariance matrix, since
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
#' @family evorates_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #get rate posterior distributions
#' example.fit%means%1:nrow(example.fit$call$tree$edge)
#' #',' behavior
#' example.fit%means%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%means%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%means%c('.','dev') #all parameters
#' 
#' 
#' @export
`%means%`<-function(fit,select){
  is.char<-try(is.character(select),silent=TRUE)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  flag<-FALSE
  if(.is.evorates.element(fit)){
    type<-.get.element.type(fit)
    if(type!='chains'){
      flag<-TRUE
    }else{
      fit<-list(chains=fit)
    }
  }else if(!inherits(fit,'evorates_fit')){
    flag<-TRUE
  }
  if(flag){
    stop("the %means% operator only accepts fitted evolving rates model fits (class 'evorates_fit') or loose chains elements on left hand side")
  }
  out<-.int.means(fit,select)
  .simplify.element(out)
}

#%sampler% subsumed into other functions, %MAPs% straight-up deleted

#' Extract posterior distribution diagnostics from a fitted evolving rates model
#'
#'
#' This is an operator for efficiently extracting diagnostic summary statistics of posterior distributions for
#' particular parameters from output of the function \code{fit.evorates}.
#'
#'
#' @param fit An object of class "\code{evorates_fit}", the output of a call to \code{fit.evorates}.
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
#' because evorates model outputs only store values for the lower triangle of the covariance matrix, since
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
#' @family evorates_fit operators
#' 
#' 
#' @examples
#' #requires example fitted model object
#' #get all rate posterior distribution diagnostics
#' example.fit%diagnostics%1:nrow(example.fit$call$tree$edge)
#' #',' behavior
#' example.fit%diagnostics%'X2,X1_evocov'
#' #'dev' behavior
#' example.fit%diagnostics%'.' #all parameters EXCEPT for rate deviation parameters
#' example.fit%diagnostics%c('.','dev') #all parameters
#' #specific diagnostics extraction
#' example.fit%diagnostics%list('R_0','inits') #initial values
#' example.fit%diagnostics%list('R_0','ess') #bulk and tail effective sample sizes
#' 
#' 
#' @export
`%diagnostics%`<-function(fit,select){
  is.char<-try(is.character(select),silent=TRUE)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  flag<-FALSE
  if(.is.evorates.element(fit)){
    type<-.get.element.type(fit)
    if(type!='chains'&type!='diagnostics'){
      flag<-TRUE
    }else{
      if(type=='chains'){
        fit<-list(chains=fit)
      }else{
        fit<-list(param.diags=fit)
      }
    }
  }else if(!inherits(fit,'evorates_fit')){
    flag<-TRUE
  }
  if(flag){
    stop("the %diagnostics% operator only accepts fitted evolving rates model fits (class 'evorates_fit') or loose chains elements on left hand side")
  }
  out<-.int.diagnostics(fit,select)
  .simplify.element(out)
}

#still think %post.probs% would be pointless, but I'll consider it in the future...
#7/27: %quantiles% for whole large_fit array is only longer than a simple apply implementation by
#~200 milliseconds...

#ADDING ability to select along first dimension would be nice, I think...done
#Allowing for name matching too? Would make things more complicated...

#' @export
`%select%`<-function(element,select){
  is.char<-try(is.character(select),silent=TRUE)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  element<-.coerce.to.3D(element)
  tmp<-.select.iterations(element,select)
  element<-tmp[[1]]
  select<-tmp[[2]]
  if(is.numeric(select)){
    select<-paste0('^',dimnames(element)[[2]][select],'$')
  }
  out<-.int.op(element,select)
  .simplify.element(out)
}