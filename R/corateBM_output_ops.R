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

#still think %post.probs% would be pointless, but I'll consider it in the future...
#7/27: %quantiles% for whole large_fit array is only longer than a simple apply implementation by
#~200 milliseconds...