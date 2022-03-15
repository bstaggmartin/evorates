####MAIN FUNCTIONS####

#' Extract posterior distribution samples from a fitted evorates model
#'
#'
#' This operator extracts samples from the posterior distributions for particular parameters from the
#' output of \code{fit.evorates} or \code{output.evorates}
#'
#'
#' @param fit An object of class "\code{evorates_fit}" or "\code{param_block}"
#' @param select A list with two elements: 1) A character or numeric vector. If of class character, the given
#' text is matched to parameter
#' names using regular expressions. If of class numeric, the numbers are matched to indices of rate parameters
#' ("\code{R_i}"), including the root ("\code{R_0}"). 2) A vector of integers specifying particular iterations
#' to extract from the chain. If unsupplied, all iterations will be extracted.
#' 
#' 
#' @return A numeric vector, matrix, or 3D array of class "\code{param_block}". The dimensions will always go
#' in the order of iterations, then parameters, then chains, collapsing any dimensions of length 1 and storing
#' associated name information as attributes.
#' 
#' 
#' @details 
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
  .proc.op('chains',fit,select,deparsed.select=deparse(substitute(select)))
}

#' Extract posterior distribution quantiles from a fitted evorates model
#' 
#' 
#' @export
`%quantiles%`<-function(fit,select){
  .proc.op('quantiles',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','quantiles'))
}

#' Extract posterior distribution means from a fitted evorates model
#' 
#' 
#' @export
`%means%`<-function(fit,select){
  .proc.op('means',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','means'))
}

#' Extract posterior distribution diagnostics from a fitted evorates model
#' 
#' 
#' @export
`%diagnostics%`<-function(fit,select){
  .proc.op('diagnostics',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','diagnostics'))
}

####SIMPLIFIED FUNCTION NAMES####

#' @export
`%c%`<-function(fit,select){
  .proc.op('chains',fit,select,deparsed.select=deparse(substitute(select)))
}

#' @export
`%q%`<-function(fit,select){
  .proc.op('quantiles',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','quantiles'))
}

#' @export
`%m%`<-function(fit,select){
  .proc.op('means',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','means'))
}

#' @export
`%d%`<-function(fit,select){
  .proc.op('diagnostics',fit,select,deparsed.select=deparse(substitute(select)),choices=c('chains','diagnostics'))
}

####SELECT####

#' Extract subsets from a parameter block
#' 
#' 
#' @export
`%select%`<-function(x,select){
  .proc.select(x,select,deparsed.select=deparse(substitute(select)))
}

#' @export
`%s%`<-function(x,select){
  .proc.select(x,select,deparsed.select=deparse(substitute(select)))
}