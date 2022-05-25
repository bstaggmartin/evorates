#' Extract branchwise rate parameters from a fitted evorates model
#' 
#' 
#' This is a convenient function for extracting branchwise rate parameters from an
#' \code{evorates_fit} object. It even works when the object has no branchwise
#' rate parameters (i.e., a Brownian Motion model).
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param select A numeric vector specifying edge indices in \code{fit$call$tree} for  which
#' to extract branchwise rate parameters. These can be negative to exclude edges as well. All
#' edges indices are extracted by default.
#' @param type A string specifying whether to extract posterior samples ("\code{chains}"),
#' quantiles ("\code{quantiles}"), means ("\code{means}"), or diagnostics ("\code{diagnostics}").
#' Can be any unambiguous abbreviation of these strings as well. Defaults to extracting
#' posterior samples.
#' @param extra.select A numeric, integer, or character vector specifying the specific samples/
#' quantiles/diagnostics to extract, depending on \code{type}. Defaults to \code{NULL},
#' which generally extracts all available samples/quantiles/diagnostics. See documentation on
#' param_block operators for details (\link{grapes-chains-grapes}, \link{grapes-quantiles-grapes},
#' \link{grapes-means-grapes}, or \link{grapes-diagnostics-grapes}). 
#' @param simplify \code{TRUE} or \code{FALSE}: should the resulting \code{param_block} array be 
#' simplified? If \code{TRUE} (the default), dimensions of length 1 in the result are
#' automatically collapsed, with corresponding information stored as attributes (this is the default
#' behavior of param_block operators).
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} set to whatever
#' \code{type} is. The dimension of these arrays will generally go in the order of
#' iterations/quantiles/diagnostics, then parameters, then chains. Any dimensions of length 1 are
#' collapsed and stored as attributes if \code{simplify} is \code{TRUE}. Note that branchwise
#' rate parameters go by the name \code{R_i}, where \code{i} is the index of the corresponding
#' edge in \code{fit$call$tree}.
#' 
#' 
#' @details The edge indices of \code{fit$call$tree} can be viewed by running:
#' \code{plot(fit$call$tree); edgelabels()}.
#' 
#' In the case of a \code{fit} with constrained rate variance (\code{R_sig2}) and trend
#' (\code{R_mu}) parameters, \code{fit} doesn't contain explicit branchwise rate parameters (\code{
#' R_i}, where \code{i} is the index of the corresponding edge). In this case, a simple Brownian
#' Motion model of trait evolution has been fitted, and the branchwise rate parameters are all
#' equivalent to the rate at the root (\code{R_0}). As such, this function just creates a \code{
#' param_block} array of \code{R_0} replicated the appropriate number of times, with parameters
#' renamed with the appropriate \code{R_i}'s.
#' 
#' 
#' @family parameter extraction functions
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #get posterior samples of branchwise rates
#' Rs <- get.R(cet_fit)
#' #let's say we want to get some of the edges towards the root
#' plot(cet_fit$call$tree); edgelabels()
#' #select particular edges
#' Rs <- get.R(cet_fit, select = c(1, 2, 9, 28, 29, 30))
#' #maybe specific samples too?
#' Rs <- get.R(cet_fit, select = c(1, 2, 9, 28, 29, 30),
#'             extra.select=c(1, 23, 47))
#' 
#' #this may be a more common way to get edge indices; say we want to look at the genus Mesoplodon
#' edges <- get.clade.edges(cet_fit$call$tree, "Mesoplodon")
#' Meso.Rs <- get.R(cet_fit, select = edges)
#' #or at everything EXCEPT Mesoplodon
#' notMeso.Rs <- get.R(cet_fit, select = -edges)
#' 
#' #could also look at quantiles, means, or diagnostics
#' med.Rs <- get.R(cet_fit, select = edges,
#'                 type = "quantiles",
#'                 extra.select = 0.5)
#' mean.Rs <- get.R(cet_fit, select = edges,
#'                  type = "means")
#' init.Rs <- get.R(cet_fit, select = edges,
#'                  type = "diagnostics",
#'                  extra.select = c("inits", "ess"))
#'                  
#' #here's an example of what happens when you don't simplify the result
#' med.Rs <- get.R(cet_fit, select = edges,
#'                 type = "quantiles",
#'                 extra.select = 0.5,
#'                 simplify = FALSE)
#'                  
#' #note that all this is pretty much equivalent to running <fit> %<type>% list(<select>, <extra.select>)
#' #at least when <fit> includes branchwise rate parameters
#' 
#' 
#' @export
get.R<-function(fit,
                select=seq_len(Nedge(fit)),
                type=c('chains','quantiles','means','diagnostics'),
                extra.select=NULL,
                simplify=TRUE){
  select<-.check.edge.indices(select,
                              deparse(substitute(select)),
                              Nedges=Nedge(fit))
  type<-.match.type(type)
  no.Rs<-fit$call$constrain.Rsig2&!fit$call$trend
  if(no.Rs){
    in.select<-select
    select<-'^R_0$'
  }
  out<-.call.op(type,fit,list(select,extra.select),FALSE)
  if(no.Rs){
    out<-.call.select(out,rep(1,length(in.select)))
    names(out)<-paste0('R_',in.select)
  }
  if(simplify){
    out<-.simplify.par(out)
  }
  out
}