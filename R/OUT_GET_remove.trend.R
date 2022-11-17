#' Extract detrended branchwise rate parameters from a fitted evorates model
#' 
#' 
#' This function extracts "detrended" branchwise rate parameters from an
#' \code{evorates_fit} object. It even works when the object has no branchwise
#' rate parameters (i.e., a Brownian Motion model) or has no trend parameter.
#' 
#' 
#' @param fit An object of class "\code{evorates_fit}".
#' @param select A numeric vector specifying edge indices in \code{fit$call$tree} for  which
#' to extract dtrended branchwise rate parameters. These can be negative to exclude edges as well.
#' All edges indices are extracted by default.
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
#' collapsed and stored as attributes if \code{simplify} is \code{TRUE}. Note that detrended
#' branchwise rate parameters go by the name \code{uncent_Rdev_i} (i.e., an "uncentered rate
#' deviation"), where \code{i} is the index of the corresponding edge in \code{fit$call$tree}.
#' 
#' 
#' @details The edge indices of \code{fit$call$tree} can be viewed by running:
#' \code{plot(fit$call$tree); edgelabels()}.
#' 
#' In the case of a \code{fit} with constrained rate variance (\code{R_sig2}) parameters, the
#' detrended branchwise rate parameters are all equivalent to the rate at the root (\code{R_0}). 
#' As such, this function just creates a \code{param_block} array of \code{R_0} replicated the
#' appropriate number of times, with parameters renamed with the appropriate \code{uncent_Rdev_i}'s 
#' (i.e., an "uncentered rate deviations"), where \code{i} is the index of the corresponding edge
#' in \code{fit$call$tree}.
#' 
#' 
#' @seealso more information on detrending is available in the documentation for \link{fit.evorates},
#' \link{input.evorates}, \link{run.evorates}, and \link{output.evorates}, particularly in the section
#' on parameter definitions
#' 
#' 
#' @family parameter extraction functions
#' 
#' 
#' @examples
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #get posterior samples of detrended branchwise rates
#' uncent.Rdevs <- remove.trend(cet_fit)
#' #let's say we want to get detrended rates for some of the edges towards the root
#' plot(cet_fit$call$tree); edgelabels()
#' #select particular edges
#' uncent.Rdevs <- remove.trend(cet_fit, select = c(1, 2, 9, 28, 29, 30))
#' #maybe specific samples too?
#' uncent.Rdevs <- remove.trend(cet_fit, select = c(1, 2, 9, 28, 29, 30),
#'                              extra.select=c(1, 23, 47))
#' 
#' #this may be a more common way to get edge indices; say we want to look at detrended rates in the genus Mesoplodon
#' edges <- get.clade.edges(cet_fit$call$tree, "Mesoplodon")
#' Meso.uncent.Rdevs <- remove.trend(cet_fit, select = edges)
#' #or at everything EXCEPT Mesoplodon
#' notMeso.uncent.Rdevs <- remove.trend(cet_fit, select = -edges)
#' 
#' #could also look at quantiles, means, or diagnostics
#' med.uncent.Rdevs <- remove.trend(cet_fit, select = edges,
#'                                  type = "quantiles",
#'                                  extra.select = 0.5)
#' mean.uncent.Rdevs <- remove.trend(cet_fit, select = edges,
#'                                   type = "means")
#' init.uncent.Rdevs <- remove.trend(cet_fit, select = edges,
#'                                   type = "diagnostics",
#'                                   extra.select = c("inits", "ess"))
#'                  
#' #here's an example of what happens when you don't simplify the result
#' med.uncent.Rdevs <- remove.trend(cet_fit, select = edges,
#'                                  type = "quantiles",
#'                                  extra.select = 0.5,
#'                                  simplify = FALSE)
#'                  
#' #note that all this equivalent to detrending <fit> %<type>% list(<select>, <extra.select>)
#' #at least when <fit> includes branchwise rate parameters
#' 
#' 
#' @export
remove.trend<-function(fit,
                       select=seq_len(Nedge(fit)),
                       type=c('chains','quantiles','means','diagnostics'),
                       extra.select=NULL,
                       simplify=TRUE){
  select<-.check.edge.indices(select,
                              deparse(substitute(select)),
                              Nedge(fit))
  type<-.match.type(type)
  if(fit$call$trend){
    if(!fit$call$constrain.Rsig2){
      Rmu<-fit$chains%select%'R_mu'
      er<-edge.ranges(fit)[select,,drop=FALSE]
      el<-fit$call$tree$edge.length[select]
      out<-get.R(fit,type='chains',select=select,simplify=FALSE)-
        (-log(abs(Rmu))-log(el)+log(abs(exp(Rmu*er[,2])-exp(Rmu*er[,1]))))
      if(type=='diagnostics'){
        inits.Rmu<-fit$diagnostics%select%list('R_mu','inits')
        inits.out<-get.R(fit,type='diagnostics',select=select,extra.select='inits',simplify=FALSE)-
          (-log(abs(inits.Rmu))-log(el)+log(abs(exp(inits.Rmu*er[,2])-exp(inits.Rmu*er[,1]))))
      }
    }else{
      fit$call$trend<-FALSE #forces get.R to extract R_0
    }
  }
  if(!fit$call$trend){
    out<-get.R(fit,
               select=select,
               type=type,
               extra.select=extra.select,
               simplify=FALSE)
  }else{
    out<-list(chains=out,sampler.control=1) #"cheating"--just so sampler.params isn't NULL
    out<-.call.op(type,out,list('.',extra.select),FALSE)
    if(type=='diagnostics'){
      inds<-dimnames(out)[[1]]=='inits'
      out[inds,,]<-rep(inits.out,each=sum(inds))
    }
  }
  names(out)<-paste0('uncent_Rdev_',select)
  if(simplify){
    out<-.simplify.par(out)
  }
  out
}