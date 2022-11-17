#' Generate a parameter block of standard normal samples
#' 
#' 
#' This function initializes an array of class "\code{param_block}" containing samples from
#' standard normal distribuion (i.e., mean 0 and standard deviation 1). There are most helpful for
#' calculating distributions of expected rates over some time interval.
#' 
#' 
#' @param niter The number of iterations/rows of the resulting \code{param_block} array.
#' @param n The number of parameters/columns of the resulting \code{param_block} array.
#' @param nchains The number of chains/slices of the reuslting \code{param_block} array.
#' @param simplify \code{TRUE} or \code{FALSE}: should the resulting \code{param_block} array be 
#' simplified? If \code{TRUE} (the default), dimensions of length 1 in the result are
#' automatically collapsed, with corresponding information stored as attributes (this is the default
#' behavior of param_block operators).
#' @param fit An object of class "\code{evorates_fit}" or "\code{param_block}". If not \code{NULL}
#' (the default), the resulting \code{param_block} array will copy the names and dimensions of
#' \code{fit} except for parameters/columns.
#' @param chain.names A character vector of names for the chains/slices of the resulting
#' \code{param_block} array.
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} of "\code{chains}".
#' The dimension of these arrays will generally go in the order of
#' iterations, then parameters, then chains. Any dimensions of length 1 are
#' collapsed and stored as attributes if \code{simplify} is \code{TRUE}. Parameters
#' are all named "\code{N(0;1)i}", where \code{i} denotes the \code{i}th column of the resulting array.
#' 
#' 
#' @details The resulting normal samples may be transformed to have different means/standard deviations
#' by multiplying with and adding constants--see examples below.
#' 
#' 
#' @seealso \link{param_block-class} for general information on \code{param_block} arrays and 
#' \code{\link[=grapes-chains-grapes]{\%chains\%}()},
#' \code{\link[=grapes-quantiles-grapes]{\%quantiles\%}()},
#' \code{\link[=grapes-means-grapes]{\%means\%}()},
#' \code{\link[=grapes-diagnostics-grapes]{\%diagnostics\%}()},
#' and \code{\link[=grapes-select-grapes]{\%select\%}()} for more information on
#' \code{param_block} operators.
#' 
#' 
#' @examples
#' rnorm.par(1500, 1, 4)
#' #compare without simplification
#' rnorm.par(1500, 1, 4, simplify = FALSE)
#' 
#' #using a fitted evorates model to automatically determine dimensions/names
#' data("cet_fit")
#' rnorm.par(n = 2, fit = cet_fit)
#' #or a param_block
#' parblock <- rnorm.par(500, 3, 1)
#' rnorm.par(n = 5, fit = parblock)
#' 
#' #transforming to non-standard normal distribution
#' means <- c(-1, 5, 2)
#' sds <- c(2, 0.5, 1)
#' prof.plot(parblock * sds + means, alpha = 0.5, smooth = TRUE)
#' 
#' #distribution of fold-changes expected after 1 million years
#' rates.per.my<-exp(setNames(sqrt(cet_fit %chains% "R_sig2") * rnorm.par(fit = cet_fit) + cet_fit %chains% "R_mu",
#'                            "rates_per_my"))
#' prof.plot(rates.per.my)
#' #after 10 million years
#' rates.per.10my<-exp(setNames(sqrt(cet_fit %chains% "R_sig2") * sqrt(10) * rnorm.par(fit = cet_fit) + 10 * cet_fit %chains% "R_mu",
#'                              "rates_per_10my"))
#' prof.plot(rates.per.10my)
#' #easier to interpret on log scale maybe
#' prof.plot(log(rates.per.10my))            
#' 
#' 
#' @export
rnorm.par<-function(niter=1,n=1,nchains=1,simplify=TRUE,fit=NULL,chain.names=paste('chain',seq_len(nchains))){
  if(!is.null(fit)){
    if(.is.par(fit)){
      fit<-.make.par.3D(fit)
      dims<-dim(fit)
      niter<-dims[1]
      nchains<-dims[length(dims)]
      tmp<-dimnames(fit)
      chain.names<-tmp[[length(tmp)]]
    }else if(inherits(fit,'evorates_fit')){
      niter<-fit$sampler.control$iter-fit$sampler.control$warmup
      nchains<-fit$sampler.control$chains
      tmp<-attr(fit$chains,'chains')
      if(is.null(tmp)){
        tmp<-dimnames(fit$chains)
        chain.names<-tmp[[length(tmp)]]
      }else{
        chain.names<-tmp
      }
    }
  }
  out<-array(rnorm(niter*n*nchains),
             c(niter,n,nchains),
             c(iterations=list(NULL),
               parameters=list(paste('N(0;1)',if(n==1) '' else 1:n,sep='')),
               chains=list(chain.names)))
  attr(out,'param_type')<-'chains'
  out<-.add.par.class(out)
  if(simplify){
    out<-.simplify.par(out)
  }
  out
}