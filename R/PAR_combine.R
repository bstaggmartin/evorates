.combine.par<-function(x){
  x<-x[lengths(x)>0]
  x<-lapply(x,.make.par.3D)
  out.type<-.get.par.type(x[[1]])
  out.dimnames<-lapply(x,dimnames)
  out.dimnames[[1]][[2]]<-unlist(lapply(out.dimnames,'[[',2))
  out.dims<-dim(x[[1]])
  out.dims[2]<-length(out.dimnames[[1]][[2]])
  out<-array(NA,out.dims,out.dimnames[[1]])
  foo<-function(x,chain){
    x[,,chain]
  }
  for(i in 1:out.dims[3]){
    out[,,i]<-unlist(lapply(x,foo,chain=i))
  }
  attr(out,'param_type')<-out.type
  .add.par.class(out)
}

#' @export
c.param_block<-function(...,fit=NULL){
  ls<-list(...)
  types<-lapply(ls,.get.par.type)
  pot.probs<-lengths(types)==0
  types[lengths(types)==0]<-'select'
  if(is.null(fit)){
    if(all(pot.probs)){
      stop('all provided param_blocks appear to be selection(s) to extract from an evorates fit object, yet no fit was supplied')
    }else if (any(pot.probs)){
      warning(paste(ls[pot.probs],collapse=', '),'appear to be selection(s) to extract from an evorates fit object, yet no fit was supplied: ignored these selections')
      ls<-ls[!pot.probs]
      types<-types[!post.probs]
    }
  }
  types<-unlist(types)
  type.nms<-c('quantiles','means','diagnostics')
  type.which<-lapply(type.nms,'==',l=types)
  type.has<-unlist(lapply(type.which,any))
  if(sum(type.has)>1){
    stop('c() can only be used to combined chains param_blocks with those of a single other type')
  }
  out.type<-type.nms[type.has]
  type.inds<-unlist(type.which[type.has])
  if(is.null(type.inds)){
    type.inds<-rep(FALSE,length(ls))
  }
  type.par<-ls[type.inds]
  #coerce param_blocks of quantiles, means, or diagnostics types to be compatible
  #(return an error if they're not)
  if(length(type.par)==1){
    type.par<-list(.expand.par(type.par[[1]]))
  }else if(length(type.par)>1){
    check<-do.call(.compatible.dims.check,type.par)
    if(!check[[4]]){
      if(!check[[5]][2]){
        stop('only param_blocks with the same chains can be combined')
      }
      #should never get to below portion with means...
      nms<-lapply(check[[3]],'[[',1)
      uni.nms<-unique(unlist(nms))
      foo<-function(x){
        all(unlist(lapply(nms,function(ii) x%in%ii)))
      }
      avail.nms<-uni.nms[unlist(lapply(uni.nms,foo))]
      if(length(avail.nms)==0){
        stop('param_blocks of ',out.type,' type contain no matching quantities')
      }
      if(out.type=='diagnostics'){
        inds<-lapply(nms,match,x=avail.nms)
      }else{
        inds<-rep(list(avail.nms),length(check[[1]]))
      }
      check[[1]]<-lapply(check[[1]],function(ii) setNames(list(ii),out.type))
      check[[1]]<-lapply(seq_along(check[[1]]),function(ii) .call.op(out.type,check[[1]][[ii]],list('.',inds[[ii]]),FALSE))
    }
    type.par<-check[[1]]
  }else{
    out.type<-'chains'
  }
  if(out.type=='chains'){
    nms<-NULL
  }else{
    nms<-dimnames(type.par[[1]])[[1]]
  }
  #deal with chain stuff
  chain.inds<-types=='chains'
  chain.par<-ls[chain.inds]
  if(length(chain.par)==1){
    chain.par<-list(.expand.par(chain.par[[1]]))
  }else if(length(chain.par)>1){
    check<-do.call(.compatible.dims.check,chain.par)
    if(!check[[4]]){
      if(!check[[5]][2]){
        stop('only param_blocks with the same chains can be combined')
      }
      #check for sampler parameters
      flag<-FALSE
      lens<-unlist(lapply(check[[2]],'[',1))
      if(length(unique(lens))>2){
        flag<-TRUE
      }else{
        pot.sampler<-which(lens==max(lens))
        param.names<-lapply(check[[3]][pot.sampler],'[[',2)
        sampler.names<-paste0(c('accept_stat','stepsize','treedepth','n_leapfrog','divergent','energy','prior','lik','post'),'__')
        is.sampler.names<-unlist(lapply(param.names,function(ii) 
          all(ii%in%sampler.names)))
        if(all(is.sampler.names)){
          target.len<-min(lens)
          diags.len<-max(lens)
          check[[1]][pot.sampler]<-lapply(
            check[[1]][pot.sampler],
            function(ii) ii[-(1:(diags.len-target.len)),,,drop=FALSE]
          )
        }else{
          flag<-TRUE
        }
      }
      if(flag){
        stop('chain param_blocks appear to contain mismatching iterations')
      }
    }
    chain.par<-check[[1]]
  }
  if(out.type!='chains'){
    chain.par<-lapply(chain.par,function(ii) setNames(list(ii),'chains'))
    chain.par<-lapply(chain.par,.call.op,type=out.type,select=list('.',nms),check.sampler=FALSE)
  }
  if(!is.null(fit)){
    which.select<-which(types=='select')
    for(i in which.select){
      ls[[i]]<-.call.op(out.type,fit,list(ls[[i]],nms),TRUE)
    }
  }
  ls[type.inds]<-type.par
  ls[chain.inds]<-chain.par
  .simplify.par(.combine.par(ls))
}

#' Combine parameter blocks
#' 
#' 
#' This function is a wrapper for the \code{c()} method for \code{param_block} arrays that allows you
#' to also specify a fitted evorates model from which to extract parameters specified via \code{...}.
#' 
#' 
#' @param ... \code{param_block} arrays and/or character/numeric parameter selections to extract from \code{fit}.
#' Additional arguments controlling how parameters are extracted from \code{fit} are determined
#' automatically based on any \code{param_block} arrays in \code{...}. If \code{...} only consist
#' of parameter selections, the function defaults to extracting posterior samples from \code{fit}. All 
#' \code{param_block} arrays in \code{...} must have compatible rows and slices (i.e., same number
#' of chains, same names).
#' \code{param_block} arrays with a \code{param_type} of "\code{chains}" have rows compatible with any other
#' \code{param_block} array (they can be coerced to another \code{param_type} on the fly), with the exception
#' of other \code{chains} \code{param_block} arrays that
#' have a different number of rows. Otherwise, \code{param_block} arrays generally have compatible rows with other arrays
#' of the same \code{param_type} provided they share least one row (i.e., quantile or diagnostic) in common.
#' @param fit An object of class "\code{evorates_fit}" from which to extract parameters given by the 
#' character/numeric elements of \code{...}. Ignored if \code{NULL} (the default).
#' 
#' 
#' @return An array of class "\code{param_block}" with a \code{param_type} determined by the elements
#' of \code{...} (tends to default "\code{chains}" when in doubt). The dimensions of these arrays will generally
#' go in the order of iterations/quantiles/diagnostics, then parameters, then chains.  Any dimensions of length 1 are
#' collapsed and stored as attributes. The resulting array contains all parameters specified by \code{...} as columns,
#' ordered as they are in \code{...}.
#' 
#' 
#' @details I had to create a separate function for this because R does not default to the \code{c()} method for
#' \code{param_block} arrays if character or numeric arguments are passed to \code{c()}. This function
#' simply ensures R calls the proper method no matter what's passed to \code{...}.
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
#' #get whale/dolphin evorates fit
#' data("cet_fit")
#' 
#' #how are average rates for some focal clades affected by R_sig2/R_mu estimates?
#' parblock <- get.bg.rate(fit = cet_fit,
#'                         node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                                c('Mesoplodon','Orca','Globicephalinae')),
#'                         )
#' parblock <- par.c("R_mu", "R_sig2", parblock, fit = cet_fit)
#' plot(parblock %chains% "Mesoplodon" ~ parblock %chains% "R_sig2"
#' plot(parblock %chains% "Orca" ~ parblock %chains% "R_mu")
#' 
#' #automatic conversion based on param_block type
#' parblock <- get.bg.rate(fit = cet_fit,
#'                         node.groups = setNames(list('Mesoplodon','Orcinus',c('Pseudorca','Feresa')),
#'                                                c('Mesoplodon','Orca','Globicephalinae')),
#'                         type = "quantiles")
#' par.c("R_mu", "R_sig2", parblock, fit = cet_fit)
#' 
#' #can mix numeric and character selections
#' par.c("R_0", "R_mu", 1, 28, fit = cet_fit)
#' 
#' 
#' @export
par.c<-function(...,fit=NULL){
  c.param_block(...,fit=fit)
}
