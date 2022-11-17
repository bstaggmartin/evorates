#add auto-conversion capabilities similar to combining functions?
#inherits param_type of l if both l and r are param_blocks
#' @export
Ops.param_block<-function(l,r=NULL){
  make.Ops.param_blocks.pw<-mget('make.Ops.param_blocks.pw',
                                 envir=parent.frame(),
                                 ifnotfound=FALSE)[[1]]
  if(is.null(r)){
    type<-.get.par.type(l)
    nms<-names(l)
    if(!is.list(nms)){
      nms<-list(nms)
    }
    names(l)<-lapply(nms,
                     function(ii) paste0('(',
                                         '%',.Generic,'%',
                                         ii,
                                         ')'))
    l<-.strip.par.class(l)
    out<-do.call(.Generic,list(l))
    out<-.add.par.class(out)
  }else{
    is.l<-.Method[1L]=='Ops.param_block'
    is.r<-.Method[2L]=='Ops.param_block'
    if(is.l&is.r){
      type<-.get.par.type(l)
      ll<-deparse(substitute(l))
      rr<-deparse(substitute(r))
      out<-.int.par.math(l,r,.Generic,ll,rr,
                             make.Ops.param_blocks.pw)
    }else{
      out<-NULL
      if(is.l){
        type<-.get.par.type(l)
        if(!is.null(dim(r))){
          stop('math between param_blocks and plain matrices/arrays not (yet) allowed')
        }
        vec<-length(r)!=1
        if(vec){
          l<-.make.par.3D(l)
        }
        nms<-names(l)
        if(!is.list(nms)){
          nms<-list(nms)
        }
        if(vec){
          lens<-c(length(nms[[1]]),length(r))
          if(lens[1]!=lens[2]){
            if(any(!lens)) max.len<-0 else max.len<-max(lens)
            l<-l[,rep(1:lens[1],length.out=max.len),,drop=FALSE]
            l<-.add.par.class(l)
            r<-rep(r,length.out=max.len)
          }
        }
        names(l)<-lapply(nms,
                         function(ii) paste0('(',ii,
                                             '%',.Generic,'%',
                                             r,')'))
        l<-.strip.par.class(l)
        if(vec){
          out<-sweep(l,2,r,function(x,y) do.call(.Generic,list(x,y)))
        }
      }else{
        type<-.get.par.type(r)
        if(!is.null(dim(l))){
          stop('math between plain matrices/arrays and param_blocks not (yet) allowed')
        }
        vec<-length(l)!=1
        if(vec){
          r<-.make.par.3D(r)
        }
        nms<-names(r)
        if(!is.list(nms)){
          nms<-list(nms)
        }
        if(vec){
          lens<-c(length(l),length(nms[[1]]))
          if(lens[1]!=lens[2]){
            if(any(!lens)) max.len<-0 else max.len<-max(lens)
            l<-rep(l,length.out=max.len)
            r<-r[,rep(1:lens[2],length.out=max.len),,drop=FALSE]
            r<-.add.par.class(r)
          }
        }
        names(r)<-lapply(nms,
                         function(ii) paste0('(',l,
                                             '%',.Generic,'%',
                                             ii,')'))
        r<-.strip.par.class(r)
        if(vec){
          out<-sweep(r,2,l,function(x,y) do.call(.Generic,list(y,x)))
        }
      }
      if(is.null(out)){
        out<-do.call(.Generic,list(l,r))
      }
      if(is.l) tmp<-l else tmp<-r
      for(i in c('quantiles','diagnostics','parameters','chains')){
        if(!is.null(attr(tmp,i))){
          attr(out,i)<-attr(tmp,i)
        }
      }
      out<-.add.par.class(out)
    }
  }
  attr(out,'param_type')<-type
  out
}

#' @export
Math.param_block<-function(x,...){
  tmp<-do.call(paste,c(list(...),collapse=list(';')))
  if(nzchar(tmp)){
    tmp<-paste0(';',tmp)
  }
  nms<-names(x)
  if(!is.list(nms)){
    nms<-list(nms)
  }
  foo<-function(x){
    parens<-grep('^\\(',x)
    x[parens]<-gsub('\\)$','',gsub('^\\(','',x[parens]))
    paste0(.Generic,
           '(',
           x,
           tmp,
           ')')
  }
  names(x)<-lapply(nms,foo)
  x<-.strip.par.class(x)
  out<-do.call(.Generic,c(list(x),list(...)))
  out<-.add.par.class(out)
  out
}

#will require combine function, but idea is simple--just apply .Generic across parameters to get sums of parameter combos, etc.
#' @export
Summary.param_block<-function(...,na.rm=TRUE){
  out<-.expand.par(do.call(c,list(...)))
  type<-.get.par.type(out)
  outnames<-dimnames(out)
  param.names<-outnames[[2]]
  nms<-paste(param.names,collapse=';')
  parens<-grep('^\\(',nms)
  nms[parens]<-gsub('\\)$','',gsub('^\\(','',nms[parens]))
  if(eval(.Generic)=='range'){
    out<-.add.par.class(aperm(apply(out,c(1,3),.Generic,na.rm=na.rm),
                              c(2,1,3)))
    outnames[[2]]<-paste0(c('min','max'),
                          '(',
                          nms,
                          ')')
    dimnames(out)<-outnames
    attr(out,'param_type')<-type
  }else{
    out<-.add.par.class(apply(out,c(1,3),.Generic,na.rm=na.rm))
    attr(out,'parameters')<-paste0(.Generic,
                                   '(',
                                   nms,
                                   ')')
    attr(out,'param_type')<-type
    out<-.expand.par(out)
  }
  out
}

.int.par.math<-function(l,r,.Generic,in.l,in.r,make.Ops.param_blocks.pw){
  lnames<-names(l)
  if(is.list(lnames)){
    lnames<-lnames[[which.max(lengths(lnames))]]
  }
  rnames<-names(r)
  if(is.list(rnames)){
    rnames<-rnames[[which.max(lengths(rnames))]]
  }
  lsing<-length(lnames)==1
  rsing<-length(rnames)==1
  check<-.compatible.dims.check(l,r)
  if(!check[[4]]){
    stop('Dimensional mismatch between ',
         in.l,
         ' and ',
         in.r,
         ': did these param_blocks come from from the same evorates_fit, and were they extracted in the same way?')
  }
  if(lsing|rsing){
    lr<-check[[1]]
    lrdims<-check[[2]]
    #in the case that both are singular, you don't necessarily need to coerce anything to 3D...
    #below is most of the way there, but you have to account for cases where number of parameter dimensions mismatches...
    # if(lsing&rsing){
    #   lr<-lapply(lr,.strip.par.class)
    #   out<-do.call(.Generic,lr)
    #   out<-.add.par.class(out)
    #   nms<-check[[3]]
    #   names(out)<-lapply(1,
    #                     )
    # }
    if(lsing) sing<-0 else sing<-1
    notsing<-(!sing)+1
    sing<-sing+1
    lr[[sing]]<-aperm(.make.par.3D(lr[[sing]]),c(1,3,2))
    dims<-length(lrdims[[notsing]])
    lr[[notsing]]<-aperm(lr[[notsing]],
                         (1:dims)[c(1,dims,2:(dims-1))])
    nms<-check[[3]][[notsing]]
    nms<-nms[-c(1,length(nms))]
    out<-lr[[notsing]]
    out[]<-do.call(.Generic,list(as.vector(lr[[1]]),as.vector(lr[[2]])))
    if(dims>3){
      out<-aperm(out,c(seq.int(3,dims),1,2))
    }else{
      out<-aperm(out,c(1,3,2))
    }
    out<-.add.par.class(out)
    if(lsing){
      names(out)<-lapply(nms,
                         function(ii) paste0('(',names(.add.par.class(l)),
                                             '%',.Generic,'%',
                                             ii,')'))
    }else{
      names(out)<-lapply(nms,
                         function(ii) paste0('(',ii,
                                             '%',.Generic,'%',
                                             names(.add.par.class(r)),')'))
    }
  }else{
    l<-.make.par.3D(check[[1]][[1]])
    r<-.make.par.3D(check[[1]][[2]])
    ldims<-dim(l)
    rdims<-dim(r)
    lnames<-dimnames(l)
    rnames<-dimnames(r)
    if(make.Ops.param_blocks.pw){
      out.param.names<-paste0('(',rep(lnames[[2]],rdims[2]),
                              '%',.Generic,'%',
                              rep(rnames[[2]],each=ldims[2]),')')
      out<-.add.par.class(array(NA,
                                c(ldims[1],length(out.param.names),ldims[3]),
                                c(lnames[1],parameters=list(out.param.names),lnames[3])))
      counter<-1
      for(i in 1:rdims[2]){
        for(j in 1:ldims[2]){
          out[,counter,]<-do.call(.Generic,
                                  c(list(l[,j,]),
                                    list(r[,i,])))
          counter<-counter+1
        }
      }
    }else{
      lens<-c(ldims[2],rdims[2])
      if(lens[1]!=lens[2]){
        max.len<-max(lens)
        l<-l[,rep(1:lens[1],length.out=max.len),,drop=FALSE]
        lnames[[2]]<-rep(lnames[[2]],length.out=max.len)
        r<-r[,rep(1:lens[2],length.out=max.len),,drop=FALSE]
        rnames[[2]]<-rep(rnames[[2]],length.out=max.len)
      }
      out<-.add.par.class(do.call(.Generic,list(.strip.par.class(l),.strip.par.class(r))))
      names(out)<-paste0('(',lnames[[2]],
                         '%',.Generic,'%',
                         rnames[[2]],')')
    }
  }
  out
}

#' Do pairwise comparison among parameter blocks
#' 
#' Wrap this function around a mathematical expression involving \code{param_block} arrays
#' to switch to "pairwise comparison mode". This changes the default behavior from recycling
#' columns of \code{param_block} arrays to ensure each array has the same number of parameters to
#' instead performing mathematical operations between \emph{each pair} of parameters in the arrays.
#' 
#' @param expr A mathematical expression generally involving two or more \code{param_block} arrays.
#' All arrays must have compatible iterations/quantiles/diagnostics (typically the first
#' dimension/rows) and chains (typically the third dimension/slices).
#' @param fit An optional object of class "\code{evorates_fit}" or "\code{param_block}". If
#' not \code{NULL} (the default), any unrecognized \emph{object names} in \code{expr} will be
#' interpreted as parameters to extract from \code{fit} \emph{using the} \code{\%chains\%}
#' \emph{operator only} (the delayed evaluation makes it difficult to automatically
#' select different operators based on \code{expr}). Note that this function will only pass
#' referenced \emph{object names}, not strings, to the \code{\%chains\%} operator. As a result, only
#' unquoted names can be used to specify selections from \code{fit}; just as with \code{data.frame}
#' objects, backticks (\code{``}) can be used to enclose object names including non-standard
#' characters (see examples below). For more details on how parameters are extracted, see
#' \link{grapes-chains-grapes}.
#' 
#' 
#' @return Typically another array of class "\code{param_block}" with a \code{param_type} set
#' to that of the inputted arrays. The dimension of these arrays will generally go in the order of
#' iterations/quantiles/diagnostics, then parameters, then chains. The array will have the same
#' iterations/quantiles/diagnostics and chains as the inputted arrays but will have a new set of
#' parameters corresponding to every pairwise combination of parameters among the inputted arrays.
#' As is typical with doing math on \code{param_block} arrays, the result is never simplified (i.e.,
#' dimension of length 1 are never collapsed),
#' such that complicated expressions are evaluated more efficiently.
#' 
#' Resulting parameters are arranged in "column-major" order, meaning that each parameter on the right hand
#' side of any binary mathematical operation (e.g., \code{+}, \code{>}, etc.) will be recycled to
#' match up with the total number of parameters on the left hand side. For example,
#' \code{cet_fit \%chains\% 3:4 - cet_fit \%chains\% 1:2} would result in parameters
#' \code{R_3\%-\%R_1}, \code{R_4\%-\%R_1}, \code{R_3\%-\%R_2}, and \code{R_4\%-\%R_2}, in that order.
#' 
#' Note that, since \code{expr} could technically be anything, you could get a completely different
#' output from this function if \code{expr} doesn't include any \code{param_block} objects (e.g.,
#' \code{pwc(1)} would just return 1).
#' 
#' 
#' @details This function can be a bit finicky at times (particularly when taking advantage of the
#' \code{fit} argument), but I find it too helpful to get rid of at this point. Basically, any time
#' you do math with \code{param_block} arrays, R checks whether a variable in the parent environment,
#' \code{make.Ops.param_blocks.pw}, is \code{TRUE} or \code{FALSE}. If \code{TRUE}, R does math with
#' \code{param_block} arrays in pairwise mode. This function is thus nothing more than a wrapper that
#' evaluates \code{expr} with \code{make.Ops.param_blocks.pw} set to \code{TRUE} in the parent 
#' environment. Any time you start messing with scoping and going across environments in R, weird
#' side effects are possible.
#' 
#' Practically speaking, just be sure to look out for ambiguous situations where you try to use an object
#' name that already exists in your current R environment to refer to something to extract from \code{fit}!
#' When in doubt about this, it is safest to keep \code{fit} as \code{NULL} and manually extract the
#' necessary \code{param_block} arrays within \code{expr} or as separate objects later referred to
#' by \code{expr}.
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
#' #here we compare branchwise rates for edges 3 and 4 to those for 1 and 2
#' pwc.result <- pwc(fit %chains% 3:4 - fit%chains% 1:2)
#' #compare to
#' norm.result <- fit %chains% 3:4 - fit%chains% 1:2
#' 
#' #notice that we cannot use numeric selections in conjunction with the fit argument
#' #such selections are just interpreted as numbers
#' pwc.result <- pwc(3:4 - 1:2, fit = cet_fit) #returns c(2, 2)
#' #note that character selections are similarly just interpreted as strings
#' \dontrun{pwc.result <- pwc(c("^R_3$", "^R_4$") - c("^R_1$", "^R_2$"), fit = cet_fit) #returns error}
#' #you have to use backticks instead to make the selections object names rather than strings
#' pwc.result <- pwc(c(`^R_3$`, `^R_4$`) - c(`^R_1$`, `^R_2$`), fit = cet_fit)
#' #object names work without backticks if they don't include special symbols like ^ or $
#' pwc.result <- pwc(c(R_sig2, R_mu) - c(`^R_1$`, `^R_2$`), fit = cet_fit)
#' #also, note that specifying fit always results in selecting posterior samples/"chains", so this won't work
#' \dontrun{pwc.result <- pwc(c(`^R_3$`, `^R_4$`) - cet_fit %quantiles% 1:2, fit = cet_fit) #returns error}
#' #in this case, it would be better to use more manual selection to get the desired result
#' pwc.result <- pwc(cet_fit %quantiles% 3:4 - cet_fit %quantiles% 1:2, fit = cet_fit)
#' #so yeah, the fit argument can be a bit clunky...
#' 
#' #as a more practical example...
#' #maybe we want to compare branchwise rates within the Globicephalinae subfamily?
#' globs <- fit %chains% get.clade.edges(cet_fit$call$tree, c("Pseudorca", "Feresa"))
#' pwc.result <- pwc(globs - globs)
#' #notice how the results here can "blow up" and be hard to interpret...
#' #there are 10 edges in Globicephalinae clade so pairwise comparison results in 10 * 10 = 100 parameters!
#' #a helpful technique here is to instead summarize the results in a matrix; here's a function for this:
#' foo <- function(x, par1, par2 = NULL){
#'    only.par1.flag <- FALSE
#'    if(is.null(par2)){
#'       par2 <- par1
#'       only.par1.flag <- TRUE
#'    }
#'    nms1 <- names(par1)
#'    nms2 <- names(par2)
#'    out <- matrix(x, length(nms1), length(nms2),
#'                  dimnames = list(nms1, nms2))
#'    if(only.par1.flag) diag(out) <- NA
#'    out
#' }
#' #now let's look at the median differences (averaged across each chain) in rates:
#' pwc.result <- foo(rowMeans(pwc.result %quantiles% list(".", 0.5)), globs)
#' #or even the posterior probability associated with the difference:
#' pwc.result <- foo(rowMeans(pwc(globs > globs) %means% "."), globs)
#' #note that parameters on the left hand side of ">" correspond to rows and those on the right to columns
#' #so pwc.result here shows the posterior probability that the row parameter is greater than the column one
#' #we could of course compare branchwise rates among different groups too using this framework
#' meso <- cet_fit %chains% get.clade.edges(cet_fit$call$tree, "Mesoplodon")
#' pwc.result <- foo(rowMeans(pwc(globs > meso) %means% "."), globs, meso)
#' #of course, this is still confusing when we don't know the numbers of each edge in the phylogeny
#' 
#' #here's perhaps a better example where we compare all "tip rates" with better parameter names
#' tip.rates <- setNames(cet_fit %chains% tip.edges(cet_fit), cet_fit$call$tree$tip.label)
#' pwc.result <- foo(rowMeans(pwc(tip.rates > tip.rates) %means% "."), tip.rates)
#' #the result is a large matrix, but it does allow you to compare rates among all tips in the phylogeny
#' #it may be helpful to plot such results, though I haven't spent too much time experimenting with this
#' image(x = 1:nrow(pwc.result), y = 1:ncol(pwc.result), z = t(pwc.result)[,nrow(pwc.result):1],
#'       col = colorRampPalette(c('blue', 'steelblue', 'gray', 'pink2', 'red'))(12),
#'       xlab = '', xaxt = 'n', ylab = '', yaxt = 'n')
#' row.genera <- rev(sapply(strsplit(rownames(pwc.result), '_'), '[', 1))
#' axis(2, at = 1:nrow(pwc.result), labels = row.genera, las = 2)
#' col.genera <- sapply(strsplit(colnames(pwc.result), '_'), '[', 1)
#' axis(1, at = 1:ncol(pwc.result), labels = col.genera, las = 2)
#' #redder colors indicate a row tip likely exhibiting a faster rate than a column tip and vice versa for bluer colors
#' 
#' 
#' @export
pwc<-function(expr,fit=NULL){
  make.Ops.param_blocks.pw<-TRUE
  expr<-substitute(expr)
  if(!is.null(fit)){
    for(i in all.vars(expr)){
      if(!exists(i)){
        #doesn't seem important enough to deal with right now, but you may want to switch to another
        #param_type depending on the param_type of any existing variables in expr
        #^the above would actually be really hard to do because of asynchronous evaluation
        #would be better to provide a type and extra.select argument if it becomes enough of a problem
        assign(i,fit%chains%i)
      }
    }
  }else{
    rm(fit)
  }
  eval(expr)
}