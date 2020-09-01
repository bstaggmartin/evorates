#' @export
combine.elements<-function(in.params,fit=NULL,element=NULL,select.extra=NULL,simplify=T){
  params<-in.params
  if(.is.corateBM.element(params)){
    params<-list(params)
  }
  if(!is.list(params)){
    params<-as.list(params)
  }
  types<-sapply(params,function(ii) if(.is.corateBM.element(ii)) 'element' else 'select')
  if(all(types=='select')&is.null(fit)){
    stop(deparse(substitute(in.params)),' appears to consist of strings/numbers specifying parameters to extract out of a corateBM_fit, but no corateBM_fit is supplied')
  }
  if(any(types=='select')&is.null(fit)){
    warning(deparse(substitute(in.params)),' contains ',paste(params[which(types=='select')],collapse=', '),', which appear to be strings/numbers specifying parameters to extract out of a corateBM_fit, but no corateBM_fit is supplied: these strings/numbers were excluded')
    params<-params[-which(types=='select')]
    types<-types[-which(types=='select')]
  }
  params[types=='element']<-lapply(params[types=='element'],function(ii) .coerce.to.3D(.expand.element(ii)))
  select.params<-NULL
  if(sum(types=='select')>0){
    select<-unlist(params[types=='select'])
    if(is.null(element)){
      if(sum(types=='element')>0){
        element.types<-sapply(params[types=='element'],.get.element.type)
        if(any(element.types%in%c('ambiguous','unrecognized'))|length(unique(element.types))>1){
          stop('element is unspecified, but element type based on provided loose elements is ambiguous: try specifying which element you wish to extract and double-check all loose elements are of the same type (chains, quantiles, means, etc.)')
        }else{
          element<-unique(element.types)
        }
      }else{
        element<-'chains'
      }
    }else{
      try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs','diagnostics','sampler')),silent=T)
      if(inherits(try.element,'try-error')){
        stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', 'MAPs', 'diagnostics', or 'sampler'")
      }
      element<-try.element
    }
    if(element=='quantiles'|element=='diagnostics'){
      if(is.null(select.extra)&sum(types=='select')>0){
        dim1.names<-lapply(params[types=='element'],function(ii) dimnames(ii)[[1]])
        if(length(unique(dim1.names))>1){
          stop('mismatching quantiles and/or parameter diagnostics in provided loose elements')
        }else if(element=='quantiles'){
          select.extra<-as.numeric(substr(dim1.names[[1]],0,nchar(dim1.names[[1]])-1))/100
        }else{
          select.extra<-dim1.names
        }
      }
    }
    if(!is.null(select.extra)){
      select<-list(select,select.extra)
    }
    select.params<-do.call(paste0('.int.',element),list(fit=fit,select=select))
  }
  if(is.null(params[types=='element'])){
    out<-select.params
  }else{
    if(!is.null(select.params)){
      out<-c(list(select.params),params[types=='element'])
    }else{
      out<-params[types=='element']
    }
    out.dim<-vector(mode='list',length=3)
    out.dimnames<-vector(mode='list',length=3)
    for(i in 1:3){
      if(i==2){
        out.dim[[i]]<-sum(sapply(out,function(ii) dim(ii)[i]))
        out.dimnames[[i]]<-unlist(lapply(out,function(ii) dimnames(ii)[[i]]))
      }else{
        out.dim[[i]]<-unique(sapply(out,function(ii) dim(ii)[i]))
        out.dimnames[[i]]<-unique(lapply(out,function(ii) dimnames(ii)[[i]]))
      }
    }
    if(all(lengths(c(out.dim[-2],out.dimnames[-2]))==1)){
      for(i in c(1,3)){
        out.dim[[i]]<-out.dim[[i]][[1]]
        out.dimnames[i]<-out.dimnames[[i]][1]
      }
      out.dim<-unlist(out.dim)
    }else{
      stop('dimensional mismatch in elements specified by ',deparse(substitute(in.params)),': did these all come from the same corateBM_fit, and were they all extracted using the same parameters?')
    }
    names(out.dimnames)<-c(names(dimnames(out[[1]]))[1],'parameters','chains')
    out.arr<-array(NA,out.dim,out.dimnames)
    for(i in 1:length(out)){
      tmp.param.names<-dimnames(out[[i]])[[2]]
      out.arr[,tmp.param.names,]<-out[[i]]
    }
    out<-out.arr
  }
  if(simplify){
    out<-.simplify.element(out)
  }
  out
}
#can get duplicates of the same parameter and probs not super efficient since it runs a separate call to %chains% each time it runs...
#could cannabalize this function to add parameters to your corateBM_fit...
#yeah, I did a minor update to begin generalizing this, but ultimately I'd like it to be a combine.elements function...
#need to figure out what to do in cases where the element is ambiguous
#also, should 4-D arrays be coerced to 3-D with .expand.element?
#now generalized to simply combine all elements it's face with, including ones specified by select. Tries its best to match
#element type, but doesn't do anything like coercing provided elements to other types on the fly (considering doing this in the
#future...)

#' @export
`%select%`<-function(element,select){
  is.char<-try(is.character(select),silent=T)
  if(inherits(is.char,'try-error')){
    select<-deparse(substitute(select))
  }
  element<-.coerce.to.3D(.expand.element(element))
  if(is.numeric(select)){
    select<-paste0('^',dimnames(element)[[2]][select],'$')
  }
  out<-.int.op(element,select)
  .simplify.element(out)
}
