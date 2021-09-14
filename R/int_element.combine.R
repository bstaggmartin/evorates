.combine.elements<-function(x){
  x<-x[lengths(x)>0]
  x<-lapply(x,.coerce.to.3D)
  out.type<-.get.element.type(x[[1]])
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
  attr(out,'element')<-out.type
  out
}

#' @export
c.loose_element<-function(...,fit=NULL){
  ls<-list(...)
  existence<-lapply(ls,function(ii) try(is.character(ii),silent=TRUE))
  for(i in 1:length(ls)){
    if(inherits(existence,'try-error')){
      ls[[i]]<-deparse(substitute(ls[[i]]))
    }
  }
  types<-lapply(ls,.get.element.type)
  pot.probs<-lengths(types)==0
  types[lengths(types)==0]<-'select'
  if(is.null(fit)){
    if(all(pot.probs)){
      stop('all provided elements appear to be selection(s) to extract from an evorates fit object, yet no fit was supplied')
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
    stop('c() can only be used to combined chains loose elements with those of a single other type')
  }
  out.type<-type.nms[type.has]
  type.inds<-unlist(type.which[type.has])
  if(is.null(type.inds)){
    type.inds<-rep(FALSE,length(ls))
  }
  type.ele<-ls[type.inds]
  #coerce elements of quantiles, means, or diagnostics types to be compatible
  #(return an error if they're not)
  if(length(type.ele)==1){
    type.ele<-list(.expand.element(type.ele[[1]]))
  }else if(length(type.ele)>1){
    check<-do.call(.check.dims.compat,type.ele)
    if(!check[[4]]){
      if(!check[[5]][2]){
        stop('only loose_elements with the same chains can be combined')
      }
      #should never get to below portion with means...
      nms<-lapply(check[[3]],'[[',1)
      uni.nms<-unique(unlist(nms))
      foo<-function(x){
        all(unlist(lapply(nms,function(ii) x%in%ii)))
      }
      avail.nms<-uni.nms[unlist(lapply(uni.nms,foo))]
      if(length(avail.nms)==0){
        stop('loose elements of ',out.type,' type contain no matching quantities')
      }
      if(out.type=='diagnostics'){
        inds<-lapply(nms,match,x=avail.nms)
      }else{
        inds<-rep(list(avail.nms),length(check[[1]]))
      }
      check[[1]]<-lapply(check[[1]],function(ii) setNames(list(ii),out.type))
      check[[1]]<-lapply(1:length(check[[1]]),function(ii) 
        do.call(paste0('.int.',out.type),
                list(check[[1]][[ii]],select=list('.',inds[[ii]]))))
    }
    type.ele<-check[[1]]
  }else{
    out.type<-'chains'
  }
  if(out.type=='chains'){
    nms<-NULL
  }else{
    nms<-dimnames(type.ele[[1]])[[1]]
  }
  #deal with chain stuff
  chain.inds<-types=='chains'
  chain.ele<-ls[chain.inds]
  if(length(chain.ele)==1){
    chain.ele<-list(.expand.element(chain.ele[[1]]))
  }else if(length( chain.ele)>1){
    check<-do.call(.check.dims.compat,chain.ele)
    if(!check[[4]]){
      if(!check[[5]][2]){
        stop('only loose_elements with the same chains can be combined')
      }
      #check for sampler parameters
      flag<-FALSE
      lens<-unlist(lapply(check[[2]],'[',1))
      if(length(unique(lens))>2){
        flag<-TRUE
      }else{
        pot.sampler<-which(lens==max(lens))
        param.names<-lapply(check[[3]][pot.sampler],'[[',2)
        is.sampler.names<-unlist(lapply(param.names,function(ii) 
          all(grep(.get.sampler.names(),ii))))
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
        stop('loose chain elements appear to contain mismatching iterations')
      }
    }
    chain.ele<-check[[1]]
  }
  if(out.type!='chains'){
    chain.ele<-lapply(chain.ele,function(ii) setNames(list(ii),'chains'))
    chain.ele<-lapply(chain.ele,function(ii) 
      do.call(paste0('.int.',out.type),list(ii,list('.',nms))))
  }
  if(!is.null(fit)){
    which.select<-which(types=='select')
    for(i in which.select){
      ls[[i]]<-do.call(paste0('.int.',out.type),list(fit,list(ls[[i]],nms)))
    }
  }
  ls[type.inds]<-type.ele
  ls[chain.inds]<-chain.ele
  .simplify.element(.combine.elements(ls))
}

#' @export
ele.c<-function(...,fit=NULL){
  evorates:::c.loose_element(...,fit=fit)
}
