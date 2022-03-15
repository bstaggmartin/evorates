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
    check<-do.call(.check.dims.compat,type.par)
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
    check<-do.call(.check.dims.compat,chain.par)
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

#' @export
par.c<-function(...,fit=NULL){
  c.param_block(...,fit=fit)
}
