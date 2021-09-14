#add auto-conversion capabilities similar to combining functions?

#' @export
Ops.loose_element<-function(l,r=NULL){
  make.Ops.loose_elements.pw<-mget('make.Ops.loose_elements.pw',
                                   envir=parent.frame(),
                                   ifnotfound=FALSE)[[1]]
  if(is.null(r)){
    type<-.get.element.type(l)
    nms<-names(l)
    if(!is.list(nms)){
      nms<-list(nms)
    }
    names(l)<-lapply(nms,
                     function(ii) paste0('(',
                                         '%',.Generic,'%',
                                         ii,
                                         ')'))
    l<-.strip.ele(l)
    out<-do.call(.Generic,list(l))
    out<-.add.ele(out)
    
  }else{
    is.l<-.Method[1L]=='Ops.loose_element'
    is.r<-.Method[2L]=='Ops.loose_element'
    if(is.l&is.r){
      type<-.get.element.type(l)
      ll<-deparse(substitute(l))
      rr<-deparse(substitute(r))
      out<-.simplify.element(.int.element.math(l,r,.Generic,ll,rr,
                                               make.Ops.loose_elements.pw))
    }else{
      out<-NULL
      if(is.l){
        type<-.get.element.type(l)
        vec<-length(r)>1
        if(vec){
          r<-as.vector(r)
          l<-.coerce.to.3D(l)
        }
        nms<-names(l)
        if(!is.list(nms)){
          nms<-list(nms)
        }
        if(vec){
          lens<-c(length(nms[[1]]),length(r))
          if(lens[1]!=lens[2]){
            max.len<-max(lens)
            l<-l[,rep(1:lens[1],length.out=max.len),,drop=FALSE]
            l<-.add.ele(l)
            r<-rep(r,length.out=max.len)
          }
        }
        names(l)<-lapply(nms,
                         function(ii) paste0('(',ii,
                                             '%',.Generic,'%',
                                             r,')'))
        l<-.strip.ele(l)
        if(vec){
          out<-.simplify.element(
            sweep(l,2,r,function(x,y) do.call(.Generic,list(x,y)))
          )
        }
      }else{
        type<-.get.element.type(r)
        vec<-length(l)>1
        if(vec){
          l<-as.vector(l)
          r<-.coerce.to.3D(r)
        }
        nms<-names(r)
        if(!is.list(nms)){
          nms<-list(nms)
        }
        if(vec){
          lens<-c(length(l),length(nms[[1]]))
          if(lens[1]!=lens[2]){
            max.len<-max(lens)
            l<-rep(l,length.out=max.len)
            r<-r[,rep(1:lens[2],length.out=max.len),,drop=FALSE]
            r<-.add.ele(r)
          }
        }
        names(r)<-lapply(nms,
                         function(ii) paste0('(',l,
                                             '%',.Generic,'%',
                                             ii,')'))
        r<-.strip.ele(r)
        if(vec){
          out<-.simplify.element(
            sweep(r,2,l,function(x,y) do.call(.Generic,list(y,x)))
          )
        }
      }
      if(is.null(out)){
        out<-do.call(.Generic,list(l,r))
      }
      out<-.add.ele(out)
    }
  }
  attr(out,'element')<-type
  out
}

#' @export
Math.loose_element<-function(x,...){
  tmp<-do.call(paste,c(list(...),collapse=list(';')))
  if(nzchar(tmp)){
    tmp<-paste0(';',tmp)
  }
  nms<-names(x)
  if(!is.list(nms)){
    nms<-list(nms)
  }
  names(x)<-lapply(nms,
                   function(ii) paste0(.Generic,
                                       '(',
                                       ii,
                                       tmp,
                                       ')'))
  x<-.strip.ele(x)
  out<-do.call(.Generic,c(list(x),list(...)))
  out<-.add.ele(out)
  out
}

#will require combine function, but idea is simple--just apply .Generic across parameters to get sums of parameter combos, etc.
#' @export
Summary.loose_element<-function(...,na.rm=TRUE){
  out<-do.call(c,list(...))
  element.type<-.get.element.type(out)
  outdims<-dim(out)
  outnames<-dimnames(out)
  if(eval(.Generic)=='range'){
    out.param.names<-paste0(c('min','max'),
                            '(',
                            paste(outnames[[2]],collapse=';'),
                            ')')
    out<-.add.ele(aperm(apply(out,c(1,3),.Generic,na.rm=na.rm),c(2,1,3)))
    outnames[[2]]<-out.param.names
    dimnames(out)<-outnames
  }else{
    out.param.names<-paste0(.Generic,
                            '(',
                            paste(outnames[[2]],collapse=';'),
                            ')')
    out<-.add.ele(apply(out,c(1,3),.Generic,na.rm=na.rm))
    attr(out,'parameters')<-out.param.names
  }
  attr(out,'element')<-element.type
  .simplify.element(out)
}

.int.element.math<-function(l,r,.Generic,in.l,in.r,make.Ops.loose_elements.pw){
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
  check<-.check.dims.compat(l,r)
  if(!check[[4]]){
    stop('Dimensional mismatch between ',
         in.l,
         ' and ',
         in.r,
         ': did these elements come from from the same evorates_fit, and were they extracted in the same way?')
  }
  if(lsing|rsing){
    lr<-check[[1]]
    lrdims<-check[[2]]
    #in the case that both are singular, you don't necessarilly need to coerce anything to 3D...
    #below is most of the way there, but you have to account for cases where number of parameter dimensions mismatches...
    # if(lsing&rsing){
    #   lr<-lapply(lr,.strip.ele)
    #   out<-do.call(.Generic,lr)
    #   out<-.add.ele(out)
    #   nms<-check[[3]]
    #   names(out)<-lapply(1,
    #                     )
    # }
    if(lsing) sing<-0 else sing<-1
    notsing<-(!sing)+1
    sing<-sing+1
    lr[[sing]]<-aperm(.coerce.to.3D(lr[[sing]]),c(1,3,2))
    dims<-length(lrdims[[notsing]])
    lr[[notsing]]<-aperm(lr[[notsing]],
                        (1:dims)[c(1,dims,2:(dims-1))])
    nms<-check[[3]][[notsing]]
    nms<-nms[-c(1,length(nms))]
    out<-lr[[notsing]]
    out[]<-do.call(.Generic,list(as.vector(lr[[1]]),as.vector(lr[[2]])))
    if(dims>3){
      out<-aperm(out,c(3:dims,1,2))
    }else{
      out<-aperm(out,c(1,3,2))
    }
    out<-.add.ele(out)
    if(lsing){
      names(out)<-lapply(nms,
                         function(ii) paste0('(',names(.add.ele(l)),
                                             '%',.Generic,'%',
                                             ii,')'))
    }else{
      names(out)<-lapply(nms,
                         function(ii) paste0('(',ii,
                                             '%',.Generic,'%',
                                             names(.add.ele(r)),')'))
    }
  }else{
    l<-.coerce.to.3D(check[[1]][[1]])
    r<-.coerce.to.3D(check[[1]][[2]])
    ldims<-dim(l)
    rdims<-dim(r)
    lnames<-dimnames(l)
    rnames<-dimnames(r)
    if(make.Ops.loose_elements.pw){
      out.param.names<-paste0('(',rep(lnames[[2]],rdims[2]),
                              '%',.Generic,'%',
                              rep(rnames[[2]],each=ldims[2]),')')
      out<-.add.ele(array(NA,
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
      out<-.add.ele(do.call(.Generic,list(.strip.ele(l),.strip.ele(r))))
      names(out)<-paste0('(',lnames[[2]],
                         '%',.Generic,'%',
                         rnames[[2]],')')
    }
  }
  out
}

#' @export
pwc<-function(expr,fit=NULL){
  make.Ops.loose_elements.pw<-TRUE
  if(inherits(expr,'formula')){
    expr<-expr[[length(expr)]]
  }else if(typeof(expr)!='language'){
    if(is.character(expr)){
      expr<-parse(text=expr)
    }else{
      expr<-substitute(expr)
    }
  }
  if(!is.null(fit)){
    for(i in all.vars(expr)){
      if(!exists(i)){
        assign(i,fit%chains%i)
      }
    }
  }
  eval(expr)
}