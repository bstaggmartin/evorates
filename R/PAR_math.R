#add auto-conversion capabilities similar to combining functions?

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
          stop('math between param_blocks and plain matrices/arrays not allowed')
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

#' @export
pwc<-function(expr,fit=NULL){
  make.Ops.param_blocks.pw<-TRUE
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
  }else{
    rm(fit)
  }
  eval(expr)
}