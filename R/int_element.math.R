#' @export
Ops.loose_element<-function(l,r=NULL){
  if(is.null(r)){
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
    ll<-deparse(substitute(l))
    rr<-deparse(substitute(r))
    is.l<-.Method[1L]=='Ops.loose_element'
    is.r<-.Method[2L]=='Ops.loose_element'
    if(is.l&is.r){
      out<-.simplify.element(.int.element.math(l,r,.Generic,ll,rr))
    }else{
      if(is.l){
        nms<-names(l)
        if(!is.list(nms)){
          nms<-list(nms)
        }
        names(l)<-lapply(nms,
                         function(ii) paste0('(',ii,
                                             '%',.Generic,'%',
                                             rr,')'))
        l<-.strip.ele(l)
      }else{
        nms<-names(r)
        if(!is.list(nms)){
          nms<-list(nms)
        }
        names(r)<-lapply(nms,
                         function(ii) paste0('(',ll,
                                             '%',.Generic,'%',
                                             ii,')'))
        r<-.strip.ele(r)
      }
      out<-do.call(.Generic,list(l,r))
      out<-.add.ele(out)
    }
  }
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

.int.element.math<-function(l,r,.Generic,in.l,in.r){
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
  if(lsing|rsing){
    check<-.check.dims.compat(l,r)
    if(is.null(check)){
      stop('Dimensional mismatch between ',
           in.l,
           ' and ',
           in.r,
           ': did these elements come from from the same evorates_fit, and were they extracted in the same way?')
    }else{
      l<-check[[1]]
      r<-check[[2]]
      ldims<-length(check[[3]])
      rdims<-length(check[[4]])
      if(lsing){
        l<-aperm(.coerce.to.3D(check[[1]]),c(1,3,2))
        dims<-length(check[[4]])
        r<-aperm(check[[2]],
                 (1:dims)[c(1,dims,2:(dims-1))])
        out<-r
        nms<-check[[6]][-c(1,length(check[[6]]))]
      }else{
        dims<-length(check[[3]])
        l<-aperm(check[[1]],
                 (1:dims)[c(1,dims,2:(dims-1))])
        r<-aperm(.coerce.to.3D(check[[2]]),c(1,3,2))
        out<-l
        nms<-check[[5]][-c(1,length(check[[5]]))]
      }
      out[]<-do.call(.Generic,list(as.vector(l),as.vector(r)))
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
    }
  }else{
    check<-.check.dims.compat(l,r)
    if(is.null(check)){
      stop('Dimensional mismatch between ',
           in.l,
           ' and ',
           in.r,
           ': did these elements come from from the same evorates_fit, and were they extracted in the same way?')
    }else{
      l<-.coerce.to.3D(check[[1]])
      r<-.coerce.to.3D(check[[2]])
      ldims<-dim(l)
      rdims<-dim(r)
      lnames<-dimnames(l)
      rnames<-dimnames(r)
      out.param.names<-paste0('(',rep(lnames[[2]],rdims[2]),
                              '%',.Generic,'%',
                              rep(rnames[[2]],each=ldims[2]),
                              ')')
      out<-array(NA,
                 c(ldims[1],length(out.param.names),ldims[3]),
                 c(lnames[1],parameters=list(out.param.names),lnames[3]))
      counter<-1
      for(i in 1:rdims[2]){
        for(j in 1:ldims[2]){
          out[,counter,]<-do.call(.Generic,
                                  c(list(l[,j,]),
                                    list(r[,i,])))
          counter<-counter+1
        }
      }
    }
  }
  out
}