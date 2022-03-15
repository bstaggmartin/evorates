#update to support joint or parameter-wise scaling--only does scaling for each parameter/chain combo at moment, which is undesirable
#' @export
scale.loose_element<-function(x,...){
  .add.ele(scale(.strip.ele(x),...))
}

#' @export
cut.loose_element<-function(x,...,joint=TRUE,simplify=TRUE){
  if(joint){
    out<-cut(.strip.ele(x),...)
    levels(x)<-levels(out)
    x[]<-out
    levels(x)<-levels(out)
    x<-.add.ele(x)
  }else{
    x<-.expand.element(x)
    element.type<-.get.element.type(x)
    foo<-function(ii,...){
      list(cut(ii,...))
    }
    outdims<-dim(x)
    outnames<-dimnames(x)
    dimlen<-length(dim(x))
    param.dims<-seq_len(dimlen)[-c(1,dimlen)]
    tmp<-apply(x,param.dims,foo,...)
    x<-aperm(x,c(1,dimlen,param.dims))
    x[]<-unlist(lapply(tmp,function(ii) as.numeric(ii[[1]])))
    x<-aperm(x,c(1,param.dims+1,2))
    attr(x,'element')<-element.type
    x<-.add.ele(x)
    levels(x)<-array(lapply(tmp,function(ii) levels(ii[[1]])),
                     outdims[param.dims],outnames[param.dims])
    if(simplify){
      x<-.simplify.element(x)
    }
  }
  nms<-names(x)
  if(!is.list(nms)){
    nms<-list(nms)
  }
  nms<-lapply(nms,function(ii) paste0('cut(',ii,')'))
  x
}