#' @export
compare.params<-function(fit=NULL,params1,params2=NULL,post.probs=T,operation=c('-','+','/','*'),
                         element=c('chains','quantiles','means','MAPs'),select.extra=NULL){
  try.operation<-try(match.arg(operation,c('-','+','/','*')),silent=T)
  if(inherits(try.operation,'try-error')){
    stop(operation," is not an available parameter operation: please specify one of the following: '-', '+', '/', or '*'")
  }
  operation<-try.operation
  if(!is.null(element)){
    try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs')),silent=T)
    if(inherits(try.element,'try-error')){
      stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', or 'MAPs' (can also set to NULL to suppress element output)")
    }
    element<-try.element
    if(!is.null(select.extra)){
      select<-list('.|dev',select.extra)
    }else{
      select<-'.|dev'
    }
  }
  if(!(operation%in%c('-','/'))&post.probs){
    warning('post.probs can only be TRUE when parameters are being compared (i.e., subtraction/division): post.probs was set to FALSE')
    post.probs<-F
  }
  if(is.null(element)&!post.probs){
    stop('both element and posterior probability output suppressed: nothing to return')
  }
  if(is.null(fit)&element=='MAPs'){
    stop('desired element is set to MAPs, but no corateBM_fit is supplied--need posterior probability info to find MAPs')
  }
  params1<-.combine.elements(params1,fit,'chains',simplify=F)
  if(is.null(params2)&operation=='-'){
    tmp<-list(chains=params1)
    class(tmp)<-'corateBM_fit'
    if(!is.null(element)){
      out<-list(do.call(paste('%',element,'%',sep=''),list(fit=tmp,select=select)))
      names(out)<-element
    }else{
      out<-list()
    }
    if(operation=='-'){
      out$post.probs<-.simplify.element(apply(params1,c(2,3),function(ii) sum(ii>0)/length(ii)))
    }else{
      out$post.probs<-.simplify.element(apply(params1,c(2,3),function(ii) sum(ii>1)/length(ii)))
    }
    if(length(out)==1){
      out[[1]]
    }else{
      out
    }
  }else{
    params2<-.combine.elements(params2,fit,'chains',simplify=F)
    out.param.names<-paste(rep(dimnames(params1)[[2]],dim(params2)[2]),
                           '%',operation,'%',
                           rep(dimnames(params2)[[2]],each=dim(params1)[2]),
                           sep='')
    out<-list(chains=array(NA,
                           c(dim(params1)[1],length(out.param.names),dim(params1)[3]),
                           c(dimnames(params1)[1],parameters=list(out.param.names),dimnames(params1)[3])))
    for(i in dimnames(params1)[[2]]){
      for(j in dimnames(params2)[[2]]){
        out$chains[,paste(i,'%',operation,'%',j,sep=''),]<-do.call(operation,
                                                                   c(list(params1[,i,]),
                                                                     list(params2[,j,])))
      }
    }
    if(post.probs){
      if(operation=='-'){
        out$post.probs<-apply(out$chains,c(2,3),function(ii) sum(ii>0)/length(ii))
      }else{
        out$post.probs<-apply(out$chains,c(2,3),function(ii) sum(ii>1)/length(ii))
      }
    }
    tmp<-list(chains=out$chains)
    if(element=='MAPs'){
      tmp$sampler.params<-fit$sampler.params
    }
    class(tmp)<-'corateBM_fit'
    out$chains<-NULL
    if(!is.null(element)){
      out[[element]]<-do.call(paste('.int.',element,sep=''),list(fit=tmp,select=select))
    }
    out<-lapply(rev(out),.simplify.element)
    if(length(out)==1){
      out[[1]]
    }else{
      out
    }
  }
}
#add an option to do marginal posterior probabilities?