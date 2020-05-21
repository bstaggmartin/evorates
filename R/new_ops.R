#' @export
check.n.proc<-function(fit,element){
  if(!inherits(fit,'corateBM_fit')){
    stop(paste("the %",element,"% operator only accepts fitted autocorrelated BM objects (class 'corateBM_fit') on left hand side",sep=""))
  }
  if(length(dim(fit[[element]]))==0){
    new.element<-array(fit[[element]],c(1,length(fit[[element]]),1))
    if(element=='quantiles'){
      dimnames(new.element)<-list(attr(fit[[element]],'quantiles'),'parameters'=names(fit[[element]]),NULL)
    }else{
      dimnames(new.element)<-list(NULL,'parameters'=names(fit[[element]]),NULL)
    }
    new.element
  }else if(length(dim(fit[[element]]))==2){
    if(names(dimnames(fit[[element]]))[1]=='parameters'){
      new.element<-array(fit[[element]],c(1,dim(fit[[element]])))
      dimnames(new.element)<-c(list(NULL),dimnames(fit[[element]]))
      new.element
    }else{
      new.element<-array(fit[[element]],c(dim(fit[[element]]),1))
      dimnames(new.element)<-c(dimnames(fit[[element]]),list(NULL))
      new.element
    }
  }else{
    fit[[element]]
  }
}

#' @export
int.op<-function(fit,element,select){
  select<-as.vector(select)
  param.names<-dimnames(fit[[element]])[[2]]
  if(is.numeric(select)){
    n.R<-max(as.numeric(substr(param.names,
                               regexpr('\\[',param.names)+1,
                               regexpr('\\]',param.names)-1)),
             na.rm=T)
    if(all(select>n.R)){
      stop("couldn't find any corresponding parameters")
    }
    if(any(select>n.R)){
      warning(paste("couldn't find parameters corresponding to:",
                    paste(select[select>n.R],collapse=', ')))
      inds<-select[select<=n.R]
    }
    inds<-4+sort(unique(select))
  }else if(is.character(select)){
    select<-unlist(lapply(select,strsplit,split=','))
    new.select<-sub('\\[','\\\\\\[',select)
    new.select<-sub('\\]','\\\\\\]',new.select)
    tmp<-lapply(new.select,function(ii) grep(ii,param.names))
    forbidden.inds<-grep('dev',param.names)
    for(i in grep('dev',select,invert=T)){
      tmp[[i]]<-tmp[[i]][!(tmp[[i]]%in%forbidden.inds)]
    }
    if(all(lengths(tmp)==0)){
      stop("couldn't find any corresponding parameters")
    }
    if(any(lengths(tmp)==0)){
      warning(paste("couldn't find parameters corresponding to:",
                    paste(select[lengths(tmp)==0],collapse=', ')))
    }
    inds<-sort(unique(unlist(tmp)))
  }else{
    stop(paste("the %",element,"% operator only accepts numeric or character vectors on right hand side",sep=""))
  }
  out<-fit[[element]][,inds,]
  if(length(inds)==1){
    attr(out,'parameters')<-param.names[inds]
  }
  out
}

#' @export
`%chains%`<-function(fit,select){
  fit[['chains']]<-check.n.proc(fit,'chains')
  int.op(fit,'chains',select)
}

#' @export
`%quantiles%`<-function(fit,select){
  if(is.null(fit[['quantiles']])){
    def.report.quantiles<-c(0.025,0.5,0.975)
  }else{
    fit[['quantiles']]<-check.n.proc(fit,'quantiles')
    def.report.quantiles<-as.numeric(substr(dimnames(fit[['quantiles']])[[1]],0,
                                            nchar(dimnames(fit[['quantiles']])[[1]])-1))/100
  }
  if(is.list(select)){
    if(length(select)<2){
      select[[2]]<-def.report.quantiles
    }
  }else{
    select<-list(select,def.report.quantiles)
  }
  fit[['chains']]<-check.n.proc(fit,'chains')
  fit[['quantiles']]<-apply(fit[['chains']],c(2,3),quantile,probs=select[[2]])
  if(is.matrix(fit[['quantiles']])){
    fit[['quantiles']]<-array(fit[['quantiles']],dim=c(1,dim(fit[['quantiles']])))
    dimnames(fit[['quantiles']])<-c('quantiles'=paste(def.report.quantiles*100,'%',sep=''),dimnames(fit[['chains']])[-1])
  }
  out<-int.op(fit,'quantiles',select[[1]])
  if(length(dim(out))==0){
    if(length(out)==length(select[[2]])){
      names(out)<-paste(select[[2]]*100,'%',sep='')
    }else{
      attr(out,'quantiles')<-paste(select[[2]]*100,'%',sep='')
    }
  }
  out
}

#' @export
`%MAPs%`<-function(fit,select){
  if(is.null(fit[['MAPs']])){
    fit[['chains']]<-check.n.proc(fit,'chains')
    fit[['MAPs']]<-apply(fit[['chains']],c(2,3),mean)
  }
  fit[['MAPs']]<-check.n.proc(fit,'MAPs')
  out<-int.op(fit,'MAPs',select)
  if(length(out)==1){
    names(out)<-attr(out,'parameters')
  }
  attr(out,'parameters')<-NULL
  out
}

#' @export
select.chains<-function(fit,chains){
  if(any(lengths(sapply(fit,dim))==0)){
    stop("Your autocorrelated BM fit appears to have already been simplified to a single chain")
  }
  fit$chains<-fit$chains[,,chains]
  if(!is.null(fit$quantiles)){
    report.quantiles<-as.numeric(substr(dimnames(fit$quantiles)[[1]],0,
                                        nchar(dimnames(fit$quantiles)[[1]])-1))/100
    fit$quantiles<-fit$quantiles[,,chains]
    if(length(dim(fit$quantiles))==0){
      attr(fit$quantiles,'quantiles')<-paste(report.quantiles*100,'%',sep='')
    }
  }
  if(!is.null(fit$MAPs)){
    fit$MAPs<-fit$MAPs[,chains]
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-fit$post.probs[,chains]
  }
  if(length(chains)==1){
    attr(fit,'chains')<-paste('chain',chains)
  }
  fit
}

#' @export
combine.chains<-function(fit){
  if(any(lengths(sapply(fit,dim))==0)){
    stop("Your autocorrelated BM fit appears to have already been simplified to a single chain")
  }
  fit$chains<-do.call(rbind,lapply(1:dim(fit$chains)[3],function(ii) fit$chains[,,ii]))
  names(dimnames(fit$chains))<-c('iterations','parameters')
  if(!is.null(fit$quantiles)){
    report.quantiles<-as.numeric(substr(dimnames(fit$quantiles)[[1]],0,
                                        nchar(dimnames(fit$quantiles)[[1]])-1))/100
    fit$quantiles<-apply(fit$chains,2,quantile,probs=report.quantiles)
    if(length(dim(fit$quantiles))==0){
      attr(fit$quantiles,'quantiles')<-paste(report.quantiles*100,'%',sep='')
    }
  }
  if(!is.null(fit$MAPs)){
    fit$MAPs<-apply(fit$chains,2,mean)
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-apply(fit%chains%'dev',2,function(ii) sum(ii>0)/length(ii))
  }
  fit
}