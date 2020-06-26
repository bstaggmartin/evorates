#check and process -- make sure fit is of appropriate class and coerce fit$element to an array
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

#internal operator -- select parameters out of fit$element based on edge numbers or regular expressions
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
    inds<-2+sort(unique(select))
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

#' @name operators
#' @rdname operators
#'
#' @title New operators
NULL

#' @rdname operators
#' @export
`%chains%`<-function(fit,select){
  fit[['chains']]<-check.n.proc(fit,'chains')
  int.op(fit,'chains',select)
}

#' @rdname operators
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

#' @rdname operators
#' @export
`%means%`<-function(fit,select){
  if(is.null(fit[['means']])){
    fit[['chains']]<-check.n.proc(fit,'chains')
    fit[['means']]<-apply(fit[['chains']],c(2,3),mean)
  }
  fit[['means']]<-check.n.proc(fit,'means')
  out<-int.op(fit,'means',select)
  if(length(out)==1){
    names(out)<-attr(out,'parameters')
    attr(out,'parameters')<-NULL
  }
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
  if(!is.null(fit$means)){
    fit$means<-fit$means[,chains]
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
  if(!is.null(fit$means)){
    fit$means<-apply(fit$chains,2,mean)
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-apply(fit%chains%'dev',2,function(ii) sum(ii>0)/length(ii))
  }
  fit
}

#a more bare-bones version of the coerce to array part of check.n.proc focused on the 'chains'
#element; for use in trace.plot and profile.plot
#' @export
chains.to.array<-function(chains){
  if(length(dim(chains))==0){
    new.chains<-as.matrix(chains)
    colnames(new.chains)<-attr(chains,'parameters')
    chains<-new.chains
  }
  if(length(dim(chains))==2){
    if(sum(grepl('chains',names(dimnames(chains))[2]))!=0){
      new.chains<-array(chains,c(dim(chains)[1],1,dim(chains)[2]))
      dimnames(new.chains)<-c(dimnames(chains)[1],list(attr(chains,'parameters')),dimnames(chains)[2])
      chains<-new.chains
    }else{
      new.chains<-array(chains,c(dim(chains),1))
      dimnames(new.chains)<-c(dimnames(chains),list(NULL))
      chains<-new.chains
    }
  }
  chains
}

#function for getting the final list of arguments lists in a call to test.plot or profile.plot
#' @export
get.args.master<-function(chains,separate.R,separate.X,separate.dev,together,...){
  args.master<-vector(mode='list',length=dim(chains)[2])
  names(args.master)<-paste('args.',dimnames(chains)[[2]],sep='')
  if(!separate.R){
    args.master<-args.master[!grepl('^args\\.R\\[\\d+\\]$',names(args.master))]
    if(is.null(args.master$args.R)){
      args.master<-c(args.master,'args.R'=list(NULL))
    }
  }
  if(!separate.X){
    args.master<-args.master[grepl(paste('^args.R0','Rsig2','X0','bg.rate','R\\[\\d+\\]',
                                         'R\\[\\d+\\] dev','Rmu','Ysig2','R','X','dev$',
                                         collapse='',sep='$|^args.'),
                                   names(args.master))]
    if(is.null(args.master$args.X)){
      args.master<-c(args.master,'args.X'=list(NULL))
    }
  }
  if(!separate.dev){
    args.master<-args.master[!grepl('^args\\.R\\[\\d+\\] dev$',names(args.master))]
    if(is.null(args.master$args.dev)){
      args.master<-c(args.master,'args.dev'=list(NULL))
    }
  }
  for(i in names(args.master)){
    if(i%in%names(list(...))){
      args.master[[i]]<-list(...)[[i]]
      if(grepl('^args\\.R\\[\\d+\\]$',i)){
        edge.ind<-paste('args.',
                        substr(i,regexpr('\\[',i)+1,regexpr('\\]',i)-1),
                        sep='')
        if(edge.ind%in%names(list(...))){
          warning(paste('multiple argument lists matched with parameter ',substr(i,6,nchar(i)),': argument list indicated by integer (',edge.ind,') not used',sep=''))
        }
      }
    }else if(grepl('^args\\.R\\[\\d+\\]$',i)){
      edge.ind<-paste('args.',
                      substr(i,regexpr('\\[',i)+1,regexpr('\\]',i)-1),
                      sep='')
      if(edge.ind%in%names(list(...))){
        args.master[[i]]<-list(...)[[edge.ind]]
      }
    }else{
      args.master[[i]]<-NULL
    }
    def.args<-list(...)
    #def.args<-def.args[!sapply(def.args,is.list)]
    def.args<-def.args[!grepl('args\\.',names(list(...)))]
    def.args<-def.args[!(names(def.args)%in%names(args.master[[i]]))]
    args.master[[i]]<-c(args.master[[i]],def.args)
  }
  args.master
}