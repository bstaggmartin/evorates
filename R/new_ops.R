#' @export
`%chains%`<-function(fit,select){
  if(!inherits(fit,'corateBM_fit')){
    stop("the %chains% operator only accepts fitted autocorrelated BM objects (class 'corateBM_fit') on left hand side")
  }
  if(is.matrix('$'(fit,chains))){
    new.fit<-array('$'(fit,chains),c(dim('$'(fit,chains)),1))
    dimnames(new.fit)<-c(dimnames('$'(fit,chains)),NULL)
    '$'(fit,chains)<-new.fit
  }
  select<-as.vector(select)
  column.names<-'[['(dimnames('$'(fit,chains)),2)
  if(is.numeric(select)){
    n.R<-max(as.numeric(substr(column.names,
                               regexpr('\\[',column.names)+1,
                               regexpr('\\]',column.names)-1)),
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
    tmp<-lapply(new.select,function(ii) grep(ii,column.names))
    forbidden.inds<-grep('dev',column.names)
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
    stop('the %chains% operator only accepts numeric or character vectors on right hand side')
  }
  if(length(inds)==1&'['(dim('$'(fit,chains)),3)==1){
    out<-'['('$'(fit,chains),,inds,)
    out<-as.matrix(out)
    colnames(out)<-column.names[inds]
    out
  }else{
    '['('$'(fit,chains),,inds,)
  }
}

#' @export
select.chains<-function(fit,chains){
  fit$chains<-fit$chains[,,chains]
  if(!is.null(fit$quantiles)){
    fit$quantiles<-fit$quantiles[,,chains]
  }
  if(!is.null(fit$MAPs)){
    fit$MAPs<-fit$MAPs[,chains]
  }
  if(!is.nulll(fit$post.probs)){
    fit$post.probs<-fit$post.probs[,chains]
  }
}

#' @export
merge.chains<-function(fit){
  fit$chains<-do.call(rbind,lapply(1:dim(fit$chains)[3],function(ii) fit$chains[,,ii]))
  if(!is.null(fit$quantiles)){
    report.quantiles<-as.numeric(substr(dimnames(fit$quantiles)[[1]],0,
                                        nchar(dimnames(fit$quantiles)[[1]])-1))/100
    fit$quantiles<-apply(fit$chains,2,quantile,probs=report.quantiles)
  }
  if(!is.null(fit$MAPs)){
    fit$MAPs<-apply(fit$chains,2,mean)
  }
  if(!is.null(fit$post.probs)){
    fit$post.probs<-apply(fit%chains%'dev',2,function(ii) sum(ii>0)/length(ii))
  }
  fit
}
