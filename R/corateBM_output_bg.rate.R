#' @export
get.bg.rate<-function(fit,node=NULL,element=c('chains','quantiles','means','MAPs'),edge.group=NULL,
                      name=NULL,select.extra=NULL){
  if(!inherits(fit,'corateBM_fit')){
    stop("fit must be a fitted correlated rates BM fit (class 'corateBM_fit')")
  }
  tree<-fit$call$tree
  try.element<-try(match.arg(element,c('chains','quantiles','means','MAPs')),silent=T)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a correlated rates BM fit: please specify one of the following: 'chains', 'quantiles', 'means', or 'MAPs'")
  }
  element<-try.element
  if(element=='quantiles'){
    int.element<-'chains'
  }else{
    int.element<-element
  }
  if(is.null(edge.group)){
    if(!is.null(node)){
      edge.group<-get.edge.des(tree,node)
    }else{
      stop('must specify a set of edges of interest by setting node equal to set of nodes defining a monophyletic clade OR setting edge.group equal to a vector of edge indices')
    }
  }
  edge.group<-sort(unique(edge.group))
  if(is.null(name)){
    breaks<-which(edge.group[-1]!=edge.group[-length(edge.group)]+1)
    if(length(breaks)==0){
      tmp<-cbind(edge.group[1],edge.group[length(edge.group)])
    }else{
      tmp<-rbind(c(edge.group[1],edge.group[breaks[1]]),
                 cbind(edge.group[breaks+1],c(edge.group[breaks[-1]],edge.group[length(edge.group)])))
    }
    tmp<-sapply(1:nrow(tmp),
                function(ii) if(tmp[ii,2]!=tmp[ii,1]) paste(tmp[ii,],collapse='-') else tmp[ii,1])
    name<-paste(paste('R_',tmp,sep=''),collapse='&')
  }
  R<-try(do.call(paste('.int.',int.element,sep=''),list(fit=fit,select=edge.group)))
  if(inherits(R,'try-error')){
    if(element=='quantiles'&!is.null(select.extra)){
      select<-list('R_0',select.extra)
    }else{
      select<-'R_0'
    }
    do.call(paste('%',element,'%',sep=''),list(fit=fit,select=select))
  }else{
    wgts<-tree$edge.length[edge.group]/sum(tree$edge.length[edge.group])
    out<-array(apply(R,c(1,3),function(ii) sum(ii*wgts)),
               dim=c(dim(R)[1],1,dim(R)[3]),
               c(dimnames(R)[1],parameters=name,dimnames(R)[3]))
    if(element=='quantiles'){
      if(!is.null(select.extra)){
        select<-list(name,select.extra)
      }else{
        select<-name
      }
      tmp<-list(chains=out)
      class(tmp)<-'corateBM_fit'
      out<-.int.quantiles(tmp,select)
    }
    out<-.simplify.element(out)
    out
  }
}