#allow for simultaneously getting multiple clades via passing a list to node? --> done
#add support for removing trend, if present?
#add diagnostics?
#' @export
get.bg.rate<-function(fit,node=length(fit$call$tree$tip.label)+1,edge.group=NULL,
                      element=c('chains','quantiles','means','diagnostics'),
                      geometric=TRUE,remove.trend=TRUE,
                      name=NULL,select.extra=NULL,simplify=TRUE){
  if(!inherits(fit,'evorates_fit')){
    stop("fit must be a fitted evolving rates model fit (class 'evorates_fit')")
  }
  tree<-fit$call$tree
  try.element<-try(match.arg(element,c('chains','quantiles','means','diagnostics')),silent=TRUE)
  if(inherits(try.element,'try-error')){
    stop(element," is not an available element to extract from a evolving rates model fit: please specify one of the following: 'chains', 'quantiles', 'means', or 'diagnostics'")
  }
  element<-try.element
  if(element=='quantiles'){
    int.element<-'chains'
  }else{
    int.element<-element
  }
  if(remove.trend){
    try.R<-try(fit%chains%1,silent=TRUE)
    try.Rmu<-try(fit%chains%'^R_mu$',silent=TRUE)
    if(!inherits(try.R,'try-error')&!inherits(try.Rmu,'try-error')){
      e<-nrow(tree$edge)
      is.simplified<-length(dim(fit$chains))<3
      if(!is.simplified){
        fit$chains[,paste('R',1:e,sep='_'),]<-remove.trend(fit,'chains')
      }else{
        fit$chains[,paste('R',1:e,sep='_')]<-remove.trend(fit,'chains')
      }
      if(int.element!='chains'){
        fit[[int.element]]<-NULL
      }
    }
  }
  if(!is.list(node)){
    node<-list(node)
  }
  if(!is.list(edge.group)){
    edge.group<-list(edge.group)
  }
  if(!is.list(name)){
    if(is.null(name)){
      name<-list(name)
    }else{
      name<-as.list(name)
    }
  }
  max.len<-max(length(node),length(edge.group))
  out.list<-rep(list(NULL),length.out=max.len)
  tmp<-out.list
  tmp[which(!sapply(node,is.null))]<-node[which(!sapply(node,is.null))]
  node<-tmp
  tmp<-out.list
  tmp[which(!sapply(edge.group,is.null))]<-edge.group[which(!sapply(edge.group,is.null))]
  edge.group<-tmp
  tmp<-out.list
  tmp[which(!sapply(name,is.null))]<-name[which(!sapply(name,is.null))]
  name<-tmp
  for(i in 1:length(out.list)){
    if(is.null(edge.group[[i]])&!is.null(node[[i]])){
      edge.group[[i]]<-try(get.clade.edges(tree,node[[i]]))
    }
  }
  errors<-sapply(edge.group,function(ii) inherits(ii,'try-error'))
  if(all(errors)){
    stop("all specified node groups included nodes that don't exist in tree")
  }
  if(any(errors)){
    warning("removed node groups including nodes that don't exist in tree")
    edge.group<-edge.group[!errors]
    name<-name[!errors]
    out.list<-out.list[!errors]
  }
  nulls<-sapply(edge.group,is.null)
  if(all(nulls)){
    stop('all edge groups are NULL: for each background rate you must specify a set of edges of interest by either setting node equal to set of nodes defining a monophyletic clade OR setting edge.group equal to a vector of edge indices')
  }
  if(any(nulls)){
    warning('removed NULL edge groups: for each background rate you must specify a set of edges of interest by either setting node equal to set of nodes defining a monophyletic clade OR setting edge.group equal to a vector of edge indices')
    edge.group<-edge.group[!nulls]
    name<-name[!nulls]
    out.list<-out.list[!nulls]
  }
  e<-nrow(tree$edge)
  bad.indices<-sapply(edge.group,function(ii) any(ii>e))
  if(all(bad.indices)){
    stop("all edge groups include edge indices that don't exist in tree: are you basing edge groups on the right tree?")
  }
  if(any(bad.indices)){
    warning("removed edge groups including edge indices that don't exist in tree: are you basing edge groups on the right tree?")
    edge.group<-edge.group[!bad.indices]
    name<-name[!bad.indices]
    out.list<-out.list[!bad.indices]
  }
  edge.group<-lapply(edge.group,function(ii) sort(unique(ii)))
  try.R<-try(fit%chains%1,silent=TRUE)
  if(inherits(try.R,'try-error')){
    warning('fit has no rate heterogeneity: background rates will be identical for all specified edge groups')
    has.R<-F
    if((element=='quantiles'|element=='chains')&!is.null(select.extra)){
      select<-list(0,select.extra)
    }else{
      select<-0
    }
    R<-do.call(paste('.int.',element,sep=''),list(fit=fit,select=select))
    tmp.dimnames<-dimnames(R)
  }else{
    has.R<-TRUE
  }
  for(i in 1:length(out.list)){
    if(is.null(name[[i]])){
      breaks<-which(edge.group[[i]][-1]!=edge.group[[i]][-length(edge.group[[i]])]+1)
      if(length(breaks)==0){
        tmp<-cbind(edge.group[[i]][1],edge.group[[i]][length(edge.group[[i]])])
      }else{
        tmp<-rbind(c(edge.group[[i]][1],edge.group[[i]][breaks[1]]),
                   cbind(edge.group[[i]][breaks+1],c(edge.group[[i]][breaks[-1]],edge.group[[i]][length(edge.group[[i]])])))
      }
      tmp<-sapply(1:nrow(tmp),
                  function(ii) if(tmp[ii,2]!=tmp[ii,1]) paste(tmp[ii,],collapse='-') else tmp[ii,1])
      name[[i]]<-paste(paste('R_',tmp,sep=''),collapse='&')
    }
    if(!has.R){
      tmp.dimnames[[2]]<-name[[i]]
      dimnames(R)<-tmp.dimnames
      out.list[[i]]<-R
    }else{
      if(element=='chains'&!is.null(select.extra)){
        select<-list(edge.groups[[i]],select.extra)
      }else{
        select<-edge.group[[i]]
      }
      R<-do.call(paste('.int.',int.element,sep=''),list(fit=fit,select=select))
      wgts<-tree$edge.length[edge.group[[i]]]/sum(tree$edge.length[edge.group[[i]]])
      if(geometric){
        out.list[[i]]<-array(apply(R,c(1,3),function(ii) sum(ii*wgts,na.rm=TRUE)),
                             dim=c(dim(R)[1],1,dim(R)[3]),
                             c(dimnames(R)[1],parameters=name[[i]],dimnames(R)[3]))
      }else{
        out.list[[i]]<-array(log(apply(R,c(1,3),function(ii) sum(exp(ii)*wgts,na.rm=TRUE))),
                             dim=c(dim(R)[1],1,dim(R)[3]),
                             c(dimnames(R)[1],parameters=name[[i]],dimnames(R)[3]))
      }
      if(element=='quantiles'){
        if(!is.null(select.extra)){
          select<-list(name[[i]],select.extra)
        }else{
          select<-name[[i]]
        }
        tmp<-list(chains=out)
        class(tmp)<-'evorates_fit'
        out.list[[i]]<-.int.quantiles(tmp,select)
      }
    }
  }
  .combine.elements(out.list,simplify=simplify)
}