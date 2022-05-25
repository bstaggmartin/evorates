#doesn't seem to return correctly coerced results with remove.trend, geometric, and keep.R options altered!
#also, may want to not return results as named list with keep.R when there's only one edge.group
#' @export
get.bg.rate<-function(fit,
                      node.groups=NULL,
                      edge.groups=NULL,
                      type=c('chains','quantiles','means','diagnostics'),
                      extra.select=NULL,
                      simplify=TRUE,
                      remove.trend=FALSE,geometric=FALSE,
                      log=TRUE,
                      partial.match=TRUE,
                      keep.R=FALSE){
  tree<-fit$call$tree
  edge.groups<-.make.edge.groups.list(node.groups,edge.groups,tree,partial.match)
  type<-.match.type(type)
  #computational speed-ups possible by not defaulting to chain selection than coercing later...but this is much cleaner code
  FUN<-if(remove.trend) evorates::remove.trend else evorates::get.R
  Rs<-lapply(edge.groups,FUN,fit=fit,type="chains",simplify=FALSE)
  foo<-function(x){
    els<-tree$edge.length[x]
    els/sum(els)
  }
  wgts<-lapply(edge.groups,foo)
  bg.Rs<-if(geometric) Rs else lapply(Rs,exp)
  nms<-names(bg.Rs[[1]])[1]
  bg.Rs<-lapply(seq_along(edge.groups),function(ii) sum(bg.Rs[[ii]]*wgts[[ii]]))
  if(type=='diagnostics'){
    inits.Rs<-lapply(edge.groups,FUN,fit=fit,type='diagnostics',extra.select='inits',simplify=FALSE)
    inits.bg.Rs<-if(geometric) inits.Rs else lapply(inits.Rs,exp)
    inits.bg.Rs<-lapply(seq_along(edge.groups),function(ii) sum(inits.bg.Rs[[ii]]*wgts[[ii]]))
  }
  
  #name making
  tmp<-regexpr('R_\\d+|Rdev_\\d+',nms)
  before<-substr(nms,1,tmp+attr(tmp,'match.length')-1)
  before<-gsub('(R_|Rdev_)\\d+$','\\1',before)
  after<-substr(nms,tmp+attr(tmp,'match.length'),nchar(nms))
  nms.vec<-names(edge.groups)
  if(is.null(nms.vec)){
    nms.vec<-rep(NA,length(edge.groups))
  }
  counter<-1
  for(i in edge.groups){
    if(is.na(nms.vec[counter])){
      len<-length(i)
      breaks<-which(diff(i)>1)
      if(length(breaks)==0){
        tmp<-cbind(i[1],i[len])
      }else{
        tmp<-rbind(c(i[1],i[breaks[1]]),
                   cbind(i[breaks+1],
                         c(i[breaks[-1]],i[len])))
      }
      tmp<-sapply(1:nrow(tmp),
                  function(ii) if(tmp[ii,2]!=tmp[ii,1]) paste(tmp[ii,],collapse='-') else tmp[ii,1])
      
      
      nms.vec[counter]<-tmp<-paste(paste0(before,tmp,after),collapse=';')
    }else{
      tmp<-nms.vec[counter]
    }
    names(bg.Rs[[counter]])<-paste0('bg_rate(',tmp,')')
    counter<-counter+1
  }
  
  #final output
  if(log&!geometric){
    final.fun<-'log'
  }else if(!log&geometric){
    final.fun<-'exp'
  }else{
    final.fun<-NULL
  }
  foo<-function(x){
    if(!is.null(final.fun)){
      x<-do.call(final.fun,list(x))
    }
    tmp<-list(chains=x,sampler.control=1) #"cheating"--just so sampler.params isn't NULL
    out<-.call.op(type,tmp,list('.',extra.select),FALSE)
  }
  out<-lapply(bg.Rs,foo)
  if(type=='diagnostics'){
    inds<-dimnames(out[[1]])[[1]]=='inits'
    n.inds<-sum(inds)
    if(n.inds){
      for(i in seq_along(out)){
        if(!is.null(final.fun)){
          tmp<-do.call(final.fun,inits.bg.Rs[i])
        }else{
          tmp<-inits.bg.Rs[[i]]
        }
        out[[i]][inds,,]<-rep(tmp,each=n.inds)
      }
    }
  }
  if(keep.R){
    if(!log){
      final.fun<-'exp'
    }else{
      final.fun<-NULL
    }
    Rs<-lapply(Rs,foo)
    if(type=='diagnostics'){
      if(n.inds){
        for(i in seq_along(Rs)){
          if(!is.null(final.fun)){
            tmp<-do.call(final.fun,inits.Rs[i])
          }else{
            tmp<-inits.Rs[[i]]
          }
          Rs[[i]][inds,,]<-rep(tmp,each=n.inds)
        }
      }
    }
    if(simplify){
      out<-lapply(out,.simplify.par)
      Rs<-lapply(Rs,.simplify.par)
    }
    out<-lapply(seq_along(edge.groups),
                function(ii) list(R=Rs[[ii]],bg_rate=out[[ii]]))
    if(length(out)==1){
      out<-out[[1]]
    }else{
      names(out)<-nms.vec
    }
  }else{
    out<-.combine.par(out)
    if(simplify){
      out<-.simplify.par(out)
    }
  }
  
  out
}