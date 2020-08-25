#' @export
ladderize.sim<-function(sim,right=T,try.all.elements=T){
  lad.tree<-ladderize(sim$tree,right)
  ord<-match(apply(lad.tree$edge,1,paste,collapse=','),apply(sim$tree$edge,1,paste,collapse=','))
  sim$tree<-lad.tree
  if(!is.null(sim$R)){
    sim$R<-sim$R[ord]
  }
  if(try.all.elements){
    for(i in names(sim)[lengths(sim)==nrow(sim$tree$edge)]){
      if(!(i%in%c('X','trait.data','X_0','R','tree'))){
        sim[[i]]<-sim[[i]][ord]
      }
    }
  }
  sim
}

#' @export
ladderize.fit<-function(fit,right=T){
  lad.tree<-ladderize(fit$call$tree,right)
  ord<-match(apply(lad.tree$edge,1,paste,collapse=','),apply(fit$call$tree$edge,1,paste,collapse=','))
  fit$call$tree<-lad.tree
  for(i in names(fit)){
    if(!(i%in%c('call','sampler.control','sampler.params'))){
      if(is.null(dim(fit[[i]]))){
        isvec<-T
        fit[[i]]<-.expand.element(fit[[i]])
      }else{
        isvec<-F
      }
      old.dimnames<-dimnames(fit[[i]])
      for(j in c('R_[1-9][0-9]*$','R_[1-9][0-9]*_dev$')){
        tmp.inds<-grep(j,old.dimnames[['parameters']])
        if(length(tmp.inds)>0){
          tmp.inds<-tmp.inds[ord]
          new.inds<-c((1:min(tmp.inds))[-min(tmp.inds)],tmp.inds,(max(tmp.inds):length(old.dimnames[['parameters']]))[-1])
          fit[[i]]<-.index.element(fit[[i]],new.inds,which(names(old.dimnames)=='parameters'),allow.reorder=T)
        }
      }
      dimnames(fit[[i]])<-old.dimnames
      if(isvec){
        fit[[i]]<-.simplify.element(fit[[i]])
      }
    }
  }
  fit
}
