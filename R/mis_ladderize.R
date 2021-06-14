#have to check if it's rude to make method out of ladderize and add ape's ladderize to it...
#change name to ladder for less naming conflicts...
#' @export
ladder<-function(phy,...){
  UseMethod('ladder')
}

#' @export
ladder.phylo<-function(phy,right=TRUE){
  ape::ladderize(phy,right)
}

#' @export
ladder.evorates<-function(phy,right=TRUE,try.all.elements=FALSE){
  lad.tree<-ape::ladderize(phy$tree,right)
  ord<-match(apply(lad.tree$edge,1,paste,collapse=','),apply(phy$tree$edge,1,paste,collapse=','))
  phy$tree<-lad.tree
  if(!is.null(phy$R)){
    phy$R<-phy$R[ord]
  }
  if(try.all.elements){
    for(i in names(phy)[lengths(phy)==nrow(phy$tree$edge)]){
      if(!(i%in%c('X','trait.data','X_0','R','tree'))){
        phy[[i]]<-phy[[i]][ord]
      }
    }
  }
  phy
}

#' @export
ladder.evorates_fit<-function(phy,right=TRUE){
  lad.tree<-ape::ladderize(phy$call$tree,right)
  ord<-match(apply(lad.tree$edge,1,paste,collapse=','),apply(phy$call$tree$edge,1,paste,collapse=','))
  phy$call$tree<-lad.tree
  for(i in names(phy)){
    if(!(i%in%c('call','sampler.control','sampler.params'))){
      if(is.null(dim(phy[[i]]))){
        isvec<-T
        phy[[i]]<-.expand.element(phy[[i]])
      }else{
        isvec<-F
      }
      old.dimnames<-dimnames(phy[[i]])
      for(j in c('R_[1-9][0-9]*$','R_[1-9][0-9]*_dev$')){
        tmp.inds<-grep(j,old.dimnames[['parameters']])
        if(length(tmp.inds)>0){
          tmp.inds<-tmp.inds[ord]
          new.inds<-c((1:min(tmp.inds))[-min(tmp.inds)],tmp.inds,(max(tmp.inds):length(old.dimnames[['parameters']]))[-1])
          phy[[i]]<-.index.element(phy[[i]],new.inds,which(names(old.dimnames)=='parameters'),allow.reorder=T)
        }
      }
      dimnames(phy[[i]])<-old.dimnames
      if(isvec){
        phy[[i]]<-.simplify.element(phy[[i]])
      }
    }
  }
  phy
}

