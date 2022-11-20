####APE METHODS####

#' @export
Nedge.evorates_fit<-function(phy){
  Nedge.phylo(phy$call$tree)
}

#' @export
Nedge.evorates<-function(phy){
  Nedge.phylo(phy$tree)
}

#' @export
Nnode.evorates_fit<-function(phy,internal.only=TRUE){
  Nnode.phylo(phy$call$tree,internal.only)
}

#' @export
Nnode.evorates<-function(phy,internal.only=TRUE){
  Nnode.phylo(phy$tree,internal.only)
}

#' @export
Ntip.evorates_fit<-function(phy){
  Ntip.phylo(phy$call$tree)
}

#' @export
Ntip.evorates<-function(phy){
  Ntip.phylo(phy$tree)
}

####CUSTOM METHOD: GETTING MATRIX OF EDGE START + END TIMES####

#' Extract start and end times for each edge in a phylogeny
#' 
#' @export
edge.ranges<-function(phy){
  UseMethod('edge.ranges')
}

#' @export
edge.ranges.phylo<-function(phy){
  tmp<-node.depth.edgelength(phy)
  matrix(tmp[phy$edge],ncol=2)
}

#' @export
edge.ranges.evorates<-function(phy){
  edge.ranges.phylo(phy$tree)
}

#' @export
edge.ranges.evorates_fit<-function(phy){
  edge.ranges.phylo(phy$call$tree)
}

####CUSTOM METHODs FOR GETTING VARIOUS EDGE INDICES####
#"tree-walkers" now avoid for loops for the most part and work with polytomies/non-cladwise trees
#they are SUPER fast and scale well to large trees!

#' Extract indices for edges corresponding to tips in a phylogeny
#' 
#' @export
tip.edges<-function(phy,include.names=TRUE){
  UseMethod('tip.edges')
}

#' @export
tip.edges.phylo<-function(phy,include.names=TRUE){
  out<-match(1:Ntip(phy),phy$edge[,2])
  if(include.names){
    names(out)<-phy$tip.label
  }
  out
}

#' @export
tip.edges.evorates_fit<-function(phy,include.names=TRUE){
  tip.edges(phy$call$tree,include.names)
}

#' @export
tip.edges.evorates<-function(phy,include.names=TRUE){
  tip.edges(phy$tree,include.names)
}

#' Extract indices for edges sister to each edge in a phylogeny
#' 
#' @export
sis.edges<-function(phy){
  UseMethod('sis.edges')
}

#' @export
sis.edges.phylo<-function(phy){
  anc<-anc.edges.phylo(phy,commonformat=FALSE)
  len<-length(anc)
  out<-rep(list(integer(0)),len)
  inds<-!is.na(anc)
  anc2<-anc[inds]
  tmp<-seq_len(len)
  des<-c(list(which(!inds)),split(tmp[inds],anc2))
  ndes<-lengths(des)
  di<-ndes==2
  di.des<-unlist(des[di])
  if(length(di.des)>0){
    odds<-seq.int(1,length(di.des),2)
    evens<-odds+1
    odds<-di.des[odds]
    evens<-di.des[evens]
    out[odds]<-evens
    out[evens]<-odds
  }
  poly.des<-des[!di]
  if(length(poly.des)>0){
    ndes<-ndes[!di]
    unlist.poly.des<-unlist(poly.des,use.names=FALSE)
    out[unlist.poly.des]<-rep(poly.des,ndes)
    foo<-function(x){
      tmp<-out[[x]]
      tmp[tmp!=x]
    }
    out[unlist.poly.des]<-lapply(unlist.poly.des,foo)
  }
  out
}

#' @export
sis.edges.evorates<-function(phy){
  sis.edges(phy$tree)
}

#' @export
sis.edges.evorates_fit<-function(phy){
  sis.edges(phy$call$tree)
}

#' Extract indices for edge ancestral to each edge in a phylogeny
#' 
#' @export
anc.edges<-function(phy){
  UseMethod('anc.edges')
}

#' @export
anc.edges.phylo<-function(phy,commonformat=TRUE){
  out<-match(phy$edge[,1],phy$edge[,2])
  if(commonformat){
    nulls<-is.na(out)
    out<-as.list(out)
    out[nulls]<-list(integer(0))
  }
  out
}

#' @export
anc.edges.evorates<-function(phy){
  anc.edges(phy$tree)
}

#' @export
anc.edges.evorates_fit<-function(phy){
  anc.edges(phy$call$tree)
}

#' Extract indices for edges descending from each edge in a phylogeny
#' 
#' @export
des.edges<-function(phy){
  UseMethod('des.edges')
}

#' @export
des.edges.phylo<-function(phy){
  anc<-anc.edges.phylo(phy,commonformat=FALSE)
  len<-length(anc)
  out<-rep(list(integer(0)),len)
  inds<-!is.na(anc)
  anc<-anc[inds]
  tmp<-split(seq_len(len)[inds],anc)
  out[as.numeric(names(tmp))]<-tmp
  out
}

#' @export
des.edges.evorates<-function(phy){
  des.edges(phy$tree)
}

#' @export
des.edges.evorates_fit<-function(phy){
  des.edges(phy$call$tree)
}

#' Extract indices for edges descending from root node
#' 
#' @export root.edges
root.edges<-function(phy){
  UseMethod('root.edges')
}

#' @export
#' @method root.edges phylo
root.edges.phylo<-function(phy){
  which(phy$edge[,1]==Ntip(phy)+1)
}

#' @export
#' @method root.edges evorates
root.edges.evorates<-function(phy){
  root.edges(phy$tree)
}

#' @export
#' @method root.edges evorates_fit
root.edges.evorates_fit<-function(phy){
  root.edges(phy$call$tree)
}

####CUSTOM METHOD: LADDERIZATION####

#have to check if it's rude to make method out of ladderize and add ape's ladderize to it...
#change name to ladder for less naming conflicts...but could still be rude?

#' Ladderize objects with phylogenetic trees and edgewise data
#' 
#' 
#' This is a generic that extends the \code{ladderize()} function from \code{ape} to objects of class
#' \code{evorates} and \code{evorates_fit}. Use it to make tree plots prettier.
#' 
#' 
#' @param phy An object of class \code{evorates}, \code{evorates_fit}, or \code{phylo}.
#' @param ... Optional arguments including \code{right} and \code{try.all.elements}. When \code{right} is \code{TRUE},
#' smaller clades are generally sorted to the right/bottom of the phylogeny, and vice versa
#' when \code{right} is \code{FALSE} (\code{TRUE} by default). \code{try.all.elements} only works with \code{evorates}
#' objects, and works assuming any vector in the object with the same number of elements as there are edges in the
#' phylogeny corresponds to edgewise data, thereby sorting it (\code{TRUE} by default). This will additionally sort the
#' \code{post.probs} element present in most \code{evorates} objects coerced from \code{evorates_fit} objects, but may cause
#' unpredictable results in cases when users manually customize an \code{evorates} objects.
#' 
#' 
#' @return An object of the same class as \code{phy}, with reordered phylogeny and edgewise data.
#' 
#' 
#' @seealso \link[ape]{ladderize}
#' 
#' 
#' @examples
#' #requires example simulation object
#' ladder(example.sim)
#' ladder(example.sim$tree)
#' #requires example fitted model object
#' ladder(example.fit)
#' ladder(example.fit$tree)
#' 
#' 
#' @export
ladder<-function(phy,...){
  UseMethod('ladder')
}

#' @export
ladder.phylo<-function(phy,right=TRUE){
  ape::ladderize(phy,right)
}

#' @export
ladder.evorates<-function(phy,right=TRUE,try.all.elements=TRUE){
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
  loop.inds<-c('^R_[1-9][0-9]*|[%\\(]R_[1-9][0-9]*','^Rdev_[1-9][0-9]*|[%\\(]Rdev_[1-9][0-9]*')
  par.inds<-which(names(phy)!='call'&names(phy)!='sampler.control'&names(phy)!='sampler.params')
  for(i in par.inds){
    phy[[i]]<-.expand.par(phy[[i]])
    param.nms<-names(phy[[i]])
    #WILL BREAK with 4D arrays...but not really a worry for now
    for(j in loop.inds){
      tmp.inds<-grep(j,param.nms)
      if(length(tmp.inds)>0){
        phy[[i]][,tmp.inds,]<-phy[[i]][,tmp.inds[ord],,drop=FALSE]
      }
    }
  }
  phy
}