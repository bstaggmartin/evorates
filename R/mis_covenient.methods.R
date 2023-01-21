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
    inds<-which(lengths(phy)==nrow(phy$tree$edge)&
                  !(names(phy)%in%c('R_0',
                                    'R_sig2',
                                    'X_0',
                                    'R',
                                    'R_mu',
                                    'tree',
                                    'X',
                                    'Xsig2',
                                    'k',
                                    'trait.data',
                                    'Ysig2',
                                    'n.obs')
                    )
                )
    for(i in inds){
      phy[[i]]<-phy[[i]][ord]
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

####CUSTOM METHOD: DROP.TIPS####

#have to check if it's rude to make method out of drop.tip and add ape's drop.tip to it...
#change name to drop.tipS for less naming conflicts...but could still be rude?

#' Drop tips from objects with phylogenetic trees and edgewise/nodewise data
#' 
#' 
#' This is a generic that extends the \code{drop.tip()} function from \code{ape} to objects of class
#' \code{evorates} and \code{evorates_fit}. Mostly intended for dropping "pseudo-tips" used for
#' fossil calibrations, but could be helpful in a number of scenarios.
#' 
#' 
#' @param phy An object of class \code{evorates}, \code{evorates_fit}, or \code{phylo}.
#' @param tip A numeric or character vector specifying which tips to drop from \code{phy}.
#' @param partial.match \code{TRUE} or \code{FALSE}: if \code{tip} is a character vector, should it be
#' partially matched to tip labels? For example, one could write a genus name to select all
#' tips belonging to that genus. Defaults to \code{FALSE} to align with normal \link[ape]{drop.tip}
#' behavior by default.
#' @param entire.clade \code{TRUE} or \code{FALSE}: should \emph{all} descendants of the
#' most recent common ancestor of tips in \code{tip} be removed? Allows for easy exclusion of
#' entire monophyletic clades. Defaults to \code{FALSE} to align with normal \link[ape]{drop.tip}
#' behavior by default.
#' @param invert \code{TRUE} or \code{FALSE}: should \code{tip} specify tips to keep, rather than
#' to remove, as in \link[ape]{keep.tip}? Defaults to \code{FALSE} to align with normal
#' \link[ape]{drop.tip} behavior by default.
#' @param ... Optional arguments. Notably, arguments from \link[ape]{drop.tip} can only be used
#' with \code{phylo} objects. In the cases of \code{evorates} and \code{evorates_fit}, arguments
#' from \link[ape]{drop.tip} are set to their defaults with the exception of \code{collapse.singles},
#' which is set to \code{FALSE} (though I hope to allow this option to be \code{TRUE} in the future).
#' Other optional arguments in the case of \code{evorates} and
#' \code{evorates_fit} objects include:
#' \itemize{
#' \item{\code{try.all.elements} (for \code{evorates} objects only), which, if \code{TRUE} (the default),
#' causes the function to assume any vector with the same number of elements as there are edges in the
#' phylogeny corresponds to edgewise data and reshuffles indices accordingly. This will additionally sort the
#' \code{post.probs} element present in most \code{evorates} objects coerced from \code{evorates_fit} objects,
#' but may cause unpredictable results in cases when users manually customize an \code{evorates} objects.}
#' \item{\code{store.traits}, which, if \code{TRUE} (the default), causes the function to infer/store trait
#' values for all nodes in \code{phy}. This is helpful because dropping tips will inevitably alter the
#' inference of trait values--particularly ancestral state reconstructions--and this is probably undesirable
#' in most cases.}
#' \item{\code{recalc.bg} (for \code{evorates_fit} objects only), which, if \code{TRUE} (defaults to \code{FALSE}),
#' causes the function to recalculate the background rate ("\code{bg_rate}") parameter based on the new
#' phylogenetic tree, which may be desirable in some cases.}
#' \item{\code{recalc.devs} (for \code{evorates_fit} objects only), which, if \code{TRUE} (defaults to \code{FALSE}),
#' causes the function to recalculate rate deviations/posterior probabilities to determine which branchwise
#' rates are "anomalously" fast or slow based on the new phylogenetic tree, which may be desirable in some cases.
#' The arguments \code{remove.trend} and \code{geometric} may also be specified to further control of the (re)calculation
#' of rate deviations (see \link{get.bg.rate} for details).}
#' }
#' 
#' 
#' @return An object of the same class as \code{phy}, with reordered phylogeny and
#' edgewise/nodewise data.
#' 
#' 
#' @seealso \link[ape]{drop.tip}, \link{get.bg.rate}
#' 
#' 
#' @examples
#' 
#' 
#' @export
drop.tips<-function(phy,tip,partial.match,entire.clade,invert,...){
  UseMethod('drop.tips')
}

.proc.tip<-function(phy,tip,partial.match,entire.clade,invert){
  if(is.character(tip)&partial.match){
    tip<-unique(unlist(lapply(tip,grep,x=phy$tip.label),use.names=FALSE))
  }
  if(entire.clade){
    tip<-phytools::getDescendants(phy,ape::getMRCA(phy,tip))
    tip<-tip[tip<Ntip(phy)]
  }
  if(invert){
    if(is.character(tip)){
      tip<-which(phy$tip.label%in%tip)
    }
    tip<-seq_along(phy$tip.label)[-tip]
  }
  tip
}

#notably, node labels might present something of a problem...
.int.drop.tip<-function(phy,tip){
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) tip <- which(phy$tip.label %in% tip)
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip)) return(phy)
  if (length(tip) > Ntip-2) stop("Can only drop down to a minimum of 2 tips!")
  #elminiating tips and internal clades from tree...
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel)) 
      break
    keep[sel] <- FALSE
  }
  phy$edge <- phy$edge[keep, ]
  phy$edge.length <- phy$edge.length[keep]
  phy
}

.relab.nodes<-function(phy,Ntip,Nnode){
  TERMS<-!(phy$edge[,2]%in%phy$edge[,1])
  oldNo.ofNewTips<-phy$edge[TERMS,2]
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  phy$tip.label<-phy$tip.label[oldNo.ofNewTips]
  phy$Nnode <- dim(phy$edge)[1]-n+1L
  newNb<-integer(Ntip+Nnode)
  newNb[Ntip+1] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  phy
}

#' @export
#' @method drop.tips phylo
drop.tips.phylo<-function(phy,tip,partial.match=FALSE,entire.clade=FALSE,invert=FALSE,...){
  tip<-.proc.tip(phy,tip,partial.match,entire.clade,invert)
  ape::drop.tip(phy,tip,...)
}

#' @export
#' @method drop.tips evorates
drop.tips.evorates<-function(phy,tip,
                             partial.match=FALSE,entire.clade=FALSE,invert=FALSE,
                             try.all.elements=TRUE,
                             store.traits=TRUE){
  if(store.traits&nrow(phy$X)<Ntip(phy)+Nnode(phy)){
    #will need to be updated to handle multivariate stuff
    elen<-phy$tree$edge.length
    if(!is.null(phy$R)){
      elen<-elen*exp(phy$R)
    }
    elen[is.na(elen)]<-0
    tmp<-NULL
    for(i in 1:ncol(phy$X)){
      XX<-array(NA,c(1,Ntip(phy)+Nnode(phy),1),
                list(NULL,c(phy$tree$tip.label,Ntip(phy)+seq_len(Nnode(phy))),NULL))
      PP<-XX
      LL<-XX
      XX[,phy$tree$tip.label,]<-phy$X[tree$tip.label,i]
      PP[,phy$tree$tip.label,]<-Inf
      LL[,phy$tree$edge[,2],]<-elen
      LL[,Ntip(phy)+1,]<-0
      tmp<-cbind(tmp,evorates:::.anc.recon(tree,XX,LL,PP,FALSE)[[1]][1,,1])
    }
    colnames(tmp)<-colnames(phy$X)
    phy$X<-tmp
  }else if(!store.traits&nrow(phy$X)==Ntip(phy)+Nnode(phy)){
    phy$X<-phy$X[phy$tree$tip.label,,drop=FALSE]
  }
  
  tip<-.proc.tip(phy$tree,tip,partial.match,entire.clade,invert)
  drop.tree<-.int.drop.tip(phy$tree,tip)
  ord<-match(apply(drop.tree$edge,1,paste,collapse=','),apply(phy$tree$edge,1,paste,collapse=','))
  if(!is.null(phy$R)){
    phy$R<-phy$R[ord]
  }
  if(try.all.elements){
    inds<-which(lengths(phy)==nrow(phy$tree$edge)&
                  !(names(phy)%in%c('R_0',
                                    'R_sig2',
                                    'X_0',
                                    'R',
                                    'R_mu',
                                    'tree',
                                    'X',
                                    'Xsig2',
                                    'k',
                                    'trait.data',
                                    'Ysig2',
                                    'n.obs')
                  )
    )
    for(i in inds){
      phy[[i]]<-phy[[i]][ord]
    }
  }
  
  #probably could be made more efficient
  old.edge<-drop.tree$edge
  drop.tree<-.relab.nodes(drop.tree,Ntip=Ntip(phy),Nnode=Nnode(phy))
  old.nodes<-as.vector(old.edge)
  inds<-which(!duplicated(old.nodes))
  old.nodes<-old.nodes[inds]
  new.nodes<-drop.tree$edge[inds]
  if(nrow(phy$X)<Ntip(phy)+Nnode(phy)){
    old.nodes<-old.nodes[old.nodes<=Ntip(phy)]
    new.nodes<-new.nodes[new.nodes<=Ntip(drop.tree)]
  }
  old.nodes<-c(phy$tree$tip.label,Ntip(phy)+seq_len(Nnode(phy)))[old.nodes]
  ord<-order(new.nodes)
  new.nodes<-c(drop.tree$tip.label,Ntip(drop.tree)+seq_len(Nnode(drop.tree)))[new.nodes]
  phy$X<-phy$X[old.nodes[ord],,drop=FALSE]
  rownames(phy$X)<-new.nodes[ord]
  
  phy$tree<-drop.tree
  
  phy
}

#add option to collapse singles down the road?
#' @export
#' @method drop.tips evorates_fit
drop.tips.evorates_fit<-function(phy,tip,
                                 partial.match=FALSE,entire.clade=FALSE,invert=FALSE,
                                 store.traits=TRUE,
                                 recalc.bg=FALSE,
                                 recalc.devs=FALSE,remove.trend=TRUE,geometric=TRUE){
  if(store.traits){
    phy<-get.post.traits(phy,store.in.fit=TRUE,simplify=FALSE)
  }else{
    phy<-.drop.stored.traits(phy)
  }
  
  #code borrowed from ladder mainly
  #first just align the edges properly
  tip<-.proc.tip(phy$call$tree,tip,partial.match,entire.clade,invert)
  drop.tree<-.int.drop.tip(phy$call$tree,tip)
  ord<-match(apply(drop.tree$edge,1,paste,collapse=','),apply(phy$call$tree$edge,1,paste,collapse=','))
  
  loop.inds<-c('^R_[1-9][0-9]*|[%\\(]R_[1-9][0-9]*','^Rdev_[1-9][0-9]*|[%\\(]Rdev_[1-9][0-9]*')
  par.inds<-which(names(phy)!='call'&names(phy)!='sampler.control'&names(phy)!='sampler.params')
  for(i in par.inds){
    param.nms<-names(phy[[i]])
    par.type<-.get.par.type(phy[[i]])
    #WILL BREAK with 4D arrays...but not really a worry for now
    for(j in loop.inds){
      tmp.inds<-grep(j,param.nms)
      if(length(tmp.inds)>0){
        ndrops<-length(tmp.inds[-ord])
        keeps<-tmp.inds[seq_len(length(tmp.inds)-ndrops)]
        drops<-length(tmp.inds)-(ndrops-1):0
        phy[[i]][,keeps,]<-phy[[i]][,tmp.inds[ord],,drop=FALSE]
        param.nms<-param.nms[-drops]
        phy[[i]]<-phy[[i]][,-rev(tmp.inds)[seq_len(ndrops)],,drop=FALSE]
        phy[[i]]<-.add.par.class(phy[[i]])
        attr(phy[[i]],"param_type")<-par.type
      }
    }
  }
  #now edges are appropriately reordered
  
  #relabel nodes in drop.tree
  if(store.traits){
    old.edge<-drop.tree$edge
  }
  drop.tree<-.relab.nodes(drop.tree,Ntip=Ntip(phy),Nnode=Nnode(phy))
  
  #relabel nodes in stored trait data
  #probably could be made more efficient
  if(store.traits){
    trait.name<-colnames(phy$call$trait.data)
    nms<-paste0(trait.name,'_',c(phy$call$tree$tip.label,seq_len(Nnode(phy))+Ntip(phy)))
    inds<-match(seq_len(max(phy$call$tree$edge)),as.vector(old.edge))
    nas<-is.na(inds)
    old.pos<-match(nms[!nas],names(phy$chains))
    new.nodes<-drop.tree$edge[inds[!nas]]
    ord<-order(new.nodes)
    new.pos<-ord-1+min(old.pos)
    new.nms<-new.nodes[ord]
    tmp<-new.nms<=Ntip(drop.tree)
    new.nms[tmp]<-drop.tree$tip.label[new.nodes[tmp]]
    new.nms<-paste0(trait.name,'_',new.nms)
    drops<-seq_len(length(nms)-length(new.pos))+max(new.pos)
    for(i in c("chains","quantiles","means","diagnostics")){
      if(!is.null(phy[[i]])){
        phy[[i]][,new.pos,]<-phy[[i]][,old.pos,,drop=FALSE]
        names(phy[[i]])[new.pos]<-new.nms
        phy[[i]]<-.add.par.class(phy[[i]][,-drops,,drop=FALSE])
        attr(phy[[i]],'param_type')<-i
      }
    }
  }
  
  phy$call$tree<-drop.tree
  
  if(recalc.bg&any(names(phy$chains)=="bg_rate")){
    #re-calculate background rate
    bg<-get.bg.rate(phy)
    phy$chains[,"bg_rate",]<-bg
    for(i in c("quantiles","means","diagnostics")){
      if(!is.null(phy[[i]])){
        extra.select<-dimanmes(phy[[i]])[[1]]
        phy[[i]][,"bg_rate",]<-.call.op(i,list(chains=bg,sampler.params=1),list(".",extra.select),FALSE)
        if(i=="diagnostics"){
          inds<-extra.select=='inits'
          n.inds<-sum(inds)
          phy[[i]][inds,"bg_rate",]<-rep(get.bg.rate(phy,type="diagnostics",extra.select="^inits$"),each=n.inds)
        }
      }
    }
  }
  if(recalc.devs){
    #re-calculate Rdevs and posterior probabilities...
    R<-get.bg.rate(phy,simplify=FALSE,keep.R=TRUE,
                   remove.trend=remove.trend,geometric=geometric)
    R<-R$R-R$bg_rate
    names(R)<-paste0("Rdev_",seq_along(names(R)))
    tmp<-R
    tmp[tmp==0]<-NA
    phy$post.probs<-.call.op('means',list(chains=tmp>0,sampler.params=1),'.',FALSE)
    pos<-grepl("^Rdev_[1-9][0-9]*",names(phy$chains))
    for(i in c("chains","quantiles","means","diagnostics")){
      if(!is.null(phy[[i]])){
        extra.select<-dimanmes(phy[[i]])[[1]]
        tmp<-.call.op(i,list(chains=R,sampler.params=1),list(".",extra.select),FALSE)
        if(i=="diagnostics"){
          inds<-extra.select=='inits'
          n.inds<-sum(inds)
          inits.R<-get.bg.rate(phy,simplify=FALSE,keep.R=TRUE,
                               remove.trend=TRUE,geometric=geometric,
                               type="diagnostics",extra.select="^inits$")
          inits.R<-inits.R$R-inits.R$bg_rate
          tmp[inds,pos,]<-rep(inits.R,each=n.inds)
        }
        if(has.stored){
          phy[[i]][,pos,]<-tmp
        }else{
          phy[[i]]<-.combine.par(list(phy[[i]],tmp))
        }
      }
    }
  }
  
  phy
}