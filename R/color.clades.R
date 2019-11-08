#'@title Produces an edge-wise color map for a given phylogenetic tree where monophyletic clades are assigned to a given set of
#' colors based on MRCAs 
#'@author Bruce
#' 
#'@description (note: if colors are assigned to nested clades, the color for the smaller clade will override that for the
#' larger clade within the smaller clade, technically allowing paraphyletic clades to be assigned unqiue colors). You already have a working version of this integrated into the 'plot.simBM' method, but I think visualizing simulations in this
#' package will go a lot smoother if you make coloring phylogeny edges a more modular process and make this a separate function;
#' combined with evolve.colors below, as well as planned 'mask.clades/lineages' and 'jitter.colors' functions, users will have
#' flexible, intuitive tools to coloring their simulations the way they see fit. 9/10 update: modified coloration system--moved this function from 'simBM' to here, made alter.cols handle more aspects of color management internally
#' 9/11 thoughts: modify the way this function (and its lineage counterpart) handles base colors such that it can take already
#' made colmaps and recolor specified clades and lineages. 11/5 thoughts: if you make these functions capable of inversion it will render the 'alter functions' obsolete, which will
#' simplify the coloring system to be more intuitive. Nah, I still like the alter.lin/clade idea, honestly--inversion doesn't make predictable sense in a scenario where you are
#' multiple colors to multiple clades/lineages. 11/5 update: completed 9/11 thoughts for color.clades, color.lineages, and alter.colmap. Basic system: generate colmaps using evolve.colors(), color.clades(), or color.lineages(), and further modify colmaps using
#' color.clades()/color.lineages() (set base.cols equal to the colmap to be modified in this case), alter.colmap, and jitter.colors
#
#'@param tree tree of class 'phylo
#'@param MRCAs most recent common ancestor
#'@param alph 
#'


color.clades<-function(tree, MRCAs, cols=palette()[2:(length(MRCAs)+1)], alph=NA, base.cols=palette()[1], base.alph=NA){
  
  if(length(cols)<length(MRCAs)){
    cols<-rep(cols,length.out=length(MRCAs))
  }
  if(length(alph)<length(MRCAs)){
    alph<-rep(alph,length.out=length(MRCAs))
  }
  cols<-alter.cols(cols,alph=alph)
  clade.edges<-lapply(MRCAs,function(nn) which(tree$edge[,2]%in%phytools::getDescendants(tree,nn)))
  inc.ord<-order(lengths(clade.edges),decreasing=T)
  clade.edges<-clade.edges[inc.ord]
  for(c in 1:length(clade.edges)){
    clade.edges[[c]]<-clade.edges[[c]][!(clade.edges[[c]]%in%unlist(clade.edges[-c]))]
  }
  clade.edges[inc.ord]<-clade.edges
  clade.score<-rep(1:length(clade.edges),lengths(clade.edges))
  if(length(base.cols)<nrow(tree$edge)){
    base.cols<-rep(base.cols,length.out=nrow(tree$edge))
  }
  if(length(base.alph)<nrow(tree$edge)){
    base.alph<-rep(base.alph,length.out=nrow(tree$edge))
  }
  colmap<-alter.cols(base.cols,base.alph)
  colmap[unlist(clade.edges)]<-cols[clade.score]
  colmap
}
